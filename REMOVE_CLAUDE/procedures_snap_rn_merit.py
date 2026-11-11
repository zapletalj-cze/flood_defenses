"""
Snap MERIT river network segments to OSM waterways and DTM local minima.

Densifies MERIT line vertices, projects them onto nearest OSM lines,
then refines positions to DTM cell centers with lowest elevation
within a focal window. Includes spike removal and over-shift reversion.

Usage of this class is optional

Author: Jakub Zapletal
Version: 1.0
Date: 04/07/2026
"""

import geopandas as gpd
import numpy as np
from gis_tool.gis import Raster
from gis_tool.helpers import Parameters
from osgeo import gdal
from shapely.geometry import LineString, MultiLineString, Point


class MERIT_snapper:
    def __init__(self):
        self.yaml_path = None
        self.general_parameters = None
        self.dtm_path = None
        self.rn_path = None
        self.rn_osm_path = None
        self.results_dir = None
        self.domain_name = "domain"
        self.area_id = "area"

    def set_general_parameters(self, path_yaml):
        self.yaml_path = path_yaml
        self.general_parameters = Parameters.load_local_parameters(
            self.yaml_path, "general_parameters"
        )

    def set_resolution_parameters(self, resolution):
        self.resolution = resolution

    def set_dtm_parameters(self, path_yaml):
        self.yaml_path = path_yaml
        self.dtm_parameters = Parameters.load_local_parameters(
            self.yaml_path, "dtm_parameters"
        )

    def set_domain_name(self, name):
        self.domain_name = name

    def set_area_id(self, area_id):
        self.area_id = area_id

    def set_mask_path(self, path):
        self.mask_path = path

    def set_results_dir(self, path):
        self.results_dir = path

    def set_dtm_path(self, dtm_paths):
        self.dtm_path = dtm_paths

    def set_rn_path(self, rn_path):
        self.rn_path = rn_path

    def set_rn_osm_path(self, rn_osm_path):
        self.rn_osm_path = rn_osm_path

    def _dtm_meta(self):
        ds = gdal.Open(self.dtm_path)
        if ds is None:
            raise FileNotFoundError(f"Could not open DTM: {self.dtm_path}")
        band = ds.GetRasterBand(1)
        gt = ds.GetGeoTransform()
        return {
            "ds": ds,
            "band": band,
            "gt": gt,
            "origin_x": float(gt[0]),
            "origin_y": float(gt[3]),
            "dx": float(gt[1]),
            "dy": float(abs(gt[5])),
            "xsize": int(ds.RasterXSize),
            "ysize": int(ds.RasterYSize),
            "nodata": band.GetNoDataValue(),
            "epsg": Raster.get_raster_info(self.dtm_path, ["projection_epsg"])[
                "projection_epsg"
            ],
        }

    @staticmethod
    def _ensure_crs(gdf, target_epsg):
        if gdf.crs is None:
            raise ValueError("Input vector has undefined CRS.")
        if gdf.crs.to_epsg() != int(target_epsg):
            return gdf.to_crs(epsg=int(target_epsg))
        return gdf

    @staticmethod
    def _densify_coords(coords, max_segment_length):
        if len(coords) < 2 or max_segment_length is None or max_segment_length <= 0:
            return coords

        dense = [coords[0]]
        for i in range(1, len(coords)):
            x0, y0 = coords[i - 1]
            x1, y1 = coords[i]
            dx = x1 - x0
            dy = y1 - y0
            dist = float(np.hypot(dx, dy))
            if dist > max_segment_length:
                n_parts = int(np.ceil(dist / max_segment_length))
                for k in range(1, n_parts):
                    t = k / n_parts
                    dense.append((x0 + dx * t, y0 + dy * t))
            dense.append((x1, y1))
        return dense

    @staticmethod
    def _explode_merit_vertices(gdf_merit, max_segment_length=200.0):
        records = []
        meta = []

        for feature_idx, geom in gdf_merit.geometry.items():
            if geom is None or geom.is_empty:
                meta.append((feature_idx, None, []))
                continue

            if geom.geom_type == "LineString":
                parts = [geom]
                geom_type = "LineString"
            elif geom.geom_type == "MultiLineString":
                parts = list(geom.geoms)
                geom_type = "MultiLineString"
            else:
                meta.append((feature_idx, None, []))
                continue

            part_sizes = []
            for part_idx, part in enumerate(parts):
                coords = list(part.coords)
                coords = MERIT_snapper._densify_coords(coords, max_segment_length)
                part_sizes.append(len(coords))
                for vertex_idx, coord in enumerate(coords):
                    records.append(
                        {
                            "feature_idx": feature_idx,
                            "part_idx": part_idx,
                            "vertex_idx": vertex_idx,
                            "geometry": Point(float(coord[0]), float(coord[1])),
                        }
                    )
            meta.append((feature_idx, geom_type, part_sizes))

        vertices = gpd.GeoDataFrame(records, geometry="geometry", crs=gdf_merit.crs)
        vertices["vertex_uid"] = np.arange(len(vertices), dtype=np.int64)
        return vertices, meta

    @staticmethod
    def _rebuild_merit_geometry(gdf_merit, snapped_vertices, geometry_meta):
        by_feature = {}
        grouped = snapped_vertices.sort_values(
            ["feature_idx", "part_idx", "vertex_idx"]
        ).groupby(["feature_idx", "part_idx"])
        for (fidx, pidx), sub in grouped:
            coords = [(pt.x, pt.y) for pt in sub.geometry.values]
            if len(coords) >= 2:
                by_feature.setdefault(fidx, {})[pidx] = LineString(coords)

        out = gdf_merit.copy()
        for feature_idx, geom_type, part_sizes in geometry_meta:
            if geom_type is None:
                continue
            parts = by_feature.get(feature_idx, {})
            if not parts:
                continue
            ordered = [parts[k] for k in sorted(parts.keys()) if k in parts]
            if not ordered:
                continue
            if geom_type == "LineString":
                out.at[feature_idx, "geometry"] = ordered[0]
            else:
                out.at[feature_idx, "geometry"] = MultiLineString(ordered)
        return out

    @staticmethod
    def _project_to_line(point, line):
        if line is None or line.is_empty:
            return point
        try:
            d = line.project(point)
            return line.interpolate(d)
        except (AttributeError, TypeError, ValueError):
            return point

    @staticmethod
    def _xy_to_rc(x, y, origin_x, origin_y, dx, dy):
        col = int(np.floor((x - origin_x) / dx))
        row = int(np.floor((origin_y - y) / dy))
        return row, col

    @staticmethod
    def _rc_to_xy(row, col, origin_x, origin_y, dx, dy):
        x = origin_x + (col + 0.5) * dx
        y = origin_y - (row + 0.5) * dy
        return x, y

    def _remove_spike_vertices(
        self,
        snapped_vertices,
        spike_angle_deg=75.0,
        spike_length_max=120.0,
    ):
        records = []
        grouped = snapped_vertices.sort_values(
            ["feature_idx", "part_idx", "vertex_idx"]
        ).groupby(["feature_idx", "part_idx"])

        for (fidx, pidx), sub in grouped:
            coords = [(pt.x, pt.y) for pt in sub.geometry.values]
            cleaned = [(float(x), float(y)) for x, y in coords]
            for _ in range(6):
                if len(cleaned) < 3:
                    break
                _chg = False
                _out = [cleaned[0]]
                for i in range(1, len(cleaned) - 1):
                    _a = np.asarray(cleaned[i - 1], dtype=np.float64)
                    _b = np.asarray(cleaned[i], dtype=np.float64)
                    _c = np.asarray(cleaned[i + 1], dtype=np.float64)
                    _ab = float(np.hypot(*(_b - _a)))
                    _bc = float(np.hypot(*(_c - _b)))
                    if _ab == 0 or _bc == 0:
                        _chg = True
                        continue
                    _v1 = _a - _b
                    _v2 = _c - _b
                    _denom = float(np.hypot(*_v1) * np.hypot(*_v2))
                    if _denom == 0:
                        _chg = True
                        continue
                    _cosang = float(np.clip(np.dot(_v1, _v2) / _denom, -1.0, 1.0))
                    _ang = float(np.degrees(np.arccos(_cosang)))
                    if _ang < float(spike_angle_deg) and min(_ab, _bc) <= float(
                        spike_length_max
                    ):
                        _chg = True
                        continue
                    if (
                        _ang >= 175.0
                        and min(_ab, _bc) <= float(spike_length_max) * 0.25
                    ):
                        _chg = True
                        continue
                    _out.append((float(_b[0]), float(_b[1])))
                _out.append(cleaned[-1])
                cleaned = _out
                if not _chg:
                    break
            if len(cleaned) < 2:
                cleaned = coords[:2]

            for vidx, (x, y) in enumerate(cleaned):
                records.append(
                    {
                        "feature_idx": int(fidx),
                        "part_idx": int(pidx),
                        "vertex_idx": int(vidx),
                        "geometry": Point(float(x), float(y)),
                    }
                )

        out = gpd.GeoDataFrame(records, geometry="geometry", crs=snapped_vertices.crs)
        out["vertex_uid"] = np.arange(len(out), dtype=np.int64)
        return out

    @staticmethod
    def _revert_over_shifted_parts(
        snapped_vertices,
        original_vertices,
        max_part_shift_p95=25.0,
        max_part_shift_max=40.0,
    ):
        if snapped_vertices.empty or original_vertices.empty:
            return snapped_vertices

        out = snapped_vertices.copy()
        dx = out.geometry.x.values - original_vertices.geometry.x.values
        dy = out.geometry.y.values - original_vertices.geometry.y.values
        out["_shift"] = np.hypot(dx, dy)

        part_stats = (
            out.groupby(["feature_idx", "part_idx"])["_shift"]
            .agg(shift_max="max", shift_p95=lambda s: float(np.percentile(s, 95)))
            .reset_index()
        )
        bad_parts = part_stats[
            (part_stats["shift_p95"] > float(max_part_shift_p95))
            | (part_stats["shift_max"] > float(max_part_shift_max))
        ][["feature_idx", "part_idx"]]

        if not bad_parts.empty:
            bad_idx = set(zip(bad_parts["feature_idx"], bad_parts["part_idx"]))
            mask = out.apply(
                lambda r: (r["feature_idx"], r["part_idx"]) in bad_idx,
                axis=1,
            )
            out.loc[mask, "geometry"] = original_vertices.loc[mask, "geometry"].values

        out = out.drop(columns=["_shift"])
        return out

    def snap_merit_network(self):
        search_radius = 300.0
        max_osm_offset = 150.0
        max_osm_guided_shift_distance = 180.0
        max_shift_distance_without_osm = 0.0
        focal_window_size = 9
        max_shift_distance = 120.0
        max_segment_length = 20.0

        if self.rn_path is None:
            raise ValueError("MERIT path not set. Use set_rn_path().")
        if self.dtm_path is None:
            raise ValueError("DTM path not set. Use set_dtm_path().")

        osm_path = self.rn_osm_path
        if osm_path is None:
            raise ValueError("OSM path not set. Use set_rn_osm_path().")

        if focal_window_size < 1 or focal_window_size % 2 == 0:
            raise ValueError("focal_window_size must be odd and >= 1.")

        dtm_meta = self._dtm_meta()
        band = dtm_meta["band"]
        half_window = focal_window_size // 2

        gdf_merit = gpd.read_file(self.rn_path)
        gdf_osm = gpd.read_file(osm_path)

        gdf_merit = self._ensure_crs(gdf_merit, dtm_meta["epsg"])
        gdf_osm = self._ensure_crs(gdf_osm, dtm_meta["epsg"])

        gdf_merit_check = gdf_merit[
            gdf_merit.geometry.notnull() & ~gdf_merit.geometry.is_empty
        ].copy()
        gdf_merit_check = gdf_merit_check.explode(index_parts=False).reset_index(
            drop=True
        )
        gdf_merit_check = gdf_merit_check[
            gdf_merit_check.geometry.geom_type.isin(["LineString", "MultiLineString"])
        ].copy()
        gdf_merit_check = gdf_merit_check.explode(index_parts=False).reset_index(
            drop=True
        )

        gdf_osm = gdf_osm[
            gdf_osm.geometry.notnull() & ~gdf_osm.geometry.is_empty
        ].copy()
        gdf_osm = gdf_osm.explode(index_parts=False).reset_index(drop=True)
        gdf_osm = gdf_osm[
            gdf_osm.geometry.geom_type.isin(["LineString", "MultiLineString"])
        ].copy()
        gdf_osm = gdf_osm.explode(index_parts=False).reset_index(drop=True)

        merit_bbox = gdf_merit_check.total_bounds  # minx, miny, maxx, maxy
        from shapely.geometry import box as shapely_box

        merit_extent = shapely_box(
            merit_bbox[0] - search_radius,
            merit_bbox[1] - search_radius,
            merit_bbox[2] + search_radius,
            merit_bbox[3] + search_radius,
        )
        gdf_osm = gdf_osm[gdf_osm.intersects(merit_extent)].copy()

        if gdf_merit_check.empty:
            raise RuntimeError("MERIT file has no line segments to snap.")
        if gdf_osm.empty:
            raise RuntimeError("OSM file has no line segments to snap.")

        merit_check = gpd.sjoin_nearest(
            gdf_merit_check[["geometry"]],
            gdf_osm[["geometry"]],
            how="left",
            max_distance=float(search_radius),
            distance_col="seg_dist",
        )

        missing_merit = int(merit_check["index_right"].isna().sum())
        if missing_merit > 0:
            raise ValueError(
                "Segment coverage check failed between MERIT and OSM files. "
                f"MERIT segments without OSM match: {missing_merit}, "
                f"search_radius={search_radius}."
            )

        vertices, geometry_meta = self._explode_merit_vertices(
            gdf_merit,
            max_segment_length=float(max_segment_length),
        )
        if vertices.empty:
            raise RuntimeError("No MERIT vertices extracted.")

        vertices_orig = vertices.copy()
        candidate = vertices.copy()
        candidate["used_osm"] = False

        if not gdf_osm.empty:
            v_join = vertices[
                ["vertex_uid", "feature_idx", "part_idx", "vertex_idx", "geometry"]
            ].copy()
            o_join = gdf_osm[["geometry"]].copy().reset_index(drop=True)
            o_join["osm_idx"] = o_join.index

            joined = gpd.sjoin_nearest(
                v_join,
                o_join,
                how="left",
                max_distance=float(search_radius),
                distance_col="osm_dist",
            )
            joined = joined.dropna(subset=["osm_idx"]).copy()
            if not joined.empty:
                joined = joined.sort_values(["vertex_uid", "osm_dist", "osm_idx"])
                joined = joined.drop_duplicates(subset=["vertex_uid"], keep="first")

                joined = joined[joined["osm_dist"] <= float(max_osm_offset)].copy()

                if not joined.empty:
                    osm_geom = o_join.set_index("osm_idx").geometry
                    proj_points = []
                    for _, row in joined.iterrows():
                        p = vertices.at[int(row["vertex_uid"]), "geometry"]
                        line = osm_geom.loc[int(row["osm_idx"])]
                        proj_points.append(self._project_to_line(p, line))
                    joined["candidate_geom"] = proj_points

                    uid_to_point = dict(
                        zip(joined["vertex_uid"].astype(int), joined["candidate_geom"])
                    )
                    mask = candidate["vertex_uid"].isin(uid_to_point.keys())
                    candidate.loc[mask, "geometry"] = candidate.loc[
                        mask, "vertex_uid"
                    ].map(uid_to_point)
                    candidate.loc[mask, "used_osm"] = True

        p_orig = list(vertices_orig.geometry.values)
        p_candidates = list(candidate.geometry.values)
        used_osm_flags = list(candidate["used_osm"].values)

        stage_points = []
        for p0, p_candidate, used_osm in zip(p_orig, p_candidates, used_osm_flags):
            if bool(used_osm):
                _dx = p_candidate.x - p0.x
                _dy = p_candidate.y - p0.y
                _d = float(np.hypot(_dx, _dy))
                if _d <= float(max_osm_guided_shift_distance) or _d == 0:
                    p_stage = p_candidate
                else:
                    _ratio = float(max_osm_guided_shift_distance) / _d
                    p_stage = Point(p0.x + _dx * _ratio, p0.y + _dy * _ratio)
            else:
                p_stage = p0
            stage_points.append(p_stage)

        snapped_geoms = []
        for p_stage, used_osm in zip(stage_points, used_osm_flags):
            _row, _col = self._xy_to_rc(
                p_stage.x,
                p_stage.y,
                dtm_meta["origin_x"],
                dtm_meta["origin_y"],
                dtm_meta["dx"],
                dtm_meta["dy"],
            )
            _row = int(np.clip(_row, 0, dtm_meta["ysize"] - 1))
            _col = int(np.clip(_col, 0, dtm_meta["xsize"] - 1))
            _x0 = max(0, _col - half_window)
            _y0 = max(0, _row - half_window)
            _x1 = min(dtm_meta["xsize"], _col + half_window + 1)
            _y1 = min(dtm_meta["ysize"], _row + half_window + 1)
            _arr = band.ReadAsArray(_x0, _y0, _x1 - _x0, _y1 - _y0)
            if _arr is None or _arr.size == 0:
                _bx, _by = self._rc_to_xy(
                    _row,
                    _col,
                    dtm_meta["origin_x"],
                    dtm_meta["origin_y"],
                    dtm_meta["dx"],
                    dtm_meta["dy"],
                )
            else:
                _arr = np.asarray(_arr)
                _nodata = dtm_meta["nodata"]
                if _nodata is not None:
                    _valid = _arr != _nodata
                    if not np.any(_valid):
                        _best_r, _best_c = _row, _col
                    else:
                        _masked = np.where(_valid, _arr, np.inf)
                        _i = int(np.argmin(_masked))
                        _rr, _cc = np.unravel_index(_i, _masked.shape)
                        _best_r = _y0 + _rr
                        _best_c = _x0 + _cc
                else:
                    _i = int(np.argmin(_arr))
                    _rr, _cc = np.unravel_index(_i, _arr.shape)
                    _best_r = _y0 + _rr
                    _best_c = _x0 + _cc
                _bx, _by = self._rc_to_xy(
                    _best_r,
                    _best_c,
                    dtm_meta["origin_x"],
                    dtm_meta["origin_y"],
                    dtm_meta["dx"],
                    dtm_meta["dy"],
                )
            p_dtm = Point(_bx, _by)
            local_max_shift = float(max_shift_distance)
            if bool(used_osm):
                local_max_shift = min(
                    local_max_shift, float(max_osm_guided_shift_distance)
                )
            else:
                local_max_shift = min(
                    local_max_shift, float(max_shift_distance_without_osm)
                )
            _dx = p_dtm.x - p_stage.x
            _dy = p_dtm.y - p_stage.y
            _d = float(np.hypot(_dx, _dy))
            if _d <= local_max_shift or _d == 0:
                p_final = p_dtm
            else:
                _ratio = local_max_shift / _d
                p_final = Point(p_stage.x + _dx * _ratio, p_stage.y + _dy * _ratio)
            snapped_geoms.append(p_final)

        snapped_vertices = vertices.copy()
        snapped_vertices.geometry = gpd.GeoSeries(snapped_geoms, crs=vertices.crs)

        min_angle_deg = 35.0
        max_passes = 2
        records = []
        grouped = snapped_vertices.sort_values(
            ["feature_idx", "part_idx", "vertex_idx"]
        ).groupby(["feature_idx", "part_idx"])

        for (fidx, pidx), sub in grouped:
            pts = [(float(pt.x), float(pt.y)) for pt in sub.geometry.values]
            if len(pts) < 3:
                for vidx, (x, y) in enumerate(pts):
                    records.append(
                        {
                            "feature_idx": int(fidx),
                            "part_idx": int(pidx),
                            "vertex_idx": int(vidx),
                            "geometry": Point(x, y),
                        }
                    )
                continue

            for _ in range(int(max_passes)):
                if len(pts) < 3:
                    break
                changed = False
                out_pts = [pts[0]]

                for i in range(1, len(pts) - 1):
                    a = np.asarray(pts[i - 1], dtype=np.float64)
                    b = np.asarray(pts[i], dtype=np.float64)
                    c = np.asarray(pts[i + 1], dtype=np.float64)
                    v1 = a - b
                    v2 = c - b
                    n1 = float(np.hypot(*v1))
                    n2 = float(np.hypot(*v2))
                    if n1 == 0 or n2 == 0:
                        changed = True
                        continue

                    cosang = float(np.dot(v1, v2) / (n1 * n2))
                    cosang = float(np.clip(cosang, -1.0, 1.0))
                    ang = float(np.degrees(np.arccos(cosang)))

                    if ang < float(min_angle_deg):
                        changed = True
                        continue

                    out_pts.append((float(b[0]), float(b[1])))

                out_pts.append(pts[-1])
                pts = out_pts
                if not changed:
                    break

            if len(pts) < 2:
                pts = [(float(pt.x), float(pt.y)) for pt in sub.geometry.values[:2]]

            for vidx, (x, y) in enumerate(pts):
                records.append(
                    {
                        "feature_idx": int(fidx),
                        "part_idx": int(pidx),
                        "vertex_idx": int(vidx),
                        "geometry": Point(x, y),
                    }
                )

        snapped_vertices = gpd.GeoDataFrame(
            records,
            geometry="geometry",
            crs=snapped_vertices.crs,
        )
        snapped_vertices["vertex_uid"] = np.arange(
            len(snapped_vertices), dtype=np.int64
        )

        snapped_vertices = self._remove_spike_vertices(snapped_vertices)

        gdf_out = self._rebuild_merit_geometry(
            gdf_merit, snapped_vertices, geometry_meta
        )

        gdf_out.to_file(self.rn_path, driver="GPKG", index=False)

        return self.rn_path
