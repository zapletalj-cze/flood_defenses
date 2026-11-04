# Standard library imports
import ctypes
import datetime
import heapq
import json
import logging
import multiprocessing
import os
import shutil
import sys
import time
import tomllib
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures.process import BrokenProcessPool
from math import acos, ceil, degrees, floor, sqrt
from os.path import exists, join

import geopandas as gpd
import numpy as np
import pandas as pd
import pyflwdir
import pyogrio
import shapely
from numba import njit
from osgeo import gdal, gdal_array, ogr, osr
from scipy import ndimage
from shapely.geometry import LineString, Point, Polygon, box
from shapely.ops import linemerge, nearest_points, unary_union
from shapely.wkt import loads

# glob_params = Parameters()
gdal.UseExceptions()
gdal.SetCacheMax(4 * pow(1024, 3))
log = logging.getLogger("gis")


class Files:
    @staticmethod
    def set_workspace(path, folders):
        if not exists(path):
            os.makedirs(path)
        for folder in folders:
            folder_path = join(path, folder)
            if not exists(folder_path):
                os.makedirs(folder_path)

    @staticmethod
    def get_value_from_dict(item, key, default=None):
        return item.get(key, default)

    def flush_temp(self):
        try:
            temp_path = os.path.join(self.model_root, "temp")
            if os.path.exists(temp_path):
                shutil.rmtree(temp_path)
                print(f"Temp folder '{temp_path}' has been flushed.")
            else:
                print(f"Temp folder '{temp_path}' does not exist.")
        except Exception as e:
            print(f"Failed to flush temp folder: {e}")

    def flush_temp_subdir(self, subdir_path):
        try:
            temp_path = os.path.join(self.model_root, "temp")
            target_path = os.path.join(temp_path, subdir_path)

            if os.path.exists(target_path) and os.path.isdir(target_path):
                shutil.rmtree(target_path)
                print(f"Subdirectory '{subdir_path}' in temp has been flushed.")
            else:
                print(f"Subdirectory '{subdir_path}' does not exist in temp.")
        except Exception as e:
            print(f"Failed to flush subdirectory '{subdir_path}' in temp: {e}")

    @staticmethod
    def get_root_path():
        current_path = os.path.dirname(os.path.abspath(__file__))
        while os.path.basename(current_path) != "IF-HydroSim":
            current_path = os.path.dirname(current_path)
        return current_path

    @staticmethod
    def split_list(input_list, num_sublists):
        """
        Split a list into a defined number of sublists.
        """
        avg = len(input_list) / float(num_sublists)
        out_list = []
        index_l = 0

        while index_l < len(input_list):
            out_list.append(input_list[int(index_l) : int(index_l + avg)])
            index_l += avg

        return out_list

    @staticmethod
    def split_dict(input_dict, num_sublists):
        """
        Split a dictionary into a defined number of sublists of (key, value) pairs.
        """
        items = list(input_dict.items())
        avg = len(items) / float(num_sublists)
        sublists = []
        index_l = 0

        while index_l < len(items):
            sublists.append(items[int(index_l) : int(index_l + avg)])
            index_l += avg

        list_dict = []
        for sublist in sublists:
            list_dict.append(dict(sublist))
        return list_dict


class Raster:
    def __init__(self):
        self.extent_epsg = None
        self.snap_raster = None
        self.snap = False

    @staticmethod
    def to_array(path, band_number=1):
        """
        Read raster as TIF to numpy array
        :param path: path to raster file
        :return: numpy array
        """
        ds = gdal.Open(path)
        if ds is None:
            raise FileNotFoundError(f"Could not open input raster: {path}")
        array = ds.GetRasterBand(band_number).ReadAsArray()
        return array

    @staticmethod
    def to_array_mem(ds_raster, band_number=1):
        """
        Read raster as MEM raster to numpy array
        :param ds_raster: MEM raster dataset
        :return: numpy array
        """
        array = ds_raster.GetRasterBand(band_number).ReadAsArray()
        return array

    @staticmethod
    def save_raster(output_path, ds):
        """
        Saves GDAL Dataset to file as TIF
        ::param output_path: path to output raster file
        :param ds: existing GDAL Dataset to save
        """
        driver = gdal.GetDriverByName("GTiff")
        new_ds = driver.Create(
            output_path,
            ds.RasterXSize,
            ds.RasterYSize,
            ds.RasterCount,
            ds.GetRasterBand(1).DataType,
            options=["COMPRESS=ZSTD", "PREDICTOR=2", "TILED=YES"],
        )
        new_ds.SetGeoTransform(ds.GetGeoTransform())
        new_ds.SetProjection(ds.GetProjection())
        for i in range(1, ds.RasterCount + 1):
            band = new_ds.GetRasterBand(i)
            band.WriteArray(ds.GetRasterBand(i).ReadAsArray())
        new_ds.FlushCache()
        new_ds = None

    def clip_by_vector(
        in_raster: str,
        out_raster: str,
        vector_file: str,
        epsg_out: int,
        cell_size_out,
        data_type,
        CPU_AVAILABLE=16,
        no_data=-9999,
        bounds=None,
        interpolation="nearest",
    ):
        """
        Clips raster by vector file and reprojects to new projection and resolution
        :param in_raster: input raster file
        :param out_raster: output raster file
        :param vector_file: vector file for clipping
        :param epsg_out: EPSG code of output projection
        :param cell_size_out: output cell size
        :param data_type: output data type in GDAL format
        """
        ds_in = gdal.Open(in_raster)
        if ds_in is None:
            raise FileNotFoundError(f"Could not open input raster: {in_raster}")
        srs_in = osr.SpatialReference()
        srs_in.ImportFromWkt(ds_in.GetProjection())
        srs_out = osr.SpatialReference()
        srs_out.ImportFromEPSG(epsg_out)
        options = gdal.WarpOptions(
            format="GTiff",
            cutlineDSName=vector_file,
            xRes=cell_size_out,
            yRes=cell_size_out,
            dstSRS=srs_out,
            srcSRS=srs_in,
            outputBounds=bounds,
            resampleAlg=interpolation,
            srcNodata=ds_in.GetRasterBand(1).GetNoDataValue(),
            dstNodata=no_data,
            outputType=data_type,
            creationOptions=[
                "COMPRESS=ZSTD",
                # "PREDICTOR=2",
                "TILED=YES",
                "BIGTIFF=YES",
            ],
            warpOptions=[
                f"NUM_THREADS={CPU_AVAILABLE}",
            ],
        )
        gdal.Warp(out_raster, ds_in, options=options)
        ds_in = None

    def clip_by_vector_mem(
        in_raster: str,
        vector_file,
        epsg_out: int,
        cell_size_out,
        data_type,
        CPU_AVAILABLE=1,
        no_data=-9999,
        bounds=None,
    ):
        """
        Clips raster by vector file and reprojects to new projection and resolution
        :param in_raster: input raster file
        :param vector_file: vector file path (str) or OGR datasource object
        :param epsg_out: EPSG code of output projection
        :param cell_size_out: output cell size
        :param data_type: output data type in GDAL format
        """
        ds_in = gdal.Open(in_raster)
        if ds_in is None:
            raise FileNotFoundError(f"Could not open input raster: {in_raster}")
        srs_in = osr.SpatialReference()
        srs_in.ImportFromWkt(ds_in.GetProjection())
        srs_out = osr.SpatialReference()
        srs_out.ImportFromEPSG(epsg_out)

        # Handle both file paths and OGR datasource objects
        if isinstance(vector_file, str):
            warp_options = gdal.WarpOptions(
                format="MEM",
                cutlineDSName=vector_file,
                xRes=cell_size_out,
                yRes=cell_size_out,
                dstSRS=srs_out,
                srcSRS=srs_in,
                outputBounds=bounds,
                srcNodata=ds_in.GetRasterBand(1).GetNoDataValue(),
                dstNodata=no_data,
                outputType=data_type,
                multithread=True,
            )
        else:
            warp_options = gdal.WarpOptions(
                format="MEM",
                cutlineLayer=vector_file,
                xRes=cell_size_out,
                yRes=cell_size_out,
                dstSRS=srs_out,
                srcSRS=srs_in,
                outputBounds=bounds,
                srcNodata=ds_in.GetRasterBand(1).GetNoDataValue(),
                dstNodata=no_data,
                outputType=data_type,
                multithread=True,
            )
        mem_ds = gdal.Warp("", ds_in, options=warp_options)
        ds_in = None
        return mem_ds

    @staticmethod
    def from_array_mem(
        array: np.array,
        output_raster: str | None = None,
        source_raster: str | None = None,
        geotransform: list | None = None,
        size_xy: list | None = None,
        epsg_out: int | None = None,
        data_type=None,
        nodata=None,
        raster_driver="GTiff",
    ):
        """
        Converts array to raster file or in memory raster
        :param array: numpy array
        :param output_raster: path to output raster file
        :param source_raster: path to source raster file for configuration
        :param geotransform: geotransform of output raster
        :param epsg_out: projection of output raster as EPSG code`
        :param data_type: data type of output raster, use gdal_dt from helpers.py
        :param no_data: no data value of output raster
        :param raster_driver: driver for output raster, supported GTiff and MEM
        """
        # Set projection for output raster
        projection = osr.SpatialReference()
        if epsg_out is None:
            try:
                epsg_out = Raster().extent_epsg
            except Exception:
                epsg_out = Raster.get_raster_info(source_raster, ["projection_epsg"])[
                    "projection_epsg"
                ]
        projection.ImportFromEPSG(epsg_out)
        projection = projection.ExportToWkt()
        # Get necessary parameters from source raster
        info_required = []
        if geotransform is None:
            info_required.append("geotransform")
        if data_type is None:
            info_required.append("data_type")
        if nodata is None:
            info_required.append("nodata")
        if size_xy is None:
            info_required.append("size")
        if len(info_required) > 0:
            if source_raster is None:
                raise ValueError(
                    "Path to output raster file must be provided if configuration is not set."
                )
        info = Raster.get_raster_info(source_raster, info_required)
        for key, value in info.items():
            if key == "geotransform":
                geotransform = value
            elif key == "data_type":
                data_type = value
            elif key == "nodata":
                nodata = value if value is not None else -9999
            elif key == "size":
                size_xy = value
        if geotransform is None:
            geotransform = info["geotransform"]
        if data_type is None:
            data_type = info["data_type"]
        if nodata is None:
            nodata = info["nodata"]
        if size_xy is None:
            size_xy = info["size"]
        # Create output raster
        driver = gdal.GetDriverByName(raster_driver)
        if raster_driver == "MEM":
            ds = driver.Create("", size_xy[0], size_xy[1], 1, data_type)
        else:
            ds = driver.Create(
                output_raster,
                size_xy[0],
                size_xy[1],
                1,
                data_type,
                options=["COMPRESS=ZSTD", "PREDICTOR=2", "TILED=YES"],
            )
        ds.SetGeoTransform(geotransform)
        ds.SetProjection(projection)
        ds.GetRasterBand(1).WriteArray(array)
        ds.GetRasterBand(1).SetNoDataValue(nodata)
        if raster_driver == "MEM":
            return ds
        ds.FlushCache()
        ds = None

    def create_snap_raster(self):
        """
        Create a snap raster for the extent, 2x2 pixels
        :return: path to snap raster
        """
        path = join(self.model_root, "config", "grid")
        Files.set_workspace(join(self.model_root, "config"), ["grid"])
        snap_raster = join(path, "snap_raster.tif")
        if not exists(snap_raster):
            xmin = self.extent[0]
            ymin = self.extent[2]
            cell_size = self.default_resolution
            data = np.ones((2, 2), dtype=np.uint8)
            driver = gdal.GetDriverByName("GTiff")
            ds = driver.Create(snap_raster, 2, 2, 1, gdal.GDT_Byte)
            ds.SetGeoTransform(
                (xmin, cell_size, 0, ymin + 2 * cell_size, 0, -cell_size)
            )
            srs = osr.SpatialReference()
            srs.ImportFromEPSG(self.extent_epsg)
            ds.SetProjection(srs.ExportToWkt())
            ds.GetRasterBand(1).WriteArray(data)
            ds.FlushCache()
            ds = None
        self.snap_raster = snap_raster

    @staticmethod
    def get_raster_nodata(in_raster: str, band_number: int = 1):
        """
        Get raster nodata value
        :param in_raster: input raster file
        :param band_number: band number to get the nodata value from
        :return: nodata value
        """
        ds = gdal.Open(in_raster)
        if ds is None:
            raise FileNotFoundError(f"Could not open input raster: {in_raster}")
        band = ds.GetRasterBand(band_number)
        no_data = band.GetNoDataValue()
        return no_data

    @staticmethod
    def raster_to_vect(
        in_raster,
        out_vector,
        layer_name,
        raster_val_name,
        epsg,
        filter_values: list | None = None,
        area_threshold=None,
    ):
        """
        Converts vecror to raster using Polygonize function, expected:
        - single band raster file
        - categorical or binary raster
        :in raster: input raster file
        :out_vector: output vector file
        :raster_val_name: name of the raster value field in output vector
        :layer_name: same as basename of vector file
        :epsg: EPSG code of output vector
        """
        # Set vector driver for final output
        fe = os.path.splitext(out_vector)[1]
        if fe == ".gpkg":
            drv = ogr.GetDriverByName("GPKG")
        elif fe == ".shp":
            if len(raster_val_name) > 8:
                raster_val_name = raster_val_name[:8]
                print(
                    f"Field name too long, truncating to 8 characters: {raster_val_name}"
                )
            drv = ogr.GetDriverByName("ESRI Shapefile")
        elif fe == ".geojson":
            drv = ogr.GetDriverByName("GeoJSON")
        else:
            raise ValueError(f"Unsupported output vector format: {fe}")

        # Set inmemory driver
        dst_layername = layer_name
        drv_mem = ogr.GetDriverByName("Memory")  # Still using Memory instead of MEM
        dst_ds = drv_mem.CreateDataSource("in_memory")

        # If input raster is INMEMORY, use datasource directly
        try:
            drv_raster = in_raster.GetDriver()
            drivername = drv_raster.ShortName
        except Exception:
            pass

        if drivername == "MEM":
            ds_raster = in_raster
        else:
            ds_raster = gdal.Open(in_raster)
        band = ds_raster.GetRasterBand(1)
        sp_ref = osr.SpatialReference()
        sp_ref.SetFromUserInput(f"EPSG:{epsg}")

        dst_layer = dst_ds.CreateLayer(dst_layername, srs=sp_ref)

        # Polygonize inmemory, way faster than direct writing
        fld = ogr.FieldDefn(raster_val_name, ogr.OFTInteger)
        dst_layer.CreateField(fld)
        dst_field = dst_layer.GetLayerDefn().GetFieldIndex(raster_val_name)
        gdal.Polygonize(band, None, dst_layer, dst_field, [], callback=None)
        if filter_values is not None and len(filter_values) == 1:
            filter_values = 2 * filter_values  # TEMPORARY FIX?
        dst_layer.SetAttributeFilter(f"{raster_val_name} IN {tuple(filter_values)}")

        # filter on area if set
        if area_threshold is not None:
            dst_layer.CreateField(ogr.FieldDefn("area", ogr.OFTReal))
            # option A using transaction per row
            dst_layer.StartTransaction()
            for feature in dst_layer:
                geom = feature.GetGeometryRef()
                area = geom.GetArea()
                feature.SetField("area", area)
                dst_layer.SetFeature(feature)
            dst_layer.CommitTransaction()
            dst_layer.SetAttributeFilter(f"area > {area_threshold}")
        gdal.VectorTranslate(out_vector, dst_ds)
        dst_ds = None
        ds_raster = None
        return True

    @staticmethod
    def raster_to_mem_vector(
        in_raster,
        raster_val_name,
        epsg,
        filter_values: list | None = None,
        area_threshold=None,
    ):
        """
        Converts vecror to raster using Polygonize function, expected:
        - single band raster file
        - categorical or binary raster
        :in raster: input raster file
        :raster_val_name: name of the raster value field in output vector
        :layer_name: same as basename of vector file
        :epsg: EPSG code of output vector
        return ogr memory vector
        """
        # Set vector driver for final output
        drv = ogr.GetDriverByName("MEM")
        # Set inmemory driver
        dst_ds = drv.CreateDataSource("in_memory")
        dst_layername = "MEM layer"
        try:
            drv_rast = in_raster.GetDriver()
            drivername = drv_rast.ShortName
        except Exception:
            pass

        if drivername == "MEM":
            ds = in_raster
        else:
            ds = gdal.Open(in_raster)
        band = ds.GetRasterBand(1)
        sp_ref = osr.SpatialReference()
        sp_ref.SetFromUserInput(f"EPSG:{epsg}")
        dst_layer = dst_ds.CreateLayer(dst_layername, srs=sp_ref)

        # Polygonize inmemory, way faster than direct writing
        fld = ogr.FieldDefn(raster_val_name, ogr.OFTInteger)
        dst_layer.CreateField(fld)
        dst_field = dst_layer.GetLayerDefn().GetFieldIndex(raster_val_name)
        gdal.Polygonize(band, None, dst_layer, dst_field, [], callback=None)
        if filter_values is not None and len(filter_values) == 1:
            filter_values = 2 * filter_values  # TEMPORARY FIX?
        dst_layer.SetAttributeFilter(f"{raster_val_name} IN {tuple(filter_values)}")

        # filter on area if set
        if area_threshold is not None:
            dst_layer.CreateField(ogr.FieldDefn("area", ogr.OFTReal))
            # option A using transaction per row
            dst_layer.StartTransaction()
            for feature in dst_layer:
                geom = feature.GetGeometryRef()
                area = geom.GetArea()
                feature.SetField("area", area)
                dst_layer.SetFeature(feature)
            dst_layer.CommitTransaction()
            dst_layer.SetAttributeFilter(f"area > {area_threshold}")
        ds = None
        return dst_ds

    ### VRT manipulation
    @staticmethod
    def raster_list_to_vrt(in_rasters: list, output_vrt, no_data, resolution="highest"):
        """
        Creates a VRT file from a list of raster files
        :param in_rasters: list of raster files
        :param output_vrt: path to output VRT file
        :param no_data: nodata value to set in the VRT file and source rasters
        :param resolution: resolution of the output VRT file, see gdal.BuildVRTOptions
        """
        options = gdal.BuildVRTOptions(
            resolution=resolution, srcNodata=no_data, VRTNodata=no_data
        )
        gdal.BuildVRT(output_vrt, in_rasters, options=options)
        return output_vrt

    @staticmethod
    def raster_to_vrt(in_raster, out_vrt_dataset, no_data=-9999):
        """
        Converts a raster to a VRT (Virtual Raster) file and sets a NoData value if specified.
        :param in_raster: Input raster file path.
        :param out_vrt_dataset: Output VRT file path.
        :param no_data: NoData value to set in the VRT file.
        """
        # Open the input raster file
        src_ds = gdal.Open(in_raster, gdal.GA_Update)
        if src_ds is None:
            raise FileNotFoundError(f"Could not open input raster: {in_raster}")

        # Create a VRT dataset
        vrt_options = gdal.BuildVRTOptions(srcNodata=no_data)
        vrt_ds = gdal.BuildVRT(out_vrt_dataset, [src_ds], options=vrt_options)
        if vrt_ds is None:
            raise RuntimeError(f"Failed to create VRT file: {out_vrt_dataset}")

        vrt_ds.FlushCache()
        vrt_ds = None

    @staticmethod
    def update_nodata_vrt(vrt_file, no_data):
        """
        Updates the NoData value in an existing VRT file.
        Assumes the VRT file has only one raster band.
        :vrt_file: Path to the VRT file.
        :no_data: New NoData value to set.
        """
        vrt_ds = gdal.Open(vrt_file, gdal.GA_Update)
        if vrt_ds is None:
            raise FileNotFoundError(f"Could not open VRT file: {vrt_file}")
        if vrt_ds is None:
            raise FileNotFoundError(f"Could not open VRT file: {vrt_file}")

        band = vrt_ds.GetRasterBand(1)
        if band is None:
            raise RuntimeError("No raster band found in the VRT file.")
        # set Nodata value
        band.SetNoDataValue(no_data)
        vrt_ds = None
        vrt_ds.FlushCache()

    @staticmethod
    def setVRT_ColorInterp(raster_path, color_interp=gdal.GCI_GrayIndex):
        """
        Set color interpretation of a VRT file/raster
        :param raster_path: path to VRT/raster file
        :param color_interp: color interpretation
        """
        ds = gdal.Open(raster_path, gdal.GA_Update)
        band = ds.GetRasterBand(1)
        band.SetColorInterpretation(color_interp)
        ds.FlushCache()
        ds = None

    @staticmethod
    def getVRT_ColorInterp(vrt_path):
        """
        Get color interpretation of a VRT file
        :param vrt_path: path to VRT file
        """
        try:
            with gdal.Open(vrt_path) as ds:
                band = ds.GetRasterBand(1)
                color_interp = band.GetColorInterpretation()
                return color_interp
        except Exception as e:
            print(f"Failed to get color interpretation: {e}")
            return None

    @staticmethod
    def get_extent_raster(in_raster: str):
        """
        Extracts extent of raster file
        :param in_raster: input raster file
        :return extent as [xmin, xmax, ymin, ymax]
        """
        ds = gdal.Open(in_raster)
        if ds is None:
            raise FileNotFoundError(f"Could not open input raster: {in_raster}")
        geotransform = ds.GetGeoTransform()
        x_min = geotransform[0]
        x_max = geotransform[0] + ds.RasterXSize * geotransform[1]
        y_max = geotransform[3]
        y_min = geotransform[3] - ds.RasterYSize * abs(geotransform[5])
        extent = [x_min, x_max, y_min, y_max]
        return extent

    @staticmethod
    def get_extent_raster_mem(in_raster):
        """
        Extracts extent of raster file
        :param in_raster: input raster file
        :return extent as [xmin, xmax, ymin, ymax]
        """
        ds = in_raster
        if ds is None:
            raise FileNotFoundError(f"Could not open input raster: {in_raster}")
        geotransform = ds.GetGeoTransform()
        x_min = geotransform[0]
        x_max = geotransform[0] + ds.RasterXSize * geotransform[1]
        y_max = geotransform[3]
        y_min = geotransform[3] - ds.RasterYSize * abs(geotransform[5])
        extent = [x_min, x_max, y_min, y_max]
        return extent

    @staticmethod
    def get_raster_info(in_raster: str, info: list):
        """
        Extracts raster info
        :param in_raster: input raster file
        :param info: raster properties to extract, supported [extent, nodata, cell_size, projection_epsg, geotransform, data_type, size]
        :return dict with info
        """
        ds = gdal.Open(in_raster)
        if ds is None:
            raise FileNotFoundError(f"Could not open input raster: {in_raster}")
        export = {}
        if isinstance(info, str):
            info = [info]
        if "extent" in info:
            extent = Raster.get_extent_raster(in_raster)
            export["extent"] = extent
        if "nodata" in info:
            nodata = ds.GetRasterBand(1).GetNoDataValue()
            export["nodata"] = nodata
        if "cell_size" in info:
            geotransform = ds.GetGeoTransform()
            cell_size = abs(geotransform[1])
            export["cell_size"] = cell_size
        if "size" in info:
            size = [ds.RasterXSize, ds.RasterYSize]
            export["size"] = size
        if "projection_epsg" in info:
            srs = osr.SpatialReference()
            srs.ImportFromWkt(ds.GetProjection())
            export["projection_epsg"] = int(srs.GetAttrValue("AUTHORITY", 1))
        if "geotransform" in info:
            geotransform = ds.GetGeoTransform()
            export["geotransform"] = geotransform
        if "data_type" in info:
            data_type = ds.GetRasterBand(1).DataType
            export["data_type"] = data_type
        return export

    def create_empty_raster(
        output_path: str,
        extent: list,
        cell_size: float,
        epsg_code: int,
        nodata_value: float,
        gdal_dtype: int = gdal.GDT_Int32,
    ) -> gdal.Dataset:
        min_x, max_x, min_y, max_y = extent
        width = int((max_x - min_x) // cell_size + 1)
        height = int((max_y - min_y) // cell_size + 1)
        gt = (min_x, cell_size, 0, max_y, 0, -cell_size)
        raster_driver = gdal.GetDriverByName("GTiff")
        options = [
            "TFW=YES",
            "BIGTIFF=YES",
            "TILED=YES",
            "NUM_THREADS=ALL_CPUS",
            "COMPRESS=ZSTD",
        ]
        dst_ds = raster_driver.Create(
            output_path,
            width,
            height,
            1,
            gdal_dtype,
            options=options,
        )
        if dst_ds is None:
            raise RuntimeError(f"Failed to create GDAL dataset for {output_path}")
        dst_ds.SetGeoTransform(gt)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(epsg_code)
        crs_wkt = srs.ExportToWkt()
        dst_ds.SetProjection(crs_wkt)
        band_dst = dst_ds.GetRasterBand(1)
        band_dst.SetNoDataValue(float(nodata_value))
        return dst_ds

    @staticmethod
    def get_raster_info_mem(in_raster, info: list):
        """
        Extracts raster info
        :param in_raster: dataset of in memory raster
        :param info: raster properties to extract, supported [extent, nodata, cell_size, projection_epsg, geotransform, data_type, size]
        :return dict with info
        """
        ds = in_raster
        if ds is None:
            raise FileNotFoundError(f"Could not open input raster: {in_raster}")
        export = {}
        if isinstance(info, str):
            info = [info]
        if "extent" in info:
            extent = Raster.get_extent_raster_mem(in_raster)
            export["extent"] = extent
        if "nodata" in info:
            nodata = ds.GetRasterBand(1).GetNoDataValue()
            export["nodata"] = nodata
        if "cell_size" in info:
            geotransform = ds.GetGeoTransform()
            cell_size = abs(geotransform[1])
            export["cell_size"] = cell_size
        if "size" in info:
            size = [ds.RasterXSize, ds.RasterYSize]
            export["size"] = size
        if "projection_epsg" in info:
            srs = osr.SpatialReference()
            srs.ImportFromWkt(ds.GetProjection())
            export["projection_epsg"] = int(srs.GetAttrValue("AUTHORITY", 1))
        if "geotransform" in info:
            geotransform = ds.GetGeoTransform()
            export["geotransform"] = geotransform
        if "data_type" in info:
            data_type = ds.GetRasterBand(1).DataType
            export["data_type"] = data_type
        return export

    def reproject_raster(
        self,
        in_raster: str,
        out_raster: str,
        epsg_out: int,
        cell_size_out=0,
        data_type=gdal.GDT_Float32,
        interpolation="nearest",
        CPU_AVAILABLE=4,
        no_data=-9999,
        snap=False,
        extent=None,
        in_memory=False,
    ):
        """
        Reprojects raster to new projection and resolution
        :param in_raster: input raster file
        :param out_raster: output raster file
        :param epsg_out: EPSG code of output projection
        :param cell_size_out: output cell size
        :param data_type: output data type in GDAL format
        :param snap: snap extent to snap raster
        :param extent: extent of output raster
        """
        ds_in = gdal.Open(in_raster)
        if in_memory:
            out_raster = None  # in memory raster
        if ds_in is None:
            raise FileNotFoundError(f"Could not open input raster: {in_raster}")
        srs_in = osr.SpatialReference()
        srs_in.ImportFromWkt(ds_in.GetProjection())
        srs_out = osr.SpatialReference()
        srs_out.ImportFromEPSG(epsg_out)

        if cell_size_out == 0:
            cell_size_out = self.default_resolution
        if snap:
            if extent is None:
                extent = Raster.get_extent_raster(in_raster)
            extent_snapped = Vector.align_extent_to_snap(
                extent,
                snap_raster=self.snap_raster,
                cell_size=cell_size_out,
                optimize=True,
            )
            extent = extent_snapped
        bounds = [extent[0], extent[2], extent[1], extent[3]]
        if in_memory:
            format = "MEM"
            creation_options = []
        else:
            format = "GTiff"
            creation_options = [
                "COMPRESS=ZSTD",
                "PREDICTOR=2",
                "TILED=YES",
                "BIGTIFF=YES",
            ]
        options = gdal.WarpOptions(
            format=format,
            outputBounds=bounds,
            xRes=cell_size_out,
            yRes=cell_size_out,
            dstSRS=srs_out,
            srcSRS=srs_in,
            srcNodata=ds_in.GetRasterBand(1).GetNoDataValue(),
            dstNodata=no_data,
            outputType=data_type,
            resampleAlg=interpolation,
            warpOptions=[
                f"NUM_THREADS={CPU_AVAILABLE}",
            ],
            creationOptions=creation_options,
        )
        if in_memory:
            ds_out = gdal.Warp("", ds_in, options=options)
            return ds_out  # return in memory raster
        else:
            gdal.Warp(out_raster, ds_in, options=options)

        ds_in = None

    def rasterize_to_new_raster(
        self,
        out_raster,
        in_vector,
        value,
        cell_size,
        extent,
        format="GTiff",
        nodata=-9999,
        data_type=gdal.GDT_Float32,
        all_touched=False,
    ):
        vector = None
        try:
            vector = ogr.Open(in_vector)
        except Exception:
            drv = in_vector.GetDriver()
            drv_name = drv.ShortName
            if drv_name == "MEM" or drv_name == "Memory":
                vector = in_vector
        if vector is None:
            raise FileNotFoundError(
                f"Input vector file {in_vector} could not be opened."
            )

        # Create the output raster
        driver = gdal.GetDriverByName(format)
        if driver is None:
            raise ValueError(f"Driver {format} not found.")

        # Calculate raster dimensions
        x_min, x_max, y_min, y_max = extent
        cols = int((x_max - x_min) / cell_size)
        rows = int((y_max - y_min) / cell_size)
        if out_raster is None:
            out_raster_ds = driver.Create("", cols, rows, 1, data_type)
        else:
            out_raster_ds = driver.Create(out_raster, cols, rows, 1, data_type)
        if out_raster_ds is None:
            raise OSError(f"Could not create output raster file {out_raster}")

        # Set the spatial reference and geotransform
        out_raster_ds.SetGeoTransform((x_min, cell_size, 0, y_max, 0, -cell_size))
        out_raster_ds.SetProjection(vector.GetLayer(0).GetSpatialRef().ExportToWkt())

        # Set nodata value
        out_band = out_raster_ds.GetRasterBand(1)
        out_band.SetNoDataValue(nodata)

        options = []
        if all_touched:
            options.append("ALL_TOUCHED=TRUE")
        extent_str = f"{x_min},{y_min},{x_max},{y_max}"
        options.append(f"OUTPUT_BOUNDS={extent_str}")

        gdal.RasterizeLayer(
            out_raster_ds, [1], vector.GetLayer(), burn_values=[value], options=options
        )

        out_band.FlushCache()
        out_raster_ds.FlushCache()
        if format == "MEM":
            return out_raster_ds

    @staticmethod
    def get_nearest_cell_center(gt, point_x, point_y):
        """Return center of the cell containing the point
        :param gt: geotransform of the raster
        :param point_x: x coordinate of the point
        :param point_y: y coordinate of the point

        """
        x_origin, y_origin, cell_size_x, cell_size_y = gt[0], gt[3], gt[1], gt[5]
        snapped_x = (
            x_origin
            + cell_size_x * floor((point_x - x_origin) / cell_size_x)
            + cell_size_x / 2
        )
        snapped_y = (
            y_origin
            + cell_size_y * floor((point_y - y_origin) / cell_size_y)
            + cell_size_y / 2
        )
        return snapped_x, snapped_y

    @staticmethod
    def extract_window(
        band, px, py, window_size=3, safe_mode=False, mask_value=-999999
    ):
        """
        Extract a window around a pixel in a raster band
        :param band: raster band
        :param px: x coordinate of the center pixel
        :param py: y coordinate of the center pixel
        :param window_size: size of the window
        :param safe_mode: if True, return None if window is out of bounds
        :return: window as numpy array

        """
        # if window size is odd, return error
        if window_size % 2 == 0:
            raise ValueError("Window size must be an odd number.")
        half_window = window_size // 2
        start_px = px - half_window
        start_py = py - half_window
        try:
            window = band.ReadAsArray(start_px, start_py, window_size, window_size)
        except Exception:
            if safe_mode:
                return None
            else:
                window = np.full((window_size, window_size), mask_value)
        return window

    def snap_point_on_array_pos(x, y, window_size):
        center_x_index = int(window_size / 2)
        center_y_index = int(window_size / 2)
        dif_x = x - center_x_index
        dif_y = y - center_y_index
        return dif_x, dif_y

    # --- Window-based IO for large rasters (WKT-based, no EPSG assumption) ---

    @staticmethod
    def inspect(in_raster: str):
        """
        Storage properties of a raster that decide how fast windows can be read
        :param in_raster: raster file
        :return: dict with size, pixel size, data type, nodata, compression, block size,
            overviews, file size and raw (uncompressed) size
        """
        ds = gdal.Open(str(in_raster))
        band = ds.GetRasterBand(1)
        gt = ds.GetGeoTransform()
        block_x, block_y = band.GetBlockSize()
        structure = ds.GetMetadata("IMAGE_STRUCTURE") or {}
        bytes_per_pixel = gdal.GetDataTypeSize(band.DataType) // 8
        files = ds.GetFileList() or [str(in_raster)]
        return {
            "path": str(in_raster),
            "driver": ds.GetDriver().ShortName,
            "width": ds.RasterXSize,
            "height": ds.RasterYSize,
            "bands": ds.RasterCount,
            "pixel_x": gt[1],
            "pixel_y": gt[5],
            "data_type": gdal.GetDataTypeName(band.DataType),
            "bytes_per_pixel": bytes_per_pixel,
            "nodata": band.GetNoDataValue(),
            "scale": band.GetScale(),
            "offset": band.GetOffset(),
            "compression": structure.get("COMPRESSION", "NONE"),
            "predictor": structure.get("PREDICTOR", ""),
            "block_x": block_x,
            "block_y": block_y,
            "strip_layout": block_y < 64 or block_x == ds.RasterXSize,
            "overviews": band.GetOverviewCount(),
            "overview_widths": [
                band.GetOverview(i).XSize for i in range(band.GetOverviewCount())
            ],
            "file_size_gb": sum(os.path.getsize(f) for f in files if exists(f)) / 1e9,
            "raw_size_gb": ds.RasterXSize * ds.RasterYSize * bytes_per_pixel / 1e9,
            "crs_wkt": ds.GetProjection(),
        }

    @staticmethod
    def build_vrt(in_rasters, out_vrt):
        """
        VRT over rasters that must share one CRS; stops with a message otherwise
        :param in_rasters: raster path or list of paths
        :param out_vrt: output VRT path
        :return: path to the VRT
        """
        in_rasters = [str(p) for p in Config.as_list(in_rasters)]
        reference = None
        for path in in_rasters:
            srs = osr.SpatialReference(wkt=gdal.Open(path).GetProjection())
            if reference is None:
                reference = srs
            elif not srs.IsSame(reference):
                raise SystemExit(
                    f"DTM {path} has a different CRS than {in_rasters[0]}. "
                    "Reproject it to a common CRS first."
                )
        os.makedirs(os.path.dirname(str(out_vrt)) or ".", exist_ok=True)
        gdal.BuildVRT(
            str(out_vrt), in_rasters, options=gdal.BuildVRTOptions(resolution="highest")
        ).FlushCache()
        return str(out_vrt)

    @staticmethod
    def read_window(
        in_raster,
        bounds,
        resolution,
        resample="average",
        warp_memory_mb=512,
        overview="AUTO",
        multithread=True,
        dst_wkt=None,
        src_nodata=None,
    ):
        """
        Reads a map window resampled to a grid aligned with bounds (gdal.Warp to MEM).
        With overview="AUTO" GDAL reads from the overview closest to the target
        resolution; use "NONE" for resampling methods such as max that overviews
        built with AVERAGE cannot reproduce.
        :param in_raster: raster path or open dataset
        :param bounds: (xmin, ymin, xmax, ymax) in the output CRS (dst_wkt, else the
                       raster CRS), multiples of resolution
        :param resolution: output cell size
        :param resample: GDAL resampling algorithm
        :param warp_memory_mb: warp working buffer in MB
        :param overview: "AUTO", "NONE" or overview index
        :param multithread: warp with a worker thread (GDAL_NUM_THREADS)
        :param dst_wkt: target CRS as WKT when it differs from the source CRS
        :param src_nodata: source nodata value overriding the band's own
        :return: (float32 array with NaN as nodata, geotransform)
        """
        ds = (
            gdal.Open(str(in_raster))
            if isinstance(in_raster, (str, os.PathLike))
            else in_raster
        )
        band = ds.GetRasterBand(1)
        kwargs = {
            "format": "MEM",
            "outputBounds": list(bounds),
            "xRes": resolution,
            "yRes": resolution,
            "resampleAlg": resample,
            "outputType": gdal.GDT_Float32,
            "dstNodata": float("nan"),
            "multithread": multithread,
            "warpMemoryLimit": warp_memory_mb,
            "overviewLevel": overview,
        }
        if src_nodata is not None:
            kwargs["srcNodata"] = src_nodata
        elif band.GetNoDataValue() is not None:
            kwargs["srcNodata"] = band.GetNoDataValue()
        if dst_wkt:
            kwargs["dstSRS"] = dst_wkt
        out = gdal.Warp("", ds, options=gdal.WarpOptions(**kwargs))
        array = out.GetRasterBand(1).ReadAsArray()
        scale, offset = band.GetScale(), band.GetOffset()
        if scale not in (None, 1.0) or offset not in (None, 0.0):
            array = array * (scale or 1.0) + (offset or 0.0)
        return array.astype(np.float32, copy=False), out.GetGeoTransform()

    @staticmethod
    def read_native(in_raster, bounds):
        """
        Reads the native pixels covering bounds (no resampling)
        :param in_raster: raster path or open dataset
        :param bounds: (xmin, ymin, xmax, ymax) in the raster CRS
        :return: (float32 array with NaN as nodata, geotransform of the window)
        """
        ds = (
            gdal.Open(str(in_raster))
            if isinstance(in_raster, (str, os.PathLike))
            else in_raster
        )
        band = ds.GetRasterBand(1)
        gt = ds.GetGeoTransform()
        col0 = max(int(floor((bounds[0] - gt[0]) / gt[1] + 1e-9)), 0)
        col1 = min(int(ceil((bounds[2] - gt[0]) / gt[1] - 1e-9)), ds.RasterXSize)
        row0 = max(int(floor((bounds[3] - gt[3]) / gt[5] + 1e-9)), 0)
        row1 = min(int(ceil((bounds[1] - gt[3]) / gt[5] - 1e-9)), ds.RasterYSize)
        if col1 <= col0 or row1 <= row0:
            return np.empty((0, 0), np.float32), gt
        array = band.ReadAsArray(col0, row0, col1 - col0, row1 - row0).astype(
            np.float32
        )
        if band.GetNoDataValue() is not None:
            array[array == np.float32(band.GetNoDataValue())] = np.nan
        scale, offset = band.GetScale(), band.GetOffset()
        if scale not in (None, 1.0) or offset not in (None, 0.0):
            array = array * (scale or 1.0) + (offset or 0.0)
        window_gt = (gt[0] + col0 * gt[1], gt[1], 0.0, gt[3] + row0 * gt[5], 0.0, gt[5])
        return array, window_gt

    @staticmethod
    def array_to_mem(array, geotransform, wkt, nodata=None):
        """
        Numpy array to MEM dataset; float arrays get NaN as nodata unless set
        :param array: 2D array (bool is stored as Byte)
        :param geotransform: GDAL geotransform
        :param wkt: CRS as WKT
        :param nodata: nodata value
        :return: gdal.Dataset
        """
        if array.dtype == bool:
            array = array.astype(np.uint8)
        data_type = gdal_array.NumericTypeCodeToGDALTypeCode(array.dtype)
        ds = gdal.GetDriverByName("MEM").Create(
            "", array.shape[1], array.shape[0], 1, data_type
        )
        ds.SetGeoTransform(geotransform)
        ds.SetProjection(wkt)
        band = ds.GetRasterBand(1)
        if nodata is not None:
            band.SetNoDataValue(float(nodata))
        elif np.issubdtype(array.dtype, np.floating):
            band.SetNoDataValue(float("nan"))
        band.WriteArray(array)
        return ds

    @staticmethod
    def best_compression(preferred="ZSTD"):
        """
        :param preferred: preferred GTiff compression
        :return: preferred if the GDAL build supports it, otherwise DEFLATE
        """
        options = gdal.GetDriverByName("GTiff").GetMetadataItem(
            "DMD_CREATIONOPTIONLIST"
        )
        return preferred if preferred in options else "DEFLATE"

    @staticmethod
    def write_array(path, array, geotransform, wkt, nodata=None, compress="ZSTD"):
        """
        Writes an array as tiled, compressed GeoTIFF (predictor by data type)
        :param path: output GeoTIFF
        :param array: 2D array
        :param geotransform: GDAL geotransform
        :param wkt: CRS as WKT
        :param nodata: nodata value (NaN for float arrays if not set)
        :param compress: preferred compression
        """
        os.makedirs(os.path.dirname(str(path)) or ".", exist_ok=True)
        mem = Raster.array_to_mem(array, geotransform, wkt, nodata)
        is_float = np.issubdtype(array.dtype, np.floating)
        options = [
            "TILED=YES",
            "BLOCKXSIZE=256",
            "BLOCKYSIZE=256",
            f"COMPRESS={Raster.best_compression(compress)}",
            f"PREDICTOR={3 if is_float else 2}",
            "BIGTIFF=IF_SAFER",
            "NUM_THREADS=ALL_CPUS",
        ]
        gdal.Translate(str(path), mem, format="GTiff", creationOptions=options)

    @staticmethod
    def rasterize_geometries(
        geometries, values, shape, geotransform, wkt, all_touched=False, dtype=np.int32
    ):
        """
        Burns shapely geometries with per-geometry values onto a grid
        :param geometries: iterable of shapely geometries
        :param values: burn value per geometry
        :param shape: (rows, cols)
        :param geotransform: GDAL geotransform
        :param wkt: CRS as WKT
        :param all_touched: burn every cell touched by the geometry
        :param dtype: output numpy dtype
        :return: array, 0 where nothing was burnt
        """
        vector_ds = gdal.GetDriverByName("MEM").Create("", 0, 0, 0, gdal.GDT_Unknown)
        layer = vector_ds.CreateLayer(
            "features", osr.SpatialReference(wkt=wkt), ogr.wkbUnknown
        )
        is_float = np.issubdtype(np.dtype(dtype), np.floating)
        layer.CreateField(
            ogr.FieldDefn("value", ogr.OFTReal if is_float else ogr.OFTInteger64)
        )
        definition = layer.GetLayerDefn()
        for geometry, value in zip(geometries, values):
            if geometry is None or geometry.is_empty:
                continue
            feature = ogr.Feature(definition)
            feature.SetGeometry(ogr.CreateGeometryFromWkb(shapely.to_wkb(geometry)))
            feature.SetField("value", float(value) if is_float else int(value))
            layer.CreateFeature(feature)
        raster = Raster.array_to_mem(
            np.zeros(shape, dtype), geotransform, wkt, nodata=0
        )
        gdal.RasterizeLayer(
            raster,
            [1],
            layer,
            options=[
                "ATTRIBUTE=value",
                f"ALL_TOUCHED={'TRUE' if all_touched else 'FALSE'}",
            ],
        )
        return raster.GetRasterBand(1).ReadAsArray()

    @staticmethod
    def polygonize_labels(labels, geotransform, wkt, connectedness=4):
        """
        Polygonizes a label raster (0 = background), one dissolved polygon per label
        :param labels: 2D integer array
        :param geotransform: GDAL geotransform
        :param wkt: CRS as WKT
        :param connectedness: 4 or 8
        :return: gpd.GeoDataFrame with columns label, geometry
        """
        vector_ds = gdal.GetDriverByName("MEM").Create("", 0, 0, 0, gdal.GDT_Unknown)
        layer = vector_ds.CreateLayer(
            "polygons", osr.SpatialReference(wkt=wkt), ogr.wkbPolygon
        )
        layer.CreateField(ogr.FieldDefn("label", ogr.OFTInteger))
        if labels.any():
            label_ds = Raster.array_to_mem(
                labels.astype(np.int32), geotransform, wkt, nodata=0
            )
            mask_ds = Raster.array_to_mem(
                (labels > 0).astype(np.uint8), geotransform, wkt
            )
            options = ["8CONNECTED=8"] if connectedness == 8 else []
            gdal.Polygonize(
                label_ds.GetRasterBand(1), mask_ds.GetRasterBand(1), layer, 0, options
            )
        features = [
            (
                f.GetField("label"),
                shapely.from_wkb(bytes(f.GetGeometryRef().ExportToWkb())),
            )
            for f in layer
        ]
        gdf = gpd.GeoDataFrame(
            {"label": pd.Series([f[0] for f in features], dtype="int64")},
            geometry=gpd.GeoSeries([f[1] for f in features], crs=wkt),
        )
        if gdf.empty:
            return gdf
        return gdf.dissolve(by="label").reset_index()

    @staticmethod
    def fill_nodata_cells(array, to_fill, geotransform, wkt, max_search_cells):
        """
        Interpolates selected cells with gdal.FillNodata (inverse distance from the
        nearest valid cells); only the bounding box of to_fill is processed
        :param array: 2D float array, NaN = nodata
        :param to_fill: boolean mask of cells to interpolate
        :param geotransform: GDAL geotransform
        :param wkt: CRS as WKT
        :param max_search_cells: search distance in cells
        :return: (filled float32 array, number of cells that could not be filled)
        """
        result = array.astype(np.float32, copy=True)
        if not to_fill.any():
            return result, 0
        margin = int(max_search_cells) + 1
        rows = np.flatnonzero(to_fill.any(axis=1))
        cols = np.flatnonzero(to_fill.any(axis=0))
        r0, r1 = max(rows[0] - margin, 0), min(rows[-1] + 1 + margin, array.shape[0])
        c0, c1 = max(cols[0] - margin, 0), min(cols[-1] + 1 + margin, array.shape[1])
        sub = result[r0:r1, c0:c1].copy()
        sub_fill = to_fill[r0:r1, c0:c1]
        sub[sub_fill] = np.nan
        sub_gt = (
            geotransform[0] + c0 * geotransform[1],
            geotransform[1],
            0.0,
            geotransform[3] + r0 * geotransform[5],
            0.0,
            geotransform[5],
        )
        target = Raster.array_to_mem(sub, sub_gt, wkt)
        mask = Raster.array_to_mem((~np.isnan(sub)).astype(np.uint8), sub_gt, wkt)
        gdal.FillNodata(
            target.GetRasterBand(1),
            mask.GetRasterBand(1),
            maxSearchDist=int(max_search_cells),
            smoothingIterations=0,
        )
        filled = target.GetRasterBand(1).ReadAsArray()
        done = sub_fill & ~np.isnan(filled)
        window = result[r0:r1, c0:c1]
        window[done] = filled[done]
        return result, int((sub_fill & ~done).sum())

    @staticmethod
    def footprint(in_raster, resolution):
        """
        Valid-data footprint of a raster from a decimated read; GDAL serves the read
        from overviews when they exist, otherwise only the rows and columns needed
        are decoded. The mask is dilated by one coarse cell so that narrow strips of
        data lost by decimation are kept.
        :param in_raster: raster path
        :param resolution: footprint cell size in raster units
        :return: gpd.GeoDataFrame with one dissolved (multi)polygon
        """
        ds = gdal.Open(str(in_raster))
        band = ds.GetRasterBand(1)
        gt = ds.GetGeoTransform()
        factor = max(1, int(round(resolution / abs(gt[1]))))
        out_w = max(1, ceil(ds.RasterXSize / factor))
        out_h = max(1, ceil(ds.RasterYSize / factor))
        array = band.ReadAsArray(
            buf_xsize=out_w,
            buf_ysize=out_h,
            resample_alg=gdal.GRIORA_NearestNeighbour,
        )
        valid = np.ones(array.shape, bool)
        if band.GetNoDataValue() is not None:
            valid &= array != band.GetNoDataValue()
        if np.issubdtype(array.dtype, np.floating):
            valid &= ~np.isnan(array)
        valid = ndimage.binary_dilation(valid, structure=np.ones((3, 3), bool))
        coarse_gt = (
            gt[0],
            gt[1] * ds.RasterXSize / out_w,
            0.0,
            gt[3],
            0.0,
            gt[5] * ds.RasterYSize / out_h,
        )
        polygons = Raster.polygonize_labels(
            valid.astype(np.int32), coarse_gt, ds.GetProjection()
        )
        return gpd.GeoDataFrame(
            geometry=[shapely.union_all(polygons.geometry.values)],
            crs=ds.GetProjection(),
        )


class Vector:
    @staticmethod
    def load_vector(path, bbox=None, layer=None):
        if path.endswith((".parquet", ".pqt")):
            try:
                gdf = gpd.read_parquet(path)
                return gdf
            except Exception as e:
                print(f"Unable to load {path} to GeoDataFrame: {e}")
                return None
        else:
            try:
                if bbox:
                    gdf = gpd.read_file(path, engine="pyogrio", bbox=tuple(bbox))
                else:
                    gdf = gpd.read_file(path, engine="pyogrio")
                return gdf
            except Exception as e:
                print(f"Unable to load {path} to GeoDataFrame: {e}")
                return None

    @staticmethod
    def save_vector(gdf, path, gpd_driver="GPKG"):
        try:
            gdf.to_file(path, driver=gpd_driver)
        except Exception as e:
            print(f"Unable to save GeoDataFrame to {path}: {e}")

    @staticmethod
    def save_vector_to_parquet(gdf, path):
        try:
            gdf.to_parquet(path, engine="pyarrow")
        except Exception as e:
            print(f"Unable to save GeoDataFrame to Parquet file: {path}: {e}")

    @staticmethod
    def reproject_vector(gdf, epsg_out):
        if gdf.crs.to_epsg() == epsg_out:
            return gdf
        try:
            gdf_reprojected = gdf.to_crs(epsg=epsg_out)
            return gdf_reprojected
        except Exception as e:
            print(f"Unable to reproject GeoDataFrame: {e}")
            return None

    @staticmethod
    def buffer_vector(gdf: gpd.GeoDataFrame, buffer_size: float, dissolve=False):
        if not isinstance(gdf, gpd.GeoDataFrame):
            raise TypeError("Input must be a GeoDataFrame.")

        try:
            gdf["geometry"] = gdf["geometry"].buffer(buffer_size)

            if dissolve:
                gdf = gdf.dissolve()

            return gdf
        except Exception as e:
            raise RuntimeError(f"Unable to buffer GeoDataFrame: {e}")

    @staticmethod
    def select_by_attribute_categorical(
        gdf: gpd.GeoDataFrame, attribute: str, values: list, remove=False
    ):
        """
        Filter a GeoDataFrame based on categorical values from list of values.

        :param gdf: gpd.GeoDataFrame, the input GeoDataFrame
        :param attribute: str, the name of the column to filter on
        :param values: list, the values to select or remove
        :param remove: bool, if True, removes rows matching the values; otherwise, keep them - default
        :return: gpd.GeoDataFrame, the filtered GeoDataFrame
        """
        if isinstance(values, str):
            values = [values]
        if not isinstance(gdf, gpd.GeoDataFrame):
            raise TypeError("Input must be a GeoDataFrame.")
        if attribute not in gdf.columns:
            raise ValueError(f"Attribute '{attribute}' not found in the GeoDataFrame.")
        if remove:
            gdf = gdf[~gdf[attribute].isin(values)]
        else:
            gdf = gdf[gdf[attribute].isin(values)]
        return gdf

    @staticmethod
    def filter_by_extent(gdf: gpd.GeoDataFrame, extent: list):
        """
        :param gdf - Geodataframe to be filtered
        :param extent - list of coordinates [xmin, xmax, ymin, ymax]
        """
        # Define the extent as a box (bounding box)
        extent_box = box(extent[0], extent[2], extent[1], extent[3])
        gdf_filtered = gdf[gdf.geometry.intersects(extent_box)]
        return gdf_filtered

    @staticmethod
    def merge_geodataframes(gdf_list: list):
        """
        :param gdf_list
        """
        gdf_merged = pd.concat(gdf_list, ignore_index=True)
        return gdf_merged

    @staticmethod
    def list_to_txt(input_list, path, name):
        output_file = join(path, name)
        with open(output_file, "w") as output:
            output.write("\n".join(map(str, input_list)))

    @staticmethod
    def convert_extent_to_polygon(extent):
        """
        :param extent - list of coordinates [xmin, xmax, ymin, ymax]
        :return shapely.Polygon object
        """
        coordinates = (
            (extent[0], extent[2]),
            (extent[1], extent[2]),
            (extent[1], extent[3]),
            (extent[0], extent[3]),
            (extent[0], extent[2]),
        )
        return Polygon(coordinates)

    @staticmethod
    def convert_bounds_to_polygon(bounds):
        """
        :param extent - list of coordinates [xmin, ymin, xmax, ymax]
        :return shapely.Polygon object
        """
        coordinates = (
            (bounds[0], bounds[1]),
            (bounds[2], bounds[1]),
            (bounds[2], bounds[3]),
            (bounds[0], bounds[3]),
            (bounds[0], bounds[1]),
        )
        return Polygon(coordinates)

    @staticmethod
    def convert_extent_to_ogr_geometry(extent):
        """
        :param extent - list of coordinates [xmin, xmax, ymin, ymax]
        :return ogr.Geometry object
        """
        # Create ring
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint_2D(extent[0], extent[2])
        ring.AddPoint_2D(extent[1], extent[2])
        ring.AddPoint_2D(extent[1], extent[3])
        ring.AddPoint_2D(extent[0], extent[3])
        ring.AddPoint_2D(extent[0], extent[2])

        # Create polygon
        polygon = ogr.Geometry(ogr.wkbPolygon)
        polygon.AddGeometry(ring)
        return polygon

    @staticmethod
    def ext1_to_ext2(source_ext, source_epsg, out_epsg):
        polygon = Vector.convert_extent_to_polygon(source_ext)
        gdf = gpd.GeoDataFrame(
            index=[0], crs=f"epsg:{source_epsg}", geometry=[polygon]
        ).to_crs(epsg=out_epsg)
        bounds = gdf.total_bounds
        return [bounds[0], bounds[2], bounds[1], bounds[3]]

    def bounds1_to_bounds2(bounds, source_epsg, out_epsg):
        polygon = Vector.convert_bounds_to_polygon(bounds)
        gdf = gpd.GeoDataFrame(
            index=[0], crs=f"epsg:{source_epsg}", geometry=[polygon]
        ).to_crs(epsg=out_epsg)
        bounds = gdf.total_bounds
        return bounds

    def progress_callback(complete, message, data):
        percent = int(complete * 100)  # round to integer percent
        data.update(percent)  # update the progressbar
        return 1

    def select_files_in_extent(list_files, extent):
        """
        select files where extents intersects
        :return list of file intersecting given extent
        """
        extent_polygon = Vector.convert_extent_to_polygon(extent)
        filtered_list = []
        for file in list_files:
            compare_polygon = Vector.convert_extent_to_polygon(
                Vector.get_extent_vector(file)
            )
            if shapely.intersects(extent_polygon, compare_polygon):
                filtered_list.append(file)
        return filtered_list

    def get_epsg_from_vector(file):
        """
        extracts epsg of polygon file
        :param file vector dataset or geopandas.GeodataFrame
        :return extent as epsg code as int
        """
        if isinstance(file, gpd.GeoDataFrame):
            return int(file.crs.to_epsg())

        elif os.path.splitext(file)[1] in [".shp", ".gpkg", ".geojson"]:
            data = ogr.Open(file)
            layer = data.GetLayer()
            spatial_ref = layer.GetSpatialRef()
            if spatial_ref.IsGeographic() == 1:  # this is a geographic srs
                proj_type = "GEOGCS"
                return int(spatial_ref.GetAuthorityCode(proj_type))
            else:  # this is a projected srs
                proj_type = "PROJCS"
                return int(spatial_ref.GetAuthorityCode(proj_type))
        else:
            print(f"File: {file} is not in supported format.")

    def clip_polygon(data, clipping_area):
        """
        :param data vector dataset
        :type data gpd.Geodataframe or string
        :param clipping_area vector dataset
        :type clipping_area gpd.Geodataframe or string
        if geometries are not in the same projection, CRS from param clipping area would be applied
        """
        if isinstance(data, gpd.GeoDataFrame):
            pass
        else:
            data = gpd.read_file(data)
        if isinstance(clipping_area, gpd.GeoDataFrame):
            pass
        else:
            clipping_area = gpd.read_file(clipping_area)

        data_epsg = data.crs.to_epsg()
        clipping_area_epsg = clipping_area.crs.to_epsg()
        if data_epsg != clipping_area_epsg:
            data = data.to_crs(clipping_area_epsg)
        clipped_data = gpd.clip(data, clipping_area)
        return clipped_data

    def fix_geometry(data, compare=False):
        """
        :param data vector dataset
        :type data gpd.Geodataframe or string
        :param compare if True, function will return number of rows where geometry was fixed
        :type compare bool
        :return fixed geometries
        Returns fixed geometries for geodataframe or vector geodataset.
        """
        if isinstance(data, gpd.GeoDataFrame):
            pass
        elif os.path.splitext(data)[1] in [".shp", ".gpkg", ".geojson"]:
            data = gpd.read_file(data)
        else:
            raise ValueError(f"File: {data} is not in supported format.")
        if compare:
            data["copy_geometry"] = data["geometry"]
        g = data["geometry"]
        g_v = g.make_valid()
        data["geometry"] = g_v
        if compare:
            data["changed"] = np.where(data["copy_geometry"] == data["geometry"], 0, 1)
            geom_change = data["changed"].sum()
            if geom_change > 0:
                print(f"Geometries were fixed in {geom_change} rows")
            data.drop(columns=["copy_geometry", "changed"], inplace=True)
        return data

    def align_extent_to_snap(extent, snap_raster, cell_size=0, optimize=False):
        """
        :param extent - list of coordinates [left_x, right_x, bottom_y, upper_y] (mandatory)
        :param snap_raster - reference grid for snapping (mandatory) or its extent
        :param cell_size - raster resolution, mandatory if optimize is True (mandatory if optimize is True)
        :param optimize - if true, bounding box will be resized in accordance with raster original extent (optional)
        :return extent snapped to snap raster
        """
        ds_snap = gdal.Open(snap_raster)
        geotransform_snap = list(ds_snap.GetGeoTransform())
        if (
            geotransform_snap[1] > 0
            and (extent[0] - geotransform_snap[0]) > 0
            or geotransform_snap[1] < 0
            and ((extent[0] - geotransform_snap[0]) < 0)
        ):
            meth_x = floor
        else:
            meth_x = ceil

        if (
            geotransform_snap[5] > 0
            and (extent[3] - geotransform_snap[3]) > 0
            or geotransform_snap[5] < 0
            and ((extent[3] - geotransform_snap[3]) < 0)
        ):
            meth_y = floor
        else:
            meth_y = ceil

        x_left = (
            geotransform_snap[0]
            + meth_x((extent[0] - geotransform_snap[0]) / geotransform_snap[1])
            * geotransform_snap[1]
        )
        y_upper = (
            geotransform_snap[3]
            + meth_y((extent[3] - geotransform_snap[3]) / geotransform_snap[5])
            * geotransform_snap[5]
        )
        x_right = (
            x_left
            + ceil((extent[1] - x_left) / geotransform_snap[1]) * geotransform_snap[1]
        )
        y_lower = (
            y_upper
            + ceil((extent[2] - y_upper) / geotransform_snap[5]) * geotransform_snap[5]
        )
        if optimize:
            x_left = floor((extent[0] - x_left) / cell_size) * cell_size + x_left
            x_right = ceil((extent[1] - x_right) / cell_size) * cell_size + x_right
            y_lower = floor((extent[2] - y_lower) / cell_size) * cell_size + y_lower
            y_upper = ceil((extent[3] - y_upper) / cell_size) * cell_size + y_upper

        extent_snapped = [x_left, x_right, y_lower, y_upper]
        return extent_snapped

    @staticmethod
    def get_extent_layer(layer):
        extent = layer.GetExtent()
        return [extent[0], extent[1], extent[2], extent[3]]

    @staticmethod
    def get_extent_vector(file):
        """
        Extracts extent of polygon file
        :param file vector dataset or geopandas.GeoDataFrame
        :return extent as [xmin, xmax, ymin, ymax]
        """
        if isinstance(file, gpd.GeoDataFrame):
            bounds = file.total_bounds
            return [bounds[0], bounds[2], bounds[1], bounds[3]]
        elif os.path.splitext(file)[1] in [".shp", ".gpkg", ".geojson"]:
            data = ogr.Open(file)
            layer = data.GetLayer()
            extent = Vector.get_extent_layer(layer)
            return extent
        else:
            raise ValueError(f"File: {file} is not in supported format.")

    def get_convex_hull(file, out_epsg, out_type="shapely_polygon"):
        """
        based on: https://pcjericks.github.io/py-gdalogr-cookbook/vector_layers.html#save-the-convex-hull-of-all-geometry-from-an-input-layer-to-an-output-layer
        extracts extent of convex hull in specified projection
        :param file vector dataset or gpd.GeodataFrame
        :type file str or object
        :param out_epsg epsg code
        :type out_epsg int
        :param out_type defines output's data_type (shapely.Polygon or ogr.Geometry)
        :return convex hull as Polygon
        """
        if isinstance(file, (gpd.GeoDataFrame, gpd.GeoSeries)):
            convex_hull = Polygon(file.convex_hull.get_coordinates())

        elif os.path.splitext(file)[1] in [".shp", ".gpkg", ".geojson"]:
            data = ogr.Open(file)
            layer = data.GetLayer()
            geomcol = ogr.Geometry(ogr.wkbGeometryCollection)
            for feature in layer:
                geomcol.AddGeometry(feature.GetGeometryRef())

            # Calculate convex hull
            convexhull = geomcol.ConvexHull()
            sourceprj = layer.GetSpatialRef()
            targetprj = osr.SpatialReference()
            targetprj.ImportFromEPSG(out_epsg)
            transform = osr.CoordinateTransformation(sourceprj, targetprj)
            convexhull.Transform(transform)
            if out_type == "shapely_polygon":
                convex_hull = loads(convexhull.ExportToWkt())
        else:
            print(f"Warning: File {file} is not in supported format.")
            return None
        return convex_hull

    @staticmethod
    def extract_intersection(gdf_lines1, gdf_lines2):
        """
        Extracts intersection of two line layers
        :param gdf_lines1: first line layer
        :param gdf_lines2: second line layer
        :return: gpd.GeoDataFrame with intersection
        """
        if gdf_lines1.crs != gdf_lines2.crs:
            Exception("CRS of both layers must be the same")
        intersections = gdf_lines1.unary_union.intersection(gdf_lines2.unary_union)
        # Filter only Point geometries
        try:
            points = [geom for geom in intersections.geoms if geom.geom_type == "Point"]
        except Exception:
            points = []
        return points

    @staticmethod
    def load_ogr_ds_to_geodataframe(ogr_ds):
        gdf = (pyogrio.read_dataframe(ogr_ds),)
        return gdf

    def sjoin_nearest(gdf1, gdf2, distance=0.1):
        """
        Spatial join of two geodataframe on nearest predicate
        :param gdf1: gpd.GeoDataFrame
        :param gdf2: gpd.GeoDataFrame with points
        :param distance: float, maximum distance to join
        :return: gpd.GeoDataFrame with points and lines
        """
        gdf_joined = gpd.sjoin_nearest(
            gdf1, gdf2, how="inner", max_distance=distance, distance_col="temp_dist"
        )
        return gdf_joined

    def merge_lines_on_pseudonodes(gdf_lines):
        """
        Merges lines on pseudonodes
        BEWARE: this function return geodataframe with geometry only
        :param gdf_lines: gpd.GeoDataFrame with Linestring as geometry
        :return: gpd.GeoDataFrame with merged lines
        """
        crs_orig = gdf_lines.crs
        if len(gdf_lines.index) > 1:
            geom_merged = linemerge(unary_union(gdf_lines.geometry))
            gdf_lines_merged = (
                gpd.GeoDataFrame(geometry=[geom_merged])
                .explode()
                .reset_index()
                .set_crs(crs_orig)
            )
        else:
            gdf_lines_merged = gdf_lines.copy()
        return gdf_lines_merged

    @staticmethod
    def check_single_layer(path):
        """
        Stops with a clear message unless the vector file contains exactly one layer
        :param path: vector file
        """
        ds = ogr.Open(str(path))
        if ds is None:
            raise SystemExit(f"Cannot open vector file: {path}")
        if ds.GetLayerCount() != 1:
            names = [ds.GetLayer(i).GetName() for i in range(ds.GetLayerCount())]
            raise SystemExit(
                f"{path} contains {ds.GetLayerCount()} layers {names}; every vector input "
                "must contain exactly one layer (convert it with gdal.VectorTranslate)."
            )

    @staticmethod
    def read_merged(paths, wkt, bbox=None, where=None, keep_columns=None):
        """
        Reads one or more single-layer vector files, reprojects them to wkt and merges
        them; the column source holds the file stem
        :param paths: path or list of paths
        :param wkt: target CRS as WKT
        :param bbox: optional (xmin, ymin, xmax, ymax) in the target CRS
        :param where: optional SQL WHERE clause
        :param keep_columns: attribute columns to keep when present; True keeps all
        :return: gpd.GeoDataFrame with keep_columns, source and geometry
        """
        frames = []
        for path in Config.as_list(paths):
            Vector.check_single_layer(path)
            kwargs = {}
            if bbox is not None:
                kwargs["bbox"] = gpd.GeoSeries([box(*bbox)], crs=wkt)
            if where:
                kwargs["where"] = where
            gdf = gpd.read_file(path, **kwargs)
            if gdf.crs is None:
                log.warning("%s has no CRS; assuming the DTM CRS.", path)
                gdf = gdf.set_crs(wkt)
            gdf = gdf.to_crs(wkt)
            gdf["source"] = os.path.splitext(os.path.basename(str(path)))[0]
            if keep_columns is True:
                columns = [c for c in gdf.columns if c not in ("source", "geometry")]
            else:
                columns = [c for c in (keep_columns or []) if c in gdf.columns]
            frames.append(gdf[columns + ["source", "geometry"]])
        gdf = gpd.GeoDataFrame(
            pd.concat(frames, ignore_index=True), geometry="geometry", crs=wkt
        )
        return gdf[gdf.geometry.notna() & ~gdf.geometry.is_empty].reset_index(drop=True)

    @staticmethod
    def write_gpkg(gdf, path):
        """
        Writes a GeoDataFrame as its own single-layer GeoPackage (replaces the file)
        :param gdf: gpd.GeoDataFrame
        :param path: output .gpkg
        """
        os.makedirs(os.path.dirname(str(path)) or ".", exist_ok=True)
        if exists(str(path)):
            os.remove(str(path))
        gdf.to_file(str(path), driver="GPKG")

    @staticmethod
    def length_near_boundary(polygons, lines, distance):
        """
        Length of lines within distance of each polygon boundary
        :param polygons: gpd.GeoDataFrame of polygons
        :param lines: gpd.GeoDataFrame of lines in the same CRS
        :param distance: search distance in CRS units
        :return: numpy array with one length per polygon
        """
        if polygons.empty or lines.empty:
            return np.zeros(len(polygons))
        zones = shapely.buffer(shapely.boundary(polygons.geometry.values), distance)
        idx_zone, idx_line = lines.sindex.query(zones, predicate="intersects")
        lengths = shapely.length(
            shapely.intersection(zones[idx_zone], lines.geometry.values[idx_line])
        )
        return np.bincount(idx_zone, weights=lengths, minlength=len(polygons))

    @staticmethod
    def layer_names(path):
        """
        :param path: vector file
        :return: list of layer names (stops with a message if the file cannot be opened)
        """
        ds = ogr.Open(str(path))
        if ds is None:
            raise SystemExit(f"Cannot open vector file: {path}")
        return [ds.GetLayer(i).GetName() for i in range(ds.GetLayerCount())]

    @staticmethod
    def single_layer_copy(path, layer, out_path):
        """
        Copies one layer of a multi-layer vector source into its own single-layer
        GeoPackage with gdal.VectorTranslate, so that the copy can be read without
        selecting a layer
        :param path: multi-layer vector file
        :param layer: name of the layer to copy
        :param out_path: output .gpkg (replaced)
        :return: out_path as string
        """
        if layer not in Vector.layer_names(path):
            raise SystemExit(
                f"{path} has no layer {layer!r}; layers: {Vector.layer_names(path)}"
            )
        os.makedirs(os.path.dirname(str(out_path)) or ".", exist_ok=True)
        if exists(str(out_path)):
            os.remove(str(out_path))
        copied = gdal.VectorTranslate(
            str(out_path),
            str(path),
            options=gdal.VectorTranslateOptions(format="GPKG", layers=[layer]),
        )
        if copied is None:
            raise SystemExit(f"Layer {layer!r} of {path} could not be copied.")
        copied = None
        return str(out_path)


class Geometry:
    # from StackOverflow to be tested
    def find_nearest_line(point, lines):
        return min(lines, key=lambda line: line.distance(point))

    def extract_points_on_line(line, point, dist: list):
        proj_point = nearest_points(line, point)[1]
        line_length = line.length
        distances = [
            line.project(proj_point) - dist[0],
            line.project(proj_point) + dist[1],
        ]

        start_dist = max(0, distances[0])
        end_dist = min(line_length, distances[1])
        start_point = line.interpolate(start_dist)
        end_point = line.interpolate(end_dist)
        return start_point, end_point

    def create_points_on_line(line, distance):
        points = []
        for i in np.arange(0, line.length, distance):
            point = line.interpolate(i)
            points.append(point)
        return points

    def construct_point(x, y):
        return Point(x, y)

    def construct_line(points):
        return LineString(points)

    def calculate_angle(p1, p2, p3):
        a = sqrt((p2.x - p1.x) ** 2 + (p2.y - p1.y) ** 2)
        b = sqrt((p3.x - p2.x) ** 2 + (p3.y - p2.y) ** 2)
        c = sqrt((p3.x - p1.x) ** 2 + (p3.y - p1.y) ** 2)
        if a == 0 or b == 0:
            return 180
        angle = acos((a**2 + b**2 - c**2) / (2 * a * b))
        return degrees(angle)

    def close_holes(poly: Polygon) -> Polygon:
        """
        Close polygon holes by limitation to the exterior ring.

        :param poly: Input shapely Polygon
        :return: Polygon geometry

            :Example:

            df.geometry.apply(lambda p: close_holes(p))

        """
        if poly.interiors:
            return Polygon(list(poly.exterior.coords))
        else:
            return poly

    def create_fishnet(extent, cell_size, block_size, overlap):
        """
        Creates regular square grid
        :param extent: list of coordinates [xmin, xmax, ymin, ymax]
        :param cell_size: size of the cell in meters
        :param block_size: size of the edge in number of cells
        :param overlap: overlap between cells in meters
        :return: list of shapely.Polygon objects representing the grid cells
        """
        xmin, xmax, ymin, ymax = extent
        block_size = block_size * cell_size
        x, y = (xmin, ymin)
        geom_array = []
        while y <= ymax:
            while x <= xmax:
                geom = Polygon(
                    [
                        (x, y),
                        (x, y + block_size),
                        (x + block_size, y + block_size),
                        (x + block_size, y),
                        (x, y),
                    ]
                )
                geom_array.append(geom)
                x += block_size - overlap
            x = xmin
            y += block_size - overlap
        return geom_array


class CulvertsSpecific:
    @staticmethod
    def create_fishnet(extent, cell_size, block_size, overlap):
        """
        Creates regular square grid
        :param extent: list of coordinates [xmin, xmax, ymin, ymax]
        :param cell_size: size of the cell in meters
        :param block_size: size of the edge in number of cells
        :param overlap: overlap between cells in meters
        :return: list of shapely.Polygon objects representing the grid cells
        """
        xmin, xmax, ymin, ymax = extent
        block_size = block_size * cell_size
        x, y = (xmin, ymin)
        geom_array = []
        while y <= ymax:
            while x <= xmax:
                geom = Polygon(
                    [
                        (x, y),
                        (x, y + block_size),
                        (x + block_size, y + block_size),
                        (x + block_size, y),
                        (x, y),
                    ]
                )
                geom_array.append(geom)
                x += block_size - overlap
            x = xmin
            y += block_size - overlap
        return geom_array


class Hydro:
    @staticmethod
    def fill_depression_wang(dtm_array):
        filled_dtm, d8 = pyflwdir.dem.fill_depressions(elevtn=dtm_array)
        return filled_dtm

    def flow_accumulation(d8):
        flow_acc = pyflwdir.flow_accumulation(d8)
        return flow_acc

    @staticmethod
    def breach_to_streams(dem, stream, burn=0.0, burn_buffer_cells=0.0, eps=1e-5):
        """
        Stream cell reached downstream by every cell. Routing runs on a conditioned
        copy of dem: tapered stream burning (cells lowered by burn at the stream,
        decreasing linearly to zero at burn_buffer_cells, chamfer distance), complete
        breaching (priority flood from the stream cells, depressions carved along the
        flood path, flats given an eps gradient), then D8 steepest descent. Without the
        burning, paths can run along channel cells parallel to the rasterized stream and
        reach it only far downstream. The input dem is not modified.
        :param dem: 2D float array, NaN = nodata
        :param stream: 2D bool array of stream cells (the only outlets)
        :param burn: lowering of the stream cells in the conditioned copy
        :param burn_buffer_cells: distance in cells over which the lowering tapers off
        :param eps: elevation step used for carving and flats
        :return: 2D int32 array of flat indices of the stream cell reached, -1 if none
        """
        if dem.size >= 2**31:
            raise ValueError(
                "Window too large for 32-bit cell indices; use smaller units."
            )
        return _breach_to_streams(
            np.ascontiguousarray(dem, dtype=np.float32),
            np.ascontiguousarray(stream, dtype=np.bool_),
            float(burn),
            float(burn_buffer_cells),
            eps,
        )

    @staticmethod
    def priority_flood(surface, seeds):
        """
        Lowest level at which water starting at the seeds reaches each cell
        (4-connected, spill = max(spill of the previous cell, surface)). Cells reached at
        the current minimal level go through a plain queue (improved priority flood).
        :param surface: 2D float array, NaN = impassable
        :param seeds: 2D bool array of seed cells (level 0)
        :return: (spill float32 with NaN where unreached,
                  seed int32 flat index of the seed each minimal path starts from,
                  order int32 in which cells were reached, -1 where unreached)
        """
        if surface.size >= 2**31:
            raise ValueError(
                "Window too large for 32-bit cell indices; use smaller units."
            )
        spill, seed, order = _priority_flood(
            np.ascontiguousarray(surface, dtype=np.float32),
            np.ascontiguousarray(seeds, dtype=np.bool_),
        )
        spill[~np.isfinite(spill)] = np.nan
        return spill, seed, order

    @staticmethod
    def distance_cells(mask, limit):
        """
        Chamfer distance (weights 1 and sqrt 2) in cells to the nearest True cell
        :param mask: 2D bool array
        :param limit: distances are capped at this value
        :return: 2D float32 array
        """
        return _chamfer_distance(
            np.ascontiguousarray(mask, dtype=np.bool_), np.float32(limit)
        )


class MosaicRasters:
    """
    Class for creating mosaic for rasters, min added in 1.0.3

    :param rasters: List of input rasters
    :param output: Output file as tif or any other
    :param cell_size: Output cell size
    :param epsg: EPSG code is integer
    :param data_type: data type as gdal.GDT_Int16 or other, gdal.GDT_Float32
    :param no_data: Output no data value
    :param progress: Print progress sof merging


    :return:   String with saved tif

    """

    def __init__(
        self,
        rasters: list,
        output: str,
        cell_size: float | None = None,
        epsg: int | None = None,
        data_type=gdal.GDT_Int16,
        no_data: float = -32768,
        compress: bool = True,
        progress: bool = True,
        extent: list | None = None,
        compress_type=None,
        round_extent: bool = False,
    ):
        if compress_type is None:
            self.compress_type = ["COMPRESS=LZW"]
        else:
            self.compress_type = compress_type
        self.rasters = rasters
        self.output = output
        self.cell_size = cell_size
        if not self.cell_size:
            self.cell_size = self.get_resolution(rasters[0])[0]
        self.data_type = data_type
        self.no_data = no_data
        self.epsg = epsg
        self.progress = progress
        self.round_extent = round_extent

        if extent is not None:
            self.extent = extent
            self.bands = 1
        else:
            self.extent, self.bands = self._get_max_extent()

        if self.epsg is None:
            self.epsg = self.get_epsg(rasters[0])
        self.compress = compress
        self.dst = self._initialize_raster()
        self.dst = None
        self.dst = gdal.Open(self.output, gdal.GA_Update)

    @staticmethod
    def get_epsg(raster: str) -> int:
        """
        Extracting epsg if not specified on init
        """
        d = gdal.Open(raster)
        proj = osr.SpatialReference(wkt=d.GetProjection())
        return proj.GetAttrValue("AUTHORITY", 1)

    @staticmethod
    def get_resolution(raster):
        """
        Method which will return restolution in x and y direction

        :param raster: Raster file as string
        :return: X, Y cell size
        """
        src = gdal.Open(raster)
        gt = src.GetGeoTransform()
        cell_size_x = gt[1]
        cell_size_y = gt[5]
        src = None
        return cell_size_x, cell_size_y

    def _get_max_extent(self):
        """
        Private method which will get extent from all rasters

        :return: [total_xmin, total_xmax, total_ymin, total_ymax]
        """
        total_xmax, total_xmin, total_ymax, total_ymin = [], [], [], []
        bands_total = []
        for file in self.rasters:
            src = gdal.Open(file, gdal.GA_ReadOnly)

            bands_total.append(src.RasterCount)
            upx, xres, xskew, upy, yskew, yres = src.GetGeoTransform()

            ## no of cols/rows
            cols = src.RasterXSize
            rows = src.RasterYSize
            ## calculation of raster extent
            xmin = upx + 0 * xres + rows * xskew
            ymin = upy + 0 * yskew + rows * yres

            xmax = upx + cols * xres + 0 * xskew
            ymax = upy + cols * yskew + 0 * yres

            total_xmax.append(xmax)
            total_xmin.append(xmin)
            total_ymax.append(ymax)
            total_ymin.append(ymin)
            src = None
        if self.round_extent:
            return [
                round(min(total_xmin), 2),
                round(max(total_xmax), 2),
                round(min(total_ymin), 2),
                round(max(total_ymax), 2),
            ], max(bands_total)
        else:
            return [
                min(total_xmin),
                max(total_xmax),
                min(total_ymin),
                max(total_ymax),
            ], max(bands_total)

    def _initialize_raster(self):
        """
        Private method which will create output file
        """
        extent = self.extent
        if self.progress:
            print(extent)
        tiff_width = int((extent[1] - extent[0]) // self.cell_size)
        tiff_height = int((extent[3] - extent[2]) // self.cell_size)
        if self.progress:
            print(f"Size of new raster is {tiff_width} x {tiff_height}")

        gt = extent[0], self.cell_size, 0, extent[3], 0, -self.cell_size

        raster = gdal.GetDriverByName("GTiff")
        if self.compress:
            options = [
                "TFW=YES",
                "BIGTIFF=YES",
                "TILED=YES",
                "NUM_THREADS=ALL_CPUS",
            ]
            options.extend(self.compress_type)
            dst_ds = raster.Create(
                self.output,
                int(tiff_width),
                int(tiff_height),
                self.bands,
                self.data_type,
                options=options,
            )
        else:
            dst_ds = raster.Create(
                self.output,
                int(tiff_width),
                int(tiff_height),
                self.bands,
                self.data_type,
                options=["TFW=YES", "BIGTIFF=YES", "TILED=YES"],
            )
        ## setting upper left corner and resolution, y has to be minus
        dst_ds.SetGeoTransform(gt)

        srs = osr.SpatialReference()
        if isinstance(self.epsg, int) or "EPSG" in self.epsg or "epsg" in self.epsg:
            srs.ImportFromEPSG(self.epsg)
        elif "ESRI" in self.epsg:
            srs.ImportFromESRI(self.epsg.split(":")[1])
        elif "proj" in self.epsg:
            srs.ImportFromProj4(self.epsg)
        else:
            srs.ImportFromWkt([self.epsg])

        srs = srs.ExportToWkt()
        dst_ds.SetProjection(srs)
        band = dst_ds.GetRasterBand(1)

        ## we can change no data value if needed, this one is UnsignedInteger16Bit
        band.SetNoDataValue(self.no_data)

        # array = np.zeros((int(tiff_height), int(tiff_width)))
        # band.WriteArray(array)

        return dst_ds

    def mosaic_maximum(self):
        """
        Calling this method on the class will start the mosaicing of maximum per all rasters.

        :return: None

        :Example:

        ::

         MosaicRasters(*arg).mosaic_maximum()
        """
        for ix, file in enumerate(self.rasters):
            try:
                if self.progress:
                    print(f"{ix + 1} out of {len(self.rasters)}")

                ## write raster into array
                src = gdal.Open(file)

                for band_number in range(src.RasterCount):
                    band_dst = self.dst.GetRasterBand(band_number + 1)
                    band = src.GetRasterBand(band_number + 1)
                    desc = band.GetDescription()
                    band_dst.SetDescription(desc)

                    no_data_src = band.GetNoDataValue()

                    cols = src.RasterXSize
                    rows = src.RasterYSize

                    ## upper left corner
                    gt = src.GetGeoTransform()
                    # print(gt)
                    # print(dst.GetGeoTransform())
                    minx = gt[0]
                    maxy = gt[3]

                    y, x = maxy - self.extent[3], self.extent[0] - minx

                    locy, locx = (
                        int(abs(y / self.cell_size)),
                        int(abs(x / self.cell_size)),
                    )

                    # print(f'\tRequesting array of size {rows} x {cols} at location {locx} {locy}')

                    ## comparing if array already written is smaller or not
                    array_src = band.ReadAsArray(0, 0, cols, rows)
                    array_src = np.where(
                        (array_src == no_data_src), self.no_data, array_src
                    )

                    array_dst = band_dst.ReadAsArray(locx, locy, cols, rows)

                    array_to_write = np.where(
                        array_dst != self.no_data,
                        np.where(
                            array_src != self.no_data,
                            np.where(array_dst < array_src, array_src, array_dst),
                            array_dst,
                        ),
                        np.where(array_src != self.no_data, array_src, self.no_data),
                    )

                    self.dst.GetRasterBand(band_number + 1).WriteArray(
                        array_to_write, locx, locy
                    )
                src = None
            except Exception as e:
                print(e)

        self.dst = None
        return self.output

    def mosaic_minimum(self):
        """
        :Example:

        ::

         MosaicRasters(*arg).mosaic_minimum()

        """

        for ix, file in enumerate(self.rasters):
            try:
                if self.progress:
                    print(f"{ix + 1} out of {len(self.rasters)}")
                ## write raster into array
                src = gdal.Open(file)

                for band_number in range(src.RasterCount):
                    band_dst = self.dst.GetRasterBand(band_number + 1)
                    band = src.GetRasterBand(band_number + 1)
                    no_data_src = band.GetNoDataValue()

                    cols = src.RasterXSize
                    rows = src.RasterYSize

                    ## upper left corner
                    gt = src.GetGeoTransform()
                    # print(gt)
                    # print(dst.GetGeoTransform())
                    minx = gt[0]
                    maxy = gt[3]

                    y, x = maxy - self.extent[3], self.extent[0] - minx

                    locy, locx = (
                        int(abs(y / self.cell_size)),
                        int(abs(x / self.cell_size)),
                    )

                    # print(f'\tRequesting array of size {rows} x {cols} at location {locx} {locy}')

                    ## comparing if array already written is smaller or not
                    array_src = band.ReadAsArray(0, 0, cols, rows)
                    array_src = np.where(
                        (array_src == no_data_src), self.no_data, array_src
                    )

                    array_dst = band_dst.ReadAsArray(locx, locy, cols, rows)
                    ## if DEST raster has data
                    array_to_write = np.where(
                        array_dst != self.no_data,
                        ## if also src has data
                        np.where(
                            array_src != self.no_data,
                            ## if output is higher than small
                            np.where(
                                array_dst > array_src,
                                ## put small
                                array_src,
                                ## leave as was
                                array_dst,
                            ),
                            ## if src is no_data, put data from big
                            array_dst,
                        ),
                        ## if no_data in dest but data in src
                        np.where(
                            array_src != self.no_data,
                            ## leave as src small
                            array_src,
                            ## else put no_data if all is no_data
                            self.no_data,
                        ),
                    )

                    self.dst.GetRasterBand(band_number + 1).WriteArray(
                        array_to_write, locx, locy
                    )
            except Exception as e:
                print(e)
        self.dst = None
        return self.output

    def mosaic_in_order(self):
        """
        Mosaic rasters in order as they are provided in list.
        First raster will be on top, last on the bottom.

        :Example:

        ::

         MosaicRasters(*arg).mosaic_in_order()

        """

        for ix, file in enumerate(reversed(self.rasters)):
            if self.progress:
                print(f"{ix + 1} out of {len(self.rasters)}")

            ## write raster into array
            src = gdal.Open(file)
            for band_number in range(src.RasterCount):
                band_dst = self.dst.GetRasterBand(band_number + 1)

                band = src.GetRasterBand(band_number + 1)
                no_data_src = band.GetNoDataValue()

                cols = src.RasterXSize
                rows = src.RasterYSize

                ## upper left corner
                gt = src.GetGeoTransform()
                # print(gt)
                # print(dst.GetGeoTransform())
                minx = gt[0]
                maxy = gt[3]

                y, x = maxy - self.extent[3], self.extent[0] - minx

                locy, locx = int(abs(y / self.cell_size)), int(abs(x / self.cell_size))

                # print(f'\tRequesting array of size {rows} x {cols} at location {locx} {locy}')

                ## comparing if array already written is smaller or not
                array_src = band.ReadAsArray(0, 0, cols, rows)
                array_src = np.where(
                    (array_src == no_data_src), self.no_data, array_src
                )

                array_dst = band_dst.ReadAsArray(locx, locy, cols, rows)
                ## if DEST raster has data
                array_to_write = np.where(
                    array_dst != self.no_data,
                    ## if also src has data
                    np.where(array_src != self.no_data, array_src, array_dst),
                    ## if no_data in dest but data in src
                    np.where(
                        array_src != self.no_data,
                        ## leave as src small
                        array_src,
                        ## else put no_data if all is no_data
                        self.no_data,
                    ),
                )

                self.dst.GetRasterBand(band_number + 1).WriteArray(
                    array_to_write, locx, locy
                )
        self.dst = None
        return self.output

    def mosaic_sum(self):
        """
        Mosaic rasters as a sum.

        :Example:

        ::

         MosaicRasters(*arg).mosaic_sum()

        """

        for ix, file in enumerate(reversed(self.rasters)):
            if self.progress:
                print(f"{ix + 1} out of {len(self.rasters)}")

            ## write raster into array
            src = gdal.Open(file)
            for band_number in range(src.RasterCount):
                band_dst = self.dst.GetRasterBand(band_number + 1)

                band = src.GetRasterBand(band_number + 1)
                no_data_src = band.GetNoDataValue()
                no_data_dst = self.dst.GetRasterBand(band_number + 1).GetNoDataValue()

                cols = src.RasterXSize
                rows = src.RasterYSize

                ## upper left corner
                gt = src.GetGeoTransform()
                # print(gt)
                # print(dst.GetGeoTransform())
                minx = gt[0]
                maxy = gt[3]

                y, x = maxy - self.extent[3], self.extent[0] - minx

                locy, locx = int(abs(y / self.cell_size)), int(abs(x / self.cell_size))

                print(
                    f"\tRequesting array of size {rows} x {cols} at location {locx} {locy}"
                )

                ## comparing if array already written is smaller or not
                array_src = band.ReadAsArray(0, 0, cols, rows)
                array_src = np.where((array_src == no_data_src), 0, array_src)

                array_dst = band_dst.ReadAsArray(locx, locy, cols, rows)
                array_dst = np.where(array_dst == no_data_dst, 0, array_dst)
                ## if DEST raster has data
                array_to_write = array_dst + array_src
                array_to_write = np.where(
                    array_to_write == 0, self.no_data, array_to_write
                )
                self.dst.GetRasterBand(band_number + 1).WriteArray(
                    array_to_write, locx, locy
                )

        self.dst = None
        return self.output


# =============================================================================
# Defenses package: numba kernels, configuration and run utilities
# =============================================================================


@njit(cache=True)
def _heap_push(keys, vals, size, key, val):
    # 4-ary min-heap on parallel arrays; ties broken by value (cell index) for determinism.
    if size == keys.shape[0]:
        grown_keys = np.empty(2 * keys.shape[0], keys.dtype)
        grown_vals = np.empty(2 * vals.shape[0], vals.dtype)
        grown_keys[:size] = keys[:size]
        grown_vals[:size] = vals[:size]
        keys, vals = grown_keys, grown_vals
    i = size
    while i > 0:
        p = (i - 1) >> 2
        if keys[p] < key or (keys[p] == key and vals[p] <= val):
            break
        keys[i] = keys[p]
        vals[i] = vals[p]
        i = p
    keys[i] = key
    vals[i] = val
    return keys, vals, size + 1


@njit(cache=True)
def _heap_pop(keys, vals, size):
    key, val = keys[0], vals[0]
    size -= 1
    last_key, last_val = keys[size], vals[size]
    i = 0
    while True:
        first = 4 * i + 1
        if first >= size:
            break
        best = first
        for child in range(first + 1, min(first + 4, size)):
            if keys[child] < keys[best] or (
                keys[child] == keys[best] and vals[child] < vals[best]
            ):
                best = child
        if keys[best] < last_key or (keys[best] == last_key and vals[best] < last_val):
            keys[i] = keys[best]
            vals[i] = vals[best]
            i = best
        else:
            break
    keys[i] = last_key
    vals[i] = last_val
    return key, val, size


@njit(cache=True)
def _fifo_push(fifo, tail, value):
    if tail == fifo.shape[0]:
        grown = np.empty(2 * fifo.shape[0], fifo.dtype)
        grown[:tail] = fifo[:tail]
        fifo = grown
    fifo[tail] = value
    return fifo, tail + 1


@njit(cache=True)
def _chamfer_distance(mask, limit):
    # Two-pass chamfer distance (weights 1 and sqrt 2, in cells) to the nearest True cell, capped at limit.
    rows, cols = mask.shape
    d = np.full((rows, cols), limit, np.float32)
    diag = np.float32(sqrt(2.0))
    for i in range(rows):
        for j in range(cols):
            if mask[i, j]:
                d[i, j] = 0.0
                continue
            v = d[i, j]
            if i > 0:
                v = min(v, d[i - 1, j] + 1.0)
                if j > 0:
                    v = min(v, d[i - 1, j - 1] + diag)
                if j < cols - 1:
                    v = min(v, d[i - 1, j + 1] + diag)
            if j > 0:
                v = min(v, d[i, j - 1] + 1.0)
            d[i, j] = v
    for i in range(rows - 1, -1, -1):
        for j in range(cols - 1, -1, -1):
            v = d[i, j]
            if i < rows - 1:
                v = min(v, d[i + 1, j] + 1.0)
                if j > 0:
                    v = min(v, d[i + 1, j - 1] + diag)
                if j < cols - 1:
                    v = min(v, d[i + 1, j + 1] + diag)
            if j < cols - 1:
                v = min(v, d[i, j + 1] + 1.0)
            d[i, j] = v
    return d


@njit(cache=True)
def _breach_to_streams(dem, stream, burn, burn_cells, eps):
    rows, cols = dem.shape
    n = rows * cols
    flat = dem.ravel()
    is_stream = stream.ravel()
    z = np.empty(n, np.float64)
    if burn > 0:
        # Tapered burning: stream cells lowered by burn, decreasing linearly to 0 at burn_cells.
        distance = _chamfer_distance(stream, np.float32(max(burn_cells, 1.0))).ravel()
        for c in range(n):
            taper = 1.0 - distance[c] / max(burn_cells, 1.0)
            z[c] = flat[c] - burn * taper if taper > 0 else flat[c]
    else:
        for c in range(n):
            z[c] = flat[c]
    parent = np.full(n, -1, np.int32)
    seen = np.zeros(n, np.bool_)
    di = np.array([-1, -1, -1, 0, 0, 1, 1, 1])
    dj = np.array([-1, 0, 1, -1, 1, -1, 0, 1])
    dist = np.array([sqrt(2.0), 1.0, sqrt(2.0), 1.0, 1.0, sqrt(2.0), 1.0, sqrt(2.0)])
    keys = np.empty(1 << 16, np.float64)
    vals = np.empty(1 << 16, np.int32)
    size = 0
    fifo = np.empty(1 << 16, np.int32)
    head = 0
    tail = 0
    for c in range(n):
        if is_stream[c] and not np.isnan(z[c]):
            seen[c] = True
            keys, vals, size = _heap_push(keys, vals, size, z[c], np.int32(c))
    while size > 0 or tail > head:
        if tail > head:
            c = fifo[head]
            head += 1
            if head == tail:
                head = 0
                tail = 0
        else:
            _, c, size = _heap_pop(keys, vals, size)
        ci = c // cols
        cj = c - ci * cols
        for k in range(8):
            ni = ci + di[k]
            nj = cj + dj[k]
            if ni < 0 or nj < 0 or ni >= rows or nj >= cols:
                continue
            nb = ni * cols + nj
            if seen[nb] or np.isnan(z[nb]):
                continue
            seen[nb] = True
            parent[nb] = c
            if z[nb] <= z[c]:
                if z[c] - z[nb] < eps:
                    # flat: a small gradient towards the outlet
                    z[nb] = z[c] + eps
                    keys, vals, size = _heap_push(keys, vals, size, z[nb], np.int32(nb))
                    continue
                # pit: carve the path back towards the outlet below the pit level; the pit lies
                # below every queued cell, so it goes to the plain queue (improved priority flood)
                target = z[nb] - eps
                p = c
                while p >= 0 and z[p] > target:
                    z[p] = target
                    target -= eps
                    p = parent[p]
                fifo, tail = _fifo_push(fifo, tail, np.int32(nb))
            else:
                keys, vals, size = _heap_push(keys, vals, size, z[nb], np.int32(nb))

    # D8 steepest descent on the conditioned surface (receivers stored in parent); every
    # receiver is lower, so paths are acyclic and resolved in O(n) with path assignment.
    for c in range(n):
        if not seen[c] or is_stream[c]:
            continue
        ci = c // cols
        cj = c - ci * cols
        best = 0.0
        receiver = parent[c]
        for k in range(8):
            ni = ci + di[k]
            nj = cj + dj[k]
            if ni < 0 or nj < 0 or ni >= rows or nj >= cols:
                continue
            nb = ni * cols + nj
            if not seen[nb]:
                continue
            slope = (z[c] - z[nb]) / dist[k]
            if slope > best:
                best = slope
                receiver = nb
        parent[c] = receiver
    reached = np.full(n, -1, np.int32)
    for c in range(n):
        if is_stream[c] and seen[c]:
            reached[c] = c
    for c in range(n):
        if not seen[c] or reached[c] >= 0:
            continue
        p = c
        while reached[p] < 0:
            p = parent[p]
        target_cell = reached[p]
        q = c
        while reached[q] < 0:
            following = parent[q]
            reached[q] = target_cell
            q = following
    return reached.reshape(rows, cols)


@njit(cache=True)
def _priority_flood(surface, seeds):
    rows, cols = surface.shape
    n = rows * cols
    flat = surface.ravel()
    is_seed = seeds.ravel()
    spill = np.full(n, np.inf, np.float32)
    seed = np.full(n, -1, np.int32)
    order = np.full(n, -1, np.int32)
    keys = np.empty(1 << 16, np.float32)
    vals = np.empty(1 << 16, np.int32)
    size = 0
    fifo = np.empty(1 << 16, np.int32)
    head = 0
    tail = 0
    for c in range(n):
        if is_seed[c] and not np.isnan(flat[c]):
            spill[c] = 0.0
            seed[c] = c
            keys, vals, size = _heap_push(
                keys, vals, size, np.float32(0.0), np.int32(c)
            )
    counter = 0
    while size > 0 or tail > head:
        if tail > head:
            c = fifo[head]
            head += 1
            if head == tail:
                head = 0
                tail = 0
            level = spill[c]
        else:
            level, c, size = _heap_pop(keys, vals, size)
            if level > spill[c]:
                continue
        if order[c] >= 0:
            continue
        order[c] = counter
        counter += 1
        ci = c // cols
        cj = c - ci * cols
        for k in range(4):
            if k == 0:
                if ci == 0:
                    continue
                nb = c - cols
            elif k == 1:
                if ci == rows - 1:
                    continue
                nb = c + cols
            elif k == 2:
                if cj == 0:
                    continue
                nb = c - 1
            else:
                if cj == cols - 1:
                    continue
                nb = c + 1
            value = flat[nb]
            if np.isnan(value) or order[nb] >= 0:
                continue
            if value <= level:
                # Reached at the current, minimal level: final at once, no heap needed.
                if level < spill[nb]:
                    spill[nb] = level
                    seed[nb] = seed[c]
                    fifo, tail = _fifo_push(fifo, tail, np.int32(nb))
            elif value < spill[nb]:
                spill[nb] = value
                seed[nb] = seed[c]
                keys, vals, size = _heap_push(keys, vals, size, value, np.int32(nb))
    return (
        spill.reshape(rows, cols),
        seed.reshape(rows, cols),
        order.reshape(rows, cols),
    )


class Config:
    """
    TOML configuration of the Defenses package. Keys are dotted names; the prefix
    "country." addresses the table [countries.<code>] of the selected country.
    Keys listed in the top-level array "placeholders" are logged as warnings when used.
    """

    # Value kinds: int, float, bool, str, dir, path, paths, opt_path, opt_paths,
    # int_list, float_list, table
    SCHEMA = {
        "general.work_dir": "dir",
        "general.output_dir": "dir",
        "general.n_workers": "int",
        "general.max_ram_gb_per_worker": "float",
        "general.gdal_cachemax_mb": "int",
        "general.warp_memory_limit_mb": "int",
        "general.gdal_threads_per_worker": "int",
        "general.log_level": "str",
        "benchmark.n_windows": "int",
        "benchmark.window_px": "int",
        "benchmark.sequential_rows": "int",
        "benchmark.seed": "int",
        "country.dtm": "paths",
        "country.levees": "paths",
        "country.rivers": "paths",
        "country.river_upstream_area_field": "str",
        "country.upstream_area_raster": "opt_path",
        "country.reference_protected_areas": "opt_paths",
        "country.reference_rp_field": "str",
        "country.overrides": "opt_path",
        "country.working_resolution_m": "float",
        "country.corridor_buffer_m": "float",
        "country.unit_core_size_m": "float",
        "country.unit_halo_m": "float",
        "country.convert_dtm": "bool",
        "country.dtm_conversion": "str",
        "country.dtm_compress": "str",
        "country.footprint_resolution_m": "float",
        "country.bytes_per_cell": "float",
        "country.crest_box_size_m": "float",
        "country.merge_overlap_ratio": "float",
        "country.drop_contradicted_truncated": "bool",
        "country.contradiction_ratio": "float",
        "country.id_rounding_m": "float",
        "country.keep_intermediate": "bool",
        "country.compartments.gap_tolerance_m": "float",
        "country.compartments.touch_tolerance_m": "float",
        "country.compartments.self_close_min_length_factor": "float",
        "country.compartments.crest_search_radius_m": "float",
        "country.compartments.barrier_removal_width_m": "float",
        "country.compartments.min_margin_m": "float",
        "country.compartments.max_level_m": "float",
        "country.compartments.min_area_m2": "float",
        "country.compartments.gated_crossings": "bool",
        "country.compartments.fill_nodata_holes": "bool",
        "country.compartments.stream_burn_m": "float",
        "country.compartments.stream_burn_buffer_m": "float",
        "country.compartments.preserve_river_cells": "bool",
        "country.compartments.cap_river_holes": "bool",
        "country.compartments.max_fill_distance_m": "float",
        "country.rp.classes": "int_list",
        "country.rp.distribution_classes": "int_list",
        "country.rp.distribution_shares": "float_list",
        "country.rp.max_rp": "int",
        "country.rp.cap_levels_m": "float_list",
        "country.rp.cap_rp": "int_list",
        "country.rp.score_weight_value": "float",
        "country.rp.score_weight_upstream": "float",
        "country.rp.levee_adjacent_distance_m": "float",
        "country.rp.storey_height_m": "float",
        "country.rp.max_damage_csv": "path",
        "country.rp.depth_damage_csv": "path",
        "country.rp.value_total_tolerance": "float",
        "country.geofabrik_region": "str",
        "country.osm_pbf": "opt_path",
        "country.endpoint_elevation_radius_m": "float",
        "sources.geofabrik_url_template": "str",
        "sources.osm_tags_hydraulic": "table",
        "sources.osm_tags_drains": "table",
        "sources.osm_tags_closing": "table",
        "sources.osm_tags_pumping": "table",
        "sources.request_pause_s": "float",
        "sources.max_retries": "int",
        "sources.ghsl_url_templates": "table",
        "sources.ghsl_tile_schema": "str",
        "sources.ghsl_tile_field": "str",
        "sources.ghsl_epoch": "int",
        "sources.ghsl_nodata_as_zero": "bool",
        "country.validation_extent": "opt_paths",
        "sources.gfm_stac_url": "str",
        "sources.gfm_collection": "str",
        "sources.gfm_assets": "table",
        "sources.gfm_date_from": "str",
        "sources.gfm_date_to": "str",
        "sources.gfm_page_limit": "int",
        "sources.gfm_token_env": "str",
        "sources.gfm_screen_resolution_m": "float",
        "sources.gfm_download_threads": "int",
        "sources.gfm_flood_values": "int_list",
        "sources.gfm_exclusion_values": "int_list",
        "sources.gfm_reference_water_values": "int_list",
        "sources.gfm_nodata_values": "int_list",
        "sources.glofas_enabled": "bool",
        "sources.glofas_dataset_id": "str",
        "sources.glofas_request": "table",
        "sources.glofas_time_keys": "table",
        "sources.glofas_area_key": "str",
        "sources.glofas_first_year": "int",
        "sources.glofas_last_year": "int",
        "sources.glofas_file_extension": "str",
        "sources.glofas_files": "opt_paths",
        "sources.glofas_discharge_variable": "str",
        "sources.glofas_upstream_area_file": "opt_path",
        "sources.glofas_upstream_area_url": "str",
        "sources.glofas_upstream_area_variable": "str",
        "sources.glofas_upstream_area_to_km2": "float",
        "sources.glofas_min_upstream_area_km2": "float",
        "sources.glofas_min_days_per_year": "int",
        "sources.glofas_valid_time_ends_day": "bool",
        "evidence.comparison_distance_m": "float",
        "evidence.min_valid_pixels": "int",
        "evidence.f_wet_comparison": "float",
        "evidence.f_dry_inside": "float",
        "evidence.f_wet_inside": "float",
        "evidence.flooded_window_days": "int",
        "evidence.glofas_search_radius_km": "float",
        "evidence.upstream_ratio_min": "float",
        "evidence.upstream_ratio_max": "float",
        "evidence.river_level_percentile": "float",
        "evidence.hint_margin_m": "float",
        "evidence.hint_dry_share": "float",
        "evidence.hint_min_scenes": "int",
        "evidence.hint_exclusion_m": "float",
        "country.lulc.enabled": "bool",
        "country.lulc.worldcover": "opt_paths",
        "country.lulc.urban_atlas": "opt_paths",
        "country.lulc.urban_atlas_layer": "str",
        "country.lulc.urban_atlas_class_field": "str",
        "country.lulc.worldcover_categories": "table",
        "country.lulc.urban_atlas_categories": "table",
        "country.lulc.max_damage_csv": "opt_path",
        "country.lulc.depth_damage_csv": "opt_path",
        "sources.worldcover_url_template": "str",
        "sources.worldcover_grid_url": "str",
        "sources.worldcover_tile_field": "str",
        "sources.worldcover_version": "str",
        "sources.worldcover_year": "int",
    }
    ALLOWED = {
        "country.dtm_conversion": {"A", "B"},
        "country.dtm_compress": {"ZSTD", "DEFLATE", "LZW"},
        "general.log_level": {"DEBUG", "INFO", "WARNING", "ERROR"},
    }

    @staticmethod
    def as_list(value):
        """
        :param value: single path, list of paths, or empty string
        :return: list of strings (empty for "" or None)
        """
        if value is None or value == "":
            return []
        if isinstance(value, (list, tuple)):
            return [str(v) for v in value]
        return [str(value)]

    @staticmethod
    def get(cfg, key):
        """
        :param cfg: configuration dict returned by load
        :param key: dotted key
        :return: value (KeyError if missing)
        """
        node = cfg
        for part in key.split("."):
            node = node[part]
        return node

    @staticmethod
    def load(path, country, keys):
        """
        Reads the TOML configuration and validates the keys a script uses: presence,
        type, file existence, allowed values and consistency between related keys.
        Stops with a list of all problems.
        :param path: TOML file
        :param country: country code, a table [countries.<code>] must exist
        :param keys: dotted keys required by the calling script
        :return: dict with the top-level tables, "country" (the country table),
            "country_code" and "used_keys"
        """
        if not exists(str(path)):
            raise SystemExit(f"Configuration file not found: {path}")
        with open(str(path), "rb") as f:
            raw = tomllib.load(f)
        countries = raw.get("countries", {})
        if country not in countries:
            raise SystemExit(
                f"Country '{country}' is not configured in {path}; "
                f"available: {sorted(countries)}"
            )
        cfg = {k: v for k, v in raw.items() if k != "countries"}
        cfg["country"] = countries[country]
        cfg["country_code"] = country
        cfg["used_keys"] = list(keys)
        errors = []
        for key in keys:
            kind = Config.SCHEMA.get(key)
            if kind is None:
                errors.append(f"{key}: not a known configuration key")
                continue
            try:
                value = Config.get(cfg, key)
            except (KeyError, TypeError):
                errors.append(f"{key}: missing")
                continue
            checked, problem = Config._check_value(kind, value)
            if problem:
                errors.append(f"{key}: {problem}")
                continue
            if key in Config.ALLOWED and checked not in Config.ALLOWED[key]:
                errors.append(
                    f"{key}: {checked!r} not in {sorted(Config.ALLOWED[key])}"
                )
            node = cfg
            parts = key.split(".")
            for part in parts[:-1]:
                node = node[part]
            node[parts[-1]] = checked
        if not errors:
            errors = Config._check_consistency(cfg, set(keys))
        if errors:
            raise SystemExit(
                f"Invalid configuration {path} (country {country}):\n  - "
                + "\n  - ".join(errors)
            )
        return cfg

    @staticmethod
    def _check_value(kind, value):
        numeric = (int, float)
        if kind == "int":
            if isinstance(value, bool) or not isinstance(value, int):
                return value, f"expected an integer, got {value!r}"
            return value, None if value >= 0 else "must not be negative"
        if kind == "float":
            if isinstance(value, bool) or not isinstance(value, numeric):
                return value, f"expected a number, got {value!r}"
            return float(value), None if value >= 0 else "must not be negative"
        if kind == "bool":
            return value, (
                None
                if isinstance(value, bool)
                else f"expected true/false, got {value!r}"
            )
        if kind == "str":
            return value, (
                None if isinstance(value, str) else f"expected a string, got {value!r}"
            )
        if kind == "dir":
            if not isinstance(value, str) or not value:
                return value, "expected a non-empty directory path"
            return value, None
        if kind in ("path", "opt_path"):
            if not isinstance(value, str):
                return value, f"expected a path string, got {value!r}"
            if not value:
                return value, None if kind == "opt_path" else "empty path"
            return value, None if exists(value) else f"file not found: {value}"
        if kind in ("paths", "opt_paths"):
            if not isinstance(value, (str, list)):
                return value, "expected a path or a list of paths"
            paths = Config.as_list(value)
            if not paths and kind == "paths":
                return paths, "at least one path is required"
            missing = [p for p in paths if not exists(p)]
            return paths, f"file(s) not found: {missing}" if missing else None
        if kind in ("int_list", "float_list"):
            if not isinstance(value, list) or not all(
                isinstance(v, numeric) and not isinstance(v, bool) for v in value
            ):
                return value, f"expected a list of numbers, got {value!r}"
            if kind == "int_list" and not all(isinstance(v, int) for v in value):
                return value, f"expected a list of integers, got {value!r}"
            return value, None
        if kind == "table":
            return value, None if isinstance(value, dict) else "expected a table"
        return value, f"unknown kind {kind}"

    @staticmethod
    def _check_consistency(cfg, keys):
        errors = []
        c = cfg["country"]
        if "country.working_resolution_m" in keys and c["working_resolution_m"] <= 0:
            errors.append("country.working_resolution_m must be positive")
        if {"country.unit_core_size_m", "country.working_resolution_m"} <= keys:
            if c["unit_core_size_m"] < 2 * c["working_resolution_m"]:
                errors.append(
                    "country.unit_core_size_m must be at least two working cells"
                )
        if {"country.unit_halo_m", "country.corridor_buffer_m"} <= keys:
            if c["unit_halo_m"] < c["corridor_buffer_m"]:
                # A core cell may lie corridor_buffer_m from its river; HAND needs that river inside the window.
                errors.append(
                    "country.unit_halo_m must be at least country.corridor_buffer_m, otherwise "
                    "cells of a core can have their river outside the unit window"
                )
        if "general.n_workers" in keys and cfg["general"]["n_workers"] < 1:
            errors.append("general.n_workers must be at least 1")
        if (
            "country.merge_overlap_ratio" in keys
            and not 0 < c["merge_overlap_ratio"] <= 1
        ):
            errors.append("country.merge_overlap_ratio must be in (0, 1]")
        if "benchmark.window_px" in keys and cfg["benchmark"]["window_px"] < 1:
            errors.append("benchmark.window_px must be positive")
        rp = c.get("rp", {})
        if {
            "country.rp.distribution_classes",
            "country.rp.distribution_shares",
        } <= keys:
            if len(rp["distribution_classes"]) != len(rp["distribution_shares"]):
                errors.append(
                    "rp.distribution_classes and rp.distribution_shares differ in length"
                )
            elif abs(sum(rp["distribution_shares"]) - 1.0) > 1e-6:
                errors.append("rp.distribution_shares must sum to 1")
        if {"country.rp.cap_levels_m", "country.rp.cap_rp"} <= keys:
            if len(rp["cap_levels_m"]) != len(rp["cap_rp"]):
                errors.append("rp.cap_levels_m and rp.cap_rp differ in length")
            elif list(rp["cap_levels_m"]) != sorted(rp["cap_levels_m"]):
                errors.append("rp.cap_levels_m must be ascending")
        if "country.rp.classes" in keys and list(rp["classes"]) != sorted(
            set(rp["classes"])
        ):
            errors.append("rp.classes must be strictly ascending")
        if {"country.rp.classes", "country.rp.distribution_classes"} <= keys:
            if not set(rp["distribution_classes"]) <= set(rp["classes"]):
                errors.append("rp.distribution_classes must be values of rp.classes")
        if {"country.rp.classes", "country.rp.cap_rp", "country.rp.max_rp"} <= keys:
            if not set(rp["cap_rp"]) | {rp["max_rp"]} <= set(rp["classes"]):
                errors.append("rp.cap_rp and rp.max_rp must be values of rp.classes")
        evidence = cfg.get("evidence", {})
        if {"evidence.f_dry_inside", "evidence.f_wet_inside"} <= keys:
            if not evidence["f_dry_inside"] < evidence["f_wet_inside"] <= 1:
                errors.append("evidence: f_dry_inside < f_wet_inside <= 1 is required")
        if {"evidence.upstream_ratio_min", "evidence.upstream_ratio_max"} <= keys:
            if (
                not 0
                < evidence["upstream_ratio_min"]
                <= 1
                <= evidence["upstream_ratio_max"]
            ):
                errors.append(
                    "evidence: upstream_ratio_min <= 1 <= upstream_ratio_max is required"
                )
        if "evidence.river_level_percentile" in keys:
            if not 0 < evidence["river_level_percentile"] <= 100:
                errors.append("evidence.river_level_percentile must be in (0, 100]")
        if (
            "evidence.hint_dry_share" in keys
            and not 0 < evidence["hint_dry_share"] <= 1
        ):
            errors.append("evidence.hint_dry_share must be in (0, 1]")
        if "evidence.hint_exclusion_m" in keys and evidence["hint_exclusion_m"] < 0:
            errors.append("evidence.hint_exclusion_m must not be negative")
        if (
            "sources.gfm_download_threads" in keys
            and cfg["sources"]["gfm_download_threads"] < 1
        ):
            errors.append("sources.gfm_download_threads must be at least 1")
        if (
            "sources.gfm_screen_resolution_m" in keys
            and cfg["sources"]["gfm_screen_resolution_m"] <= 0
        ):
            errors.append("sources.gfm_screen_resolution_m must be positive")
        lulc = c.get("lulc", {})
        for table in ["worldcover_categories", "urban_atlas_categories"]:
            if f"country.lulc.{table}" not in keys:
                continue
            for code, category in lulc[table].items():
                if table == "worldcover_categories" and not (
                    str(code).isdigit() and 0 < int(code) < 256
                ):
                    errors.append(
                        f"lulc.worldcover_categories: class {code!r} is not a raster value 1-255"
                    )
                if not (
                    isinstance(category, str)
                    and category.replace("_", "").isalnum()
                    and category.isascii()
                ):
                    errors.append(
                        f"lulc.{table}: category {category!r} must be a name of letters, digits and _"
                    )
        if {"country.lulc.urban_atlas", "country.lulc.urban_atlas_class_field"} <= keys:
            if lulc["urban_atlas"] and not lulc["urban_atlas_class_field"]:
                errors.append(
                    "lulc.urban_atlas is set but lulc.urban_atlas_class_field is empty"
                )
        for date_key in ["sources.gfm_date_from", "sources.gfm_date_to"]:
            value = cfg.get("sources", {}).get(date_key.split(".")[1], "")
            if date_key in keys and value:
                try:
                    datetime.date.fromisoformat(value)
                except ValueError:
                    errors.append(f"{date_key} must be YYYY-MM-DD or empty")
        if {"sources.glofas_first_year", "sources.glofas_last_year"} <= keys:
            if cfg["sources"]["glofas_first_year"] > cfg["sources"]["glofas_last_year"]:
                errors.append(
                    "sources.glofas_first_year must not exceed glofas_last_year"
                )
        sources = cfg.get("sources", {})
        for subset in ["sources.osm_tags_closing", "sources.osm_tags_pumping"]:
            if {subset, "sources.osm_tags_hydraulic"} <= keys:
                for key, values in sources[subset.split(".")[1]].items():
                    known = sources["osm_tags_hydraulic"].get(key, [])
                    if not set(values) <= set(known):
                        errors.append(
                            f"{subset}.{key} values {sorted(set(values) - set(known))} "
                            "are not in sources.osm_tags_hydraulic"
                        )
        return errors

    @staticmethod
    def warn_placeholders(cfg, keys=None):
        """
        Logs a warning for every used key listed in the top-level array placeholders
        :param cfg: configuration dict returned by load
        :param keys: keys to check (default: the keys validated by load)
        """
        placeholders = set(cfg.get("placeholders", []))
        for key in keys or cfg["used_keys"]:
            if key in placeholders:
                log.warning(
                    "Placeholder value in use: %s = %r", key, Config.get(cfg, key)
                )


class Run:
    @staticmethod
    def env_overrides(**variables):
        """
        Run variables set at the top of a script, overridden by environment variables
        DEFENSES_<NAME> (used by tests and batch runs). Types follow the defaults:
        bool accepts 1/0/true/false, list accepts comma-separated values.
        :return: dict name -> value in the given order
        """
        values = {}
        for name, default in variables.items():
            raw = os.environ.get(f"DEFENSES_{name}")
            if raw is None:
                values[name] = default
            elif isinstance(default, bool):
                values[name] = raw.strip().lower() in ("1", "true", "yes")
            elif isinstance(default, list):
                values[name] = [v.strip() for v in raw.split(",") if v.strip()]
            else:
                values[name] = type(default)(raw)
        return values

    @staticmethod
    def logger(name, log_dir=None, level="INFO"):
        """
        Console logging plus an optional run log file <log_dir>/<name>_<timestamp>.log
        :return: logging.Logger
        """
        root = logging.getLogger()
        root.setLevel(level)
        formatter = logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s")
        if not any(type(h) is logging.StreamHandler for h in root.handlers):
            console = logging.StreamHandler()
            console.setFormatter(formatter)
            root.addHandler(console)
        if log_dir is not None:
            os.makedirs(str(log_dir), exist_ok=True)
            path = join(str(log_dir), f"{name}_{time.strftime('%Y%m%d_%H%M%S')}.log")
            file_handler = logging.FileHandler(path, encoding="utf-8")
            file_handler.setFormatter(formatter)
            root.addHandler(file_handler)
        logging.captureWarnings(True)
        return logging.getLogger(name)

    @staticmethod
    def peak_rss_mb():
        """
        Peak resident memory of the current process in MB (Linux, macOS, Windows)
        """
        if sys.platform == "win32":

            class Counters(ctypes.Structure):
                _fields_ = [
                    ("cb", ctypes.c_ulong),
                    ("PageFaultCount", ctypes.c_ulong),
                    ("PeakWorkingSetSize", ctypes.c_size_t),
                    ("WorkingSetSize", ctypes.c_size_t),
                    ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                    ("PagefileUsage", ctypes.c_size_t),
                    ("PeakPagefileUsage", ctypes.c_size_t),
                ]

            counters = Counters()
            counters.cb = ctypes.sizeof(Counters)
            kernel32 = ctypes.windll.kernel32
            kernel32.GetCurrentProcess.restype = ctypes.c_void_p
            psapi = ctypes.windll.psapi
            psapi.GetProcessMemoryInfo.argtypes = [
                ctypes.c_void_p,
                ctypes.POINTER(Counters),
                ctypes.c_ulong,
            ]
            psapi.GetProcessMemoryInfo(
                kernel32.GetCurrentProcess(), ctypes.byref(counters), counters.cb
            )
            return counters.PeakWorkingSetSize / 2**20
        import resource

        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return peak / 2**20 if sys.platform == "darwin" else peak / 1024

    @staticmethod
    def configure_gdal(cache_mb, threads):
        """
        GDAL block cache and thread count of the current process
        :param cache_mb: GDAL_CACHEMAX in MB
        :param threads: GDAL_NUM_THREADS (int or "ALL_CPUS")
        """
        gdal.SetCacheMax(int(cache_mb) * 2**20)
        gdal.SetConfigOption("GDAL_NUM_THREADS", str(threads))

    @staticmethod
    def format_seconds(seconds):
        """
        :return: seconds as h:mm:ss
        """
        if seconds is None or not np.isfinite(seconds):
            return "?"
        seconds = int(round(seconds))
        return f"{seconds // 3600}:{seconds % 3600 // 60:02d}:{seconds % 60:02d}"

    @staticmethod
    def run_units(func, tasks, n_workers, gdal_cache_mb, gdal_threads, weights, logger):
        """
        Runs func(task) for every task in a pool of fresh processes (one process per
        task, so peak memory is measured per task and freed after it). Tasks are
        submitted in the given order; func returns a dict with a "summary" string.
        :param tasks: list of picklable dicts with key unit_id
        :param weights: relative cost per task, used for the ETA
        :return: (list of results, list of failed unit ids)
        """
        results, failed = [], []
        if not tasks:
            return results, failed
        total = float(sum(weights)) or 1.0
        finished = 0.0
        start = time.perf_counter()
        context = multiprocessing.get_context("spawn")
        with ProcessPoolExecutor(
            max_workers=n_workers,
            mp_context=context,
            max_tasks_per_child=1,
            initializer=Run.configure_gdal,
            initargs=(gdal_cache_mb, gdal_threads),
        ) as pool:
            futures = {
                pool.submit(func, task): (task, w) for task, w in zip(tasks, weights)
            }
            for future in as_completed(futures):
                task, weight = futures[future]
                finished += weight
                try:
                    result = future.result()
                    results.append(result)
                    message = result.get("summary", "")
                except BrokenProcessPool as error:
                    failed.append(task["unit_id"])
                    message = f"worker process died ({error}); possibly out of memory"
                    logger.error("Unit %s: %s", task["unit_id"], message)
                except Exception:
                    failed.append(task["unit_id"])
                    logger.exception("Unit %s failed", task["unit_id"])
                    message = "failed"
                elapsed = time.perf_counter() - start
                eta = (
                    elapsed / finished * (total - finished)
                    if finished
                    else float("nan")
                )
                logger.info(
                    "[%d/%d] %s | elapsed %s, ETA %s",
                    len(results) + len(failed),
                    len(tasks),
                    message,
                    Run.format_seconds(elapsed),
                    Run.format_seconds(eta),
                )
        return results, failed


class Http:
    @staticmethod
    def get_json(
        url, params=None, body=None, retries=4, pause_s=1.0, headers=None, timeout_s=120
    ):
        """
        JSON request with retries and exponential backoff: GET with query parameters, or
        POST of a JSON body (e.g. STAC item search and its "next" links)
        :param url: endpoint
        :param params: query parameters for GET (lists are joined with commas)
        :param body: JSON body; sent with POST when given
        :return: decoded JSON
        """
        headers = dict(headers or {})
        if params:
            query = urllib.parse.urlencode(
                {
                    k: ",".join(map(str, v)) if isinstance(v, (list, tuple)) else v
                    for k, v in params.items()
                }
            )
            url = f"{url}{'&' if '?' in url else '?'}{query}"
        data = None
        if body is not None:
            data = json.dumps(body).encode()
            headers.setdefault("Content-Type", "application/json")
        for attempt in range(retries + 1):
            try:
                request = urllib.request.Request(url, data=data, headers=headers)
                with urllib.request.urlopen(request, timeout=timeout_s) as response:
                    return json.loads(response.read().decode())
            except (urllib.error.URLError, OSError, ValueError) as error:
                if attempt == retries:
                    raise RuntimeError(f"Request to {url} failed: {error}") from error
                wait = pause_s * 2**attempt
                log.warning(
                    "Request to %s failed (%s); retry in %.0f s.", url, error, wait
                )
                time.sleep(wait)
        return None

    @staticmethod
    def download(url, path, retries=4, pause_s=1.0, headers=None, timeout_s=120):
        """
        Downloads url to path through a temporary file. Skipped when path exists with
        the size the server reports (Content-Length). Retries with exponential backoff.
        file:// URLs are supported (tests, local mirrors).
        :param url: source URL
        :param path: destination file
        :param retries: number of retries after the first attempt
        :param pause_s: base pause; attempt k waits pause_s * 2**k seconds
        :param headers: optional request headers (e.g. authorisation from an env variable)
        :return: path
        """
        headers = headers or {}
        expected = None
        try:
            request = urllib.request.Request(url, method="HEAD", headers=headers)
            with urllib.request.urlopen(request, timeout=timeout_s) as response:
                length = response.headers.get("Content-Length")
                expected = int(length) if length else None
        except (urllib.error.URLError, OSError, ValueError) as error:
            log.info(
                "No size from %s (%s); downloading without size check.", url, error
            )
        if (
            exists(str(path))
            and expected is not None
            and os.path.getsize(path) == expected
        ):
            log.info("%s exists with the expected size; download skipped.", path)
            return str(path)
        os.makedirs(os.path.dirname(str(path)) or ".", exist_ok=True)
        temporary = f"{path}.part"
        for attempt in range(retries + 1):
            try:
                request = urllib.request.Request(url, headers=headers)
                with urllib.request.urlopen(request, timeout=timeout_s) as response:
                    with open(temporary, "wb") as target:
                        shutil.copyfileobj(response, target, 2**22)
                size = os.path.getsize(temporary)
                if expected is not None and size != expected:
                    raise OSError(f"incomplete download ({size} of {expected} bytes)")
                os.replace(temporary, str(path))
                log.info("Downloaded %s (%.1f MB).", url, size / 1e6)
                return str(path)
            except (urllib.error.URLError, OSError) as error:
                if attempt == retries:
                    raise RuntimeError(f"Download of {url} failed: {error}") from error
                wait = pause_s * 2**attempt
                log.warning(
                    "Download of %s failed (%s); retry in %.0f s.", url, error, wait
                )
                time.sleep(wait)
        return str(path)


class Stats:
    @staticmethod
    def gumbel_lmoments(sample):
        """
        Gumbel (EV1) parameters from L-moments: alpha = l2 / ln 2, xi = l1 - gamma * alpha
        (gamma = Euler-Mascheroni constant, 0.5772...)
        :param sample: annual maxima
        :return: (xi, alpha); NaN when fewer than two finite values
        """
        x = np.sort(np.asarray(sample, dtype=float)[np.isfinite(sample)])
        n = len(x)
        if n < 2:
            return float("nan"), float("nan")
        b0 = x.mean()
        b1 = np.sum(np.arange(n) / (n - 1) * x) / n
        l2 = 2 * b1 - b0
        alpha = l2 / np.log(2.0)
        return float(b0 - np.euler_gamma * alpha), float(alpha)

    @staticmethod
    def gumbel_rp(q, xi, alpha):
        """
        Return period of q under a Gumbel distribution: 1 / (1 - exp(-exp(-(q - xi) / alpha)))
        :return: float or array (inf where the non-exceedance probability rounds to 1)
        """
        with np.errstate(over="ignore", divide="ignore"):
            exceedance = -np.expm1(-np.exp(-(np.asarray(q, dtype=float) - xi) / alpha))
            return 1.0 / exceedance

    @staticmethod
    def percentile_rank(values):
        """
        Percentile rank in (0, 1] (average rank for ties); NaN stays NaN, -inf ranks lowest
        :param values: array-like
        :return: numpy array
        """
        return (
            pd.Series(np.asarray(values, dtype=float))
            .rank(method="average", pct=True)
            .to_numpy()
        )

    @staticmethod
    def quantile_classes(scores, classes, shares):
        """
        Assigns classes by score so that the class shares match shares (ascending scores
        get ascending classes); ties get the same class. NaN scores get NaN.
        :param scores: array-like
        :param classes: class values, ascending
        :param shares: share per class, summing to 1
        :return: float numpy array of class values
        """
        rank = Stats.percentile_rank(scores)
        bounds = np.cumsum(shares)
        bounds[-1] = 1.0
        index = np.searchsorted(bounds, rank - 1e-9, side="left")
        result = np.full(len(rank), np.nan)
        known = ~np.isnan(rank)
        result[known] = np.asarray(classes, dtype=float)[
            np.minimum(index[known], len(classes) - 1)
        ]
        return result
