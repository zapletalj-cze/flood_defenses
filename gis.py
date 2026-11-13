# Standard library imports
import ctypes
import datetime
import heapq
import json
import logging
import multiprocessing
import os
import re
import shutil
import sys
import time
import tomllib
import urllib.error
import urllib.parse
import urllib.request
import warnings
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
        list_file = ds.GetFileList() or [str(in_raster)]
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
            "file_size_gb": sum(os.path.getsize(f) for f in list_file if exists(f))
            / 1e9,
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
        list_opt = [
            "TILED=YES",
            "BLOCKXSIZE=256",
            "BLOCKYSIZE=256",
            f"COMPRESS={Raster.best_compression(compress)}",
            f"PREDICTOR={3 if is_float else 2}",
            "BIGTIFF=IF_SAFER",
            "NUM_THREADS=ALL_CPUS",
        ]
        gdal.Translate(str(path), mem, format="GTiff", creationOptions=list_opt)

    @staticmethod
    def sample_points(path, x, y, block=2048):
        """
        Values of a raster at points (the cell containing each point) at native resolution. The raster is read
        block by block, only where points lie, so scattered points over a large raster need little memory.
        :param path: raster path
        :param x: x coordinates in the raster CRS
        :param y: y coordinates in the raster CRS
        :param block: block size in pixels
        :return: float64 array, NaN outside the raster or at nodata
        """
        ds = gdal.Open(str(path))
        band = ds.GetRasterBand(1)
        gt, nodata = ds.GetGeoTransform(), band.GetNoDataValue()
        n_cols, n_rows = ds.RasterXSize, ds.RasterYSize
        col = np.floor((np.asarray(x, float) - gt[0]) / gt[1]).astype(np.int64)
        row = np.floor((np.asarray(y, float) - gt[3]) / gt[5]).astype(np.int64)
        out = np.full(len(col), np.nan)
        inside = np.flatnonzero(
            (col >= 0) & (col < n_cols) & (row >= 0) & (row < n_rows)
        )
        if len(inside):
            key = (row[inside] // block) * (n_cols // block + 1) + col[inside] // block
            order = np.argsort(key, kind="stable")
            inside, key = inside[order], key[order]
            first = np.flatnonzero(np.r_[True, key[1:] != key[:-1]])
            for idx in np.split(inside, first[1:]):
                r0 = int(row[idx[0]] // block) * block
                c0 = int(col[idx[0]] // block) * block
                arr = band.ReadAsArray(
                    c0, r0, min(block, n_cols - c0), min(block, n_rows - r0)
                ).astype(np.float64)
                values = arr[row[idx] - r0, col[idx] - c0]
                if nodata is not None:
                    values[values == nodata] = np.nan
                out[idx] = values
        ds = None
        return out

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
            list_opt = ["8CONNECTED=8"] if connectedness == 8 else []
            gdal.Polygonize(
                label_ds.GetRasterBand(1), mask_ds.GetRasterBand(1), layer, 0, list_opt
            )
        list_feat = [
            (
                f.GetField("label"),
                shapely.from_wkb(bytes(f.GetGeometryRef().ExportToWkb())),
            )
            for f in layer
        ]
        gdf = gpd.GeoDataFrame(
            {"label": pd.Series([f[0] for f in list_feat], dtype="int64")},
            geometry=gpd.GeoSeries([f[1] for f in list_feat], crs=wkt),
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
        gdf_poly = Raster.polygonize_labels(
            valid.astype(np.int32), coarse_gt, ds.GetProjection()
        )
        return gpd.GeoDataFrame(
            geometry=[shapely.union_all(gdf_poly.geometry.values)],
            crs=ds.GetProjection(),
        )

    @staticmethod
    def benchmark_reads(
        info,
        polygon,
        n_windows,
        window_px,
        resolution,
        warp_memory_mb,
        sequential_rows,
        compress,
        rng,
    ):
        """
        Read benchmark of one raster (DTM or DSM): random native windows inside the
        valid-data polygon, the same ground windows warped to the working resolution,
        one sequential full-width read, and the compression ratio and speed of a tiled
        copy of one window. The first warp is repeated single-threaded with GDAL debug
        output; the source window it reports shows whether an overview was used.
        :param info: dict from Raster.inspect
        :param polygon: shapely polygon of valid data in the raster CRS
        :param n_windows: number of windows
        :param window_px: window size in native pixels
        :param resolution: working resolution of the warped reads
        :param warp_memory_mb: warp working buffer in MB
        :param sequential_rows: rows of the sequential read
        :param compress: compression of the tiled copy (ZSTD, DEFLATE, LZW)
        :param rng: numpy random Generator
        :return: (list of per-window records, summary dict, list of source window
            widths from the debug output)
        """
        ds = gdal.Open(info["path"])
        band = ds.GetRasterBand(1)
        gt = ds.GetGeoTransform()
        name = os.path.basename(info["path"])
        win_x, win_y = min(window_px, info["width"]), min(window_px, info["height"])
        bytes_window = win_x * win_y * info["bytes_per_pixel"]
        # Window centres drawn inside the valid footprint, so that the windows read real data.
        xs = rng.uniform(
            gt[0] + win_x / 2 * gt[1],
            gt[0] + (info["width"] - win_x / 2) * gt[1],
            200 * n_windows,
        )
        ys = rng.uniform(
            gt[3] + (info["height"] - win_y / 2) * gt[5],
            gt[3] + win_y / 2 * gt[5],
            200 * n_windows,
        )
        inside = shapely.contains_xy(polygon, xs, ys)
        list_ctr = list(zip(xs[inside], ys[inside]))[:n_windows]
        list_msg, list_rec = [], []
        sample_window = (0, 0)
        for k, (x, y) in enumerate(list_ctr):
            col0 = int(
                np.clip((x - gt[0]) / gt[1] - win_x / 2, 0, info["width"] - win_x)
            )
            row0 = int(
                np.clip((y - gt[3]) / gt[5] - win_y / 2, 0, info["height"] - win_y)
            )
            sample_window = sample_window if k else (col0, row0)
            start = time.perf_counter()
            native = band.ReadAsArray(col0, row0, win_x, win_y)
            t_native = time.perf_counter() - start
            valid_share = (
                float(np.mean(native != info["nodata"]))
                if info["nodata"] is not None
                else 1.0
            )
            # Same ground window at the working resolution, snapped outward to the working grid.
            xmin, ymax = gt[0] + col0 * gt[1], gt[3] + row0 * gt[5]
            bounds = (
                xmin,
                ymax - np.ceil(win_y * abs(gt[5]) / resolution) * resolution,
                xmin + np.ceil(win_x * gt[1] / resolution) * resolution,
                ymax,
            )
            if k == 0:
                gdal.SetConfigOption("CPL_DEBUG", "ON")
                gdal.PushErrorHandler(
                    lambda err_class, err_no, msg: list_msg.append(msg)
                )
                try:
                    Raster.read_window(
                        info["path"],
                        bounds,
                        resolution,
                        "average",
                        warp_memory_mb,
                        multithread=False,
                    )
                finally:
                    gdal.PopErrorHandler()
                    gdal.SetConfigOption("CPL_DEBUG", None)
            start = time.perf_counter()
            warped, _ = Raster.read_window(
                info["path"], bounds, resolution, "average", warp_memory_mb
            )
            t_warp = time.perf_counter() - start
            list_rec.append(
                {
                    "raster": name,
                    "window": k,
                    "col0": col0,
                    "row0": row0,
                    "valid_share": valid_share,
                    "native_s": t_native,
                    "native_mb_s": bytes_window / 1e6 / t_native,
                    "warp_s": t_warp,
                    "warp_native_equivalent_mb_s": bytes_window / 1e6 / t_warp,
                    "warp_output_cells": warped.size,
                }
            )
        list_sw = [
            int(m.group(1))
            for m in map(re.compile(r"Src=\d+,\d+,(\d+)x(\d+)").search, list_msg)
            if m
        ]
        # Sequential full-width read (pattern of a conversion pass).
        n_rows = min(sequential_rows, info["height"])
        row0 = int(rng.integers(0, info["height"] - n_rows + 1))
        start = time.perf_counter()
        band.ReadAsArray(0, row0, info["width"], n_rows)
        t_seq = time.perf_counter() - start
        seq_mb_s = info["width"] * n_rows * info["bytes_per_pixel"] / 1e6 / t_seq
        # Compression ratio and speed of the target layout, measured on one native window.
        sample = band.ReadAsArray(sample_window[0], sample_window[1], win_x, win_y)
        is_float = np.issubdtype(sample.dtype, np.floating)
        mem = gdal.GetDriverByName("MEM").Create("", win_x, win_y, 1, band.DataType)
        mem.GetRasterBand(1).WriteArray(sample)
        if info["nodata"] is not None:
            mem.GetRasterBand(1).SetNoDataValue(info["nodata"])
        start = time.perf_counter()
        gdal.Translate(
            "/vsimem/compression_sample.tif",
            mem,
            format="GTiff",
            creationOptions=[
                "TILED=YES",
                "BLOCKXSIZE=512",
                "BLOCKYSIZE=512",
                f"COMPRESS={compress}",
                f"PREDICTOR={3 if is_float else 2}",
                "NUM_THREADS=ALL_CPUS",
            ],
        )
        t_compress = time.perf_counter() - start
        ratio = gdal.VSIStatL("/vsimem/compression_sample.tif").size / bytes_window
        gdal.Unlink("/vsimem/compression_sample.tif")
        df_rec = pd.DataFrame(list_rec)
        summary = {
            "path": info["path"],
            "windows": len(list_rec),
            "native_s_median": (
                float(df_rec["native_s"].median()) if len(df_rec) else np.nan
            ),
            "native_s_p90": (
                float(df_rec["native_s"].quantile(0.9)) if len(df_rec) else np.nan
            ),
            "native_mb_s_median": (
                float(df_rec["native_mb_s"].median()) if len(df_rec) else np.nan
            ),
            "warp_s_median": (
                float(df_rec["warp_s"].median()) if len(df_rec) else np.nan
            ),
            "warp_s_p90": (
                float(df_rec["warp_s"].quantile(0.9)) if len(df_rec) else np.nan
            ),
            "warp_native_equivalent_mb_s_median": (
                float(df_rec["warp_native_equivalent_mb_s"].median())
                if len(df_rec)
                else np.nan
            ),
            "window_ground_km2": win_x * win_y * abs(gt[1] * gt[5]) / 1e6,
            "sequential_mb_s": seq_mb_s,
            "compression_ratio": ratio,
            "compression_mb_s": bytes_window / 1e6 / t_compress,
            "raw_gb": info["raw_size_gb"],
            "window_px": win_x,
        }
        return list_rec, summary, list_sw

    @staticmethod
    def map_codes(values, mapping, default=0):
        """
        Replaces raster values by codes (e.g. CORINE class index -> three-digit code)
        :param values: integer or float array (NaN or values outside the mapping give default)
        :param mapping: dict raster value -> code
        :param default: code for values without a mapping
        :return: int32 array of codes
        """
        values = np.asarray(values)
        keys = np.array([int(k) for k in mapping], dtype=np.int64)
        codes = np.array([int(v) for v in mapping.values()], dtype=np.int32)
        out = np.full(values.shape, default, np.int32)
        known = (
            np.isfinite(values)
            if values.dtype.kind == "f"
            else np.ones(values.shape, bool)
        )
        if len(keys):
            ints = np.where(known, values, -1).astype(np.int64)
            lookup = np.full(int(max(keys.max(), 0)) + 2, default, np.int32)
            lookup[keys[keys >= 0]] = codes[keys >= 0]
            inside = known & (ints >= 0) & (ints < len(lookup))
            out[inside] = lookup[ints[inside]]
        return out

    @staticmethod
    def hillshade(array, geotransform, wkt, z_factor=1.0, azimuth=315.0, altitude=45.0):
        """
        Hillshade of an elevation array with gdal.DEMProcessing to MEM
        :param array: 2D elevation array, NaN = nodata
        :param geotransform: GDAL geotransform
        :param wkt: CRS as WKT
        :param z_factor: vertical exaggeration
        :param azimuth: light direction in degrees from north
        :param altitude: light elevation in degrees
        :return: uint8 array, 0 = nodata
        """
        mem = Raster.array_to_mem(
            np.asarray(array, dtype=np.float32), geotransform, wkt
        )
        out = gdal.DEMProcessing(
            "",
            mem,
            "hillshade",
            format="MEM",
            zFactor=z_factor,
            azimuth=azimuth,
            altitude=altitude,
            computeEdges=True,
        )
        return out.GetRasterBand(1).ReadAsArray()


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
            list_name = [ds.GetLayer(i).GetName() for i in range(ds.GetLayerCount())]
            raise SystemExit(
                f"{path} contains {ds.GetLayerCount()} layers {list_name}; every vector input "
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
        list_gdf = []
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
                list_col = [c for c in gdf.columns if c not in ("source", "geometry")]
            else:
                list_col = [c for c in (keep_columns or []) if c in gdf.columns]
            list_gdf.append(gdf[list_col + ["source", "geometry"]])
        gdf = gpd.GeoDataFrame(
            pd.concat(list_gdf, ignore_index=True), geometry="geometry", crs=wkt
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
    def zones_to_grid(zones, order, shape, geotransform, wkt, field="context"):
        """
        Context of every cell of a grid from non-overlapping zone polygons
        :param zones: gpd.GeoDataFrame with the context name in field
        :param order: context names; cell value = position in order + 1, 0 = no zone
        :param shape: (rows, cols)
        :param geotransform: GDAL geotransform of the grid
        :param wkt: CRS of the grid as WKT
        :return: uint8 array
        """
        if zones.empty:
            return np.zeros(shape, np.uint8)
        zones = zones.to_crs(wkt)
        values = zones[field].map({name: k for k, name in enumerate(order, start=1)})
        if values.isna().any():
            raise SystemExit(
                f"Zones with unknown contexts {sorted(set(zones.loc[values.isna(), field]))}; "
                f"expected {list(order)}."
            )
        # Later contexts are burnt first, so earlier ones win if polygons ever overlap.
        gdf_rank = zones.assign(value=values.astype(int)).sort_values(
            "value", ascending=False
        )
        return Raster.rasterize_geometries(
            gdf_rank.geometry,
            gdf_rank["value"],
            shape,
            geotransform,
            wkt,
            False,
            np.uint8,
        )

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

    @staticmethod
    def chaikin(line, iterations=2, max_shift=float("inf")):
        """
        Chaikin corner cutting of a line; the end points stay. A vertex whose corner
        would move by more than max_shift is kept unsmoothed in every iteration. The
        shift of a corner is 1/8 of |previous - 2 * vertex + next|: the distance from the
        vertex to the midpoint of its cut chord, which the limit curve also passes.
        :param line: shapely LineString
        :param iterations: number of corner-cutting passes
        :param max_shift: largest allowed movement of a vertex
        :return: smoothed LineString
        """
        xy = shapely.get_coordinates(line)
        fixed = np.zeros(len(xy), bool)
        fixed[[0, -1]] = True
        for _ in range(iterations):
            if len(xy) < 3:
                break
            before, vertex, after = xy[:-2], xy[1:-1], xy[2:]
            shift = 0.125 * np.hypot(*(before - 2 * vertex + after).T)
            keep = fixed[1:-1] | (shift > max_shift)
            counts = np.where(keep, 1, 2)
            out = np.empty((counts.sum() + 2, 2))
            out_fixed = np.zeros(len(out), bool)
            out[0], out[-1] = xy[0], xy[-1]
            out_fixed[[0, -1]] = True
            first = 1 + np.concatenate([[0], np.cumsum(counts)[:-1]])
            out[first[keep]] = vertex[keep]
            out_fixed[first[keep]] = True
            cut = ~keep
            out[first[cut]] = 0.75 * vertex[cut] + 0.25 * before[cut]
            out[first[cut] + 1] = 0.75 * vertex[cut] + 0.25 * after[cut]
            xy, fixed = out, out_fixed
        return shapely.LineString(xy)

    @staticmethod
    def river_groups(
        geoms, upstream_km2, max_join_ratio, max_length_m, min_upstream_km2, touch
    ):
        """
        Groups river reaches that are connected and of similar size. Reaches meet at shared
        end points (within touch); the downstream neighbour of a reach is the touching reach
        with the largest upstream area, if larger than its own. Going upstream from the
        largest reaches, a reach joins the group of its downstream neighbour when the
        upstream area grows by at most max_join_ratio across the junction and the group stays
        within max_length_m; otherwise it starts a group. Reaches below min_upstream_km2 form
        group 0 (minor streams, handled together).
        :param geoms: array of LineString or MultiLineString reaches
        :param upstream_km2: upstream area per reach
        :param max_join_ratio: largest upstream area ratio (downstream / upstream) to join
        :param max_length_m: largest group length
        :param min_upstream_km2: smallest upstream area of a reach with its own group
        :param touch: distance at which end points meet
        :return: int array of group ids (0 = minor streams; groups numbered from 1)
        """
        geoms = np.asarray(geoms, dtype=object)
        area = np.nan_to_num(np.asarray(upstream_km2, float), nan=0.0)
        n = len(geoms)
        group = np.zeros(n, np.int64)
        if not n:
            return group
        parts, owner = shapely.get_parts(geoms, return_index=True)
        ends = np.concatenate(
            [shapely.get_point(parts, 0), shapely.get_point(parts, -1)]
        )
        end_owner = np.concatenate([owner, owner])
        i, j = shapely.STRtree(ends).query(ends, predicate="dwithin", distance=touch)
        a, b = end_owner[i], end_owner[j]
        pair = a != b
        a, b = a[pair], b[pair]
        # downstream neighbour: the touching reach with the largest upstream area, if larger than the reach
        down = np.full(n, -1, np.int64)
        best = np.full(n, -np.inf)
        larger = area[b] > area[a]
        for r, d in zip(a[larger], b[larger]):
            if area[d] > best[r]:
                best[r], down[r] = area[d], d
        length = shapely.length(geoms)
        group_len = {}
        next_id = 1
        for r in np.argsort(-area, kind="stable"):
            if area[r] < min_upstream_km2:
                continue
            d = down[r]
            g = group[d] if d >= 0 else 0
            if (
                g > 0
                and area[d] <= max_join_ratio * max(area[r], 1e-9)
                and group_len[g] + length[r] <= max_length_m
            ):
                group[r] = g
                group_len[g] += length[r]
            else:
                group[r] = next_id
                group_len[next_id] = length[r]
                next_id += 1
        return group

    @staticmethod
    def line_samples(geoms, spacing, step):
        """
        Points spread evenly along lines, about spacing apart (at least one per part), with the unit normal of the
        line at each point (left of the direction of digitising)
        :param geoms: array of LineString or MultiLineString
        :param spacing: approximate distance between samples
        :param step: distance before and after a sample of the two points that give the direction
        :return: (owner index per sample, xy (n, 2), normal (n, 2), length of line represented by each sample)
        """
        parts, owner = shapely.get_parts(
            np.asarray(geoms, dtype=object), return_index=True
        )
        keep = (shapely.get_type_id(parts) == 1) & (shapely.length(parts) > 0)
        parts, owner = parts[keep], owner[keep]
        length = shapely.length(parts)
        count = np.maximum(1, np.floor(length / spacing + 0.5)).astype(np.int64)
        idx = np.repeat(np.arange(len(parts)), count)
        k = np.arange(len(idx)) - np.repeat(np.cumsum(count) - count, count)
        weight = length[idx] / count[idx]
        dist = (k + 0.5) * weight
        xy = shapely.get_coordinates(shapely.line_interpolate_point(parts[idx], dist))
        d = shapely.get_coordinates(
            shapely.line_interpolate_point(parts[idx], dist + step)
        ) - shapely.get_coordinates(
            shapely.line_interpolate_point(parts[idx], dist - step)
        )
        norm = np.maximum(np.hypot(d[:, 0], d[:, 1]), 1e-12)
        normal = np.column_stack([-d[:, 1] / norm, d[:, 0] / norm])
        return owner[idx], xy, normal, weight

    @staticmethod
    def dominant_nearby(points, geoms, weights, tolerance):
        """
        Per point, the feature with the largest weight among those within the distance of
        the nearest feature plus tolerance (e.g. the main river at a confluence, where the
        tributary may be a little closer to a river cell)
        :param points: array of points
        :param geoms: array of geometries
        :param weights: weight per geometry (NaN counts as the smallest)
        :param tolerance: added distance, e.g. one working cell
        :return: int array of geometry indices (-1 where there is none)
        """
        tree = shapely.STRtree(geoms)
        i_pt, i_g = tree.query_nearest(points, all_matches=False)
        near = np.full(len(points), np.inf)
        near[i_pt] = shapely.distance(points[i_pt], geoms[i_g])
        reach = np.where(np.isfinite(near), near + tolerance, 0.0)
        p_idx, g_idx = tree.query(points, predicate="dwithin", distance=reach)
        w = np.nan_to_num(np.asarray(weights, float)[g_idx], nan=-np.inf)
        df = pd.DataFrame({"p": p_idx, "g": g_idx, "w": w}).sort_values(
            ["p", "w"], ascending=[True, False], kind="stable"
        )
        best = df.drop_duplicates("p")
        out = np.full(len(points), -1, np.int64)
        out[best["p"].to_numpy()] = best["g"].to_numpy()
        return out

    @staticmethod
    def interpolate_xy(line, distances):
        """
        Points at distances along a line, vectorized with numpy (one search in the
        cumulative segment lengths instead of one GEOS walk along the line per point)
        :param line: shapely LineString
        :param distances: distances from the start; clipped to [0, length]
        :return: (n, 2) float array of coordinates
        """
        xy = shapely.get_coordinates(line)
        seg = np.hypot(*np.diff(xy, axis=0).T)
        cum = np.r_[0.0, np.cumsum(seg)]
        d = np.clip(np.asarray(distances, dtype=float), 0.0, cum[-1])
        i = np.clip(np.searchsorted(cum, d, side="right") - 1, 0, len(seg) - 1)
        t = np.divide(d - cum[i], seg[i], out=np.zeros(len(d)), where=seg[i] > 0)
        return xy[i] + t[:, None] * (xy[i + 1] - xy[i])


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
    def hand_above_streams(dem, stream, burn, burn_buffer_cells, valid=None):
        """
        Height above the stream cell each cell drains to (HAND), the code path of
        protected_areas.py: breach_to_streams routing; cells that reach no stream cell take
        the nearest one
        :param dem: 2D float array, NaN = nodata (for HAND: the undefended DTM)
        :param stream: 2D bool array of stream cells
        :param burn: stream burning depth (see breach_to_streams)
        :param burn_buffer_cells: taper distance of the burning in cells
        :param valid: cells to compute (default: finite dem)
        :return: (hand float32, river_ref float32 elevation of the stream cell reached,
                  number of cells that took the nearest stream cell)
        """
        shape = dem.shape
        valid = np.isfinite(dem) if valid is None else valid
        reached = Hydro.breach_to_streams(dem, stream, burn, burn_buffer_cells)
        unreached = valid & (reached < 0)
        n_unreached = int(unreached.sum())
        if n_unreached:
            nearest = ndimage.distance_transform_edt(
                ~stream, return_distances=False, return_indices=True
            )
            reached[unreached] = (nearest[0].astype(np.int64) * shape[1] + nearest[1])[
                unreached
            ]
            del nearest
        river_ref = np.where(
            valid, dem.ravel()[np.maximum(reached, 0)].reshape(shape), np.nan
        ).astype(np.float32)
        return dem - river_ref, river_ref, n_unreached

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

    @staticmethod
    def least_cost_path(cost, start, end):
        """
        Least-cost path between two cells (Dijkstra, 8-connected). A step costs its
        length in cells (1, diagonal sqrt 2) times the mean cost of the two cells.
        :param cost: 2D float array of non-negative costs, NaN = impassable
        :param start: (row, col) of the start cell
        :param end: (row, col) of the end cell
        :return: (int64 array of flat cell indices from start to end, empty if end is
                  unreachable; total cost, inf if unreachable)
        """
        cost = np.ascontiguousarray(cost, dtype=np.float64)
        rows, cols = cost.shape
        for r, c in (start, end):
            if not (0 <= r < rows and 0 <= c < cols) or np.isnan(cost[r, c]):
                return np.empty(0, np.int64), float("inf")
        path, total = _least_cost_path(
            cost, np.int64(start[0] * cols + start[1]), np.int64(end[0] * cols + end[1])
        )
        return path, float(total)

    @staticmethod
    def _drop_turns(xy, part, test, passes):
        """
        Removes interior vertices whose turn fails a test, all of one pass at once, until
        nothing changes or passes are used up; the ends of every part stay
        :param xy: (n, 2) coordinates ordered by part
        :param part: part index per vertex
        :param test: function (angle_deg, shorter_segment_length, distance to the chord of the
                     two neighbours) -> bool array, True = remove
        :param passes: maximum number of passes
        :return: (xy, part) after removal
        """
        for _ in range(int(passes)):
            if len(xy) < 3:
                break
            inner = np.zeros(len(xy), bool)
            inner[1:-1] = (part[:-2] == part[1:-1]) & (part[2:] == part[1:-1])
            idx = np.flatnonzero(inner)
            v1, v2 = xy[idx - 1] - xy[idx], xy[idx + 1] - xy[idx]
            n1, n2 = np.hypot(v1[:, 0], v1[:, 1]), np.hypot(v2[:, 0], v2[:, 1])
            with np.errstate(invalid="ignore", divide="ignore"):
                cos = np.clip(np.sum(v1 * v2, axis=1) / (n1 * n2), -1.0, 1.0)
            angle = np.degrees(np.arccos(cos))
            chord = np.hypot(*(v2 - v1).T)
            with np.errstate(invalid="ignore", divide="ignore"):
                dev = np.abs(v1[:, 0] * v2[:, 1] - v1[:, 1] * v2[:, 0]) / chord
            drop = (n1 == 0) | (n2 == 0) | test(angle, np.minimum(n1, n2), dev)
            if not drop.any():
                break
            keep = np.ones(len(xy), bool)
            keep[idx[drop]] = False
            xy, part = xy[keep], part[keep]
        return xy, part

    @staticmethod
    def snap_to_channel(lines, dtm_path, guides, params):
        """
        Moves river lines into their channel. Vertices every densify_m are projected onto
        the nearest guide line within guide_search_m (shift capped at guide_max_shift_m),
        then moved to the lowest DTM cell within window_m around that point (by at most
        max_shift_m, or max_shift_without_guide_m where no guide is near). Sharp turns
        (below min_angle_deg) and short spikes (below spike_angle_deg, or nearly straight
        above straight_angle_deg, on segments up to spike_length_m, the straight case up to
        straight_length_factor of it and only within straight_max_deviation_m of the chord)
        are removed afterwards. Parts in which no vertex moves stay unchanged. Vertices at
        the same place move the same way, so connected lines stay connected. The DTM is
        read once in bands of band_rows rows.
        :param lines: array of LineString or MultiLineString in the DTM CRS
        :param dtm_path: DTM at native resolution
        :param guides: array of guide lines in the DTM CRS, or None
        :param params: dict with the keys named above plus angle_passes and spike_passes
        :return: (array of geometries, dict of statistics)
        """
        p = params
        lines = np.asarray(lines, dtype=object)
        valid = ~shapely.is_missing(lines) & ~shapely.is_empty(lines)
        parts, owner = shapely.get_parts(lines[valid], return_index=True)
        owner = np.flatnonzero(valid)[owner]
        is_line = shapely.get_type_id(parts) == 1
        parts, owner = parts[is_line], owner[is_line]
        xy, part = shapely.get_coordinates(
            shapely.segmentize(parts, p["densify_m"]), return_index=True
        )
        n = len(xy)
        stats = {"lines": int(valid.sum()), "vertices": n}
        if not n:
            return lines.copy(), stats
        stage = xy.copy()
        cap = np.full(n, float(p["max_shift_without_guide_m"]))
        # --- Guide lines: projection onto the nearest one within guide_search_m ---
        guided = np.zeros(n, bool)
        if guides is not None and len(guides):
            guides = np.asarray(guides, dtype=object)
            pts = shapely.points(xy)
            i_pt, i_g = shapely.STRtree(guides).query_nearest(
                pts, max_distance=p["guide_search_m"], all_matches=False
            )
            target = shapely.get_coordinates(
                shapely.line_interpolate_point(
                    guides[i_g], shapely.line_locate_point(guides[i_g], pts[i_pt])
                )
            )
            d = target - xy[i_pt]
            dist = np.hypot(d[:, 0], d[:, 1])
            scale = np.minimum(1.0, p["guide_max_shift_m"] / np.maximum(dist, 1e-12))
            stage[i_pt] = xy[i_pt] + d * scale[:, None]
            guided[i_pt] = True
            cap[i_pt] = min(p["max_shift_m"], p["guide_max_shift_m"])
        stats["guided_share"] = float(guided.mean())
        # --- DTM: lowest cell within the window, read in bands of rows ---
        ds = gdal.Open(str(dtm_path))
        band = ds.GetRasterBand(1)
        gt, nodata = ds.GetGeoTransform(), band.GetNoDataValue()
        half = max(int(round(p["window_m"] / abs(gt[1]) / 2)), 0)
        col = np.floor((stage[:, 0] - gt[0]) / gt[1]).astype(np.int64)
        row = np.floor((stage[:, 1] - gt[3]) / gt[5]).astype(np.int64)
        todo = (
            (cap > 0)
            & (row >= 0)
            & (row < ds.RasterYSize)
            & (col >= 0)
            & (col < ds.RasterXSize)
        )
        best = stage.copy()
        off = np.arange(-half, half + 1)
        size = len(off) ** 2
        chunk = max(1, int(2e7 // size))
        band_rows = int(p["band_rows"])
        for r0 in range(0, ds.RasterYSize, band_rows):
            sel = np.flatnonzero(todo & (row >= r0) & (row < r0 + band_rows))
            if not len(sel):
                continue
            y0, y1 = max(r0 - half, 0), min(r0 + band_rows + half, ds.RasterYSize)
            x0 = max(int(col[sel].min()) - half, 0)
            x1 = min(int(col[sel].max()) + half + 1, ds.RasterXSize)
            arr = band.ReadAsArray(x0, y0, x1 - x0, y1 - y0).astype(np.float64)
            if nodata is not None:
                arr[arr == nodata] = np.inf
            arr[~np.isfinite(arr)] = np.inf
            for k in range(0, len(sel), chunk):
                s = sel[k : k + chunk]
                rr = np.clip(
                    row[s, None, None] - y0 + off[None, :, None], 0, y1 - y0 - 1
                )
                cc = np.clip(
                    col[s, None, None] - x0 + off[None, None, :], 0, x1 - x0 - 1
                )
                rr, cc = np.broadcast_arrays(rr, cc)
                win = arr[rr, cc].reshape(len(s), -1)
                arg = np.argmin(win, axis=1)
                ok = np.isfinite(win[np.arange(len(s)), arg])
                br = rr.reshape(len(s), -1)[np.arange(len(s)), arg] + y0
                bc = cc.reshape(len(s), -1)[np.arange(len(s)), arg] + x0
                bx = gt[0] + (bc + 0.5) * gt[1]
                by = gt[3] + (br + 0.5) * gt[5]
                best[s[ok], 0], best[s[ok], 1] = bx[ok], by[ok]
        ds = None
        d = best - stage
        dist = np.hypot(d[:, 0], d[:, 1])
        scale = np.minimum(1.0, cap / np.maximum(dist, 1e-12))
        moved = stage + d * scale[:, None]
        shift = np.hypot(*(moved - xy).T)
        stats.update(
            {
                "shift_mean_m": float(shift.mean()),
                "shift_p95_m": float(np.percentile(shift, 95)),
            }
        )
        # --- Clean-up: sharp turns, then short spikes and nearly straight short steps ---
        xy2, part2 = Hydro._drop_turns(
            moved,
            part,
            lambda a, s, e: a < p["min_angle_deg"],
            p["angle_passes"],
        )
        xy2, part2 = Hydro._drop_turns(
            xy2,
            part2,
            lambda a, s, e: ((a < p["spike_angle_deg"]) & (s <= p["spike_length_m"]))
            | (
                (a >= p["straight_angle_deg"])
                & (s <= p["spike_length_m"] * p["straight_length_factor"])
                & (e <= p["straight_max_deviation_m"])
            ),
            p["spike_passes"],
        )
        stats["vertices_out"] = len(xy2)
        # every part keeps its two ends, so the new parts align with the old ones
        new_parts = shapely.linestrings(xy2, indices=part2)
        # a part collapsed into one cell, or with no vertex moved, keeps its original geometry
        still = np.bincount(part, weights=shift > 0, minlength=len(parts)) == 0
        bad = ~(shapely.length(new_parts) > 0) | still
        new_parts[bad] = parts[bad]
        stats["parts_unchanged"] = int(bad.sum())
        out = lines.copy()
        multi = shapely.get_type_id(lines) == 5
        order = np.argsort(owner, kind="stable")
        first = np.flatnonzero(np.r_[True, owner[order][1:] != owner[order][:-1]])
        for i, sub in zip(owner[order][first], np.split(new_parts[order], first[1:])):
            out[i] = shapely.multilinestrings(sub) if multi[i] else sub[0]
        return out, stats

    @staticmethod
    def crest_profiles(xy, normal, dtm_path, half_width, block=2048):
        """
        DTM profiles across lines: native DTM values every native pixel along the normal of each point, within
        half_width on both sides
        :param xy: (n, 2) points
        :param normal: (n, 2) unit normals
        :param dtm_path: DTM at native resolution
        :param half_width: profile half-width
        :param block: block size of the DTM reads
        :return: (offsets along the normal, float64 values (n, offsets), NaN = nodata or outside the DTM)
        """
        ds = gdal.Open(str(dtm_path))
        step = abs(ds.GetGeoTransform()[1])
        ds = None
        k = int(np.ceil(half_width / step - 1e-9))
        t = np.arange(-k, k + 1) * step
        px = xy[:, 0, None] + t[None, :] * normal[:, 0, None]
        py = xy[:, 1, None] + t[None, :] * normal[:, 1, None]
        z = Raster.sample_points(dtm_path, px.ravel(), py.ravel(), block)
        return t, z.reshape(len(xy), len(t))

    @staticmethod
    def _block_order(xy, dtm_path, block=2048):
        # order of points by DTM block, so that chunks of points read few blocks
        ds = gdal.Open(str(dtm_path))
        gt, n_cols = ds.GetGeoTransform(), ds.RasterXSize
        ds = None
        key = (np.floor((xy[:, 1] - gt[3]) / gt[5]) // block) * (
            n_cols // block + 1
        ) + np.floor((xy[:, 0] - gt[0]) / gt[1]) // block
        return np.argsort(key, kind="stable")

    @staticmethod
    def snap_to_crest(
        lines, dtm_path, params, fix_ends=False, touch=0.01, chunk=100000
    ):
        """
        Moves levee lines onto the crest seen in the DTM. Every vertex (densified to densify_m) looks along its
        normal within max_shift_m for the highest DTM cell; of the cells within level_tolerance_m of the highest,
        the nearest to the line is taken, so a flat crest does not move the line. It is a ridge when it stands at
        least min_height_m above the lowest cell within side_m on both sides. Vertices without a ridge take the
        offset interpolated from their neighbours; offsets are smoothed by a running median of smooth_vertices.
        Ends shared by several lines move together (connected lines stay connected); with fix_ends every end
        stays. Lines in which nothing moves keep their geometry.
        :param lines: array of LineString in the DTM CRS
        :param dtm_path: DTM at native resolution
        :param params: dict densify_m, max_shift_m, side_m, min_height_m, level_tolerance_m, smooth_vertices
        :param fix_ends: keep the end points of every line
        :param touch: distance within which end points are one shared end
        :param chunk: vertices per DTM sampling pass
        :return: (array of geometries, DataFrame per line with ridge_share, crest_height_m (median height of the
                  ridge above the higher side) and shift_median_m, dict of statistics)
        """
        p = params
        lines = np.asarray(lines, dtype=object)
        n_lines = len(lines)
        df = pd.DataFrame(
            {
                "ridge_share": np.zeros(n_lines),
                "crest_height_m": np.full(n_lines, np.nan),
                "shift_median_m": np.zeros(n_lines),
            }
        )
        stats = {"lines": n_lines, "vertices": 0}
        if not n_lines:
            return lines.copy(), df, stats
        xy, part = shapely.get_coordinates(
            shapely.segmentize(lines, p["densify_m"]), return_index=True
        )
        n = len(xy)
        stats["vertices"] = n
        first = np.r_[True, part[1:] != part[:-1]]
        last = np.r_[part[1:] != part[:-1], True]
        prev_i, next_i = np.arange(n) - 1, np.arange(n) + 1
        prev_i[first], next_i[last] = np.flatnonzero(first), np.flatnonzero(last)
        d = xy[next_i] - xy[prev_i]
        norm = np.maximum(np.hypot(d[:, 0], d[:, 1]), 1e-12)
        normal = np.column_stack([-d[:, 1] / norm, d[:, 0] / norm])
        # --- profiles across the line: highest cell near the line, lowest cells on both sides of it ---
        offset, height, ridge = np.zeros(n), np.full(n, np.nan), np.zeros(n, bool)
        order = Hydro._block_order(xy, dtm_path)
        for s in range(0, n, chunk):
            sel = order[s : s + chunk]
            t, z = Hydro.crest_profiles(
                xy[sel], normal[sel], dtm_path, p["max_shift_m"] + p["side_m"]
            )
            z_near = np.where(
                np.isnan(z) | (np.abs(t) > p["max_shift_m"] + 1e-9)[None, :], -np.inf, z
            )
            top = z_near.max(axis=1)
            dist = np.where(
                z_near >= top[:, None] - p["level_tolerance_m"],
                np.abs(t)[None, :],
                np.inf,
            )
            k = np.argmin(dist, axis=1)
            off = t[k]
            rel = t[None, :] - off[:, None]
            z_low = np.where(np.isnan(z), np.inf, z)
            left = np.where((rel < 0) & (rel >= -p["side_m"]), z_low, np.inf).min(1)
            right = np.where((rel > 0) & (rel <= p["side_m"]), z_low, np.inf).min(1)
            with np.errstate(invalid="ignore"):
                h = z_near[np.arange(len(sel)), k] - np.maximum(left, right)
            ok = np.isfinite(h) & (h >= p["min_height_m"])
            offset[sel], ridge[sel] = off, ok
            height[sel] = np.where(np.isfinite(h), h, np.nan)
        # --- offsets per line: interpolated over vertices without a ridge, running median ---
        shift = np.zeros(n)
        starts = np.flatnonzero(first)
        size = max(1, int(p["smooth_vertices"]))
        for s0, e0 in zip(starts, np.r_[starts[1:], n]):
            good = ridge[s0:e0]
            if not good.any():
                continue
            idx = np.arange(e0 - s0)
            o = np.interp(idx, idx[good], offset[s0:e0][good])
            if size > 1 and len(o) > 2:
                o = ndimage.median_filter(o, size=size, mode="nearest")
            if fix_ends:
                o[0] = o[-1] = 0.0
            shift[s0:e0] = o
        moved = xy + shift[:, None] * normal
        if not fix_ends:
            # ends at the same place move to the mean of their positions
            end_idx = np.flatnonzero(first | last)
            grid = np.round(xy[end_idx] / max(touch, 1e-6)).astype(np.int64)
            inv = np.unique(grid, axis=0, return_inverse=True)[1].ravel()
            cnt = np.bincount(inv)
            shared = cnt[inv] > 1
            for axis in (0, 1):
                mean = np.bincount(inv, moved[end_idx, axis]) / cnt
                moved[end_idx[shared], axis] = mean[inv[shared]]
        moved_dist = np.hypot(*(moved - xy).T)
        out = lines.copy()
        shapely.linestrings(moved, indices=part, out=out)
        still = np.bincount(part, weights=moved_dist > 1e-9, minlength=n_lines) == 0
        bad = still | ~(shapely.length(out) > 0)
        out[bad] = lines[bad]
        # --- statistics per line ---
        cnt = np.bincount(part, minlength=n_lines)
        df["ridge_share"] = np.bincount(part, weights=ridge, minlength=n_lines) / (
            np.maximum(cnt, 1)
        )
        df_v = pd.DataFrame(
            {"part": part, "h": np.where(ridge, height, np.nan), "s": moved_dist}
        ).groupby("part")
        df["crest_height_m"] = df_v["h"].median().reindex(range(n_lines)).to_numpy()
        df["shift_median_m"] = (
            df_v["s"].median().reindex(range(n_lines)).fillna(0.0).to_numpy()
        )
        stats.update(
            {
                "ridge_share": float(ridge.mean()),
                "shift_mean_m": float(moved_dist.mean()),
                "shift_p95_m": float(np.percentile(moved_dist, 95)),
                "lines_moved": int((~bad).sum()),
                "lines_without_ridge": int((df["ridge_share"] == 0).sum()),
            }
        )
        return out, df, stats

    @staticmethod
    def bridge_profile(segments, dtm_path, half_width, side_m=0.0, chunk=2000):
        """
        Crest the terrain holds in straight gaps, every native pixel along the segment: the highest DTM cell across
        it (within half_width on both sides) and how far that crest stands above the lowest cell within side_m
        beyond it on both sides (the ridge height: a levee body, not merely high ground)
        :param segments: array of LineStrings from one end of the gap to the other
        :param dtm_path: DTM at native resolution
        :param half_width: half-width of the band across the segment searched for the crest
        :param side_m: width of the terrain looked at beyond the crest on both sides (0: no ridge height)
        :param chunk: segments per DTM sampling pass
        :return: (lower of the two end values, lowest value along the segment, lowest ridge height along the
                  segment); -inf where the DTM is missing
        """
        segments = np.asarray(segments, dtype=object)
        n = len(segments)
        end_low, along_min = np.full(n, -np.inf), np.full(n, -np.inf)
        ridge_min = np.full(n, -np.inf)
        if not n:
            return end_low, along_min, ridge_min
        ds = gdal.Open(str(dtm_path))
        step = abs(ds.GetGeoTransform()[1])
        ds = None
        a = shapely.get_coordinates(shapely.get_point(segments, 0))
        d = shapely.get_coordinates(shapely.get_point(segments, -1)) - a
        norm = np.maximum(np.hypot(d[:, 0], d[:, 1]), 1e-12)
        normal = np.column_stack([-d[:, 1] / norm, d[:, 0] / norm])
        order = Hydro._block_order(a, dtm_path)
        for s in range(0, n, chunk):
            sel = order[s : s + chunk]
            xy, idx = shapely.get_coordinates(
                shapely.segmentize(segments[sel], step), return_index=True
            )
            t, z = Hydro.crest_profiles(
                xy, normal[sel][idx], dtm_path, half_width + side_m
            )
            z_band = np.where(
                np.isnan(z) | (np.abs(t) > half_width + 1e-9)[None, :], -np.inf, z
            )
            top = z_band.max(axis=1)
            k = np.argmax(z_band, axis=1)
            # ridge height: crest above the higher of the lowest cells within side_m on both sides of it
            rel = t[None, :] - t[k][:, None]
            z_low = np.where(np.isnan(z), np.inf, z)
            left = np.where((rel < 0) & (rel >= -side_m), z_low, np.inf).min(axis=1)
            right = np.where((rel > 0) & (rel <= side_m), z_low, np.inf).min(axis=1)
            with np.errstate(invalid="ignore"):
                ridge = np.where(
                    np.isfinite(left) & np.isfinite(right) & np.isfinite(top),
                    top - np.maximum(left, right),
                    -np.inf,
                )
            first = np.r_[True, idx[1:] != idx[:-1]]
            last = np.r_[idx[1:] != idx[:-1], True]
            end_low[sel[idx[first]]] = np.minimum(top[first], top[last])
            along_min[sel] = (
                pd.Series(top).groupby(idx).min().reindex(range(len(sel))).to_numpy()
            )
            ridge_min[sel] = (
                pd.Series(ridge).groupby(idx).min().reindex(range(len(sel))).to_numpy()
            )
        return end_low, along_min, ridge_min
        ds = gdal.Open(str(dtm_path))
        step = abs(ds.GetGeoTransform()[1])
        ds = None
        a = shapely.get_coordinates(shapely.get_point(segments, 0))
        d = shapely.get_coordinates(shapely.get_point(segments, -1)) - a
        norm = np.maximum(np.hypot(d[:, 0], d[:, 1]), 1e-12)
        normal = np.column_stack([-d[:, 1] / norm, d[:, 0] / norm])
        order = Hydro._block_order(a, dtm_path)
        for s in range(0, n, chunk):
            sel = order[s : s + chunk]
            xy, idx = shapely.get_coordinates(
                shapely.segmentize(segments[sel], step), return_index=True
            )
            _, z = Hydro.crest_profiles(xy, normal[sel][idx], dtm_path, half_width)
            top = np.where(np.isnan(z), -np.inf, z).max(axis=1)
            first = np.r_[True, idx[1:] != idx[:-1]]
            last = np.r_[idx[1:] != idx[:-1], True]
            end_low[sel[idx[first]]] = np.minimum(top[first], top[last])
            along_min[sel] = (
                pd.Series(top).groupby(idx).min().reindex(range(len(sel))).to_numpy()
            )
        return end_low, along_min

    @staticmethod
    def embankment_stations(t, z, params):
        """
        Cross-sections of a linear structure (road, railway) evaluated per station: crest, terrain on both sides,
        heights, footprint and class. Classes: 1 embankment (terrain drops on both sides by at least min_height_m),
        2 one-sided (only one side drops), 3 cutting (the line lies min_height_m below the terrain on both sides),
        0 at grade, -1 no data.
        :param t: offsets along the normal (negative = left of the line direction)
        :param z: (n, len(t)) DTM values of the profiles, NaN = nodata
        :param params: dict crest_search_m, side_near_m, side_far_m, min_height_m, footprint_fraction
        :return: dict of arrays crest_z, line_z, terrain_left, terrain_right, height_left, height_right,
                 footprint_m, cls
        """
        p = params
        a = np.abs(t)
        with np.errstate(invalid="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            crest = np.nanmax(np.where(a <= p["crest_search_m"], z, np.nan), axis=1)
            line_z = np.nanmedian(
                np.where(
                    a <= max(p["crest_search_m"] / 4, abs(t[1] - t[0])), z, np.nan
                ),
                axis=1,
            )
            band_l = (t <= -p["side_near_m"]) & (t >= -p["side_far_m"])
            band_r = (t >= p["side_near_m"]) & (t <= p["side_far_m"])
            terrain_l = np.nanmedian(np.where(band_l, z, np.nan), axis=1)
            terrain_r = np.nanmedian(np.where(band_r, z, np.nan), axis=1)
            height_l, height_r = crest - terrain_l, crest - terrain_r
            # footprint: from the crest outwards to the first cell within footprint_fraction of the height above
            # the terrain of that side
            foot = np.zeros(len(z))
            for side, terrain, height in (
                (-1, terrain_l, height_l),
                (1, terrain_r, height_r),
            ):
                limit = terrain + p["footprint_fraction"] * np.maximum(height, 0)
                reach = np.where(
                    (np.sign(t) == side) & (a <= p["side_near_m"]), a, np.inf
                )[None, :] + np.where(z <= limit[:, None], 0.0, np.inf)
                foot += np.where(
                    np.isfinite(reach.min(axis=1)), reach.min(axis=1), p["side_near_m"]
                )
        h = p["min_height_m"]
        up_l, up_r = height_l >= h, height_r >= h
        cls = np.zeros(len(z), np.int8)
        cls[up_l ^ up_r] = 2
        cls[up_l & up_r] = 1
        cut = (line_z <= terrain_l - h) & (line_z <= terrain_r - h)
        cls[cut & ~(up_l | up_r)] = 3
        cls[~(np.isfinite(crest) & np.isfinite(terrain_l) & np.isfinite(terrain_r))] = (
            -1
        )
        return {
            "crest_z": crest,
            "line_z": line_z,
            "terrain_left": terrain_l,
            "terrain_right": terrain_r,
            "height_left": height_l,
            "height_right": height_r,
            "footprint_m": foot,
            "cls": cls,
        }

    @staticmethod
    def first_event_on_path(order, seeds, event, starts):
        """
        Follows the paths of a priority flood back from start cells towards their seed and reports the first cell
        with an event code on the way. The predecessor of a cell is its 4-neighbour taken first from the queue
        (the cell that reached it), recovered from the order array of priority_flood.
        :param order: order array of priority_flood (-1 = not reached)
        :param seeds: seed cells of the flood
        :param event: int8 array, 0 = no event, > 0 = event code
        :param starts: flat indices of the start cells
        :return: (int8 code per start: the event code, 0 = start not reached, -1 = seed reached without an event;
                  int64 flat index of the event cell or of the seed, -1 if none)
        """
        return _first_event(
            np.ascontiguousarray(order, dtype=np.int32),
            np.ascontiguousarray(seeds, dtype=np.bool_),
            np.ascontiguousarray(event, dtype=np.int8),
            np.asarray(starts, dtype=np.int64),
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


@njit(cache=True)
def _first_event(order, seeds, event, starts):
    # Walks back along the flood path: the predecessor is the 4-neighbour with the smallest order.
    rows, cols = order.shape
    o = order.ravel()
    s = seeds.ravel()
    ev = event.ravel()
    m = len(starts)
    code = np.zeros(m, np.int8)
    cell = np.full(m, -1, np.int64)
    for k in range(m):
        c = starts[k]
        if o[c] < 0:
            continue
        while True:
            if ev[c] > 0:
                code[k] = ev[c]
                cell[k] = c
                break
            if s[c]:
                code[k] = -1
                cell[k] = c
                break
            ci = c // cols
            cj = c - ci * cols
            best = -1
            best_o = o[c]
            if ci > 0 and 0 <= o[c - cols] < best_o:
                best, best_o = c - cols, o[c - cols]
            if ci < rows - 1 and 0 <= o[c + cols] < best_o:
                best, best_o = c + cols, o[c + cols]
            if cj > 0 and 0 <= o[c - 1] < best_o:
                best, best_o = c - 1, o[c - 1]
            if cj < cols - 1 and 0 <= o[c + 1] < best_o:
                best, best_o = c + 1, o[c + 1]
            if best < 0:
                break
            c = best
    return code, cell


@njit(cache=True)
def _least_cost_path(cost, start, end):
    # Dijkstra on an 8-connected grid; a step costs its length (1 or sqrt 2 cells) times the mean cost of its cells.
    rows, cols = cost.shape
    flat = cost.ravel()
    n = rows * cols
    dist = np.full(n, np.inf)
    prev = np.full(n, -1, np.int64)
    done = np.zeros(n, np.bool_)
    keys = np.empty(4096, np.float64)
    vals = np.empty(4096, np.int64)
    size = 0
    dist[start] = 0.0
    keys, vals, size = _heap_push(keys, vals, size, 0.0, np.int64(start))
    d_row = np.array([-1, 1, 0, 0, -1, -1, 1, 1])
    d_col = np.array([0, 0, -1, 1, -1, 1, -1, 1])
    diag = sqrt(2.0)
    step = np.array([1.0, 1.0, 1.0, 1.0, diag, diag, diag, diag])
    while size > 0:
        d, c, size = _heap_pop(keys, vals, size)
        if done[c]:
            continue
        done[c] = True
        if c == end:
            break
        ci = c // cols
        cj = c - ci * cols
        for k in range(8):
            ni = ci + d_row[k]
            nj = cj + d_col[k]
            if ni < 0 or ni >= rows or nj < 0 or nj >= cols:
                continue
            nb = ni * cols + nj
            if done[nb] or np.isnan(flat[nb]):
                continue
            nd = d + step[k] * 0.5 * (flat[c] + flat[nb])
            if nd < dist[nb]:
                dist[nb] = nd
                prev[nb] = c
                keys, vals, size = _heap_push(keys, vals, size, nd, np.int64(nb))
    if not done[end]:
        return np.empty(0, np.int64), np.inf
    count = 0
    c = end
    while c >= 0:
        count += 1
        c = prev[c]
    path = np.empty(count, np.int64)
    c = end
    for i in range(count - 1, -1, -1):
        path[i] = c
        c = prev[c]
    return path, dist[end]


class Discharge:
    """
    River discharge of gridded hydrological products (GloFAS, EFAS), used only to turn
    observed outcomes into return periods (Stage F)
    """

    @staticmethod
    def band_times(ds):
        """
        Valid time of every band of a GRIB or NetCDF dataset
        :param ds: gdal.Dataset
        :return: int64 array of seconds since 1970-01-01 UTC
        """
        metadata = ds.GetMetadata()
        list_t = []
        for b in range(1, ds.RasterCount + 1):
            meta = ds.GetRasterBand(b).GetMetadata()
            if "GRIB_VALID_TIME" in meta:
                list_t.append(int(meta["GRIB_VALID_TIME"].split()[0]))
                continue
            # NetCDF: a valid_time dimension, else the only time dimension; a step dimension or several time
            # dimensions leave the valid time open and stop the run instead of guessing
            dims = [k for k in meta if k.startswith("NETCDF_DIM_")]
            times = [k for k in dims if "time" in k.lower()]
            if "NETCDF_DIM_valid_time" in times:
                dim = "NETCDF_DIM_valid_time"
            elif len(times) == 1 and not any(k.lower().endswith("step") for k in dims):
                dim = times[0]
            else:
                raise ValueError(
                    f"Band {b}: no unambiguous valid time among the dimensions {dims}."
                )
            units = metadata.get(f"{dim[len('NETCDF_DIM_'):]}#units", "")
            unit, _, reference = units.partition(" since ")
            step = {
                "days": 86400,
                "day": 86400,
                "hours": 3600,
                "hour": 3600,
                "minutes": 60,
                "minute": 60,
                "seconds": 1,
                "second": 1,
            }.get(unit.strip())
            if step is None or not reference:
                raise ValueError(f"Unsupported time units {units!r}.")
            t0 = pd.Timestamp(reference.strip())
            t0 = t0.tz_localize("UTC") if t0.tzinfo is None else t0.tz_convert("UTC")
            list_t.append(int(round(t0.timestamp() + float(meta[dim]) * step)))
        return np.array(list_t, np.int64)

    @staticmethod
    def period_days(valid_seconds, marks):
        """
        Day of the averaging period of each value. With marks = "end" the valid time
        closes the period (a 24 h value stamped 00 UTC belongs to the day before, a 6 h
        value stamped 06 UTC to its own day); with "start" it opens it.
        :param valid_seconds: valid times in seconds since 1970-01-01 UTC
        :param marks: "end" or "start"
        :return: datetime64[D] array
        """
        t = np.asarray(valid_seconds, np.int64) - (1 if marks == "end" else 0)
        return (t // 86400).astype("datetime64[D]")

    @staticmethod
    def aggregate_daily(days, values, how, steps_per_day):
        """
        Daily values from the steps of each day: mean (daily_mean) or maximum
        (daily_max); a day with fewer than steps_per_day finite steps is NaN
        :param days: day of each step (datetime64[D])
        :param values: (n_steps, n_cells) array
        :param how: daily_mean or daily_max
        :param steps_per_day: expected steps per day (24 / step hours)
        :return: (sorted unique days, (n_days, n_cells) float32)
        """
        order = np.argsort(days, kind="stable")
        days, values = np.asarray(days)[order], np.asarray(values, np.float64)[order]
        uniq, first = np.unique(days, return_index=True)
        finite = np.isfinite(values)
        n = np.add.reduceat(finite.astype(np.int64), first, axis=0)
        if how == "daily_max":
            agg = np.maximum.reduceat(np.where(finite, values, -np.inf), first, axis=0)
        else:
            agg = np.add.reduceat(
                np.where(finite, values, 0.0), first, axis=0
            ) / np.maximum(n, 1)
        agg[n < steps_per_day] = np.nan
        return uniq, agg.astype(np.float32)

    @staticmethod
    def annual_maxima(days, values, start_month):
        """
        Maxima per hydrological year starting on the first day of start_month (labelled
        by the calendar year it starts in). A year counts only when every day is present;
        a cell's year also needs every value finite.
        :param days: sorted unique days (datetime64[D])
        :param values: (n_days, n_cells) daily values
        :param start_month: 1-12
        :return: (years, (n_years, n_cells) maxima with NaN for incomplete years)
        """
        index = pd.DatetimeIndex(np.asarray(days, "datetime64[D]"))
        hyear = np.asarray(index.year - (index.month < start_month))
        years = np.unique(hyear)
        maxima = np.full((len(years), values.shape[1]), np.nan, np.float32)
        for i, y in enumerate(years):
            sel = hyear == y
            begin = pd.Timestamp(year=int(y), month=start_month, day=1)
            if sel.sum() < (begin + pd.DateOffset(years=1) - begin).days:
                continue
            block = values[sel]
            complete = np.isfinite(block).all(axis=0)
            maxima[i, complete] = block[:, complete].max(axis=0)
        return years, maxima

    @staticmethod
    def gumbel_fits(maxima, min_years):
        """
        Gumbel fit by L-moments per cell (column) of its complete years, and the RMSE
        between the sorted maxima and the fitted quantiles at Gringorten plotting
        positions (i - 0.44) / (n + 0.12)
        :param maxima: (n_years, n_cells), NaN for incomplete years
        :param min_years: fewer complete years give NaN
        :return: (n_years per cell, xi, alpha, fit_rmse)
        """
        n_cells = maxima.shape[1]
        out = [np.zeros(n_cells, np.int64)] + [
            np.full(n_cells, np.nan) for _ in range(3)
        ]
        for j in range(n_cells):
            x = np.sort(maxima[:, j][np.isfinite(maxima[:, j])])
            out[0][j] = len(x)
            if len(x) < max(min_years, 2):
                continue
            xi, alpha = Stats.gumbel_lmoments(x)
            p = (np.arange(1, len(x) + 1) - 0.44) / (len(x) + 0.12)
            fitted = xi - alpha * np.log(-np.log(p))
            out[1][j], out[2][j] = xi, alpha
            out[3][j] = float(np.sqrt(np.mean((x - fitted) ** 2)))
        return tuple(out)

    @staticmethod
    def detect_events(days, values, rp, min_rp, days_before, days_after):
        """
        Flood events of one daily series: runs of days with a return period of at least min_rp give one peak
        each (the day of the highest value); the window of a peak runs from days_before before it to days_after
        after it; peaks whose windows overlap are merged into one event with the highest peak
        :param days: sorted days (datetime64[D])
        :param values: daily discharge, NaN = missing
        :param rp: return period of every day
        :param min_rp: smallest return period of a peak
        :param days_before: window start before the peak, in days
        :param days_after: window end after the peak, in days
        :return: list of dicts peak_date, peak_value, peak_rp, start, end (datetime64[D]), n_peaks, in time order
        """
        values = np.asarray(values, dtype=float)
        above = np.nan_to_num(np.asarray(rp, dtype=float), nan=0.0) >= min_rp
        list_event = []
        for s, e, on in Stats.runs(above.tolist(), 1):
            if not on:
                continue
            k = s + int(np.nanargmax(values[s:e]))
            peak = {
                "peak_date": days[k],
                "peak_value": float(values[k]),
                "peak_rp": float(rp[k]),
                "start": days[k] - np.timedelta64(int(days_before), "D"),
                "end": days[k] + np.timedelta64(int(days_after), "D"),
                "n_peaks": 1,
            }
            last = list_event[-1] if list_event else None
            if last is not None and peak["start"] <= last["end"]:
                if peak["peak_value"] > last["peak_value"]:
                    for key in ("peak_date", "peak_value", "peak_rp"):
                        last[key] = peak[key]
                last["end"] = max(last["end"], peak["end"])
                last["n_peaks"] += 1
            else:
                list_event.append(peak)
        return list_event

    @staticmethod
    def match(
        lon,
        lat,
        river_area_km2,
        cell_lon,
        cell_lat,
        cell_area_km2,
        radius_km,
        ratio_min,
        ratio_max,
        same_river=None,
    ):
        """
        Discharge cell of a river location: candidates within radius_km (great circle)
        whose upstream area ratio (cell / river) lies in [ratio_min, ratio_max]. Cells
        flagged in same_river (near the river feature of the location) are preferred;
        then the smallest |log ratio|, then the smallest distance.
        :param same_river: bool per cell, or None
        :return: dict index (-1 = none), distance_km, ratio, same_river, quality
                 (good: same river; weak: other candidates only; none)
        """
        cell_lon, cell_lat = np.asarray(cell_lon, float), np.asarray(cell_lat, float)
        p1, p2 = np.radians(lat), np.radians(cell_lat)
        a = (
            np.sin((p2 - p1) / 2) ** 2
            + np.cos(p1) * np.cos(p2) * np.sin(np.radians(cell_lon - lon) / 2) ** 2
        )
        dist = 2 * 6371.0088 * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0)))
        ratio = (
            np.asarray(cell_area_km2, float) / river_area_km2
            if river_area_km2 > 0
            else np.full(len(dist), np.nan)
        )
        ok = (dist <= radius_km) & (ratio >= ratio_min) & (ratio <= ratio_max)
        same = (
            np.zeros(len(dist), bool)
            if same_river is None
            else np.asarray(same_river, bool)
        )
        result = {
            "index": -1,
            "distance_km": np.nan,
            "ratio": np.nan,
            "same_river": False,
            "quality": "none",
        }
        if not ok.any():
            return result
        pool = np.flatnonzero(ok & same) if (ok & same).any() else np.flatnonzero(ok)
        best = pool[np.lexsort((dist[pool], np.abs(np.log(ratio[pool]))))[0]]
        result.update(
            {
                "index": int(best),
                "distance_km": float(dist[best]),
                "ratio": float(ratio[best]),
                "same_river": bool(same[best]),
                "quality": "good" if same[best] else "weak",
            }
        )
        return result


class Evidence:
    """
    Tests of protected regions against flood scenes (the logic of gfm_evidence.py,
    shared with Stage E): a region counts as held when it is dry while its comparison
    zone, the undefended floodplain around it at the same height, is wet.
    """

    @staticmethod
    def scene_masks(scene, bounds, resolution, wkt, hand, params, warp_memory_mb):
        """
        Valid and wet cells of one scene on a grid
        :param scene: dict with flood_path, exclusion_path, reference_water_path and, optionally, advisory_path
        :param bounds: grid bounds in wkt
        :param resolution: cell size
        :param wkt: grid CRS
        :param hand: HAND on the grid (nodata cells are not valid)
        :param params: dict with nodata_values, exclusion_values, reference_water_values,
                       flood_values; advisory_values (optional, non-empty = cells with these
                       advisory flags are not valid, where the scene has an advisory layer)
        :param warp_memory_mb: warp memory in MB
        :return: (valid, wet, reference water) bool arrays
        """
        nodata = params["nodata_values"][0] if params["nodata_values"] else None
        layers = {}
        for role in ["flood", "exclusion", "reference_water"]:
            layers[role], _ = Raster.read_window(
                scene[f"{role}_path"],
                bounds,
                resolution,
                "near",
                warp_memory_mb,
                dst_wkt=wkt,
                src_nodata=nodata,
            )
        flood = layers["flood"]
        reference = np.isin(layers["reference_water"], params["reference_water_values"])
        valid = (
            np.isfinite(flood)
            & ~np.isin(flood, params["nodata_values"])
            & ~np.isin(layers["exclusion"], params["exclusion_values"])
            & ~reference
            & np.isfinite(hand)
        )
        if params.get("advisory_values") and scene.get("advisory_path"):
            advisory, _ = Raster.read_window(
                scene["advisory_path"],
                bounds,
                resolution,
                "near",
                warp_memory_mb,
                dst_wkt=wkt,
                src_nodata=nodata,
            )
            valid &= ~np.isin(advisory, params["advisory_values"])
        return valid, valid & np.isin(flood, params["flood_values"]), reference

    @staticmethod
    def comparison_zones(labels, hand, regions, reach, max_level_m, outside=None):
        """
        Comparison zone per region: cells outside all regions (labels == 0, and inside
        outside where given) within reach cells of the region, with HAND not above the
        median HAND inside it and not above max_level_m
        :param labels: region labels, 0 = none
        :param hand: HAND on the same grid
        :param regions: dict label -> region id
        :param reach: comparison distance in cells
        :param max_level_m: highest HAND of the zone
        :param outside: optional bool mask of cells allowed in zones
        :return: dict label -> (region id, (row slice, column slice), inside mask of the
                 slice, zone mask of the slice)
        """
        free = labels == 0 if outside is None else (labels == 0) & outside
        zones = {}
        for label, region in regions.items():
            inside = labels == label
            if not inside.any():
                continue
            rows_any = np.flatnonzero(inside.any(axis=1))
            cols_any = np.flatnonzero(inside.any(axis=0))
            part = (
                slice(
                    max(rows_any[0] - reach, 0),
                    min(rows_any[-1] + 1 + reach, labels.shape[0]),
                ),
                slice(
                    max(cols_any[0] - reach, 0),
                    min(cols_any[-1] + 1 + reach, labels.shape[1]),
                ),
            )
            near = Hydro.distance_cells(inside[part], reach + 1) <= reach
            level = np.nanmedian(hand[part][inside[part]])
            zone = (
                near & free[part] & (hand[part] <= level) & (hand[part] <= max_level_m)
            )
            zones[label] = (region, part, inside[part], zone)
        return zones

    @staticmethod
    def outcome(inside_valid, inside_wet, comp_valid, comp_wet, params):
        """
        Outcome of one region in one scene from cell counts
        :param params: dict with min_valid_pixels, f_wet_comparison, f_dry_inside,
                       f_wet_inside
        :return: held, flooded, inconclusive or no_test
        """
        inside_frac = inside_wet / inside_valid if inside_valid else np.nan
        comp_frac = comp_wet / comp_valid if comp_valid else np.nan
        # A dry region counts only if the undefended floodplain around it is wet in the same scene.
        if (
            inside_valid < params["min_valid_pixels"]
            or comp_valid < params["min_valid_pixels"]
            or not comp_frac >= params["f_wet_comparison"]
        ):
            return "no_test"
        if inside_frac <= params["f_dry_inside"]:
            return "held"
        if inside_frac >= params["f_wet_inside"]:
            return "flooded"
        return "inconclusive"

    @staticmethod
    def test_regions(zones, valid, wet, params):
        """
        Tests every region of comparison_zones in one scene
        :return: list of dicts: region, inside_valid, inside_wet, inside_frac,
                 comp_valid, comp_wet, comp_frac, outcome
        """
        out = []
        for label, (region, part, inside, zone) in zones.items():
            counts = [
                int((inside & valid[part]).sum()),
                int((inside & wet[part]).sum()),
                int((zone & valid[part]).sum()),
                int((zone & wet[part]).sum()),
            ]
            out.append(
                {
                    "region": region,
                    "inside_valid": counts[0],
                    "inside_wet": counts[1],
                    "inside_frac": counts[1] / counts[0] if counts[0] else np.nan,
                    "comp_valid": counts[2],
                    "comp_wet": counts[3],
                    "comp_frac": counts[3] / counts[2] if counts[2] else np.nan,
                    "outcome": Evidence.outcome(*counts, params),
                }
            )
        return out


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
        "country.river_guides": "opt_paths",
        "river_groups.max_join_ratio": "float",
        "river_groups.max_length_m": "float",
        "river_groups.min_upstream_area_km2": "float",
        "river_groups.crop_margin_m": "float",
        "levee_snap.densify_m": "float",
        "levee_snap.max_shift_m": "float",
        "levee_snap.side_m": "float",
        "levee_snap.min_height_m": "float",
        "levee_snap.level_tolerance_m": "float",
        "levee_snap.smooth_vertices": "int",
        "levee_snap.min_ridge_share": "float",
        "levee_gaps.max_gap_m": "float",
        "levee_gaps.crest_tolerance_m": "float",
        "leaks.sample_spacing_m": "float",
        "country.transport_lines": "opt_paths",
        "transport.station_spacing_m": "float",
        "transport.half_width_m": "float",
        "transport.crest_search_m": "float",
        "transport.side_near_m": "float",
        "transport.side_far_m": "float",
        "transport.min_height_m": "float",
        "transport.footprint_fraction": "float",
        "transport.min_segment_m": "float",
        "transport.break_window_m": "float",
        "transport.break_depth_m": "float",
        "transport.max_opening_m": "float",
        "transport.culvert_window_m": "float",
        "transport.culvert_min_depth_m": "float",
        "transport.stream_distance_m": "float",
        "transport.cross_sections_per_class": "int",
        "river_snap.densify_m": "float",
        "river_snap.guide_search_m": "float",
        "river_snap.guide_max_shift_m": "float",
        "river_snap.window_m": "float",
        "river_snap.max_shift_m": "float",
        "river_snap.max_shift_without_guide_m": "float",
        "river_snap.min_angle_deg": "float",
        "river_snap.angle_passes": "int",
        "river_snap.spike_angle_deg": "float",
        "river_snap.spike_length_m": "float",
        "river_snap.straight_angle_deg": "float",
        "river_snap.straight_length_factor": "float",
        "river_snap.spike_passes": "int",
        "river_snap.straight_max_deviation_m": "float",
        "river_snap.band_rows": "int",
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
        "country.protected_areas.gap_tolerance_m": "float",
        "country.protected_areas.touch_tolerance_m": "float",
        "country.protected_areas.self_close_min_length_factor": "float",
        "country.protected_areas.crest_search_radius_m": "float",
        "country.protected_areas.barrier_removal_width_m": "float",
        "country.protected_areas.min_margin_m": "float",
        "country.protected_areas.max_level_m": "float",
        "country.protected_areas.min_area_m2": "float",
        "country.protected_areas.gated_crossings": "bool",
        "country.protected_areas.fill_nodata_holes": "bool",
        "country.protected_areas.stream_burn_m": "float",
        "country.protected_areas.stream_burn_buffer_m": "float",
        "country.protected_areas.preserve_river_cells": "bool",
        "country.protected_areas.cap_river_holes": "bool",
        "country.protected_areas.max_fill_distance_m": "float",
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
        "sources.gfm_advisory_values": "int_list",
        "evidence.use_advisory_flags": "bool",
        "evidence.max_extent_files": "opt_paths",
        "events.min_peak_rp": "float",
        "events.days_before": "int",
        "events.days_after": "int",
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
        "country.dsm": "opt_paths",
        "country.building_footprints": "opt_paths",
        "country.manual_barriers": "opt_paths",
        "country.industrial_registry": "opt_paths",
        "country.ems_flood_extents": "opt_paths",
        "country.water_polygons": "opt_paths",
        "country.corine": "opt_paths",
        "country.corine_legend": "opt_path",
        "country.corine_legend_index_field": "str",
        "country.corine_legend_code_field": "str",
        "country.corine_check": "list",
        "country.context.order": "str_list",
        "country.context.corine_industrial_codes": "int_list",
        "country.context.corine_dense_urban_codes": "int_list",
        "country.context.urban_atlas_industrial_codes": "int_list",
        "country.context.urban_atlas_dense_urban_codes": "int_list",
        "country.context.ghsl_smod_dense_urban": "int_list",
        "country.context.ghsl_smod_urban": "int_list",
        "country.context.min_zone_area_m2": "float",
        **{
            f"country.context.params.{context}.{key}": kind
            for context in ("industrial", "dense_urban", "urban", "rural")
            for key, kind in (
                ("working_resolution_m", "float"),
                ("unit_core_size_m", "float"),
                ("gap_tolerance_m", "float"),
                ("closure_gap_max_m", "float"),
                ("crest_lines", "bool"),
                ("rp_distribution_classes", "int_list"),
                ("rp_distribution_shares", "float_list"),
                ("score_weight_value", "float"),
                ("score_weight_population", "float"),
                ("score_weight_upstream", "float"),
            )
        },
        **{
            f"crest_lines.{key}": "float"
            for key in (
                "transect_spacing_m",
                "transect_length_m",
                "max_bank_distance_m",
                "min_drop_m",
                "landward_window_m",
                "width_drop_m",
                "max_wall_width_m",
                "min_wall_height_m",
                "building_buffer_m",
                "max_crest_std_m",
                "link_max_m",
                "link_max_dz_m",
                "min_line_length_m",
            )
        },
        "barriers.source_priority": "str_list",
        "barriers.dedupe_distance_m": "float",
        "barriers.dedupe_max_angle_deg": "float",
        "barriers.dedupe_min_overlap": "float",
        "coarse.coarse_resolution_m": "float",
        "coarse.coarse_min_area_km2": "float",
        "coarse.coarse_match_distance_m": "float",
        "review.top_n": "int",
        "review.flag_weights": "table",
        # Stage E: implied protection (experimental)
        "implied.units": "str_list",
        "implied.recompute_hand": "bool",
        "implied.envelope_hand_m": "float",
        "implied.subsegment_length_m": "float",
        "implied.subsegment_length_grid": "float_list",
        "implied.built_source": "str",
        "implied.exclude_building_tags": "table",
        "implied.built_share_min": "float",
        "implied.built_share_min_grid": "float_list",
        "implied.built_area_min_m2": "float",
        "implied.hand_bin_m": "float",
        "implied.min_bin_area_m2": "float",
        "implied.smooth_bins": "int",
        "implied.onset_fraction": "float",
        "implied.onset_persist_bins": "int",
        "implied.quantiles": "float_list",
        "implied.primary_level_method": "str",
        "implied.sharp_window_m": "float",
        "implied.sharp_high": "float",
        "implied.component_built_min_m2": "float",
        "implied.merge_max_dlevel_m": "float",
        "implied.merge_max_gap": "int",
        "implied.corroboration_distance_m": "float",
        "implied.corroboration_min": "float",
        "implied.explicit_overlap_max": "float",
        "implied.regional_sop_rp": "int",
        "implied.implied_corroborated_max_rp": "int",
        "implied.historic_low_core_uncorroborated_cap": "bool",
        "implied.ghsl_epochs": "int_list",
        "implied.epoch_cutoff": "int",
        "implied.growth_low_share_min": "float",
        "implied.min_cells_temporal": "int",
        "implied.rp_distribution_classes": "int_list",
        "implied.rp_distribution_shares": "float_list",
        "implied.include_implied": "bool",
        "implied.include_in_review": "bool",
        **{
            f"implied.snap.{key}": "float"
            for key in (
                "step_m",
                "resolution_m",
                "max_offset_m",
                "min_window_m",
                "cost_epsilon",
                "chunk_length_m",
                "chunk_overlap_m",
                "max_length_ratio",
                "wide_river_m",
                "direction_check_m",
                "max_nodata_share",
            )
        },
        "implied.snap.chainage_origin": "str",
        "implied.plots.enabled": "bool",
        "implied.plots.dpi": "int",
        "implied.plots.max_per_type": "int",
        "implied.plots.sample_seed": "int",
        "implied.plots.show_max": "int",
        "implied.plots.plot_all_in_units": "str_list",
        # Stage F: discharge sources GloFAS and EFAS
        "discharge.sources": "str_list",
        "discharge.source_mode": "str",
        "discharge.ewds_api_url": "str",
        "discharge.raw_dir": "str",
        "discharge.aggregation": "str",
        "discharge.hydrological_year_start_month": "int",
        "discharge.min_years": "int",
        "discharge.max_rp_extrapolation": "float",
        "discharge.flooded_window_days": "int",
        "discharge.max_parallel_requests": "int",
        "discharge.bbox_buffer_deg": "float",
        "discharge.rp_disagreement_ratio": "float",
        **{
            f"discharge.{src}.{key}": kind
            for src in ("glofas", "efas")
            for key, kind in (
                ("dataset_id", "str"),
                ("request", "table"),
                ("data_format", "str"),
                ("download_format", "str"),
                ("first_year", "int"),
                ("last_year", "int"),
                ("request_unit", "str"),
                ("time_keys", "table"),
                ("area_key", "str"),
                ("step_hours", "int"),
                ("valid_time_marks", "str"),
                ("netcdf_variable", "str"),
                ("upstream_area_file", "str"),
                ("upstream_area_variable", "str"),
                ("upstream_area_to_km2", "float"),
                ("min_upstream_area_km2", "float"),
                ("search_radius_km", "float"),
                ("upstream_ratio_min", "float"),
                ("upstream_ratio_max", "float"),
                ("match_river_distance_m", "float"),
            )
        },
    }
    # Keys added after a stage was done: a missing key takes this default (logged), so older configurations stay valid.
    DEFAULTS = {
        # river groups of the protected areas per group (Stage G; defaults only, not in the toml files)
        "river_groups.max_join_ratio": 2.0,
        "river_groups.max_length_m": 30000.0,
        "river_groups.min_upstream_area_km2": 100.0,
        "river_groups.crop_margin_m": 5000.0,
        # levee lines onto the DTM crest, DTM-verified gaps, leak diagnostics (Stage G; defaults only)
        "levee_snap.densify_m": 10.0,
        "levee_snap.max_shift_m": 15.0,
        "levee_snap.side_m": 30.0,
        "levee_snap.min_height_m": 0.5,
        "levee_snap.level_tolerance_m": 0.1,
        "levee_snap.smooth_vertices": 5,
        "levee_snap.min_ridge_share": 0.5,
        "levee_gaps.max_gap_m": 200.0,
        "levee_gaps.crest_tolerance_m": 0.5,
        "leaks.sample_spacing_m": 50.0,
        # GFM advisory flags (Stage H): switch and values, also in config/advanced.toml
        "sources.gfm_advisory_values": [1, 2, 3],
        "evidence.use_advisory_flags": False,
        "evidence.max_extent_files": [],
        # flood events from discharge (Stage H), also in config/advanced.toml
        "events.min_peak_rp": 2.0,
        "events.days_before": 5,
        "events.days_after": 21,
        # transport embankments from the DTM (Stage G step 3; also in config/advanced.toml)
        "country.transport_lines": [],
        "transport.station_spacing_m": 10.0,
        "transport.half_width_m": 40.0,
        "transport.crest_search_m": 6.0,
        "transport.side_near_m": 12.0,
        "transport.side_far_m": 30.0,
        "transport.min_height_m": 0.5,
        "transport.footprint_fraction": 0.2,
        "transport.min_segment_m": 50.0,
        "transport.break_window_m": 150.0,
        "transport.break_depth_m": 1.0,
        "transport.max_opening_m": 150.0,
        "transport.culvert_window_m": 100.0,
        "transport.culvert_min_depth_m": 0.3,
        "transport.stream_distance_m": 10.0,
        "transport.cross_sections_per_class": 6,
        # comparison zone of the GFM tests, also the margin of the group rasters of file B
        "evidence.comparison_distance_m": 2000.0,
        # river snapping (Stage G): values of the snapping procedure the method follows
        "country.river_guides": [],
        "country.river_upstream_area_field": "",
        "river_snap.densify_m": 20.0,
        "river_snap.guide_search_m": 150.0,
        "river_snap.guide_max_shift_m": 180.0,
        "river_snap.window_m": 90.0,
        "river_snap.max_shift_m": 120.0,
        "river_snap.max_shift_without_guide_m": 0.0,
        "river_snap.min_angle_deg": 35.0,
        "river_snap.angle_passes": 2,
        "river_snap.spike_angle_deg": 75.0,
        "river_snap.spike_length_m": 120.0,
        "river_snap.straight_angle_deg": 175.0,
        "river_snap.straight_length_factor": 0.25,
        "river_snap.spike_passes": 6,
        "river_snap.straight_max_deviation_m": 2.0,
        "river_snap.band_rows": 512,
        # plot settings of Stage E, also used by Stage F scripts
        "implied.plots.enabled": True,
        "implied.plots.dpi": 120,
        "implied.plots.max_per_type": 200,
        "implied.plots.sample_seed": 42,
        "implied.plots.show_max": 10,
        "implied.plots.plot_all_in_units": [],
    }
    # Keys replaced by a later stage: a missing new key takes the value of its old key, with a warning (Stage F
    # replaces the GloFAS keys of [sources] and [evidence]).
    MIGRATIONS = {
        "discharge.glofas.dataset_id": "sources.glofas_dataset_id",
        "discharge.glofas.request": "sources.glofas_request",
        "discharge.glofas.time_keys": "sources.glofas_time_keys",
        "discharge.glofas.area_key": "sources.glofas_area_key",
        "discharge.glofas.first_year": "sources.glofas_first_year",
        "discharge.glofas.last_year": "sources.glofas_last_year",
        "discharge.glofas.netcdf_variable": "sources.glofas_discharge_variable",
        "discharge.glofas.upstream_area_file": "sources.glofas_upstream_area_file",
        "discharge.glofas.upstream_area_variable": "sources.glofas_upstream_area_variable",
        "discharge.glofas.upstream_area_to_km2": "sources.glofas_upstream_area_to_km2",
        "discharge.glofas.min_upstream_area_km2": "sources.glofas_min_upstream_area_km2",
        "discharge.glofas.search_radius_km": "evidence.glofas_search_radius_km",
        "discharge.glofas.upstream_ratio_min": "evidence.upstream_ratio_min",
        "discharge.glofas.upstream_ratio_max": "evidence.upstream_ratio_max",
        "discharge.flooded_window_days": "evidence.flooded_window_days",
        # one value for both sources until 2026-11-05; now per source (GloFAS cells are much larger)
        "discharge.glofas.match_river_distance_m": "discharge.match_river_distance_m",
        "discharge.efas.match_river_distance_m": "discharge.match_river_distance_m",
    }
    CONTEXTS = ("industrial", "dense_urban", "urban", "rural")
    ALLOWED = {
        "country.dtm_conversion": {"A", "B"},
        "country.dtm_compress": {"ZSTD", "DEFLATE", "LZW"},
        "general.log_level": {"DEBUG", "INFO", "WARNING", "ERROR"},
        "implied.built_source": {"footprints", "ghsl", "max"},
        "discharge.source_mode": {"prefer_efas", "prefer_glofas", "conservative_both"},
        "discharge.aggregation": {"daily_mean", "daily_max"},
        **{
            f"discharge.{s}.request_unit": {"year", "month"} for s in ("glofas", "efas")
        },
        **{
            f"discharge.{s}.valid_time_marks": {"end", "start"}
            for s in ("glofas", "efas")
        },
        "implied.snap.chainage_origin": {"downstream", "upstream"},
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
        raw = Config.read_merged(path)
        countries = raw.get("countries", {})
        if country not in countries:
            raise SystemExit(
                f"Country '{country}' is not configured in {path}; "
                f"available: {sorted(countries)}"
            )
        cfg = {
            k: v for k, v in raw.items() if k not in ("countries", "country_defaults")
        }
        # country values: [country_defaults] of advanced.toml, overridden by [countries.<code>]
        cfg["country"] = Config.merge(
            raw.get("country_defaults", {}), countries[country], "country"
        )
        cfg["country_code"] = country
        cfg["used_keys"] = list(keys)
        list_err = []
        for key in keys:
            kind = Config.SCHEMA.get(key)
            if kind is None:
                list_err.append(f"{key}: not a known configuration key")
                continue
            try:
                value = Config.get(cfg, key)
            except (KeyError, TypeError):
                old = Config.MIGRATIONS.get(key)
                try:
                    value = Config.get(cfg, old) if old else None
                except (KeyError, TypeError):
                    value = None
                if value is not None:
                    log.warning(
                        "%s missing: taken from the replaced key %s; update the configuration.",
                        key,
                        old,
                    )
                elif key in Config.DEFAULTS:
                    value = Config.DEFAULTS[key]
                    log.info(
                        "%s not configured: default %r.", key, Config.DEFAULTS[key]
                    )
                else:
                    list_err.append(f"{key}: missing")
                    continue
                node = cfg
                for part in key.split(".")[:-1]:
                    node = node.setdefault(part, {})
                node[key.split(".")[-1]] = value
            checked, problem = Config._check_value(kind, value)
            if problem:
                list_err.append(f"{key}: {problem}")
                continue
            if key in Config.ALLOWED and checked not in Config.ALLOWED[key]:
                list_err.append(
                    f"{key}: {checked!r} not in {sorted(Config.ALLOWED[key])}"
                )
            node = cfg
            list_part = key.split(".")
            for part in list_part[:-1]:
                node = node[part]
            node[list_part[-1]] = checked
        if not list_err:
            list_err = Config._check_consistency(cfg, set(keys))
        if list_err:
            raise SystemExit(
                f"Invalid configuration {path} (country {country}):\n  - "
                + "\n  - ".join(list_err)
            )
        return cfg

    @staticmethod
    def merge(base, over, prefix=""):
        """
        Values of over on top of base: tables are merged key by key, except keys whose
        schema type is "table" (a request, a template set), which over replaces whole
        :param base: dict of defaults
        :param over: dict of user values
        :param prefix: dotted schema path of both dicts ("" for the top level)
        :return: new dict
        """
        out = dict(base)
        for k, v in over.items():
            path = f"{prefix}.{k}" if prefix else k
            if (
                isinstance(v, dict)
                and isinstance(out.get(k), dict)
                and Config.SCHEMA.get(path) != "table"
            ):
                out[k] = Config.merge(out[k], v, path)
            else:
                out[k] = v
        return out

    # country sub-tables renamed since a stage was done: old name -> new name (read with a warning)
    RENAMED_TABLES = {"compartments": "protected_areas"}

    @staticmethod
    def rename_tables(raw, path):
        """
        Moves renamed country sub-tables (RENAMED_TABLES) to their new name, in every
        [countries.<code>] table and in [country_defaults]; logs a warning for each
        :param raw: configuration dict as read
        :param path: file name for the warning
        """
        tables = list(raw.get("countries", {}).values())
        if isinstance(raw.get("country_defaults"), dict):
            tables.append(raw["country_defaults"])
        for table in tables:
            for old, new in Config.RENAMED_TABLES.items():
                if isinstance(table, dict) and old in table and new not in table:
                    table[new] = table.pop(old)
                    log.warning(
                        "%s: table %s renamed to %s; update the configuration.",
                        path,
                        old,
                        new,
                    )

    @staticmethod
    def toml_literal(value):
        """
        TOML text of a value: string, bool, int, float or a list of them
        :param value: Python value
        :return: str
        """
        if isinstance(value, bool):
            return "true" if value else "false"
        if isinstance(value, (int, float)):
            return repr(value)
        if isinstance(value, (list, tuple)):
            return "[" + ", ".join(Config.toml_literal(v) for v in value) + "]"
        return json.dumps(str(value), ensure_ascii=False)

    @staticmethod
    def set_toml_value(path, table, key, value):
        """
        Sets one key of a TOML file in place: the line of the key keeps its indentation and
        trailing comment, every other line stays as it is. A missing key is added at the end
        of its table, a missing table at the end of the file. The result must parse, otherwise
        the file is left unchanged.
        :param path: TOML file
        :param table: dotted table name ("general", "countries.POL.rp"); "" = top level
        :param key: key name
        :param value: new value (see toml_literal)
        :return: None
        """
        with open(str(path), encoding="utf-8") as f:
            text = f.read()
        lines = text.splitlines()
        header = re.compile(r"^\s*(\[\[?)\s*([^\[\]]+?)\s*\]\]?\s*(#.*)?$")
        heads = [
            (i, m.group(1), m.group(2))
            for i, line in enumerate(lines)
            if (m := header.match(line))
        ]
        if table == "":
            start, end = 0, (heads[0][0] if heads else len(lines))
        else:
            hit = [k for k, h in enumerate(heads) if h[1] == "[" and h[2] == table]
            start = heads[hit[0]][0] + 1 if hit else None
            end = (
                heads[hit[0] + 1][0] if hit and hit[0] + 1 < len(heads) else len(lines)
            )
        literal = Config.toml_literal(value)
        pattern = re.compile(rf"^(\s*){re.escape(key)}\s*=\s*(.*)$")
        if start is None:
            lines += ["", f"[{table}]", f"{key} = {literal}"]
        else:
            done = False
            for i in range(start, end):
                m = pattern.match(lines[i])
                if not m:
                    continue
                # the value may continue on the next lines (multi-line arrays)
                last = i
                while last < end - 1:
                    try:
                        tomllib.loads(
                            "v = " + "\n".join([m.group(2)] + lines[i + 1 : last + 1])
                        )
                        break
                    except tomllib.TOMLDecodeError:
                        last += 1
                comment = ""
                if last == i:
                    # the comment starts at the first '#' after a value that parses on its own
                    rest = m.group(2)
                    for pos in [k for k, ch in enumerate(rest) if ch == "#"]:
                        try:
                            tomllib.loads(f"v = {rest[:pos]}")
                        except tomllib.TOMLDecodeError:
                            continue
                        gap = len(rest[:pos]) - len(rest[:pos].rstrip())
                        comment = " " * max(gap, 1) + rest[pos:]
                        break
                lines[i : last + 1] = [f"{m.group(1)}{key} = {literal}{comment}"]
                done = True
                break
            if not done:
                last = start - 1
                for i in range(start, end):
                    if lines[i].strip() and not lines[i].lstrip().startswith("#"):
                        last = i
                lines.insert(last + 1, f"{key} = {literal}")
        new_text = "\n".join(lines) + ("\n" if text.endswith("\n") or not text else "")
        tomllib.loads(new_text)
        with open(str(path), "w", encoding="utf-8") as f:
            f.write(new_text)

    @staticmethod
    def read_merged(path):
        """
        Reads a configuration file on top of its defaults: advanced.toml in the same
        folder (or the file named by the top-level key "advanced"; "" = none). Without
        that file the configuration is read alone. Placeholder lists are joined.
        :param path: TOML file of the user
        :return: dict
        """
        if not exists(str(path)):
            raise SystemExit(f"Configuration file not found: {path}")
        with open(str(path), "rb") as f:
            raw = tomllib.load(f)
        Config.rename_tables(raw, str(path))
        name = raw.pop("advanced", "advanced.toml")
        base_path = (
            join(os.path.dirname(os.path.abspath(str(path))), name) if name else ""
        )
        if (
            not base_path
            or not exists(base_path)
            or os.path.samefile(base_path, str(path))
        ):
            return raw
        with open(base_path, "rb") as f:
            base = tomllib.load(f)
        Config.rename_tables(base, base_path)
        placeholders = list(
            dict.fromkeys(base.get("placeholders", []) + raw.get("placeholders", []))
        )
        merged = Config.merge(base, raw)
        merged["placeholders"] = placeholders
        return merged

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
            list_path = Config.as_list(value)
            if not list_path and kind == "paths":
                return list_path, "at least one path is required"
            list_miss = [p for p in list_path if not exists(p)]
            return list_path, f"file(s) not found: {list_miss}" if list_miss else None
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
        if kind == "str_list":
            if not isinstance(value, list) or not all(
                isinstance(v, str) for v in value
            ):
                return value, f"expected a list of strings, got {value!r}"
            return value, None
        if kind == "list":
            return value, None if isinstance(value, list) else "expected a list"
        return value, f"unknown kind {kind}"

    @staticmethod
    def _check_consistency(cfg, keys):
        list_err = []
        c = cfg["country"]
        if "country.working_resolution_m" in keys and c["working_resolution_m"] <= 0:
            list_err.append("country.working_resolution_m must be positive")
        if {"country.unit_core_size_m", "country.working_resolution_m"} <= keys:
            if c["unit_core_size_m"] < 2 * c["working_resolution_m"]:
                list_err.append(
                    "country.unit_core_size_m must be at least two working cells"
                )
        if {"country.unit_halo_m", "country.corridor_buffer_m"} <= keys:
            if c["unit_halo_m"] < c["corridor_buffer_m"]:
                # A core cell may lie corridor_buffer_m from its river; HAND needs that river inside the window.
                list_err.append(
                    "country.unit_halo_m must be at least country.corridor_buffer_m, otherwise "
                    "cells of a core can have their river outside the unit window"
                )
        if "general.n_workers" in keys and cfg["general"]["n_workers"] < 1:
            list_err.append("general.n_workers must be at least 1")
        if (
            "country.merge_overlap_ratio" in keys
            and not 0 < c["merge_overlap_ratio"] <= 1
        ):
            list_err.append("country.merge_overlap_ratio must be in (0, 1]")
        if "benchmark.window_px" in keys and cfg["benchmark"]["window_px"] < 1:
            list_err.append("benchmark.window_px must be positive")
        rp = c.get("rp", {})
        if {
            "country.rp.distribution_classes",
            "country.rp.distribution_shares",
        } <= keys:
            if len(rp["distribution_classes"]) != len(rp["distribution_shares"]):
                list_err.append(
                    "rp.distribution_classes and rp.distribution_shares differ in length"
                )
            elif abs(sum(rp["distribution_shares"]) - 1.0) > 1e-6:
                list_err.append("rp.distribution_shares must sum to 1")
        if {"country.rp.cap_levels_m", "country.rp.cap_rp"} <= keys:
            if len(rp["cap_levels_m"]) != len(rp["cap_rp"]):
                list_err.append("rp.cap_levels_m and rp.cap_rp differ in length")
            elif list(rp["cap_levels_m"]) != sorted(rp["cap_levels_m"]):
                list_err.append("rp.cap_levels_m must be ascending")
        if "country.rp.classes" in keys and list(rp["classes"]) != sorted(
            set(rp["classes"])
        ):
            list_err.append("rp.classes must be strictly ascending")
        if {"country.rp.classes", "country.rp.distribution_classes"} <= keys:
            if not set(rp["distribution_classes"]) <= set(rp["classes"]):
                list_err.append("rp.distribution_classes must be values of rp.classes")
        if {"country.rp.classes", "country.rp.cap_rp", "country.rp.max_rp"} <= keys:
            if not set(rp["cap_rp"]) | {rp["max_rp"]} <= set(rp["classes"]):
                list_err.append("rp.cap_rp and rp.max_rp must be values of rp.classes")
        evidence = cfg.get("evidence", {})
        if {"evidence.f_dry_inside", "evidence.f_wet_inside"} <= keys:
            if not evidence["f_dry_inside"] < evidence["f_wet_inside"] <= 1:
                list_err.append(
                    "evidence: f_dry_inside < f_wet_inside <= 1 is required"
                )
        if {"evidence.upstream_ratio_min", "evidence.upstream_ratio_max"} <= keys:
            if (
                not 0
                < evidence["upstream_ratio_min"]
                <= 1
                <= evidence["upstream_ratio_max"]
            ):
                list_err.append(
                    "evidence: upstream_ratio_min <= 1 <= upstream_ratio_max is required"
                )
        if "evidence.river_level_percentile" in keys:
            if not 0 < evidence["river_level_percentile"] <= 100:
                list_err.append("evidence.river_level_percentile must be in (0, 100]")
        if (
            "evidence.hint_dry_share" in keys
            and not 0 < evidence["hint_dry_share"] <= 1
        ):
            list_err.append("evidence.hint_dry_share must be in (0, 1]")
        if "evidence.hint_exclusion_m" in keys and evidence["hint_exclusion_m"] < 0:
            list_err.append("evidence.hint_exclusion_m must not be negative")
        if (
            "sources.gfm_download_threads" in keys
            and cfg["sources"]["gfm_download_threads"] < 1
        ):
            list_err.append("sources.gfm_download_threads must be at least 1")
        if (
            "sources.gfm_screen_resolution_m" in keys
            and cfg["sources"]["gfm_screen_resolution_m"] <= 0
        ):
            list_err.append("sources.gfm_screen_resolution_m must be positive")
        if "country.context.order" in keys:
            list_order = c["context"]["order"]
            if sorted(list_order) != sorted(Config.CONTEXTS):
                list_err.append(
                    f"context.order must list each of {list(Config.CONTEXTS)} exactly once"
                )
            elif list_order[-1] != "rural":
                list_err.append(
                    "context.order must end with rural (rural is the rest of the corridor)"
                )
        lulc = c.get("lulc", {})
        for table in ["worldcover_categories", "urban_atlas_categories"]:
            if f"country.lulc.{table}" not in keys:
                continue
            for code, category in lulc[table].items():
                if table == "worldcover_categories" and not (
                    str(code).isdigit() and 0 < int(code) < 256
                ):
                    list_err.append(
                        f"lulc.worldcover_categories: class {code!r} is not a raster value 1-255"
                    )
                if not (
                    isinstance(category, str)
                    and category.replace("_", "").isalnum()
                    and category.isascii()
                ):
                    list_err.append(
                        f"lulc.{table}: category {category!r} must be a name of letters, digits and _"
                    )
        if {"country.lulc.urban_atlas", "country.lulc.urban_atlas_class_field"} <= keys:
            if lulc["urban_atlas"] and not lulc["urban_atlas_class_field"]:
                list_err.append(
                    "lulc.urban_atlas is set but lulc.urban_atlas_class_field is empty"
                )
        for date_key in ["sources.gfm_date_from", "sources.gfm_date_to"]:
            value = cfg.get("sources", {}).get(date_key.split(".")[1], "")
            if date_key in keys and value:
                try:
                    datetime.date.fromisoformat(value)
                except ValueError:
                    list_err.append(f"{date_key} must be YYYY-MM-DD or empty")
        if {"sources.glofas_first_year", "sources.glofas_last_year"} <= keys:
            if cfg["sources"]["glofas_first_year"] > cfg["sources"]["glofas_last_year"]:
                list_err.append(
                    "sources.glofas_first_year must not exceed glofas_last_year"
                )
        implied = cfg.get("implied", {})
        snap = implied.get("snap", {})
        for key in ["resolution_m", "step_m", "chunk_length_m", "cost_epsilon"]:
            if f"implied.snap.{key}" in keys and snap[key] <= 0:
                list_err.append(f"implied.snap.{key} must be positive")
        if {"implied.snap.chunk_length_m", "implied.snap.chunk_overlap_m"} <= keys:
            if snap["chunk_overlap_m"] >= snap["chunk_length_m"]:
                list_err.append(
                    "implied.snap.chunk_overlap_m must be shorter than chunk_length_m"
                )
        if "implied.snap.max_length_ratio" in keys and snap["max_length_ratio"] < 1:
            list_err.append("implied.snap.max_length_ratio must be at least 1")
        if "implied.snap.max_nodata_share" in keys and snap["max_nodata_share"] > 1:
            list_err.append("implied.snap.max_nodata_share must not exceed 1")
        if {"implied.snap.min_window_m", "implied.snap.resolution_m"} <= keys:
            if snap["min_window_m"] < snap["resolution_m"]:
                list_err.append(
                    "implied.snap.min_window_m must be at least one cell (resolution_m)"
                )
        if {"implied.primary_level_method", "implied.quantiles"} <= keys:
            methods = {"onset"} | {f"q{q:.2f}" for q in implied["quantiles"]}
            if implied["primary_level_method"] not in methods:
                list_err.append(
                    f"implied.primary_level_method must be one of {sorted(methods)}"
                )
        if "implied.plots.dpi" in keys and implied["plots"]["dpi"] < 1:
            list_err.append("implied.plots.dpi must be at least 1")
        dis = cfg.get("discharge", {})
        if "discharge.sources" in keys and not set(dis["sources"]) <= {
            "glofas",
            "efas",
        }:
            list_err.append("discharge.sources may contain glofas and efas only")
        if (
            "discharge.hydrological_year_start_month" in keys
            and not 1 <= dis["hydrological_year_start_month"] <= 12
        ):
            list_err.append("discharge.hydrological_year_start_month must be 1-12")
        for src in ("glofas", "efas"):
            if f"discharge.{src}.step_hours" in keys and dis[src]["step_hours"] not in (
                1,
                3,
                6,
                12,
                24,
            ):
                list_err.append(f"discharge.{src}.step_hours must divide 24 hours")
            if {
                f"discharge.{src}.upstream_ratio_min",
                f"discharge.{src}.upstream_ratio_max",
            } <= keys:
                if dis[src]["upstream_ratio_min"] > dis[src]["upstream_ratio_max"]:
                    list_err.append(
                        f"discharge.{src}.upstream_ratio_min exceeds upstream_ratio_max"
                    )
        return list_err

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
        # COUNTRY left at the script default "XX": the country named in [general] country of the configuration
        if (
            values.get("COUNTRY") == "XX"
            and os.environ.get("DEFENSES_COUNTRY") is None
            and "CONFIG_PATH" in values
        ):
            try:
                general = Config.read_merged(values["CONFIG_PATH"]).get("general", {})
            except (SystemExit, OSError, tomllib.TOMLDecodeError):
                general = {}
            if general.get("country"):
                values["COUNTRY"] = str(general["country"])
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
    def layer_files(layer):
        """
        Files of a layer of protected areas: "A" with all water sources together, "B" per river group
        :param layer: "A" or "B"
        :return: dict layer, areas (protected areas), rp (with return periods), suffix (of derived outputs)
        """
        layer = str(layer).strip().upper()
        if layer == "A":
            return {
                "layer": "A",
                "areas": "protected_areas.gpkg",
                "rp": "protected_areas_rp.gpkg",
                "suffix": "",
            }
        if layer == "B":
            return {
                "layer": "B",
                "areas": "protected_areas_by_group.gpkg",
                "rp": "protected_areas_by_group_rp.gpkg",
                "suffix": "_b",
            }
        raise SystemExit(f"LAYER must be A or B, not {layer!r}.")

    @staticmethod
    def unit_layer_path(unit_dir, name, group_id=None):
        """
        Path of a unit file of a layer: the file itself for file A (group_id None), groups/<stem>_g<id><ext> for
        a river group of file B
        :param unit_dir: unit folder
        :param name: file name, e.g. "labels.tif"
        :param group_id: river group, None or NaN for file A
        :return: str
        """
        if group_id is None or (isinstance(group_id, float) and np.isnan(group_id)):
            return join(str(unit_dir), name)
        stem, ext = os.path.splitext(name)
        os.makedirs(join(str(unit_dir), "groups"), exist_ok=True)
        return join(str(unit_dir), "groups", f"{stem}_g{int(group_id)}{ext}")

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
        list_res, list_fail = [], []
        if not tasks:
            return list_res, list_fail
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
                    list_res.append(result)
                    message = result.get("summary", "")
                except BrokenProcessPool as error:
                    list_fail.append(task["unit_id"])
                    message = f"worker process died ({error}); possibly out of memory"
                    logger.error("Unit %s: %s", task["unit_id"], message)
                except Exception:
                    list_fail.append(task["unit_id"])
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
                    len(list_res) + len(list_fail),
                    len(tasks),
                    message,
                    Run.format_seconds(elapsed),
                    Run.format_seconds(eta),
                )
        return list_res, list_fail


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
    def runs(values, min_length=1):
        """
        Runs of equal values in a sequence; runs shorter than min_length take the value of the run before them
        (the first run: of the run after it), then neighbours with the same value are joined
        :param values: 1D sequence
        :param min_length: shortest run kept
        :return: list of (start index, end index exclusive, value)
        """
        values = list(values)
        list_run = []
        for k, v in enumerate(values):
            if list_run and list_run[-1][2] == v:
                list_run[-1][1] = k + 1
            else:
                list_run.append([k, k + 1, v])
        changed = True
        while changed and len(list_run) > 1:
            changed = False
            for k, (s, e, v) in enumerate(list_run):
                if e - s < min_length:
                    other = list_run[k - 1] if k > 0 else list_run[k + 1]
                    list_run[k][2] = other[2]
                    changed = True
                    break
            joined = [list_run[0]]
            for s, e, v in list_run[1:]:
                if joined[-1][2] == v:
                    joined[-1][1] = e
                else:
                    joined.append([s, e, v])
            list_run = joined
        return [tuple(r) for r in list_run]

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
    def gumbel_q(rp, xi, alpha):
        """
        Discharge of a return period under a Gumbel distribution (inverse of gumbel_rp):
        xi - alpha * ln(-ln(1 - 1 / rp))
        :param rp: return period(s) in years, > 1
        :return: float or array
        """
        return xi - alpha * np.log(-np.log1p(-1.0 / np.asarray(rp, dtype=float)))

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

    @staticmethod
    def binned_quantiles(area, edges, quantiles):
        """
        Area-weighted quantiles of a variable given as area per bin, interpolated
        linearly within the bin that reaches the quantile
        :param area: (n, n_bins) area per bin and row
        :param edges: (n_bins + 1,) bin edges
        :param quantiles: quantiles in [0, 1]
        :return: (n, len(quantiles)) values; NaN for rows without area
        """
        area = np.atleast_2d(np.asarray(area, dtype=float))
        edges = np.asarray(edges, dtype=float)
        rows = np.arange(len(area))
        total = area.sum(axis=1)
        cum = np.divide(
            np.cumsum(area, axis=1),
            total[:, None],
            out=np.zeros_like(area),
            where=total[:, None] > 0,
        )
        out = np.full((len(area), len(quantiles)), np.nan)
        for j, q in enumerate(quantiles):
            b = np.argmax(cum >= q - 1e-12, axis=1)
            before = np.where(b > 0, cum[rows, np.maximum(b - 1, 0)], 0.0)
            here = np.divide(
                area[rows, b], total, out=np.zeros(len(area)), where=total > 0
            )
            t = np.divide(q - before, here, out=np.zeros(len(area)), where=here > 0)
            out[:, j] = edges[b] + np.clip(t, 0.0, 1.0) * (edges[b + 1] - edges[b])
        out[total <= 0] = np.nan
        return out

    @staticmethod
    def onset(share, smooth_bins, fraction, persist_bins):
        """
        Onset of a share curve per row: the lowest bin where the smoothed share (centred
        moving average over smooth_bins, NaN ignored) is at least fraction times its
        maximum and stays so for persist_bins consecutive bins
        :param share: (n, n_bins) share per bin, NaN for bins that are too small
        :param smooth_bins: width of the moving average in bins
        :param fraction: share of the maximum that starts the onset
        :param persist_bins: consecutive bins that must stay above it, the onset bin included
        :return: (onset bin per row, -1 where none; smoothed share (n, n_bins), NaN where
                  the bin itself is NaN)
        """
        share = np.atleast_2d(np.asarray(share, dtype=float))
        n, m = share.shape
        valid = np.isfinite(share)
        before = smooth_bins // 2
        after = smooth_bins - 1 - before
        pad = ((0, 0), (before + 1, after))
        sums = np.cumsum(np.pad(np.where(valid, share, 0.0), pad), axis=1)
        counts = np.cumsum(np.pad(valid.astype(float), pad), axis=1)
        win_sum = sums[:, smooth_bins:] - sums[:, :m]
        win_count = counts[:, smooth_bins:] - counts[:, :m]
        smoothed = np.divide(
            win_sum, win_count, out=np.full((n, m), np.nan), where=win_count > 0
        )
        smoothed[~valid] = np.nan
        peak = np.nanmax(np.where(np.isfinite(smoothed), smoothed, -np.inf), axis=1)
        above = np.isfinite(smoothed) & (smoothed >= fraction * peak[:, None])
        above &= peak[:, None] > 0
        run = np.zeros((n, m + 1), np.int64)
        for i in range(m - 1, -1, -1):
            run[:, i] = np.where(above[:, i], run[:, i + 1] + 1, 0)
        ok = run[:, :m] >= persist_bins
        return np.where(ok.any(axis=1), np.argmax(ok, axis=1), -1), smoothed

    @staticmethod
    def sharpness(smoothed, onset_bin, window_bins):
        """
        Sharpness of an onset: (mean smoothed share in [onset, onset + window) - mean in
        [onset - window, onset)) / maximum smoothed share; NaN bins are ignored, an empty
        window below the onset counts as 0
        :param smoothed: (n, n_bins) smoothed share (Stats.onset)
        :param onset_bin: onset bin per row, -1 = none (NaN result)
        :param window_bins: window width in bins
        :return: (n,) sharpness
        """
        smoothed = np.atleast_2d(np.asarray(smoothed, dtype=float))
        n, m = smoothed.shape
        valid = np.isfinite(smoothed)
        sums = np.concatenate(
            [np.zeros((n, 1)), np.cumsum(np.where(valid, smoothed, 0.0), axis=1)],
            axis=1,
        )
        counts = np.concatenate([np.zeros((n, 1)), np.cumsum(valid, axis=1)], axis=1)
        rows = np.arange(n)
        o = np.clip(np.asarray(onset_bin), 0, m)
        hi = np.minimum(o + window_bins, m)
        lo = np.maximum(o - window_bins, 0)
        up = np.divide(
            sums[rows, hi] - sums[rows, o],
            counts[rows, hi] - counts[rows, o],
            out=np.full(n, np.nan),
            where=counts[rows, hi] - counts[rows, o] > 0,
        )
        down = np.divide(
            sums[rows, o] - sums[rows, lo],
            counts[rows, o] - counts[rows, lo],
            out=np.zeros(n),
            where=counts[rows, o] - counts[rows, lo] > 0,
        )
        peak = np.nanmax(np.where(valid, smoothed, -np.inf), axis=1)
        result = np.divide(up - down, peak, out=np.full(n, np.nan), where=peak > 0)
        result[np.asarray(onset_bin) < 0] = np.nan
        return result

    @staticmethod
    def merge_levels(k, eligible, level, max_dlevel, max_gap):
        """
        Groups of subsegments along one river side. Consecutive eligible subsegments
        whose levels differ by at most max_dlevel form a group; a run of at most max_gap
        positions without an eligible subsegment between two such neighbours is bridged
        with the mean of their levels.
        :param k: subsegment positions (unique)
        :param eligible: bool per subsegment
        :param level: level per subsegment (used where eligible)
        :param max_dlevel: largest level difference within a group
        :param max_gap: longest bridged run of positions
        :return: (group per subsegment, -1 = none; level used, NaN = none; bridged flag)
        """
        k = np.asarray(k)
        order = np.argsort(k)
        ks = k[order]
        el = np.asarray(eligible, bool)[order]
        lv = np.asarray(level, float)[order]
        n = len(ks)
        group = np.full(n, -1, np.int64)
        used = np.full(n, np.nan)
        bridged = np.zeros(n, bool)
        g, last = -1, -1
        for i in np.flatnonzero(el):
            if last >= 0:
                gap = ks[i] - ks[last] - 1
                if gap <= max_gap and abs(lv[i] - lv[last]) <= max_dlevel:
                    between = np.arange(last + 1, i)
                    group[between] = g
                    used[between] = (lv[i] + lv[last]) / 2
                    bridged[between] = True
                    group[i], used[i], last = g, lv[i], i
                    continue
            g += 1
            group[i], used[i], last = g, lv[i], i
        out = [np.empty_like(a) for a in (group, used, bridged)]
        for a, b in zip(out, (group, used, bridged)):
            a[order] = b
        return out[0], out[1], out[2]
