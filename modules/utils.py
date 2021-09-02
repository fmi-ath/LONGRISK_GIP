from pathlib import Path
import numpy as np
import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.mask import mask
from rasterio.merge import merge
from shapely.geometry import Polygon
import geopandas as gpd
from fiona.crs import from_epsg
import pycrs
from pyproj import Proj, transform
from osgeo import gdal, osr
import os, glob
from math import floor

def ascii_to_geotiff(ascii_files_path, GTiff_files_path, xmin, ymax, search_criteria = '*.txt', CRS_code = 3067,
    xmax = None, xrotation = 0, xres = None, yrotation = 0, ymin = None, yres = None):

    """A given ascii file(s) is(are) transformed into GeoTiff file.

        Parameters
        ----------
        ascii_files_path : str
            Path to ascii file(s).
        GTiff_files_path : str
            Path where GTiff files will be located.
        xmin : float
            X coordinate of upper-left corner of the image
        ymax : float
            Y coordinate of upper-left corner of the image
        search_criteria : list or str
            List of criteria for searhing the ascii files to be converted.
        CRS_code : into
            EPSG code for which the coordinates are stablished
        xmax : float
            X coordinate of lower-right corner of the image
        ymin : float
            Y coordinate of lower-right corner of the image
        xres : float
            X resolution of the image
        yres : float
            Y resolution of the image
        xrotation : float
            X rotation of the image
        yrotation : float
            Y rotation of the image


        Returns
        -------
        out : file
            GTiff converted file
    """

    if os.path.splitext(os.path.split(ascii_files_path)[1])[1] == '':

        if not isinstance(search_criteria, list):

            search_criteria = [search_criteria]

        path = []

        for hints in search_criteria:

            # Make a search criteria to select the DEM files
            q = os.path.join(ascii_files_path, hints)
            #print(q)

            # glob function can be used to list files from a directory with specific criteria
            path += glob.glob(q)

    else:

        path = glob.glob(raster_to_crop_path)

    for file in path:

        array = np.genfromtxt(file, dtype=np.float32)

        nrows, ncols = np.shape(array)

        if (xres is None) and (xmax is not None):

            xres = (xmax - xmin) / float(ncols)

        elif (xres is None) and (xmax is None):

            raise AttributeError('If xres not given, xmax has to be given')

        if (yres is None) and (ymin is not None):

            yres = (ymax - ymin) / float(nrows)

        elif (yres is None) and (ymin is None):

            raise AttributeError('If yres not given, ymin has to be given')

        # Upper-left corner coordinates
        geotransform = (xmin, xres, xrotation, ymax, yrotation, -yres)

        GTiff_destination_file = os.path.join(GTiff_files_path, os.path.splitext(os.path.split(file)[1])[0] + '.tif')

        output_raster = gdal.GetDriverByName('GTiff').Create(GTiff_destination_file, ncols, nrows, 1 ,gdal.GDT_Float32)  # Open the file
        output_raster.SetGeoTransform(geotransform)  # Specify its coordinates
        srs = osr.SpatialReference()                 # Establish its coordinate encoding
        srs.ImportFromEPSG(CRS_code)                     # This one specifies WGS84 lat long.
                                                     # Anyone know how to specify the
                                                     # IAU2000:49900 Mars encoding?
        output_raster.SetProjection( srs.ExportToWkt() )   # Exports the coordinate system
                                                           # to the file
        output_raster.GetRasterBand(1).WriteArray(array)   # Writes my array to the raster

        output_raster.FlushCache()

    return

def rain_relocation(GTiff_files_path, xmin, ymax, X_target, Y_target, X_radar_rain,
                    Y_radar_rain, search_criteria = '*.txt', CRS_code = 3067, xmax = None, xrotation = 0,
                    xres = None, yrotation = 0, ymin = None, yres = None):

    """Relocates a desired pixel region of the image (X_radar_rain, Y_radar_rain) to a desired region (X_target, Y_target).
        Produces a new GTiff file.

        Parameters
        ----------
        GTiff_files_path : str
            Path where GTiff files will be located.
        xmin : float
            X coordinate of upper-left corner of the image
        ymax : float
            Y coordinate of upper-left corner of the image
        X_target : float
            X coordinate of region over which rain pixels will be positioned
        Y_target : float
            Y coordinate of region over which rain pixels will be positioned
        X_radar_rain : float
            X distance between rain pixels and radar (center of the image)
        Y_radar_rain : float
            Y distance between rain pixels and radar (center of the image)
        search_criteria : list or str
            List of criteria for searhing the ascii files to be converted.
        CRS_code : into
            EPSG code for which the coordinates are stablished
        xmax : float
            X coordinate of lower-right corner of the image
        ymin : float
            Y coordinate of lower-right corner of the image
        xres : float
            X resolution of the image
        yres : float
            Y resolution of the image
        xrotation : float
            X rotation of the image
        yrotation : float
            Y rotation of the image


        Returns
        -------
        out : file
            GTiff relocated file
    """

    if not Path(GTiff_files_path).suffix:

        if not isinstance(search_criteria, list):

            search_criteria = [search_criteria]

        path = []

        for hints in search_criteria:

            # Make a search criteria to select the DEM files
            q = os.path.join(GTiff_files_path, hints)
            #print(q)

            # glob function can be used to list files from a directory with specific criteria
            path += glob.glob(q)

    else:

        path = glob.glob(GTiff_files_path)

    for file in path:

        dataset = gdal.Open(file, gdal.GA_ReadOnly) # Note GetRasterBand() takes band no. starting from 1 not 0
        proj = osr.SpatialReference(wkt = dataset.GetProjection())
        geotransform = dataset.GetGeoTransform()
        band = dataset.GetRasterBand(1)
        array = band.ReadAsArray()

        nrows, ncols = np.shape(array)

        if (xres is None) and (xmax is not None):

            xres = (xmax - xmin) / float(ncols)

        elif (xres is None) and (xmax is None):

            raise AttributeError('If xres not given, xmax has to be given')

        if (yres is None) and (ymin is not None):

            yres = (ymax - ymin) / float(nrows)

        elif (yres is None) and (ymin is None):

            raise AttributeError('If yres not given, ymin has to be given')

        X = xmin
        Y = ymax

        x_rel = X_radar_rain
        y_rel = Y_radar_rain

        # Distance between radar (center of the image) and upper left corner of the image
        x = X + ncols / 2 * xres
        y = Y - nrows / 2 * yres

        # Distance between Copenhagen and upper-left corner of the radar image (ncols = nrows = 480)
        x_C = x + x_rel
        y_C = y + y_rel

        x_new = 2 * X - x_C
        y_new = 2 * Y - y_C

        Xmin = x_new - (X - X_target)
        Ymax = y_new - (Y - Y_target)

        # Upper-left corner coordinates
        geotransform = (Xmin, xres, xrotation, Ymax, yrotation, -yres)

        GTiff_destination_file = os.path.join(GTiff_files_path, os.path.splitext(os.path.split(file)[1])[0] + '_relocated.tif')

        output_raster = gdal.GetDriverByName('GTiff').Create(GTiff_destination_file, ncols, nrows, 1 ,gdal.GDT_Float32)  # Open the file
        output_raster.SetGeoTransform(geotransform)  # Specify its coordinates
        srs = osr.SpatialReference()                 # Establish its coordinate encoding
        srs.ImportFromEPSG(CRS_code)                     # This one specifies WGS84 lat long.
                                                     # Anyone know how to specify the
                                                     # IAU2000:49900 Mars encoding?
        output_raster.SetProjection( srs.ExportToWkt() )   # Exports the coordinate system
                                                           # to the file
        output_raster.GetRasterBand(1).WriteArray(array)   # Writes my array to the raster

        output_raster.FlushCache()

    return

def raster_check_projection(dst_fp, ref_fp = None, optional_crs = None):

    """Checks the projection of a given raster with respect to a reference raster or reference CRS.

        Parameters
        ----------
        dst_fp : str
            Path of the source file.
        ref_fp : str
            Path of the reference file.
        optional_crs : str
            String indicating the CRS of reference. Ex: '4326' indicating EPSG:4326


        Returns
        -------
        out : bool
            Boolean indicating if projection is different (True) or not (False)
    """

    message = '\n--> Files have same projection (no reprojection needed)'
    reproj = False

    if ref_fp is not None:

        m = ' and (reference file) ' + os.path.splitext(os.path.split(ref_fp)[1])[0]

    elif optional_crs is not None:

        m = ' and (reference CRS) ' + optional_crs

    print('\nChecking projection consistency between: (source file) ' + os.path.splitext(os.path.split(dst_fp)[1])[0] + m)

    with rasterio.open(dst_fp) as src:

        if ref_fp is not None:

            with rasterio.open(ref_fp) as data:

                if data.crs != src.crs:

                    reproj = True
                    message = '\n--> Files have different projection:\n' + '\n  Reference CRS: ' + str(data.crs) + '\n  Source CRS: ' + str(src.crs)

        elif (ref_fp is None) and (optional_crs is not None):

            if optional_crs != src.crs:

                reproj = True
                message = '\n--> File has different projection to the one specified:\n' + '\n   Reference CRS: ' + optional_crs + '\n   Source CRS: ' + str(src.crs)

        else:

            raise AttributeError('Indicate either a reference file or optional crs to compare projection')

    print(message)

    return reproj

def raster_reproject(dst_fp, reproj_fp = None, ref_fp = None, optional_crs = None, search_criteria = '*.tif'):

    """
    Based on:
    https://rasterio.readthedocs.io/en/latest/topics/reproject.html

        Reprojects a given raster or set of rasters with respect to a reference raster or reference CRS.

        Parameters
        ----------
        dst_fp : str
            Path of the source file.
        reproj_fp : str
            Path where the reprojected file will be located.
        ref_fp : str
            Path of the reference file.
        optional_crs : str
            String indicating the CRS of reference. Ex: '4326' indicating EPSG:4326
        search_criteria : list or str
            List of criteria for searhing the ascii files to be converted.


        Returns
        -------
        out : file
            Reprojected file
    """

    if ref_fp is not None:

        with rasterio.open(ref_fp) as data:

            ref_crs = data.crs

    elif (ref_fp is None) and (optional_crs is not None):

        ref_crs = optional_crs

    else:

        raise AttributeError('Indicate either a reference file or optional crs reproject')

    if os.path.splitext(os.path.split(dst_fp)[1])[1] == '':

        if not isinstance(search_criteria, list):

            search_criteria = [search_criteria]

        path = []

        for hints in search_criteria:

            # Make a search criteria to select the DEM files
            q = os.path.join(dst_fp, hints)
            #print(q)

            # glob function can be used to list files from a directory with specific criteria
            path += glob.glob(q)

    else:

        path = glob.glob(dst_fp)

    for input_file in path:

        with rasterio.open(input_file) as src:

            if ref_crs != src.crs:

                transform, width, height = calculate_default_transform(
                    src.crs, ref_crs, src.width, src.height, *src.bounds)
                kwargs = src.meta.copy()
                kwargs.update({
                    'crs': ref_crs,
                    'transform': transform,
                    'width': width,
                    'height': height
                })

                if reproj_fp is not None:

                    if os.path.splitext(os.path.split(reproj_fp)[1])[1] == '':

                        output_fp = os.path.join(reproj_fp, os.path.splitext(os.path.split(input_file)[1])[0] + '_reproj.tif')

                    else:

                        output_fp = reproj_fp

                else:

                    output_fp = os.path.join(os.path.split(input_file)[0], os.path.splitext(os.path.split(input_file)[1])[0] + '_reproj.tif')

                with rasterio.open(output_fp, 'w', **kwargs) as dst:
                    for i in range(1, src.count + 1):
                        reproject(
                            source=rasterio.band(src, i),
                            destination=rasterio.band(dst, i),
                            src_transform=src.transform,
                            src_crs=src.crs,
                            dst_transform=transform,
                            dst_crs=ref_crs,
                            resampling=Resampling.nearest)

    return

def raster_merge(rasters_folder_path, merged_file_path, search_criteria = "L*.tif"):

    """
    Based on:
    https://autogis-site.readthedocs.io/en/latest/notebooks/Raster/raster-mosaic.html

        Merges a set of rasters into one raster.

        Parameters
        ----------
        rasters_folder_path : str
            Path of raster files to merge.
        reproj_fp : str
            Path where the merged file will be located.
        search_criteria : list or str
            List of criteria for searhing the ascii files to be converted.


        Returns
        -------
        out : file
            Merged file
    """

    # File and folder paths
    out_fp = os.path.join(merged_file_path)

    if not isinstance(search_criteria, list):

        search_criteria = [search_criteria]

    dem_fps = []

    for hints in search_criteria:

        # Make a search criteria to select the DEM files
        q = os.path.join(rasters_folder_path, hints)
        #print(q)

        # glob function can be used to list files from a directory with specific criteria
        dem_fps += glob.glob(q)

    # Files that were found:
    #print(dem_fps)

    # List for the source files
    src_files_to_mosaic = []

    # Iterate over raster files and add them to source -list in 'read mode'
    for fp in sorted(dem_fps):
        src = rasterio.open(fp)
        src_files_to_mosaic.append(src)

    #print(src_files_to_mosaic)

    """
    # Create 4 plots next to each other
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(ncols=4, nrows=1, figsize=(12, 4))

    # Plot first four files
    show(src_files_to_mosaic[0], ax=ax1)
    show(src_files_to_mosaic[1], ax=ax2)
    show(src_files_to_mosaic[2], ax=ax3)
    show(src_files_to_mosaic[3], ax=ax4)

    # Do not show y-ticks values in last three axis
    for ax in [ax2, ax3, ax4]:
        ax.yaxis.set_visible(False)
    """
    # Merge function returns a single mosaic array and the transformation info
    mosaic, out_trans = merge(src_files_to_mosaic)

    # Plot the result
    #show(mosaic, cmap='terrain')

    # Copy the metadata
    out_meta = src.meta.copy()

    #print(src.meta['crs'])

    # Update the metadata
    out_meta.update({"driver": "GTiff",
                     "height": mosaic.shape[1],
                     "width": mosaic.shape[2],
                     "transform": out_trans,
                     #"crs": "+proj=utm +zone=35 +ellps=GRS80 +units=m +no_defs "
                     "crs": src.meta['crs']
                     }
                    )

    # Write the mosaic raster to disk
    with rasterio.open(out_fp, "w", **out_meta) as dest:
        dest.write(mosaic)

    return

def getFeatures(gdf):

    """
    Based on:
    https://autogis-site.readthedocs.io/en/latest/notebooks/Raster/clipping-raster.html?highlight=crop

    Function to parse features from GeoDataFrame in such a manner that rasterio wants them
    """

    import json
    return [json.loads(gdf.to_json())['features'][0]['geometry']]

def raster_crop(raster_to_crop_path, cropped_file_path = None, search_criteria = "*.tif",
                vector_file = None, polygon_coords = None, polygon_crs = None, mask_value = -9999.,
                mask_all_touched = False):

    """
    Based on:
    https://autogis-site.readthedocs.io/en/latest/notebooks/Raster/clipping-raster.html?highlight=crop

        Crops a raster or set of rasters according to shapefile or Polygon.

        Parameters
        ----------
        raster_to_crop_path : str
            Path of raster file(s) to crop.
        cropped_file_path : str
            Path where the cropped file(s) will be located.
        search_criteria : list or str
            List of criteria for searhing the ascii files to be converted.
        vector_file : str
            Path of the vector file.
        polygon_coords : arr
            Array of tuples indicating the coordinates or points of the polygon
        polygon_crs : int
            Reference system for the polygon coordinates
        mask_value : float
            Value for pixels falling out of the polygon or shapefile


        Returns
        -------
        out : file
            Cropped file(s)

        Example
        -------
        polygon_coords = [(x0, y0), (x1, y1), ..., (xn, yn), (x0, y0)]
        crs reference system of the form "EPSG:XXXX"
        polygon_crs = XXXX
    """

    if vector_file is not None:

        geo = gpd.read_file(vector_file)

    elif (polygon_coords is not None) and (polygon_crs is not None):

        geo = gpd.GeoDataFrame({'geometry': Polygon(polygon_coords)}, index=[0], crs = from_epsg(polygon_crs))

    else:

        raise AttributeError('Vector or Polygon needed for cropping')

    if os.path.splitext(os.path.split(raster_to_crop_path)[1])[1] == '':

        if not isinstance(search_criteria, list):

            search_criteria = [search_criteria]

        path = []

        for hints in search_criteria:

            # Make a search criteria to select the DEM files
            q = os.path.join(raster_to_crop_path, hints)
            #print(q)

            # glob function can be used to list files from a directory with specific criteria
            path += glob.glob(q)

    else:

        path = glob.glob(raster_to_crop_path)

    for input_file in path:

        # Read the data
        data = rasterio.open(input_file)

        # Project the Polygon into same CRS as the grid
        geo = geo.to_crs(crs = data.crs.data)

        # Print crs
        #print(geo.crs)

        coords = getFeatures(geo)
        #print(coords)

        # Clip the raster with Polygon
        out_img, out_transform = mask(dataset = data, shapes = coords, crop = True, all_touched = mask_all_touched, nodata = mask_value)

        # Copy the metadata
        out_meta = data.meta.copy()
        #print(out_meta)

        # Parse EPSG code
        #print(data.crs.data)
        epsg_code = int(data.crs.data['init'][5:])
        #print(epsg_code)

        out_meta.update({"driver": "GTiff",
                         "height": out_img.shape[1],
                         "width": out_img.shape[2],
                         "transform": out_transform,
                         "crs": data.meta['crs']}
                                 )

        if cropped_file_path is not None:

            if os.path.splitext(os.path.split(cropped_file_path)[1])[1] == '':

                output_fp = os.path.join(cropped_file_path, os.path.splitext(os.path.split(input_file)[1])[0] + '_cropped.tif')

            else:

                output_fp = cropped_file_path

        else:

            output_fp = os.path.join(os.path.split(input_file)[0], os.path.splitext(os.path.split(input_file)[1])[0] + '_cropped.tif')

        with rasterio.open(output_fp, "w", **out_meta) as dest:
                dest.write(out_img)

    return

def vector_reproject(input_vector_file, reference_vector_file, output_file):

    """Reprojects a given vector or set of vetor with respect to a reference vector file.

        Parameters
        ----------
        input_vector_file : str
            Path of the vector file to reproject.
        reference_vector_file : str
            Path of the reference vector file.
        output_file : str
            Path of the new vector file reprojected.


        Returns
        -------
        out : file
            Reprojected file
    """

    geo1 = gpd.read_file(input_vector_file)
    geo2 = gpd.read_file(reference_vector_file)

    geo1_reproj  = geo1.to_crs(geo2.crs)

    geo1_reproj.to_file(output_file)

    return

def vector_intersection(vector_file_1, vector_file_2, output_file):

    """
    Based on:
    https://autogis-site.readthedocs.io/en/latest/notebooks/L4/geometric-operations.html

        Intersects two vector files for a desired shape of the region.

        Parameters
        ----------
        vector_file_1 : str
            Path of the vector file 1.
        vector_file_2 : str
            Path of the vector file 2.
        output_file : str
            Path of the new vector file reprojected.


        Returns
        -------
        out : file
            Intersected file
    """

    geo1 = gpd.read_file(vector_file_1)
    geo2 = gpd.read_file(vector_file_2)

    try:

        assert geo2.crs == geo1.crs

    except AssertionError:

        reproj_file = os.path.join(os.path.split(vector_file_2)[0], os.path.splitext(os.path.split(vector_file_2)[1])[0] + '_reproj.shp')

        vector_reproject(vector_file_2, vector_file_1, reproj_file)

        geo2 = gpd.read_file(reproj_file)

    intersection = gpd.overlay(geo1, geo2, how = 'intersection')

    intersection.to_file(output_file)

    return

def set_raster_resolution(input_file, reference_file, output_file = None, binary = False, mask_value = None):

    """Checks the size and resolution of a given raster with respect to a reference raster.

        Parameters
        ----------
        input_file : str
            Path of the source file.
        reference_file : str
            Path of the reference file.
        binary : bool
            Boolean indicating if data in the input file are binary or continuous
        mask_value : float
            Float indicating if there is a value that should not be taking into account when averaging or expanding


        Returns
        -------
        out : file
            File with same path as input file but resolution changed according to reference file
    """

    dataset = gdal.Open(reference_file, gdal.GA_ReadOnly) # Note GetRasterBand() takes band no. starting from 1 not 0
    proj = osr.SpatialReference(wkt = dataset.GetProjection())
    ref_geotransform = dataset.GetGeoTransform()
    band = dataset.GetRasterBand(1)
    ref_array = band.ReadAsArray()

    CRS_code = int(proj.GetAttrValue('AUTHORITY',1))

    ref_shape = np.shape(ref_array)

    ref_rows, ref_cols = ref_shape

    dataset = gdal.Open(input_file, gdal.GA_ReadOnly) # Note GetRasterBand() takes band no. starting from 1 not 0
    proj = osr.SpatialReference(wkt = dataset.GetProjection())
    geotransform = dataset.GetGeoTransform()
    band = dataset.GetRasterBand(1)
    in_array = band.ReadAsArray()

    in_shape = np.shape(in_array)

    in_rows, in_cols = in_shape

    new_array = np.empty((in_rows, in_cols))

    if in_shape == ref_shape:

        return

    if in_rows > ref_rows:

        ratio = floor(in_rows / ref_rows)

        j = 0

        if binary:

            for i in range(ref_rows):

                data_rows = in_array[j:j + ratio]

                if mask_value is not None:

                    data_rows = np.where(data_rows != mask_value, data_rows, 0)

                for k in range(np.shape(data_rows)[1]):

                    new_array[i, k] = np.bincount(data_cols[:, k]).argmax()

                j = i * ratio

        else:

            for i in range(ref_rows):

                data_rows = in_array[j:j + ratio]

                if mask_value is not None:

                    data_rows = np.where(data_rows != mask_value, data_rows, 0)

                new_array[i] = np.mean(data_rows, axis = 0)

                j = i * ratio

    elif in_rows < ref_rows:

        ratio = floor(ref_rows / in_rows)

        j = 0

        for i in range(ref_rows):

            new_array[j:j + ratio] = in_array[i]

            j = i * ratio

    if in_cols > ref_cols:

        ratio = floor(in_cols / ref_cols)

        j = 0

        if binary:

            for i in range(ref_cols):

                data_cols = in_array[:, j:j + ratio]

                if mask_value is not None:

                    data_cols = np.where(data_cols != mask_value, data_cols, 0)

                for k in range(np.shape(data_cols)[0]):

                    new_array[k, i] = np.bincount(data_cols[k, :]).argmax()

                j = i * ratio

        else:

            for i in range(ref_cols):

                data_cols = in_array[:, j:j + ratio]

                if mask_value is not None:

                    data_cols = np.where(data_cols != mask_value, data_cols, 0)

                new_array[:, i] = np.mean(data_cols, axis = 1)

                j = i * ratio

    elif in_cols < ref_cols:

        ratio = floor(ref_cols / in_cols)

        j = 0

        for i in range(ref_cols):

            at = np.transpose(in_array)
            bt = np.transpose(new_array)

            #new_array[:, j:j + ratio] = in_array[:, i]
            bt[j:j + ratio] = at[i]

            new_array = np.transpose(bt)

            j = i * ratio

    new_array_2 = np.empty(ref_shape)
    new_array_2 = new_array[0:ref_rows, 0:ref_cols]

    if output_file is None:

        output_file = input_file

    output_raster = gdal.GetDriverByName('GTiff').Create(output_file, ref_cols, ref_rows, 1 ,gdal.GDT_Float32)  # Open the file
    output_raster.SetGeoTransform(ref_geotransform)  # Specify its coordinates
    srs = osr.SpatialReference()                 # Establish its coordinate encoding
    srs.ImportFromEPSG(CRS_code)                     # This one specifies WGS84 lat long.
                                                 # Anyone know how to specify the
                                                 # IAU2000:49900 Mars encoding?
    output_raster.SetProjection( srs.ExportToWkt() )   # Exports the coordinate system
                                                       # to the file
    output_raster.GetRasterBand(1).WriteArray(new_array_2)   # Writes my array to the raster

    output_raster.FlushCache()

    return
