import os, glob, json
import concurrent.futures
from math import floor
from pathlib import Path
import fnmatch
import numpy as np
import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.mask import mask
from rasterio.merge import merge
from shapely.geometry import Polygon
import geopandas as gpd
from osgeo import gdal, osr
from modules.common import save_GTiff_raster

def _ensure_list(list_or_any):
    """Check that list_or_any is a list, if not then put it into one.
    For example: 'string' -> ['string']"""
    if not isinstance(list_or_any, list):
        return [list_or_any]
    return list_or_any

def _glob_path(search_path, search_criteria):
    """Find files that match search criteria within the folder

    Args:
        search_path (str or pathlib.Path): Path to-be-searched
        search_criteria (str or list): Match criteria

    Returns:
        list: list of paths matching the search criteria pattern(s)
    """
    if Path(search_path).suffix == '':
        search_criteria = _ensure_list(search_criteria)
        path = []
        for hints in search_criteria:
            # Make a search criteria to select the DEM files
            q = os.path.join(search_path, hints)
            # glob function can be used to list files from a directory with specific criteria
            path.extend(glob.glob(q))
    else:
        path = glob.glob(search_path)

    return path

def ascii_to_geotiff(GTiff_files_path, config):
    """A given ascii file(s) is(are) transformed into GeoTiff file.

        Parameters
        ----------
        config : dictionary containing rain section from the config file

        Returns
        -------
        out : file
            GTiff converted file
    """
    xmax = config['x_max']
    xres = config['x_res']
    if xres is None and xmax is None:
        raise ValueError("If 'xres' not given, 'xmax' has to be given")

    yres=config['y_res']
    ymin=config['y_min']
    if yres is None and ymin is None:
        raise ValueError("If 'yres' not given, 'ymin' has to be given")

    srs = osr.SpatialReference()
    CRS_code = config['EPSG_code'] or 3067
    srs.ImportFromEPSG(CRS_code)
    projection_wkt = srs.ExportToWkt()
    
    search_criteria = config['ascii_search_criteria'] or '*.txt'
    path = _glob_path(config['ascii_files_path'], search_criteria)

    # Figure out geotransform
    array = np.loadtxt(path[0])
    nrows, ncols = np.shape(array)

    xres = xres or ((xmax - config['x_min']) / float(ncols))
    yres = yres or ((config['y_max'] - ymin) / float(nrows))
    xrotation = config['x_rotation'] or 0
    yrotation = config['y_rotation'] or 0
    geotransform = (config['x_min'], xres, xrotation, config['y_max'], yrotation, -1 * yres)

    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = []
        for filename in path:
            p = executor.submit(_process_ascii, filename, GTiff_files_path, geotransform, projection_wkt)
            futures.append(p)

def _process_ascii(fname, output_path, geotransform, projection_wkt):
    """Convert ascii file to geotiff
    Args:
        fname (file, str, pathlib.Path): ASCII file
        output_path (str or pathlib.Path): Output folder
        geotransform (tuple): GDAL GeoTransform tuple
        projection_wkt (str): Projection string
    """
    array = np.loadtxt(fname)
    nrows, ncols = array.shape

    destination_file = os.path.join(output_path, Path(fname).with_suffix('.tif').name)
    # Tätä ei voinut koodata uudelleen siten että käytti alla olevalle apufunktiota. 
    # Syy on varmaan että funktio ei suoritu heti vaan menee suoritus pinoon tjsp.
    driver = gdal.GetDriverByName('GTiff')
    output_raster = driver.Create(destination_file, ncols, nrows, 1, gdal.GDT_Float32)
    output_raster.SetGeoTransform(geotransform)
    output_raster.SetProjection(projection_wkt)
    output_raster.GetRasterBand(1).WriteArray(array)

    output_raster.FlushCache()
    output_raster = None

def rain_relocation(GTiff_files_path, config):
    """Relocates a desired pixel region of the image (X_radar_rain, Y_radar_rain) to a desired
    region (X_target, Y_target). Produces a new GTiff file.

        Parameters
        ----------
        GTiff_files_path : str
            Path where GTiff files will be located.
        config : dictionary of the rain parameters from the config file.

        Returns
        -------
        out : file
            GTiff relocated file
    """
    CRS_code = config['EPSG_code'] or 3067
    xrotation = config['x_rotation'] or 0
    yrotation = config['y_rotation'] or 0
    
    xres=config['x_res']
    xmax=config['x_max']
    if xres is None and xmax is None:
        raise ValueError('If xres not given, xmax has to be given')
    yres=config['y_res']
    ymin=config['y_min']
    if yres is None and ymin is None:
        raise ValueError('If yres not given, ymin has to be given')

    # GDAL setup
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(CRS_code)
    projection_wkt = srs.ExportToWkt()

    search_criteria = config['relocation_search_criteria'] if config['relocation_search_criteria'] else '*.txt'
    path = _glob_path(GTiff_files_path, search_criteria)
    
    for filename in path:
        # This prevents extra relocated files.
        if fnmatch.fnmatch(filename, '*relocated*'):
            print(filename + " skipped as it was already relocated according its filename")
            continue
        # Note GetRasterBand() takes band no. starting from 1 not 0
        dataset = gdal.Open(filename, gdal.GA_ReadOnly)
        # proj = osr.SpatialReference(wkt = dataset.GetProjection())
        if not dataset:
            print(filename)
            dataset = None
            continue
    
        geotransform = dataset.GetGeoTransform()
        band = dataset.GetRasterBand(1)
        array = band.ReadAsArray()
        # Close opened file
        band = None
        dataset = None

        nrows, ncols = np.shape(array)
        xres = xres or ((xmax - config['x_min']) / float(ncols))
        yres = yres or ((config['y_max'] - ymin) / float(nrows))
        Xmin = config['X_target'] - (ncols / 2 * xres + config['X_radar_rain'])
        Ymax = config['Y_target'] + (nrows / 2 * yres - config['Y_radar_rain'])
        geotransform = (Xmin, xres, xrotation, Ymax, yrotation, -1*yres)
        GTiff_destination_file = os.path.join(GTiff_files_path,f'{Path(filename).stem}_relocated.tif')
        save_GTiff_raster(projection_wkt, geotransform, array, GTiff_destination_file, driver=gdal.GetDriverByName('GTiff'))

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

    if ref_fp:
        m = f' and (reference file) {Path(ref_fp).name}'
    elif optional_crs:
        m = f' and (reference CRS) {optional_crs}'
    else: # ref_fp is None and optional_crs is None
        raise ValueError('Indicate either a reference file or optional crs to compare projection')

    print(f'\nChecking projection consistency between: (source file) {Path(dst_fp).name} {m}')

    with rasterio.open(dst_fp) as src:
        if ref_fp:
            with rasterio.open(ref_fp) as data:
                if data.crs != src.crs:
                    reproj = True
                    message = (f'\n--> Files have different projection:\n\n  '
                               f'Reference CRS: {data.crs}\n  Source CRS: {src.crs}')

        elif optional_crs:
            if optional_crs != src.crs:
                reproj = True
                message = (f'\n--> File has different projection to the one specified:\n\n   '
                           f'Reference CRS: {optional_crs}\n   Source CRS: {src.crs}')

    print(message)

    return reproj

def raster_reproject(dst_fp, reproj_fp=None, ref_fp=None, optional_crs=None,
                     search_criteria='*.tif'):
    """ Reprojects a given raster or set of rasters with respect to a reference raster or
    reference CRS.

    Based on:
    https://rasterio.readthedocs.io/en/latest/topics/reproject.html

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
    if ref_fp:
        with rasterio.open(ref_fp) as data:
            ref_crs = data.crs
    elif optional_crs:
        ref_crs = optional_crs
    else:
        raise ValueError('Indicate either a reference file or optional crs reproject')

    path = _glob_path(dst_fp, search_criteria)

    for input_file in path:
        input_file = Path(input_file)

        if reproj_fp is None:
            output_fp = input_file.with_name(f'{input_file.stem}_reproj.tif')
        elif Path(reproj_fp).suffix == '':
            output_fp = os.path.join(reproj_fp, f'{input_file.stem}_reproj.tif')
        else:
            output_fp = reproj_fp

        with rasterio.open(input_file) as src:
            if ref_crs == src.crs:
                continue

            transform, width, height = calculate_default_transform(
                src.crs, ref_crs, src.width, src.height, *src.bounds
            )

            kwargs = src.meta.copy()
            kwargs.update({
                'crs': ref_crs,
                'transform': transform,
                'width': width,
                'height': height
            })

            with rasterio.open(output_fp, 'w', **kwargs) as dst:
                for i in range(1, src.count + 1):
                    reproject(
                        source=rasterio.band(src, i),
                        destination=rasterio.band(dst, i),
                        src_transform=src.transform,
                        src_crs=src.crs,
                        dst_transform=transform,
                        dst_crs=ref_crs,
                        resampling=Resampling.nearest
                    )

def raster_merge(rasters_folder_path, merge_fp, search_criteria = "L*.tif"):
    """ Merges a set of rasters into one raster.

    Based on:
    https://autogis-site.readthedocs.io/en/latest/notebooks/Raster/raster-mosaic.html

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
    out_fp = os.path.join(merge_fp)

    # Remove any suffix to replicate pre-refactoring behaviour here
    path_no_suffix = Path(rasters_folder_path).with_suffix('')
    dem_fps = _glob_path(path_no_suffix, search_criteria)

    # List for the source files
    src_files_to_mosaic = []

    # Iterate over raster files and add them to source -list in 'read mode'
    for fp in sorted(dem_fps):
        src = rasterio.open(fp)
        src_files_to_mosaic.append(src)

    # Merge function returns a single mosaic array and the transformation info
    mosaic, out_trans = merge(src_files_to_mosaic)

    # Copy the metadata
    out_meta = src.meta.copy()

    # Update the metadata
    out_meta.update({
        "driver": "GTiff",
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": out_trans,
        "crs": src.meta['crs']
    })

    # Write the mosaic raster to disk
    with rasterio.open(out_fp, "w", **out_meta) as dest:
        dest.write(mosaic)

def getFeatures(gdf):
    """Function to parse features from GeoDataFrame in such a manner that rasterio wants them

    Based on: https://autogis-site.readthedocs.io/en/latest/notebooks/Raster/clipping-raster.html
    """
    return [json.loads(gdf.to_json())['features'][0]['geometry']]

def raster_crop(raster_to_crop_path, config, cropped_fp = None, search_criteria = "*.tif"):
    """Crop a raster or set of rasters according to shapefile or Polygon.

    Based on:
    https://autogis-site.readthedocs.io/en/latest/notebooks/Raster/clipping-raster.html?highlight=crop

        Parameters
        ----------
        raster_to_crop_path : str
            Path of raster file(s) to crop.
        config: dict
            Dictionary of the values related to the cropping section of the config file.
        cropped_fp : str
            Path where the cropped file(s) will be located.
        search_criteria : list or str
            List of criteria for searhing the ascii files to be converted.

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
    x_coords=config['x_polygon_coords']
    y_coords=config['y_polygon_coords']
    polygon_crs=config['polygon_crs']
    mask_value=config['mask_value'] or -9999.

    if x_coords and y_coords:
        polygon_coords = list(zip(x_coords, y_coords))
    else:
        polygon_coords = None

    if config['shapefile_file_path']:
        geo = gpd.read_file(config['shapefile_file_path'])
    elif polygon_coords and polygon_crs:
        geo = gpd.GeoDataFrame({'geometry': Polygon(polygon_coords)}, index=[0],
                               crs=f"EPSG:{polygon_crs}")
    else:
        raise ValueError('Vector or Polygon needed for cropping')

    path = _glob_path(raster_to_crop_path, search_criteria)

    if cropped_fp:
        cropped_fp = Path(cropped_fp)

    for input_file in path:
        input_file = Path(input_file)

        # Read the data
        data = rasterio.open(input_file)

        # Project the Polygon into same CRS as the grid
        geo = geo.to_crs(crs = data.crs.data)
        coords = getFeatures(geo)

        # Clip the raster with Polygon
        out_img, out_transform = mask(dataset=data, shapes=coords, crop=True,
                                      all_touched=True, nodata=mask_value)

        # Copy the metadata
        out_meta = data.meta.copy()
        out_meta.update({"driver": "GTiff",
                         "height": out_img.shape[1],
                         "width": out_img.shape[2],
                         "transform": out_transform,
                         "crs": data.meta['crs']}
                        )

        if cropped_fp:
            if not cropped_fp.suffix:
                output_fp = cropped_fp / f'{input_file.stem}_cropped.tif'
            else:
                output_fp = cropped_fp
        else:
            output_fp = input_file.with_name(f'{input_file.stem}_cropped.tif')

        with rasterio.open(output_fp, "w", **out_meta) as dest:
            dest.write(out_img)

def vector_reproject(input_vector_file, reference_vector_file, output_file):
    """Reproject a given vector or set of vectors with respect to a reference vector file.

    Parameters
    ----------
    input_vector_file : str | Path
        Path of the vector file to reproject.
    reference_vector_file : str | Path
        Path of the reference vector file.
    output_file : str | Path
        Path of the new vector file reprojected.
    """
    geo1 = gpd.read_file(input_vector_file)
    geo2 = gpd.read_file(reference_vector_file)
    geo1_reproj = geo1.to_crs(geo2.crs)
    geo1_reproj.to_file(output_file)

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

    if geo2.crs != geo1.crs:
        vector_file_2 = Path(vector_file_2)
        reproj_file = vector_file_2.with_name(f"{vector_file_2.stem}_reproj.shp")

        vector_reproject(vector_file_2, vector_file_1, reproj_file)

        geo2 = gpd.read_file(reproj_file)

    intersection = gpd.overlay(geo1, geo2, how = 'intersection')
    intersection.to_file(output_file)

def _load_gdal_raster(filename):
    """[summary]

    Args:
        filename (str): Raster file

    Returns:
        tuple of (numpy.array, dict): raster data and its spatial reference and geotransform
    """
    dataset = gdal.Open(filename, gdal.GA_ReadOnly)
    band = dataset.GetRasterBand(1)
    srs = osr.SpatialReference(wkt=dataset.GetProjection())
    geotransform = dataset.GetGeoTransform()
    data = band.ReadAsArray()

    band = None
    dataset = None

    metadata = {
        'srs': srs,
        'geotransform': geotransform,
    }
    return (data, metadata)

def write_modified_raster_to_file(input_output_file, ref_meta, new_array_2):
    ref_proj = ref_meta['srs']
    CRS_code = int(ref_proj.GetAttrValue('AUTHORITY',1))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(CRS_code)
    save_GTiff_raster(srs.ExportToWkt(), ref_meta['geotransform'], new_array_2, input_output_file)

def process_raster_widths(in_array, ref_cols, new_array):
    _, in_cols = in_array.shape
    if in_cols > ref_cols:
        ratio = floor(in_cols / ref_cols)
        for i in range(ref_cols):
            data_cols = in_array[:, i*ratio:(i+1)*ratio]
            data_cols = np.where(data_cols != 241, data_cols, 0)
            new_array[:, i] = np.mean(data_cols, axis = 1)
    else:
        ratio = floor(ref_cols / in_cols)
        for i in range(ref_cols):
            at = np.transpose(in_array)
            bt = np.transpose(new_array)
            bt[(i*ratio):((i+1)*ratio)] = at[i]
            new_array = np.transpose(bt)
    return new_array

def process_raster_heights(in_array, ref_rows, new_array):
    in_rows, _ = in_array.shape
    if in_rows > ref_rows:
        ratio = floor(in_rows / ref_rows)
        for i in range(ref_rows):
            data_rows = in_array[i*ratio:(i+1)*ratio]
            data_rows = np.where(data_rows != 241, data_rows, 0)
            new_array[i] = np.mean(data_rows, axis = 0)
    else:
        ratio = floor(ref_rows / in_rows)
        for i in range(ref_rows):
            new_array[(i * ratio):(i+1) * ratio] = in_array[i]
    return new_array

def set_raster_resolution(input_output_file, reference_file):
    """Checks the size and resolution of a given raster with respect to a reference raster.

        Parameters
        ----------
        input_file : str
            Path of the source file.
        reference_file : str
            Path of the reference file.

        Returns
        -------
        out : file
            File with same path as input file but resolution changed according to reference file
    """
    in_array, _ = _load_gdal_raster(input_output_file)
    ref_array, ref_meta = _load_gdal_raster(reference_file)

    if in_array.shape == ref_array.shape:
        return

    ref_rows, ref_cols = ref_array.shape
    in_rows, in_cols = in_array.shape

    new_array = np.empty((in_rows, in_cols))
    new_array = process_raster_heights(in_array, ref_rows, new_array)
    new_array = process_raster_widths(in_array, ref_cols, new_array)
    
    write_modified_raster_to_file(input_output_file, 
                                  ref_meta, 
                                  new_array[0:ref_rows, 0:ref_cols])