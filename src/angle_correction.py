'''
This script is used to perform incidence angle correction on Sentinel-1 data

usage: angle_correction.py [-h] 
                           [-o OUTDIR] 
                           [-i LOAD_SAR] 
                           [--ref_angle REF_ANGLE]
                           [--linear LINEAR] 
                           [--sqr_cosine]

optional arguments:
  -h, --help            show this help message and exit
  -o OUTDIR, --outdir OUTDIR
                        Specify an output directory
  -i LOAD_SAR, --load_sar LOAD_SAR
                        Provide a path to the sar data
  --ref_angle REF_ANGLE
                        Provide a value for incidence angle normalziation
  --linear LINEAR       Performs angle correction based on linear regression. It requires a shapefile (points) as an input. 
                        Point data should be colected from near to far range
  --sqr_cosine          Performs angle correction based on the square cosine
                        
                        
Example:

Incidence angel correction using linear regression with a default value (33)
-----------------------------------------------------------------------
python3 angle_correction.py -i raw/subset_S1A_IW_GRDH_20210105_ice.tif 
                            -o output/sar_normalized_sandard.tif  
                            --linear vector/water_subset.shp
                            
                            
Incidence angle correction using linear regression at 40 degrees 
------------------------------------------------------------
python3 angle_correction.py -i raw/subset_S1A_IW_GRDH_20210105_ice.tif 
                            -o output/sar_normalized_sandard.tif  
                            --ref_angle 40
                            --linear vector/water_subset.shp
                            
Incidence angle correction using square cosine at 40 degrees 
------------------------------------------------------------
python3 angle_correction.py -i raw/subset_S1A_IW_GRDH_20210105_ice.tif 
                            -o output/sar_normalized_sandard.tif  
                            --sqr_cosine
                            --ref_angle 40
         
                            
'''


import argparse
from pathlib import Path
import numpy as np
from scipy.stats import linregress

import geopandas as gpd
import pandas as pd

from rasterstats import zonal_stats
import rasterio



def write(raster_input, profile, raster_output):
    profile.update(
        dtype=raster_input.dtype,
        count=raster_input.shape[0],
        nodata=0,
        compress='lzw')

    with rasterio.open(raster_output, 'w', **profile) as out:
        out.write(raster_input)

        
def read(path):
    with rasterio.open(path) as src:
        rasters = src.read()
        sigma_vv = rasters[0]
        sigma_vh = rasters[1]
        angle = rasters[2]
        profile = src.profile
        if profile['count'] != 3:
            print (f"It seems that the number of bands in the data is {profile['count']}" '\n' 
            "The SAR product must contain 3 bands in the following order. 1) co-polarized band, 2)cross-polarized band, 3)incidence angle," '\n'
            "Bands in the SAR product not in that order might produces unexpected resuls"
    )
    return sigma_vv, sigma_vh, angle, profile


def regress(raster_extract,angle_extract):
    poly_fit = np.polyfit(angle_extract, raster_extract, 1)
    poly_val = np.polyval(poly_fit, angle_extract)
    slope, intercept, r_value, p_value, std_err = linregress(angle_extract, raster_extract)
    # print ('slope:', np.round(slope,2))
    # print ('intercept:', np.round(intercept,2))
    # print ("r-squared:", np.round(r_value**2,4))
    return poly_val, slope


def sigma_stats(gdf, raster, profile):
    sigma_allvals = []
    dict_stats_sigma = {}
    for indx, row in gdf.iterrows():
        geom = row['geometry'] 
        idd = row['id']
        v = pd.DataFrame({'id':idd,'geometry':geom}, index=[0])

        stats = zonal_stats(v['geometry'],
                            raster,
                            affine=profile['transform'],
                            categorical=True,
                            all_touched=True,
                            nodata=-999,
                       )
        dict_stats_sigma.setdefault(indx,[]).append(stats)
    
    # convert dictionary into a list
    for k,v in dict_stats_sigma.items():
        value_sigma = v[0][0]
        for sig in value_sigma:
            sigma_allvals.append(sig)     

    return sigma_allvals 
    

def angle_stats(gdf, raster, profile):  
    angle_allvals = []
    dict_stats_angle = {}
    for indx, row in gdf.iterrows():
        geom = row['geometry'] 
        # get the first column which is the id
        cols = gdf.columns
        idd = cols[0]
        v = pd.DataFrame({'id':idd,'geometry':geom}, index=[0])

        stats = zonal_stats(v['geometry'],
                            raster,
                            affine=profile['transform'],
                            categorical=True,
                            all_touched=True,
                            nodata=-999,
                       )
        dict_stats_angle.setdefault(indx,[]).append(stats)
        
    # convert dictionary into a list
    for k,v in dict_stats_angle.items():
        value_angle = v[0][0]
        for a in value_angle:
            angle_allvals.append(a)
        
    return angle_allvals 


def linear_angle_correction(gdf, co_db, angle, profile, ref_angle=33):
    stats_sigma = sigma_stats(gdf, co_db, profile)
    stats_angle = angle_stats(gdf, angle, profile)
    polyval, slope = regress(stats_sigma, stats_angle)
    
    co_db_corrected = co_db - slope * (angle - ref_angle)
    return co_db_corrected


def square_cosine_angle_correction(sigma_vv, angle, ref_angle=33):

    sigma_nought_vv = sigma_vv * np.cos(np.deg2rad(angle))**2
    sigma_corr_vv = (sigma_nought_vv * np.cos(np.deg2rad(ref_angle))**2) / np.cos(np.deg2rad(angle))**2
    return sigma_corr_vv

def stack(coPol, crossPol, angle):
    stacked = np.stack([coPol, crossPol, angle])
    return stacked


if __name__ == "__main__":
    print('Incidence angle normalization...', '\n',)
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument(
        '-o',
        '--outdir',
        required=False,
        help='Specify an output directory'
    )

    parser.add_argument(
        '-i',
        '--load_sar',
        type=str,
        required=False,
        help='Provide a path to the sar data'
    )
    

    parser.add_argument(
        '--ref_angle',
        type=int,
        required=False,
        help='Provide a value to be used for incidence angle normalziation'
    )
    
    parser.add_argument(
        '--linear',
        type=str,
        required=False,
        help='Performs angle correction based on linear regression. \
            It requires a shapefile (points) as an input. Point data should be colected \
            from near to far range'
    )
    
    parser.add_argument(
        '--sqr_cosine',
        action='store_true',
        required=False,
        help='Performs angle correction based on the square cosine'
    )
    
    args = parser.parse_args()
    
    
    if not args.outdir:
        parser.error('Please provide an output path (use option -o)')
        
    # Load SAR data
    coPol_db, crossPol_db, angle, profile = read(args.load_sar)
    
    if args.linear:
        # load shapefile
        gdf = gpd.read_file(args.linear)
        
        if args.ref_angle:
            # perform angle correction specifying the angle
            coPol_corrected = linear_angle_correction(gdf, coPol_db, angle, profile, args.ref_angle)
            crossPol_corrected = linear_angle_correction(gdf, crossPol_db, angle, profile, args.ref_angle)
            stacked = stack(coPol_corrected, crossPol_corrected, angle)
            
        elif not args.ref_angle:        
            # perform angle correction using the defoult value of 33 degrees
            coPol_corrected = linear_angle_correction(gdf, coPol_db, angle, profile)
            crossPol_corrected = linear_angle_correction(gdf, crossPol_db, angle, profile)
            stacked = stack(coPol_corrected, crossPol_corrected, angle)
            
    if args.sqr_cosine:
        if args.ref_angle:
            # perform angle correction specifying the angle
            coPol_corrected = square_cosine_angle_correction(coPol_db, angle, args.ref_angle)
            crossPol_corrected = square_cosine_angle_correction(crossPol_db, angle, args.ref_angle)
            stacked = stack(coPol_corrected, crossPol_corrected, angle)
            
        elif not args.ref_angle:
            # perform angle correction specifying the angle
            coPol_corrected = square_cosine_angle_correction(coPol_db, angle)
            crossPol_corrected = square_cosine_angle_correction(crossPol_db, angle)
            stacked = stack(coPol_corrected, crossPol_corrected, angle)
            
    
    # write 
    write(stacked, profile, args.outdir)
    print ('DONE')