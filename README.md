# Sentinel-1 incidence angle normalization

Sentinel-1 data acquisition in the mode of Interferometric Wide Swath(IW) and Extra Wide Swath (EW) acquire data over wide areas which results in progressive reduction in brightness from near to far range. The backscatter coefficient values depend to a great extent on the incident angle. This can be problematic and affect the detection and classification of sea surface features. Therefore, incidence angle normalization is required to reduce the variation of backscatter energy over the SAR scene

## How to perform incidence angle correction
- **Provide calibrated Sentinel-1 data**

  The **scr/angle_correction.py** script requires a Sentinel-1 product that contains 3 bands
  - Co-polarized band (calibrated in db units)
  - Cross-polarized band (calibrated in db units)
  - Incidence angle

  If the Sentinel-1 data is not in the format described above, the user can run the bash script `S1_preprocessing.sh`. It takes the raw GRD Sentinel-1 data and follows a series of steps (based on the graphs created in SNAP) to pre-process the data. It outputs the pre-processed SAR image which contains the 3 bands mentioned above. There are two graphs in the **graphs** folder. The one named named **s1_preprocessing_land.xml** is used to pre-process Sentinel-1 data over land. The one named named **s1_preprocessing_ocean.xml** is used to pre-process Sentinel-1 data over ocean
  
- **Collect points over Sentinel-1 swath**

  Digitize points (e.g QGIS or ArcGIS) over the SAR scene from near to far range. Ideally, the points should be closely spaced. The more points are collected (very close to each other) the more accurate the angle correction will be. Ideally, the points should be digitized over the same landcover type to avoid any bias

- **Run angle_correction script**

Examples on how to run the script are given within the script

```
usage: angle_correction.py [-h] 
                           [-o OUTDIR] 
                           [-i LOAD_SAR] 
                           [--vector VECTOR] 
                           [--ref_angle REF_ANGLE]

optional arguments:
  -h, --help            show this help message and exit
  -o OUTDIR, --outdir OUTDIR
                        Specify an output directory
  -i LOAD_SAR, --load_sar LOAD_SAR
                        Provide a path to the sar data
  --vector VECTOR       Provide a path to the vector file (shapefile)
  --ref_angle REF_ANGLE
                        Provide a value for incidence angle normalziation
```
