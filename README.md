# optint_wod_o2
Optimal Interpolation of historical dissolved oxygen data from the World Ocean Database. 
This code set requires python (jupyter notebook) and MATLAB. 
See requirements.txt for the list of required python packages

## Objective: 
- Assess global dissolved oxygen inventory trend

## Procedure:
- Obtain raw data from WOD select
- Organize bottle and CTD O2 data into monthly 1deg x 1deg bin from 0-1000m
- Apply optimal interpolation to determine monthly climatology
- Subtract monthly climatology to determine monthly O2 anomaly
- Aggregate monthly O2 anomaly into 5 year running mean
- Apply optimal interpolation to determine gridded yearly O2 anomaly

### 1. Obtain the ragged netCDF array from WOD Select: 
[WOD Select](https://www.ncei.noaa.gov/access/world-ocean-database-select/dbsearch.html)
- Under "Selection Criteria", choose "Dataset", "Measured Variables" and "Data Exclusion using WOD QC Flags"
- Under Dataset, select "Ocean Station Data" for bottle O2 data
- Under Measured Variables, select "Oxygen" for 1 and 2
- Under Data Exclusion, select "Oxygen" and "ALL" for both profile and sample
- Proceed to download the data
- Repeat the above process for CTD oxygen (separate download)
- Pre-downloaded OSD and CTD data as of May 2022 used for the 2022 Liege Colloquium is available [here](https://www.dropbox.com/sh/ivfo1yicivwaea7/AADiYhXFX8NROevV3yucOrLca?dl=0)

### 2. Binning
- Execute the following two codes: 
- bin_wod_OSD.ipynb
- bin_wod_CTD.ipynb
- Execute two scripts to generate annually binned data

### 3. Monthly climatology
- Execute the following one code: 
- stat_mean_bin.ipynb
- Generate statistical mean climatology data
- Next perform optimal interpolation to establish monthly climatology
- This part is in MATLAB, in two functions: 
- use the function objmap_o2clim(month,zlevel) and loop over 12 months and 47 zlevels
- Execite the following MATLAB command:
> for month=1:12; for zlevel=1:47; objmap_o2clim(month,zlevel); end; end;
- This function requires two supporting files (basin_name.txt and basin_mask_01.nc) available from the Dropbox link
- use the function gen_netcdf_clim to generate a single netCDF file
- Execute the following MATLAB command:
> gen_netcdf_clim;

### 4. Generate anomaly fields
- Execute the following one code: 
- pentadal_o2anom.ipynb
- Subtract the monthly climatology to calculate the anomaly
- Aggregate monthly into yearly anomaly
- Combine OSD and CTD data into a single field after 1987
- Apply 5 year running mean

### 5. Perform optimal interpolation
- This part is done by the MATLAB script
- use the objmap_o2 function and loop over the years
- Execute the following MATLAB command:
> for t=1:52; objmap_o2(t); end;
- use the function gen_netcdf to generate a single netCDF file
- Execute the following MATLAB command: 
> gen_netcdf;
