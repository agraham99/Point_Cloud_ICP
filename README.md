# Iterative Closest Point (ICP) matching as a warping correction tool between two point clouds

## Introduction

Point clouds generated from Structure from Motion (SfM) are known to exhibit systematic distortions relative to absolute coordinate systems, especially in the vertical direction, largely due to radial distortion effects caused by imperfections in consumer grade cameras (James and Robson 2014).  Where proper camera calibration is challenging or infeasible, this process serves as a point cloud post processing procedure which corrects offsets between an imprecise SfM point cloud known to contain distortion, and a precisely georeferenced reference point cloud.  

James, Mike R., and Stuart Robson. 2014. “Mitigating Systematic Error in Topographic Models Derived from UAV and Ground-Based Image Networks.” Earth Surface Processes and Landforms 39 (10): 1413–1420. doi:10.1002/esp.3609.

## Prerequisites 

### R packages 

```
R version 3.6.0 (2019-04-26)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows >= 8 x64 (build 9200)

Matrix products: default

locale:
[1] LC_COLLATE=English_Canada.1252  LC_CTYPE=English_Canada.1252    LC_MONETARY=English_Canada.1252
[4] LC_NUMERIC=C                    LC_TIME=English_Canada.1252    

attached base packages:
[1] grid      parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] tictoc_1.0          forcats_0.4.0       stringr_1.4.0       dplyr_0.8.0.1       purrr_0.3.2        
 [6] readr_1.3.1         tidyr_0.8.3         tibble_2.1.1        tidyverse_1.2.1     rlas_1.3.2         
[11] gtable_0.3.0        gridExtra_2.3       lattice_0.20-38     polynom_1.4-0       spatial_7.3-11     
[16] ggplot2_3.1.1       rgdal_1.4-3         mapview_2.6.3       concaveman_1.0.0    spatstat_1.59-0    
[21] rpart_4.1-15        nlme_3.1-139        spatstat.data_1.4-0 rgeos_0.4-3         doParallel_1.0.14  
[26] iterators_1.0.10    foreach_1.4.4       lidR_2.0.2          raster_2.8-19       sp_1.3-1           

loaded via a namespace (and not attached):
 [1] httr_1.4.0            jsonlite_1.6          viridisLite_0.3.0     splines_3.6.0        
 [5] modelr_0.1.4          shiny_1.3.2           assertthat_0.2.1      stats4_3.6.0         
 [9] cellranger_1.1.0      yaml_2.2.0            pillar_1.3.1          backports_1.1.4      
[13] glue_1.3.1            digest_0.6.18         promises_1.0.1        polyclip_1.10-0      
[17] rvest_0.3.3           colorspace_1.4-1      htmltools_0.3.6       httpuv_1.5.1         
[21] Matrix_1.2-17         plyr_1.8.4            pkgconfig_2.0.2       broom_0.5.2          
[25] haven_2.1.0           xtable_1.8-4          scales_1.0.0          webshot_0.5.1        
[29] tensor_1.5            satellite_1.0.1       later_0.8.0           spatstat.utils_1.13-0
[33] mgcv_1.8-28           generics_0.0.2        withr_2.1.2           lazyeval_0.2.2       
[37] cli_1.1.0             readxl_1.3.1          magrittr_1.5          crayon_1.3.4         
[41] mime_0.6              deldir_0.1-16         xml2_1.2.0            class_7.3-15         
[45] tools_3.6.0           data.table_1.12.2     hms_0.4.2             munsell_0.5.0        
[49] compiler_3.6.0        e1071_1.7-1           rlang_0.3.4           classInt_0.3-3       
[53] units_0.6-2           rstudioapi_0.10       htmlwidgets_1.3       goftest_1.1-1        
[57] crosstalk_1.0.0       base64enc_0.1-3       codetools_0.2-16      abind_1.4-5          
[61] DBI_1.0.0             R6_2.4.0              lubridate_1.7.4       KernSmooth_2.23-15   
[65] stringi_1.4.3         Rcpp_1.0.1            sf_0.7-4              png_0.1-7            
[69] leaflet_2.0.2         tidyselect_0.2.5     
```        

- [CloudCompare](https://www.danielgm.net/cc/) version 2.10 or 2.11 alpha


# Usage

## ICP_moving_window.R

This script serves to create the observation points which are location estimations of the X, Y and Z offset values across the area.

The reference and data (to-be-aligned) point-cloud are clipped to the extent of each window.  The point density of the clipped point-clouds are reduced to the match one another using the [subsample tool](https://www.cloudcompare.org/doc/wiki/index.php?title=Edit%5CSubsample) found in CloudCompare with the 'Spatial' option.  The  ICP algorithm is then performed on the resulting subsampled point clouds.

![Image of ICP Moving Window](https://github.com/spireaero/ICP/blob/master/images/moving_window.png)  

### Parameters

**CloudCompare**: path to the locally installed CloudCompare.exe for interfacing the commandline version

**crs**: Coordinate System in Proj4txt form


**STEP**: desired distance in meters between moving window centers (ICP observation points)

**WIN_SIZE**: width/height of each square moving window


**CANOPY_ONLY**: Boolean (T/F) whether or not to use 'top-points only' AKA Digital Surface Model (DSM) points

**SMALL_TILE_CORES**: number of corse to use (LAScatalog) for picking out canopy points

**SMALL_TILE_WINDOW**: window size or gridcell size in meters to designate DSM points within the larger moving window

**SMALL_TILE_BUFFER**: optional buffer for small window size

**SUBSAMPLE_DISTANCE**: minimum distance between points of both data and reference point clouds to subsample to before running the ICP

**RANDOM_SAMPLE_LIMIT**: maximum number of points randomly sampled and used in the RMS calculation for the ICP procedure

**ICP_OVERLAP**: estimated final 2-dimensional top down view overlap between the two clouds after running ICP





