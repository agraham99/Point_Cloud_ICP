# Iterative Closest Point (ICP) matching as a warping correction tool between two point clouds

## Introduction

Point clouds generated from Structure from Motion (SfM) are known to exhibit systematic distortions relative to absolute coordinate systems, especially in the vertical direction, largely due to radial distortion effects caused by imperfections in consumer grade cameras (James and Robson 2014).  Where proper camera calibration is challenging or infeasible, this process serves as a point cloud post processing procedure which corrects offsets between an imprecise SfM point cloud known to contain distortion, and a precisely georeferenced reference point cloud.  

James, Mike R., and Stuart Robson. 2014. “Mitigating Systematic Error in Topographic Models Derived from UAV and Ground-Based Image Networks.” Earth Surface Processes and Landforms 39 (10): 1413–1420. doi:10.1002/esp.3609.

## Prerequisites 

- R version 3.4.2 
### R packages 

```
R version 3.4.2 (2017-09-28)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows >= 8 x64 (build 9200)

Matrix products: default

locale:
[1] LC_COLLATE=English_Canada.1252  LC_CTYPE=English_Canada.1252   
[3] LC_MONETARY=English_Canada.1252 LC_NUMERIC=C                   
[5] LC_TIME=English_Canada.1252    

attached base packages:
[1] parallel  grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] UsefulRFunctions_0.1.0 tictoc_1.0             forcats_0.3.0         
 [4] stringr_1.3.1          dplyr_0.7.8            purrr_0.2.5           
 [7] readr_1.1.1            tidyr_0.8.1            tibble_1.4.2          
[10] tidyverse_1.2.1        rlas_1.3.0             doParallel_1.0.14     
[13] iterators_1.0.8        foreach_1.4.4          gtable_0.2.0          
[16] gridExtra_2.3          lattice_0.20-38        lidR_2.0.1            
[19] polynom_1.3-9          spatial_7.3-11         ggplot2_3.1.0         
[22] rgeos_0.3-26           raster_2.6-7           rgdal_1.2-15          
[25] sp_1.2-5              

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.0        lubridate_1.7.4   assertthat_0.2.0  R6_2.3.0         
 [5] cellranger_1.1.0  plyr_1.8.4        backports_1.1.1   httr_1.3.1       
 [9] pillar_1.3.0      rlang_0.3.0.1     lazyeval_0.2.1    readxl_1.1.0     
[13] rstudioapi_0.7    data.table_1.11.4 Matrix_1.2-11     labeling_0.3     
[17] munsell_0.5.0     broom_0.5.0       compiler_3.4.2    modelr_0.1.2     
[21] pkgconfig_2.0.2   mgcv_1.8-20       tidyselect_0.2.5  codetools_0.2-15 
[25] crayon_1.3.4      withr_2.1.2       ggpubr_0.1.7      nlme_3.1-131     
[29] jsonlite_1.5      magrittr_1.5      scales_1.0.0      cli_1.0.1        
[33] stringi_1.2.4     bindrcpp_0.2.2    xml2_1.1.1        tools_3.4.2      
[37] glue_1.3.0        hms_0.4.2         yaml_2.1.14       colorspace_1.3-2 
[41] rvest_0.3.2       bindr_0.1.1       haven_2.0.0
```        

- CloudCompare version 2.10 or 2.11 alpha
https://www.danielgm.net/cc/

## Usage

### ICP_moving_window.R



