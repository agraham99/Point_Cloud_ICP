# Iterative Closest Point (ICP) matching as a warping correction tool between two point clouds

### Introduction

Point clouds generated from Structure from Motion (SfM) are known to exhibit systematic distortions relative to absolute coordinate systems, especially in the vertical direction, largely due to radial distortion effects caused by imperfections in consumer grade cameras (James and Robson 2014).  Where proper camera calibration is challenging or infeasible, this process serves as a point cloud post processing procedure which corrects offsets between an imprecise SfM point cloud known to contain distortion, and a precisely georeferenced reference point cloud.  

James, Mike R., and Stuart Robson. 2014. “Mitigating Systematic Error in Topographic Models Derived from UAV and Ground-Based Image Networks.” Earth Surface Processes and Landforms 39 (10): 1413–1420. doi:10.1002/esp.3609.

### Prerequisites 

#### R packages 

```
> sessionInfo()
R version 3.6.0 (2019-04-26)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows >= 8 x64 (build 9200)

Matrix products: default

locale:
[1] LC_COLLATE=English_Canada.1252  LC_CTYPE=English_Canada.1252    LC_MONETARY=English_Canada.1252
[4] LC_NUMERIC=C                    LC_TIME=English_Canada.1252    

attached base packages:
[1] parallel  grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] mapview_2.6.3          concaveman_1.0.0       spatstat_1.59-0        rpart_4.1-15          
 [5] nlme_3.1-139           spatstat.data_1.4-0    UsefulRFunctions_0.1.0 lidR_2.0.3            
 [9] sf_0.7-4               tictoc_1.0             forcats_0.4.0          stringr_1.4.0         
[13] dplyr_0.8.0.1          purrr_0.3.2            readr_1.3.1            tidyr_0.8.3           
[17] tibble_2.1.1           tidyverse_1.2.1        rlas_1.3.2             doParallel_1.0.14     
[21] iterators_1.0.10       foreach_1.4.4          gtable_0.3.0           gridExtra_2.3         
[25] lattice_0.20-38        polynom_1.4-0          spatial_7.3-11         ggplot2_3.1.1         
[29] rgeos_0.4-3            raster_2.8-19          rgdal_1.4-3            sp_1.3-1              

loaded via a namespace (and not attached):
 [1] colorspace_1.4-1      deldir_0.1-16         class_7.3-15          leaflet_2.0.2        
 [5] rprojroot_1.3-2       satellite_1.0.1       base64enc_0.1-3       fs_1.2.7             
 [9] rstudioapi_0.10       listenv_0.7.0         remotes_2.0.4         lubridate_1.7.4      
[13] xml2_1.2.0            codetools_0.2-16      splines_3.6.0         knitr_1.22           
[17] polyclip_1.10-0       pkgload_1.0.2         jsonlite_1.6          broom_0.5.2          
[21] png_0.1-7             shiny_1.3.2           compiler_3.6.0        httr_1.4.0           
[25] backports_1.1.4       assertthat_0.2.1      Matrix_1.2-17         lazyeval_0.2.2       
[29] cli_1.1.0             later_0.8.0           htmltools_0.3.6       prettyunits_1.0.2    
[33] tools_3.6.0           glue_1.3.1            Rcpp_1.0.1            cellranger_1.1.0     
[37] crosstalk_1.0.0       xfun_0.6              globals_0.12.4        ps_1.3.0             
[41] rvest_0.3.3           mime_0.6              devtools_2.0.2        goftest_1.1-1        
[45] future_1.12.0         scales_1.0.0          promises_1.0.1        hms_0.4.2            
[49] spatstat.utils_1.13-0 yaml_2.2.0            curl_3.3              memoise_1.1.0        
[53] stringi_1.4.3         desc_1.2.0            e1071_1.7-1           pkgbuild_1.0.3       
[57] rlang_0.3.4           pkgconfig_2.0.2       evaluate_0.13         tensor_1.5           
[61] htmlwidgets_1.3       labeling_0.3          processx_3.3.0        tidyselect_0.2.5     
[65] plyr_1.8.4            magrittr_1.5          R6_2.4.0              generics_0.0.2       
[69] DBI_1.0.0             pillar_1.3.1          haven_2.1.0           withr_2.1.2          
[73] mgcv_1.8-28           units_0.6-2           abind_1.4-5           modelr_0.1.4         
[77] crayon_1.3.4          KernSmooth_2.23-15    rmarkdown_1.13        progress_1.2.0       
[81] usethis_1.5.0         readxl_1.3.1          data.table_1.12.2     callr_3.2.0          
[85] webshot_0.5.1         digest_0.6.18         classInt_0.3-3        xtable_1.8-4         
[89] httpuv_1.5.1          stats4_3.6.0          munsell_0.5.0         viridisLite_0.3.0    
[93] sessioninfo_1.1.1    
```        

- [CloudCompare](https://www.danielgm.net/cc/) version 2.10 or 2.11 alpha

----------------------------------------------------------------------------------------------------------------------------------


### Usage


#### Overview

1. Generate shift observations using [ICP_moving_window.R](https://github.com/spireaero/ICP/blob/master/ICP_moving_window.md)
2. Filter, visualize and generate predictive models of the shifts across space using [ICP_visualization.R](https://github.com/spireaero/ICP/blob/master/ICP_visualization.md)
3. Apply the non-constant transformations to the original point clouds according to the predictive models using [ICP_unwarp.R](https://github.com/spireaero/ICP/blob/master/ICP_unwarp.md)

#### [ICP Moving Window](https://github.com/spireaero/ICP/blob/master/ICP_moving_window.md)

[ICP_moving_window.R](https://github.com/spireaero/ICP/blob/master/ICP_moving_window.R) serves to create the observation points which are location estimations of the X, Y and Z offset values across the area.

The reference and data (to-be-aligned) point-cloud are clipped to the extent of each window.  The point density of the clipped point-clouds are reduced to the match one another using the [subsample tool](https://www.cloudcompare.org/doc/wiki/index.php?title=Edit%5CSubsample) found in CloudCompare with the 'Spatial' option.  The  ICP algorithm is then performed on the resulting subsampled point clouds.

![Image of ICP Moving Window](https://github.com/spireaero/ICP/blob/master/images/README_Figure_1.png)  

**Figure 1.** Schematic of moving window operation where ICP estimated results are stored at each point.

See the Markdown file for [ICP Moving Window](https://github.com/spireaero/ICP/blob/master/ICP_moving_window.md)

----------------------------------------------------------------------------------------------------------------------------------

#### [ICP Visualization](https://github.com/spireaero/ICP/blob/master/ICP_visualization.md)

[ICP_visualization.R](https://github.com/spireaero/ICP/blob/master/ICP_visualization.R) is used for setting filtering criteria to remove noise from the ICP observation points.  It is also used to generate and visualize the results of predictive modelling using the ICP observation points as input.

See the Markdown file for [ICP visualization](https://github.com/spireaero/ICP/blob/master/ICP_visualization.md)

----------------------------------------------------------------------------------------------------------------------------------

#### [ICP Unwarp](https://github.com/spireaero/ICP/blob/master/ICP_unwarp.md)

[ICP_unwarp.R](https://github.com/spireaero/ICP/blob/master/ICP_unwarp.R) is used to apply the shifts according to the models developed in *ICP Visualization*.

See the Markdown file for [ICP Unwarp](https://github.com/spireaero/ICP/blob/master/ICP_unwarp.md)



