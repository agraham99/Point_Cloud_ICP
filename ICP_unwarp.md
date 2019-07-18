ICP Unwarp
================
Alex Graham
July 15, 2019

------------------------------------------------------------------------

#### Overview

This script uses polynomial model objects created in *ICP\_visualization.R* and applies a translational shift to a point cloud.

#### Inputs

Read the raw point cloud to be shifted and the *lm* objects generated using *ICP\_visualization.R*

``` r
raw_cloud = "<path to las/laz file(s)>"

icp.lm.x = readRDS(file = "data/lm/icp_lm_x.rds")
icp.lm.y = readRDS(file = "data/lm/icp_lm_y.rds")
icp.lm.z = readRDS(file = "data/lm/icp_lm_z.rds")
```

#### Apply shifts using the LAScatalog processing engine

For more info see:

-   <https://cran.r-project.org/web/packages/lidR/vignettes/lidR-catalog-apply-examples.html>

-   <https://gis.stackexchange.com/questions/311150/non-constant-transformation-of-photogrammetric-point-clouds-in-r>

First create a generic function usable for mulitple classes

``` r
lasshift = function(las)
{
  UseMethod("lasshift", las)
}
```

Create a function which works on class LAS and performs the coordinate shifting

``` r
lasshift.LAS = function(las)
{
  # dataframe containing the x and y coords of the original cloud
  g = data.frame(list(las$X, las$Y))
  print(names(g))
  names(g) = c('x', 'y')
  
  # predict the shifts
  xshift = stats::predict(icp.lm.x, newdata = g)
  yshift = stats::predict(icp.lm.y, newdata = g)
  zshift = stats::predict(icp.lm.z, newdata = g)
  
  # apply the shift values 
  las$X = las$X + xshift
  las$Y = las$Y + yshift
  las$Z = las$Z + zshift
  
  return(las)
}
```

Create a catalog\_apply compatible function

``` r
lasshift.LAScluster = function(las)
{
  las <- readLAS(las)                        
  if (is.empty(las)) return(NULL)           
  
  las <- lasshift(las)
  return(las)
}
```

Create a method for LAScatalog object

``` r
lasshift.LAScatalog = function(las)
{
  # Force some options
  # opt_select(las) <-  "*" 
  opt_chunk_buffer(las) <- 0
  
  # options <- list(need_output_file = TRUE)
  
  output <- catalog_apply(las, lasshift)
  output <- unlist(output)
  output <- catalog(output)
  return(output)
}
```

Specify the input catalog object and some options, then run the lasshift method on the catalog

``` r
ctg = catalog(raw_cloud)
# opt_select(ctg) <- "xyz"         # read only the coordinates.
opt_chunk_size(ctg) <- 100       # process in tiles of ___ meters
opt_output_files(ctg) <- "<path to newly shifted point clouds>/{ID}_shifted"

# Objects in the global environnement are not automatically exported in each worker. 
# Set opt_cores to 1 is the solution in version 2.0.y. 
# In version 2.1.0 you can export manually some object in each worker using the package future. 
# See also https://github.com/Jean-Romain/lidR/blob/master/NEWS.md. 
opt_cores(ctg) <- 1L            

# add a buffer to the tiles
opt_chunk_buffer(ctg) <- 0       

# run the lasshift catalog_apply method. 
new_ctg = lasshift(ctg)
```
