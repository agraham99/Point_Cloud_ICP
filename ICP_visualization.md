ICP Visualization
================
Alex Graham
July 15, 2019

 

------------------------------------------------------------------------

#### Overview

This script serves to visualize and filter the ICP observation points which will be input to the model development (ICP\_modelling). A series of filtering criteria can be set here. The points will be written out to be used in ICP\_modelling.

``` r
# install.packages("rgdal")
# install.packages("sp")
# install.packages("raster")
# install.packages("rgeos")
# install.packages("ggplot2")
# install.packages("spatial")
# install.packages("polynom")
# install.packages("lidR")
# install.packages("grid")
# install.packages("lattice")
# install.packages("gridExtra")
# install.packages("ggplot2")
# install.packages("gtable")
# install.packages("doParallel")
# install.packages("foreach")
# install.packages("rlas")
# install.packages("lidR")
# install.packages("tidyverse")
# install.packages("dplyr")
# install.packages("tictoc")

library(rgdal)
library(sp)
library(raster)
library(rgeos)
library(ggplot2)
library(spatial)
library(polynom)
library(lidR)
library(grid)
library(lattice)
library(gridExtra)
library(ggplot2)
library(gtable)
library(doParallel)
library(foreach)
library(rlas)
library(lidR)
library(tidyverse)
library(dplyr)
library(tictoc)
library(sf)
```

Project coordinate system (proj4text)

``` r
crs = "+proj=utm +zone=11 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "
```

 
-

#### Parse ICP observation points data

 

Read the csv created as a result of ICP\_moving\_window.R

The csv should contain columns:

*i, isnull, coverage, ICP, x, y, RMS, r1c1, r2c1, r3c1, r1c2, r2c2, r3c2, r1c3, r2c3, r3c3, r1c4, r2c4, r3c4*

Where 'r' and 'c' refer to row and column of the ICP transformation matrix

``` r
# points.f = paste0("<path to the csv file of ICP observation points created using ICP_moving_window.R>")
points.f = paste0("D:/JOE_RAKOFSKY/ICP_points/step_150_win_30_canopy_TRUE_icp_obs.csv")
points.df = read.csv(points.f, header=T, stringsAsFactors = F)
```

 

Gather all points where ICP was run and therefore RMS is not NA

``` r
# store as a seperate variable (pts) incase we want to see these points
pts = subset(points.df, ICP == T & !is.na(RMS))
# assign pts as p.  p becomes the the points of interest for remaining analysis
p = pts
```

 

Make sure that certain columns are treated as numeric and rename translational columns

``` r
# these are columns with numeric values from the ICP operation
numcols = c(3, 5:19)
p[,numcols] = apply(p[,numcols], 2, function(x) as.numeric(as.character(x)));
# change the names of the columns associated with the translational transfoamtion values (r1c4, r2c4, and r3c4)
names(p)[17:19] = c('x_trans', 'y_trans', 'z_trans')
```

 

Make the points table into SpatialPointsDataFrame

``` r
p = SpatialPointsDataFrame(p[,5:6], p)
proj4string(p) = crs
# optional writing to shapefile
# writeOGR(p, "D:/JOE_RAKOFSKY/all_points.shp",
#          layer = basename(points.f),
#          driver = 'ESRI Shapefile',
#          overwrite_layer = T)
```

 

Optionally, points can be masked at this stage by an external shapefile

``` r
mask_layer.f = "D:/ICP/mask_layers/mask.shp"
mask_layer = readOGR(mask_layer.f)

crs(mask_layer) = crs(p)
p_masked = p[mask_layer,]

p = p_masked

# option to write out masked points
# writeOGR(p_masked, "D:/ICP/points/all_points_masked.shp",
#          layer = 'p_masked',
#           driver = 'ESRI Shapefile',
#          overwrite_layer = T)
```

 

------------------------------------------------------------------------

#### Filter points

 

Quantile threhsolds are used below to filter points based on ICP estimated rotation, RMS and translation. Trial and error may be needed to determine optimal threshold values. The values used below were derived from UAS-DAP point cloud matching with ALS cloud over AFRF Gavin Lake (near Williams Lake).

*thr.rotation* is used to filter values in the diagonal of the 3x3 rotation matrix (r1c1, r2c2, r3c3) and is a quantile value where only points greater than the threshold are selected. Values nearest 1 on the diagonal have minimal rotation so we use this to filter the results. There is a tradeoff between having sufficient number of points to have wall-to-wall coverage for the modelling stage and incluiding accurate matching observations.

*thr.RMS* is used for isolating RMS values below the quantile threhsold

*thr.xtrans*, *thr.ytrans*, and *thr.ztrans* are filtered in the same manner using a two-tailed cutoff with the *alpha* value as the quantile. Filtering should not be performed on the translation vectors that are to be modelled. For example, if we want to produce a Z translation model, we can comment out the filtering for *thr.ztrans* as is done below.

``` r
thr.rotation = 0.25
thr.RMS = 0.35
alpha = 0.01

tail_min = alpha
tail_max = 1 - alpha

thr.xtrans.min = tail_min
thr.ytrans.min = tail_min
thr.ztrans.min = tail_min

thr.xtrans.max = tail_max
thr.ytrans.max = tail_max
thr.ztrans.max = tail_max

# perform the subsetting of valid icp points
p = subset(p, r1c1 > quantile(r1c1, thr.rotation) &
             r2c2 > quantile(r2c2, thr.rotation) &
             r3c3 > quantile(r3c3, thr.rotation) &
             # filter for RMS values
             RMS < quantile(RMS, thr.RMS) &
             # filter the x and y translation low values
             x_trans > quantile(x_trans, thr.xtrans.min) &
             y_trans > quantile(y_trans, thr.ytrans.min) &
             # z_trans > quantile(z_trans, thr.ztrans.min) &
             # filter the x and y translation high values
             x_trans < quantile(x_trans, thr.xtrans.max) &
             y_trans < quantile(y_trans, thr.ytrans.max))
             # z_trans < quantile(z_trans, thr.ztrans.max))
```

 

------------------------------------------------------------------------

#### Histogram generation

 

Produce histograms of the transforamtion vlaues created from ICP. Arrange histograms 4 x 4 matrix mimicking the format of the original ICP output matrix

``` r
par(mfrow = c(3,4))
# number of bins
b = 100
# row 1 of matrix
hist(p$r1c1, b)
hist(p$r1c2, b)
hist(p$r1c3, b)
hist(p$x_trans, b)
# row 2
hist(p$r2c1, b)
hist(p$r2c2, b)
hist(p$r2c3, b)
hist(p$y_trans, b)
# row 3
hist(p$r3c1, b)
hist(p$r3c2, b)
hist(p$r3c3, b)
hist(p$z_trans, b)
```

 

------------------------------------------------------------------------

#### Spatial plots of points

 

``` r
# spatial plot of matrix results
# use trans = 'x_trans' or 'y_trans' for other directions
trans = 'z_trans'
n = abs(round((max(p[[trans]]) - min(p[[trans]]))/9, 4))
cuts.p = seq(min(p[[trans]]), max(p[[trans]]), n)
cuts.p = round(cuts.p,4)
trellis.par.set(axis.line=list(col=NA)) 
spplot(p, zcol = trans, cuts = cuts.p, key.space = 'right', cex = 0.5)

# option to write the shp file of the subset of points
f = gsub('.csv', '.shp', basename(points.f))
writeOGR(p, paste("D:/JOE_RAKOFSKY/shp/", f, sep = ''), layer = f, driver = "ESRI Shapefile", overwrite_layer = T)
```

#### Make spatial plots for each direction (x, y, z) in a loop

 

``` r
# empty list to populate with spplot objects
plots = list()
i = 0
for (trans in c('x_trans', 'y_trans', 'z_trans')){
  i = i+1

  # run this block alone for
  n = abs(round((max(p[[trans]]) - min(p[[trans]]))/7, 2))
  cuts.p = seq(min(p[[trans]]), max(p[[trans]]), n)
  cuts.p = round(cuts.p,1)
  plot = spplot(p, zcol = trans, cuts = cuts.p, key.space = 'right', cex = 0.5, xlab = trans)

  plots[[i]] = plot
}
# arrange the plots for display
do.call(grid.arrange, plots)
```