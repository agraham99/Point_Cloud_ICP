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


# coordinate system - proj4text
crs = "+proj=utm +zone=11 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "

# ------------------------------------------------------------------
# PARSE DATA

# read the csv created as a result of ICP_moving_window.R
# the csv should contain columns:  
# n, isnull, coverage,   ICP,        x,          y, RMS, r1c1, r2c1, r3c1, r1c2, r2c2, r3c2, r1c3, r2c3, r3c3, r1c4, r2c4, r3c4
# where 'r' and 'c' refer to row and column of the ICP transformation matrix
# points.f = paste0("<path to the csv file of ICP observation points created using ICP_moving_window.R>")
points.f = paste0("D:/JOE_RAKOFSKY/ICP_points/step_150_win_30_canopy_TRUE_icp_obs.csv")
points.df = read.csv(points.f, header=T, stringsAsFactors = F)

# all points where ICP was run and therefore RMS is not NA
# store as a seperate variable (pts) incase we want to see these points
pts = subset(points.df, ICP == T & !is.na(RMS))
# assign pts as p.  p becomes the the points of interest for remaining analysis
p = pts
# make sure that the columns below are numeric type
# these are columns with numeric values from the ICP operation
numcols = c(3, 5:19)
p[,numcols] = apply(p[,numcols], 2, function(x) as.numeric(as.character(x)));
# change the names of the columns associated with the translational transfoamtion values (r1c4, r2c4, and r3c4)
names(p)[17:19] = c('x_trans', 'y_trans', 'z_trans')





# ------------------------------------------------------------------
# DEFINE POINTS AS SPATIAL
p = SpatialPointsDataFrame(p[,5:6], p)
proj4string(p) = crs
# for writing to shapefile
# writeOGR(p, "D:/JOE_RAKOFSKY/mask_layers/cutblocks2011_2014_plus_roads/ICP_pts_cutblocks2011_2014_plus_roads_sl_subset.shp",
#          layer = basename(points.f),
#          driver = 'ESRI Shapefile',
#          overwrite_layer = T)


# ------------------------------------------------------------------
# MASK THE SPATIAL POINTS
mask_layer.f = "D:/JOE_RAKOFSKY/mask_layers/cutblocks2011_2014_plus_roads/cutblocks2011_2014_plus_roads_sl_subset.shp"
mask_layer = readOGR(mask_layer.f)

crs(mask_layer) = crs(p)
p_masked = p[mask_layer,]

p = p_masked

# writeOGR(p_masked, "D:/JOE_RAKOFSKY/mask_layers/cutblocks2011_2014_plus_roads/ICP_masked_points_cutblocks2011_2014_plus_roads.shp",
#          layer = 'p_masked',
#           driver = 'ESRI Shapefile',
#          overwrite_layer = T)
# 


# ------------------------------------------------------------------
# FILTER DATA
# these are the diagonal values of the 3x3 rotation matrix (first three columns are rows of the 4x4 ICP matrix)
# values nearest 1 on the diagonal have minimal rotation so we use this to filter the results. This was found thru trial and error
# We want a high number of points to use for modelling however we want these points to have accurate true ICP generated
# transforamtion values

# threhsolds to isolate "valid" icp results
# thr.rotation a quantile value where only points greater than the threshold are selected
# thr.RMS is a quantile value where only points below the threhsold are selected
thr.rotation = 0.65
thr.RMS = 0.35

# thr.xtrans has min and max to cut off tails of the distribution 
# specify the quantile threhsold for two tail exclusion
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

# ------------------------------------------------------------------
# HISTOGRAMS of matrix values
# produce histograms of the transforamtion vlaues created from ICP

# arrange histograms 4 x 4 matrix mimicking the format of the original ICP output matrix
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

# SPATIAL PLOTS of the points
# ----------------------------------------------------------------
# spatial plot of matrix results
# use trans = 'x_trans' or 'y_trans' for other directions
trans = 'z_trans'
n = abs(round((max(p[[trans]]) - min(p[[trans]]))/9, 4))
cuts.p = seq(min(p[[trans]]), max(p[[trans]]), n)
cuts.p = round(cuts.p,4)
trellis.par.set(axis.line=list(col=NA)) 
spplot(p, zcol = trans, cuts = cuts.p, key.space = 'right', cex = 0.5)

# write the shp file of the subset of points
f = gsub('.csv', '.shp', basename(points.f))
writeOGR(p, paste("D:/JOE_RAKOFSKY/shp/", f, sep = ''), layer = f, driver = "ESRI Shapefile", overwrite_layer = T)

# ------------------------------------------------------------------
# MAKE SPATIAL PLOTS of the observed shifts IN A LOOP
# empty list to populate with spplot objects
# plots = list()
# i = 0
# for (trans in c('x_trans', 'y_trans', 'z_trans')){
#   i = i+1
# 
#   # run this block alone for
#   n = abs(round((max(p[[trans]]) - min(p[[trans]]))/7, 2))
#   cuts.p = seq(min(p[[trans]]), max(p[[trans]]), n)
#   cuts.p = round(cuts.p,1)
#   plot = spplot(p, zcol = trans, cuts = cuts.p, key.space = 'right', cex = 0.5, xlab = trans)
# 
#   plots[[i]] = plot
# }
# # arrange the plots for display
# do.call(grid.arrange, plots)




# ------------------------------------------------------------------
# POLYNOMIAL MODELLING OF THE SHIFT VALUES

# the 2d coordinates of each point
x = p$x
y = p$y
# get the ICP translation columns
x.offset = p$x_trans
y.offset = p$y_trans
z.offset = p$z_trans

# ------------------------------------------------------------------
# 2D POINT PLOTTING
# change the y argument and model for X, Y or Z offsets

# plot the shifts as a function of x or y direction
ggp = ggplot()
ggp = ggp + geom_point(aes(x = y, y = z.offset), alpha = 0.1)
ggp = ggp + geom_smooth(aes(x = y, y = z.offset))
ggp


# ------------------------------------------------------------------
# generate a series of candidate polynomial models to predict the X Y and Z offsets

# define a range of polynomial factors to test for modelling across space
powersX = (seq(1,20,1))
powersY = (seq(1,20,1))

# the following creates a table to be populated with model Rsq and adjRsq values for each polymonial factor tested
# the models with highest adjRsq value for each of X Y and Z models will be chosen as final
# table will n rows equal to the product of x and y powers tested, times 3 (X, y and Z)
n = length(powersX)*length(powersY)*3
# build empty table
models.table = data.frame(Shift = character(n),
                          pwrX=numeric(n), 
                          pwrY=numeric(n), 
                          Rsq = numeric(n),
                          adjRsq = numeric(n),
                          AIC = numeric(n),
                          stringsAsFactors=FALSE) 

# i is used to index the models.table and populate it
i = 0
# switch for 'raw' argument in poly - see ?poly
raw.sw = F
for (pwrX in powersX){
  for (pwrY in powersY){
    i = i+1
    # create models with the current powers to predict the offset values from ICP
    icp.lm.x = lm(x.offset ~ poly(x, pwrX, raw= raw.sw) + poly(y, pwrY, raw= raw.sw))
    icp.lm.y = lm(y.offset ~ poly(x, pwrX, raw= raw.sw) + poly(y, pwrY, raw= raw.sw))
    icp.lm.z = lm(z.offset ~ poly(x, pwrX, raw= raw.sw) + poly(y, pwrY, raw= raw.sw))
    # model summaries
    sx = summary(icp.lm.x)
    sy = summary(icp.lm.y)
    sz = summary(icp.lm.z)
    # ANOVA of each model
    anova(anova.x <- icp.lm.x)
    anova(anova.y <- icp.lm.y)
    anova(anova.z <- icp.lm.z)
    # get AIC of each ANOVA
    AIC.x = AIC(anova.x)
    AIC.y = AIC(anova.y)
    AIC.z = AIC(anova.z)
    # make row for each model 
    rowX = list('X', pwrX, pwrY, sx$r.squared, sx$adj.r.squared, AIC.x)
    rowY = list('Y', pwrX, pwrY, sy$r.squared, sy$adj.r.squared, AIC.y)
    rowZ = list('Z', pwrX, pwrY, sz$r.squared, sz$adj.r.squared, AIC.z)
    # put rows in to the table
    models.table[i,] = rowX
    models.table[i+1,] = rowY
    models.table[i+2,] = rowZ
    i = i+2
  }
}

# subset the models for each 3D direction
models.X = subset(models.table, Shift == 'X')
models.Y = subset(models.table, Shift == 'Y')
models.Z = subset(models.table, Shift == 'Z')
# get the model with highest adjRsq
best.model.X = subset(models.X, adjRsq == max(adjRsq))
best.model.X = subset(best.model.X, pwrX == min(pwrX) & pwrY == min(pwrY))
best.model.Y = subset(models.Y, adjRsq == max(adjRsq))
best.model.Y = subset(best.model.Y, pwrX == min(pwrX) & pwrY == min(pwrY))
best.model.Z = subset(models.Z, adjRsq == max(adjRsq))
best.model.Z = subset(best.model.Z, pwrX == min(pwrX) & pwrY == min(pwrY))

# define the final modelbased on best adjRsq
icp.lm.x = lm(x.offset ~ poly(x, best.model.X$pwrX, raw= raw.sw) + poly(y, best.model.X$pwrY, raw= raw.sw))
icp.lm.y = lm(y.offset ~ poly(x, best.model.Y$pwrX, raw= raw.sw) + poly(y, best.model.Y$pwrY, raw= raw.sw))
icp.lm.z = lm(z.offset ~ poly(x, best.model.Z$pwrX, raw= raw.sw) + poly(y, best.model.Z$pwrY, raw= raw.sw))
# model summaries
summary(icp.lm.x)
summary(icp.lm.y)
summary(icp.lm.z)
# sd of residuals
sd(residuals(icp.lm.x))
sd(residuals(icp.lm.y))
sd(residuals(icp.lm.z))

# check the residual plots
plot(fitted(icp.lm.x), residuals(icp.lm.x))
plot(fitted(icp.lm.y), residuals(icp.lm.y))
plot(fitted(icp.lm.z), residuals(icp.lm.z))

# confidence intervals for variables in the models
confint(icp.lm.x, level=0.999)
confint(icp.lm.y, level=0.999)
confint(icp.lm.z, level=0.999)

# add the RMSE
RSS <- c(crossprod(icp.lm.z$residuals))
MSE <- RSS / length(icp.lm.z$residuals)
RMSE <- sqrt(MSE)

# bind the residuals from the model to tthe original points table
p = cbind(p, icp.lm.x$residuals)
names(p)[length(names(p))] = 'Xresidual'
p = cbind(p, icp.lm.y$residuals)
names(p)[length(names(p))] = 'Yresidual'
p = cbind(p, icp.lm.z$residuals)
names(p)[length(names(p))] = 'Zresidual'

# SPATIAL PLOTS of the points
# ----------------------------------------------------------------
# spatial plot of matrix results
# use trans = 'x_trans' or 'y_trans' for other directions
trans = 'Zresidual'
n = abs(round((max(p[[trans]]) - min(p[[trans]]))/11, 4))
cuts.p = seq(min(p[[trans]]), max(p[[trans]]), n)
cuts.p = round(cuts.p,4)
spplot(p, zcol = trans, cuts = cuts.p, key.space = 'right', cex = 0.5)

f = gsub('.csv', '.shp', basename(points.f))
writeOGR(p, paste("H:/AFRF_ICP/ICP_points/shp/", 'residuals_', f, sep = ''), layer = f, driver = "ESRI Shapefile", overwrite_layer = T)
# write.csv(p@data,paste("H:/AFRF_ICP/ICP_points/shp/", 'residuals_', basename(points.f), sep = ''))



# PLOTTING 1 to 1 for observed vs predicted manually
# --------------------------------------------------------------------
# library developed by P. Tompalski
library(UsefulRFunctions)

r_square <- function(obs,pred) {
  resid <- obs - pred
  1 - (var(resid, na.rm = TRUE) / var(obs, na.rm = TRUE))
}

bw<-UsefulRFunctions:::bw
scatter2 <- function (x, y, R2=T, by = NULL, axisorder = "PO", xlab = "Observed", 
                      ylab = "Predicted", title = NULL, info = T, position = 0, 
                      positionauto = T, lowerlimit = NA, upperlimit = NA, alpha = 1, normality=T,
                      add.reg.line = F, rug = F, label_text = c("n", "bias", "bias%", 
                                                                "RMSE", "RMSE%","p-value")) 
{
  if (!is.null(by)) {
    data <- data.frame(x = x, y = y, by = by)
    pts <- ggplot2::geom_point(shape = 1, size = 2, alpha = alpha, 
                               ggplot2::aes(colour = by))
  }
  else {
    data <- data.frame(x = x, y = y)
    pts <- ggplot2::geom_point(shape = 1, size = 2, alpha = alpha)
  }
  
  if(normality == T) {
    d <- UsefulRFunctions::calc.error(reference = data$x, estimate = data$y)
  } else {
    d <- UsefulRFunctions::calc.error(reference = data$x, estimate = data$y, dist.normal = F)
  }
  
  label <- paste(#label_text[1], " = ", d$count, "\n", 
    label_text[2], " = ", round(d$bias, 100), "\n", 
    label_text[3], " = ", round(d$bias_perc, 2), "\n", 
    label_text[4], " = ", round(d$RMSE, 3), "\n", 
    label_text[5], " = ", round(d$RMSE_perc, 2), "\n", 
    label_text[6], " = ", round(d$p_val, 3),
    sep = "")
  
  
  if(R2 == T) {
    R2 <- r_square(obs = x, pred = y)
    label <- paste0("R2 = ", round(R2,2),
                    "\n",
                    label
    )
  } else if (R2=="cor") {
    R2 <- cor(x,y)^2
    label <- paste0("R2 = ", round(R2,2),
                    "\n",
                    label
    )
  } else if (is.numeric(R2)) {
    label <- paste0("R2 = ", round(R2,2),
                    "\n",
                    label
    )
  }
  
  # if (!is.null(R2)) {
  #   label <- paste0("R2 = ", round(R2,2),
  #                  "\n",
  #                  label
  #                  )
  # }
  # 
  
  
  if (axisorder == "OP") {
    data <- data.frame(x = data$y, y = data$x)
    x_lab_copy <- xlab
    xlab = ylab
    ylab = x_lab_copy
  }
  if (is.na(lowerlimit)) 
    lowerlimit <- min(data[c("x", "y")], na.rm = T)
  if (is.na(upperlimit)) 
    upperlimit <- max(data[c("x", "y")], na.rm = T)
  if (position != 0) 
    positionauto <- F
  if (positionauto == T) {
    if (is.finite(d$bias_perc) & d$bias_perc < -20) 
      position <- 1
  }
  if (position == 0) {
    ann_x <- upperlimit
    ann_y <- -Inf
    ann_hjust <- 1
    ann_vjust <- -0.2
  }
  if (position == 1) {
    ann_x <- lowerlimit
    ann_y <- upperlimit
    ann_hjust <- 0
    ann_vjust <- 0.9
  }
  if (info == T) {
    ann <- ggplot2::annotate("text", x = ann_x, y = ann_y, 
                             label = label, hjust = ann_hjust, vjust = ann_vjust)
  }
  else {
    ann <- bw
  }
  if (rug == T) {
    addrug <- ggplot2::geom_rug(alpha = 0.2)
  }
  else {
    addrug <- bw
  }
  if (add.reg.line == T) {
    reg.line <- ggplot2::geom_smooth(se = FALSE, method = "lm", 
                                     colour = "red")
  }
  else {
    reg.line <- bw
  }
  plot <- ggplot2::ggplot(data = data, ggplot2::aes(x = x, 
                                                    y = y)) + pts + ggplot2::xlab(xlab) + ggplot2::ylab(ylab) + 
    ggplot2::ggtitle(title) + ggplot2::xlim(lowerlimit, upperlimit) + 
    ggplot2::ylim(lowerlimit, upperlimit) + ggplot2::geom_abline(intercept = 0, 
                                                                 slope = 1) + ann + ggplot2::theme(legend.position = "bottom") + 
    ggplot2::coord_equal(ratio = 1) + addrug + reg.line + 
    UsefulRFunctions:::bw
  return(plot)
}

# example of using the 1 to 1 scatter plot
scatter2(x = z.offset, y = icp.lm.z$fitted.values)



# # arrange plots of the offsets in each direction (6 plots)
# library(ggpubr)
# 
# ggp.list = list()
# 
# ggp.offset = function(offset, d, label){
#   ggp = ggplot()
#   ggp = ggp + geom_point(aes(x = d, y = offset), alpha = 0.1)
#   ggp = ggp + geom_smooth(aes(x = d, y = offset))
#   ggp = ggp + ggtitle(label)
#   ggp
# }
# 
# ggp.list[[1]] = ggp.offset(offset = p$x_trans, d = p$x, label = 'X_shift vs UTM X')
# ggp.list[[2]] = ggp.offset(offset = p$y_trans, d = p$x, label = 'Y_shift vs UTM X')
# ggp.list[[3]] = ggp.offset(offset = p$z_trans, d = p$x, label = 'Z_shift vs UTM X')
# 
# ggp.list[[4]] = ggp.offset(offset = p$x_trans, d = p$y, label = 'X_shift vs UTM Y')
# ggp.list[[5]] = ggp.offset(offset = p$y_trans, d = p$y, label = 'Y_shift vs UTM Y')
# ggp.list[[6]] = ggp.offset(offset = p$z_trans, d = p$y, label = 'Z_shift vs UTM Y')
# 
# ggpubr::ggarrange(plotlist = ggp.list, ncol = 3, nrow = 2)

# ----------------------------------------------------------------
# IMPLEMENT THE MODELS TO UNWARP ORIGINAL CLOUDS
# this sections uses the LASCatalog processing engine to apply the models to every point of the original point cloud



# path to original clouds
# laslist = list.files("D:/JOE_RAKOFSKY/ICP_tempdir_canopy", pattern = glob2rx('DAP_?????.las'), full.names = T)
# laslist = sample(laslist, 200)

catalog_apply_shift = function(cluster, icp.lm.z){
  las = readLAS(cluster)
  
  if (is.empty(las)) return(NULL)

  # dataframe containing the x and y coords of the original cloud
  g = data.frame(list(las$X, las$Y))
  names(g) = c('x', 'y')
  
  
  # predict the shifts
  # xshift = stats::predict(icp.lm.x, newdata = g)
  # yshift = stats::predict(icp.lm.y, newdata = g)
  zshift = stats::predict(icp.lm.z, newdata = g)
  
  # apply the shift values 
  # las$X = las$X + xshift
  # las$Y = las$Y + yshift
  las$Z = las$Z + zshift
  
  # riteLAS(las)
  return(las)

}

test.list = list.files("D:/JOE_RAKOFSKY/DAP_raw_rembuf_sl_subset", full.names = T)
test.list = test.list[1:50]

# ctg = catalog("D:/JOE_RAKOFSKY/DAP_raw_rembuf_sl_subset")
ctg = catalog(test.list)

# use laz output option
opt_laz_compression(ctg) <- TRUE

# number of cores
# opt_cores(ctg) = 8

# path and name for unwarped outputs
opt_output_files(ctg) <- paste0("D:/JOE_RAKOFSKY/mask_layers/masked_unwarp_results/", "masked_shifted_{ID}")
# run the unwarping process

# 3. Set some catalog options.
# For this dummy example, the chunk size is 80 m and the buffer is 10 m using a single core.
opt_chunk_buffer(ctg) <- 10
opt_cores(ctg)        <- 1L
opt_chunk_size(ctg)   <- 250            # small because this is a dummy example.
opt_select(ctg)       <- "xyz"         # read only the coordinates.
# opt_filter(ctg)       <- "-keep_first" # read only first returns.

# 4. Apply a user-defined function to take advantage of the internal engine
opt <- list(need_buffer = TRUE)   # catalog_apply will throw an error if buffer = 0


tic('catalog_apply unwarping cloud')

new_ctg = catalog_apply(ctg, catalog_apply_shift, icp.lm.z = icp.lm.z)

# measure the time of completion for the process
toc()



# -----------------------------------------------------------
# an additional way to apply the model to shift the original point cloud
# without using LAScatalog.  LASCatalog will fail if the input point cloud
# is not continuous coverage, ie tiles with large spaces between tiles.


model_shift = function(i)
{
  lasfile = laslist[i]
  lasfilename = paste0("D:/JOE_RAKOFSKY/mask_layers/masked_unwarp_results/", gsub(".las", '', basename(lasfile)), '_shifted.laz')
  if (!file.exists(lasfilename)){
    print('loading las...')
    las = readLAS(lasfile)
    print('finished loading las.')
    
    # dataframe containing the x and y coords of the original cloud
    g = data.frame(list(las$X, las$Y))
    names(g) = c('x', 'y')
    
  
    # predict the shifts
    # xshift = stats::predict(icp.lm.x, newdata = g)
    # yshift = stats::predict(icp.lm.y, newdata = g)
    zshift = stats::predict(icp.lm.z, newdata = g)
    
    # apply the shift values 
    # las$X = las$X + xshift
    # las$Y = las$Y + yshift
    las$Z = las$Z + zshift
    print('writing shifted las...')
    writeLAS(las, lasfilename)
  }else{
    print('file exists...')
  }
}

# path to original clouds
laslist = list.files("D:/JOE_RAKOFSKY/DAP_raw_rembuf_sl_subset", pattern = '.las', full.names = T)
# laslist = laslist[1:100]

# -----------------------------------------------
# yet another way to apply the shifting according to the model using a simple for loop

# simple loop to carry out model applied shifting
# laslist = list.files("D:/JOE_RAKOFSKY/DAP_raw_rembuf_sl_subset", pattern = '.las', full.names = T)

for (i in 1:length(laslist)){
  print(i / length(laslist))
  model_shift(i)
}



# -----------------------------------------------------------------
# this remaining portion is a work in progress


# path to original clouds
laslist = list.files("D:/JOE_RAKOFSKY/DAP_raw", pattern = ".laz", full.names = T)
laslist = sample(laslist, 20)


als = lidR::catalog("D:/JOE_RAKOFSKY/ALS_raw")
als.e = extent(als)
als.e = as(als.e, 'SpatialPolygons')
for (lasfile in laslist){
  las = readLAS(lasfile)
  ext = las@bbox
  las.e = extent(las)
  las.e = as(las.e, 'SpatialPolygons')
 
  if (length(rgeos::gIntersection(las.e, als.e)) != 0){
    xleft = ext[1,1]
    ybot = ext[2,1]
    xright = ext[1,2]
    ytop = ext[2,2]
    als.clip = lasclipRectangle(als, xleft = xleft, ybottom = ybot, xright = xright, ytop = ytop)
    
    als.clip.f = paste0("D:/JOE_RAKOFSKY/DAP_unwarp_test/step_500_win_30_canopyonly_TRUE_icp_obs_full_checks", '/', 'ALS_', basename(lasfile))
    
    writeLAS(als.clip, als.clip.f)
  }
  
}





# MAKING CHECK TRANSECTS
# This section serves to clip the reference cloud, original cloud and unwarped cloud
# to make check transects where we can visualize the results of the unwarping process

# transect parameters
t.count = 100
t.length = 200
t.width = 2

# random sample t.count number of sample points across the area
# the points will become transect centers
sample_pts = spsample(as.spatial(als), n = t.count, type = "random")
sample_pts = SpatialPointsDataFrame(sample_pts@coords, data = data.frame(sample_pts@coords))

# polygons list
t.polys = c()

# loop thru sample points and get the coords to make a transect line in a 
# random direction
for (i in 1:nrow(sample_pts)){
  
  # define the first node of transect as
  # the sample point
  node1_x = sample_pts[i,]@coords[1]
  node1_y = sample_pts[i,]@coords[2]
  node1 = c(node1_x, node1_y)
  
  # randomly pick an azimuth to develop the transect line
  az = sample(0:180, 1)
  az_rad = az*pi/180
  az_orth = az_rad + (pi/2)
  
  # define the second node of the transect
  node2_x = node1_x + t.length*sin(az_rad)
  node2_y = node1_y + t.length*cos(az_rad)
  node2 = c(node2_x, node2_y)
  
  node3_x = node1_x + t.width*sin(az_orth)
  node3_y = node1_y + t.width*cos(az_orth)
  node3 = c(node3_x, node3_y)
  
  node4_x = node2_x + t.width*sin(az_orth)
  node4_y = node2_y + t.width*cos(az_orth)
  node4 = c(node4_x, node4_y)
  
  # matrix of line coordinates
  t.rect.coords = matrix(c(node1, node2, node4, node3), nrow = 4, ncol =2, byrow=T)
  t.rect.ply = Polygon(t.rect.coords)
  
  t.polys = c(t.polys, t.rect.ply)
  
}

#list object of class Polygons
Polygons_list = c()
df = data.frame(obj_id = character(length(t.polys)),
                stringsAsFactors = F)

for(i in 1:length(t.polys)){
  df[i,] = c(as.character(i))
  ply = t.polys[[i]]
  P = Polygons(list(ply), ID = i)
  Polygons_list = c(Polygons_list, P)
}


transects.sppoly = SpatialPolygons(Polygons_list)

transects.spdf = SpatialPolygonsDataFrame(transects.sppoly, data = df)
crs(transects.spdf) = "+init=epsg:32611 +proj=utm +zone=11 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 "

f = "D:/JOE_RAKOFSKY/check_sample_shp/all_check_transects.shp"
writeOGR(transects.spdf, f, layer = 'all_transects', driver = 'ESRI Shapefile', overwrite_layer = T)


# make las trasnects
library(lidR)

dap = lidR::catalog("D:/JOE_RAKOFSKY/DAP_unwarp_test/step_500_win_30_canopyonly_TRUE_icp_obs_full")
als = lidR::catalog("D:/JOE_RAKOFSKY/ALS_raw")
dap.raw = lidR::catalog("D:/JOE_RAKOFSKY/DAP_raw")

for (i in 1:nrow(transects.spdf)){
  
  tr = transects.spdf[i,]
  
  int = rgeos::gIntersection(as.spatial(dap), tr)
  
  if (is.null(int) == F){
    tr.dap = lasclip(dap, tr)
    tr.als = lasclip(als, tr)
    tr.dap.raw = lasclip(dap.raw, tr)
    
    writeLAS(tr.dap, paste0('D:/JOE_RAKOFSKY/check_las/dap_', i, '.las'))
    writeLAS(tr.als, paste0('D:/JOE_RAKOFSKY/check_las/als_', i, '.las'))
    writeLAS(tr.dap.raw, paste0('D:/JOE_RAKOFSKY/check_las/dap_raw_', i, '.las'))
  }
 
}








