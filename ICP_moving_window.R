# R Project: ICP
# version control with Git

# This script performs a moving window operation where ICP point matching occurs in each 
# moving window
# In general, the output is a series of points with the 4x4 rotation/translation matrix of the ICP matching
# process as a row.  

# writing some bs for testing git purposes 

### install the necessary packages ###
# install.packages("lidR")
# install.packages("sp")
# install.packages("raster")
# install.packages("spatstat")
# install.packages("doParallel")
# install.packages("foreach")
# install.packages("concaveman")
# install.packages("mapview")
# install.packages("rgeos")

library(lidR)
library(sp)
library(raster)
library(spatstat)
library(doParallel)
library(foreach)
library(concaveman)
library(mapview)
library(rgeos)
library(rgdal)

# location of the CloudCompare.exe
CloudCompare = "C:/PROGRA~1/CloudCompare/CloudCompare.exe"

# project coordinate system (Proj4txt)
crs = "+proj=utm +zone=11 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "

# ------------------------------------------------------------------
# PARAMETERS (units = meters)
# STEP: distance between the center of 2 adjacent moving windows
# STEP parameter does not apply when using external points - see below
STEP = 150
# WIN_SIZE: width of square moving window
WIN_SIZE = 30

# canopy parameters
# CANOPY_ONLY: Boolean of whether to use only top of canopy points
CANOPY_ONLY = T
# SMALL_TILE_CORES: number of corse to use (LAScatalog) for picking out canopy points
SMALL_TILE_CORES = 4
# SMALL_TILE_WINDOW: window size or gridcell size to designate canopy points
SMALL_TILE_WINDOW = 2
# SMALL_TILE_BUFFER: optional buffer for small window size
SMALL_TILE_BUFFER = 0

# subsample and ICP parameters
# SUBSAMPLE_DISTANCE: minimum distance between points for point cloud to be processed in ICP 
SUBSAMPLE_DISTANCE = 0.5
# RANDOM_SAMPLE_LIMIT: limit of number of points to randomly sample for RMS calcuation for ICP at each iteration
RANDOM_SAMPLE_LIMIT = 25000
# ICP_OVERLAP: estimated overlap (%) between the 2 point clouds in each window after ICP
# should be less than 100 since some small shifting in X and Y is anticipated, large X and Y shifts should not generally be anticipated
ICP_OVERLAP = 97

# global shift
# GLOBAL_SHIFT_AUTO: if true, let cloudcompare define local coordinates for processing
# cloudcompare will preserve the absolute coordinates for any point cloud results exported
GLOBAL_SHIFT_AUTO = T
# Alternatively, set a manual global shift for local coordinates
# cloudcomapre does not like having large UTM coordinates with 5+ digits to work with - use caution if specifying manual shift
# for example: if the rough coordinates are X: 500,100 and y: 6,000,100, use GLOBAL_SHIFT_X = 500000 and GLOBAL_SHIFT_Y = 6000000
GLOBAL_SHIFT_X = 0
GLOBAL_SHIFT_Y = 0
GLOBAL_SHIFT_Z = 0

# parallel clusters to run
# N_CLUSTERS: number of processes to run in parallel using the doParallel and foreach packages
N_CLUSTERS = 8

# start row for points
# START_ROW: if processing from scratch, should be = 1
# otherwise if you are picking up a partially processed dataset, can specify a row to begin from
# Existing rows in the csv will skipped
START_ROW = 26670

# INPUTS
# DAP cloud, or cloud 'to be moved'
# DAP_FULL.path = "H:/MKRF_ICP/DAP_Raw/GC"
# DAP_FULL.path = as.vector(list.files("D:/JOE_RAKOFSKY/DAP/raw", pattern = '.laz', full.names = T)[1:1000])
DAP_FULL.path = "D:/JOE_RAKOFSKY/DAP_raw_rembuf"
DAP_FULL = lidR::catalog(DAP_FULL.path)

# ALS cloud, or reference cloud
# ALS_FULL.path = "H:/MKRF_ICP/ALS_Raw/GC"
# ALS_FULL.path = as.vector(list.files("D:/JOE_RAKOFSKY/ALS/raw", pattern = '.laz', full.names = T)[1:1000])
ALS_FULL.path = "D:/JOE_RAKOFSKY/ALS/raw"
ALS_FULL = lidR::catalog(ALS_FULL.path)

# OUTPUTS
# this the folder in which a number of ICP resultant files are written for each iteration of the moving window
# potentially written files are:
# 1. Original aligned and reference clouds clipped to each moving window bounding box
# 2. The new aligned clouds for each moving window box
# 3. Registeration matrix text file produced from the ICP procedure
ICP_OUTPUT_DIR = paste0('D:/JOE_RAKOFSKY/ICP_temp', "_step_", STEP, '_win_', WIN_SIZE, '_canopy_', CANOPY_ONLY)
if (!dir.exists(ICP_OUTPUT_DIR)){
  dir.create(ICP_OUTPUT_DIR)
}

# this is a predetermined set of points of interest where we want to perform the ICP estimations
# For example: points along known roads, or areas where harvest has occurred 
# use the T/F switch here to specify whether we are using an external points file
EXTERNAL_POINTS = F
EXTERNAL_POINTS_FILE = 'D:/JOE_RAKOFSKY/ArcGIS/Slave_Lake_Roads_Harvest/roads_and_harvest/rand_pts_roads.shp'

EXTERNAL_BOUNDARY = T
EXTERNAL_BOUNDARY_FILE = "D:/JOE_RAKOFSKY/SL_subset/SL_subset.shp"

# set a location for ICP observaion points output (csv file)
points_dir = "D:/JOE_RAKOFSKY/ICP_points"
# define the file to which icp observations are written
if (EXTERNAL_POINTS == F){
  icp_points = paste(points_dir, '/', "step_", STEP, "_win_", WIN_SIZE, "_canopy_", CANOPY_ONLY, "_icp_obs.csv", sep = '')
  
  # define the full extent objects for the las objects
  # potential here to use las boundary or convex hull of points.... 
  EXTENT_DAP = raster::extent(DAP_FULL)
  EXTENT_ALS = raster::extent(ALS_FULL)
  
  # Start the process at the top left corner of extent of DAP or 'to be moved' point cloud
  START_COORD = c(EXTENT_DAP@xmin, EXTENT_DAP@ymax)
  
  # create the grid of points, start with the bounding box of the 'to be moved' cloud
  # this is essentially the same as just using the EXTENT_DAP object
  Proj_bbox = owin(xrange = c(EXTENT_DAP@xmin, EXTENT_DAP@xmax),
                   yrange = c(EXTENT_DAP@ymin, EXTENT_DAP@ymax))
  
  # define the number of points in the grid - for grid_centers
  # divide by the STEP parameter to create the number of points needed
  N_PTS_x = (EXTENT_DAP@xmax - EXTENT_DAP@xmin) / STEP
  N_PTS_y = (EXTENT_DAP@ymax - EXTENT_DAP@ymin) / STEP
  
  # create the grid centers object
  Proj_gridpts = data.frame(gridcenters(Proj_bbox, nx = round(N_PTS_x), ny = round(N_PTS_y)))
  Proj_gridpts = SpatialPointsDataFrame(Proj_gridpts[,1:2], data = Proj_gridpts)
  crs(Proj_gridpts) = crs
}
if (EXTERNAL_POINTS == T){
  ext_fstub = gsub('.shp', '', basename(EXTERNAL_POINTS_FILE))
  icp_points = paste(points_dir, '/', ext_fstub, "_win_", WIN_SIZE, "_canopy_", CANOPY_ONLY, "_icp_obs.csv", sep = '')
  
  Proj_gridpts = readOGR(EXTERNAL_POINTS_FILE)
  crs(Proj_gridpts) = crs
  Proj_gridpts$CID = seq(1, nrow(Proj_gridpts), 1)
}

# Write the header line
# r1c1, r1c2, etc. are the values of the 3 by 4 transforamtion matrix resulting from each ICP instance in the moving window operation
# where r refers to row and c refers to column
# check if the file exists first as to not overwrite
if (!file.exists(icp_points)){
  writeLines(c("n, isnull, coverage, ICP, x, y, RMS, r1c1, r2c1, r3c1, r1c2, r2c2, r3c2, r1c3, r2c3, r3c3, r1c4, r2c4, r3c4,"), icp_points)
}

if(EXTERNAL_BOUNDARY == T){
  ext_bound = readOGR(EXTERNAL_BOUNDARY_FILE)
  ext_bound = rgeos::gBuffer(ext_bound, width = 1000)
  crs(ext_bound) = crs
  Proj_gridpts = rgeos::gIntersection(Proj_gridpts, ext_bound)
  Proj_gridpts = Proj_gridpts@coords
}
if(EXTERNAL_BOUNDARY == F){
  # rough estimate of convex hull of the DAP points
  # this reduces the number of potential ICP runs by limiting the initial extent of the project closer to the boundary of 
  # the DAP points as opposed to the bounding box
  # take systematic sample (regular grid) or 10000 points
  sample_pts = spsample(as.spatial(DAP_FULL), n = 10000, type = "regular")
  # get the concave boundary of the points
  # higher concavity argument value means more complex polygon
  c = concaveman(sample_pts@coords, concavity = 10000, length_threshold = 5)
  c.poly = as.data.frame(c)
  # give all points an ID vlalue of 1 for single polygon
  c.poly$ID = 1
  # convert the outer points to polygon
  c.poly = coords2Polygons(c[,1:2], ID = 'ID')
  # buffer the polygon by the step width to ensure full coverage of the DAP
  c.poly = buffer(c.poly, width = STEP)
  # assign the coordinate system
  crs(c.poly) = crs
  
  # return only the points which fall inside the convex hull estimate polygon
  Proj_gridpts = rgeos::gIntersection(Proj_gridpts, c.poly)
  Proj_gridpts = Proj_gridpts@coords
}

# make function for writing output to log file in for each parallel
# this allows us to keep adding to the csv while running in parallel with the 'foreach' package
catf <- function(..., file=icp_points, append=TRUE){
  cat(..., file=file, append=append)
}

# This function performs the moving window operation and checks for 
# extent matches between the two clipped point clouds at each window and decides whether or not
# to use the Run_ICP function
Moving_Window_ICP = function(i){
  # FOR TESTING
  # i = sample(1:nrow(Proj_gridpts), 1)
  # i = 215
  
  icp.dt = read.csv(icp_points, header=T, stringsAsFactors = F)
  rowcheck = subset(icp.dt, n == i)
  skiprow = F
  if (nrow(rowcheck) > 0 ){
    skiprow = T
  }
  
  
  if (skiprow == F){
    # define the corners of the current ICP box
    # this box is used to clip the 2 point clouds to be matched
    
    # get the current point as the center of the square
    pxy = as.data.frame(Proj_gridpts)[i,]
    
    # get the corners in reference to the center
    xmin = pxy$x - (WIN_SIZE / 2)
    xmax = pxy$x + (WIN_SIZE / 2)
    ymin = pxy$y - (WIN_SIZE / 2)
    ymax = pxy$y + (WIN_SIZE / 2)
    
    # define the box which we will carry out the ICP operation as spatial object
    ICP_box = as(raster::extent(xmin, xmax, ymin, ymax), "SpatialPolygons")
    proj4string(ICP_box) = crs
    
    # clip both clouds to the ICP box
    DAP_TILE = lidR::lasclipRectangle(DAP_FULL, xmin, ymin, xmax, ymax)
    ALS_TILE = lidR::lasclipRectangle(ALS_FULL, xmin, ymin, xmax, ymax)
    
    DAP_TILE.f = paste(ICP_OUTPUT_DIR, '/DAP_', i, '.las', sep = '')
    ALS_TILE.f = paste(ICP_OUTPUT_DIR, '/ALS_', i, '.las', sep = '')
    # write the files to our temp directory
    if (nrow(DAP_TILE@data) > 0 & nrow(ALS_TILE@data) > 0){
      writeLAS(DAP_TILE, DAP_TILE.f)
      writeLAS(ALS_TILE, ALS_TILE.f)
    }
    
    # start with assumption we have a valid tile
    NULL_TILE = F
    # start with assumption that we will run ICP
    ICP = T
    # Root mean square value
    RMS = NA
    # start with blank char obj which will hold the matrix converted to a row...
    # recall the header with r1c1, r2c2 etc...
    f.row.char = ''
    
    # switch for if the clip object returns NULL
    if (is.null(DAP_TILE) | is.null(ALS_TILE)){
      NULL_TILE = T
      ICP = F
    }
    
    # switch for point count thresholds are met
    # we assume that the DAP cloud will need to have at least 1000 points 
    # and the ALS cloud to have  100 points in order to have a valid ICP process
    # These may be subject to change according to the window size...
    if (nrow(DAP_TILE@data) < 100 | nrow(ALS_TILE@data) < 100){
      NULL_TILE = T
      ICP = F
    }
    
    # switch for coverage of DAP and ALS
    # this variable will hold the fractional area overlap between the two clouds
    DAP.ALS.frac_area = NA
    if (NULL_TILE == F){
      DAP_xy = cbind(DAP_TILE@data$X, DAP_TILE@data$Y)
      ALS_xy = cbind(ALS_TILE@data$X, ALS_TILE@data$Y)
      
      # 2D convex hulls (lasboundary) get the perimeter of each clipped cloud
      DAP.chull = chull(DAP_xy)
      ALS.chull = chull(ALS_xy)
      
      # corrdinate object of the conex hulls
      DAP.chull.coords = DAP_xy[c(DAP.chull, DAP.chull[1]), ]
      ALS.chull.coords = ALS_xy[c(ALS.chull, ALS.chull[1]), ]
      # polygon object of convex hulls
      DAP.chull.sp_poly = SpatialPolygons(list(Polygons(list(Polygon(DAP.chull.coords)), ID=1)))
      ALS.chull.sp_poly = SpatialPolygons(list(Polygons(list(Polygon(ALS.chull.coords)), ID=1)))
      # fractional area of DAP/ALS convex hulls
      DAP.ALS.frac_area = raster::area(DAP.chull.sp_poly) / raster::area(ALS.chull.sp_poly)
      
      # if the areas of the convex hulls are not the very similar
      # ie. there is too much discrepancy between the 2D coverage between two clipped clouds..
      if (DAP.ALS.frac_area < 0.95 | DAP.ALS.frac_area > 1.05){
        NULL_TILE = T
        ICP = F
      } 
    } # end switch for coverages
    
    # start ICP switch
    # 
    if (NULL_TILE == F){
      # set the ICP boolean to true because we are now cleared to run ICP
      ICP = T
      # define file names for the clipped point clouds
      DAP_TILE.f = paste(ICP_OUTPUT_DIR, '/DAP_', i, '.las', sep = '')
      ALS_TILE.f = paste(ICP_OUTPUT_DIR, '/ALS_', i, '.las', sep = '')
      # write the files to our temp directory
      # writeLAS(DAP_TILE, DAP_TILE.f)
      # writeLAS(ALS_TILE, ALS_TILE.f)
      # get the matrix as a row from the Run_ICP function
      f.row.char = Run_ICP(DAP_TILE.f, ALS_TILE.f, i)
      
      # remove the raw tiles
      # this can be used to save space
      # file.remove(DAP_TILE.f)
      # file.remove(ALS_TILE.f)
      
    }
    # this writes the full row defined in the header including the matrix values and the location of the ICP moving window center
    catf(log_row = paste(i,
                         NULL_TILE,
                         round(DAP.ALS.frac_area,4),
                         ICP,
                         # X and Y location of the center of the ICP curent moving window box
                         pxy[1,1],
                         pxy[1,2],
                         # ICP matrix as row of data
                         f.row.char,
                         # skip to a new line for next row
                         "\n",
                         sep = ','))
    
  } # end of skiprow switch
} # end function Moving_window_ICP

# This functions creates and fires the call to CloudComapre ICP using system()
# This function is nested in the Moving_Window_ICP function defined below
# The iterator i is inhereted from the Moving_Window_ICP function
Run_ICP = function(DAP_TILE.f, ALS_TILE.f, i){
  
  # global shift values for local coords in CloudCompare
  # global shifting can also be left as 'AUTO' where cloudcomapre will estimate the appropriate local 3D coordinate origins
  # seee the CloudComapre commmand line wiki for more details https://www.cloudcompare.org/doc/wiki/index.php?title=Command_line_mode
  if (GLOBAL_SHIFT_AUTO == F){
    GS_X = GLOBAL_SHIFT_X
    GS_Y = GLOBAL_SHIFT_Y
    GS_Z = GLOBAL_SHIFT_Z
  }
  if (GLOBAL_SHIFT_AUTO == T){
    GS_X = 'AUTO'
    GS_Y = ''
    GS_Z = ''
  }
  
  # spatial distance for subsampling
  # this is the minimum distance (m) between points after spatial subsampling
  SP = SUBSAMPLE_DISTANCE
  
  # random sample limit for ICP
  # number of points used in RMS calcualtion from the 'to be moved' cloud
  RSL = RANDOM_SAMPLE_LIMIT
  
  # estimated ICP 2D overlap of the 2 point clouds in %
  overlap = ICP_OVERLAP
  
  # build the command call first as a list object
  # the initial call 'CloudCompare' should work if it is set in Evnironment Variables 
  # This command is to spatially subsample the 'to be moved' cloud
  ss_command = c(CloudCompare,
                 '-SILENT',
                 '-C_EXPORT_FMT',
                 'LAS',
                 '-o',
                 '-GLOBAL_SHIFT',
                 GS_X,
                 GS_Y,
                 GS_Z,
                 DAP_TILE.f,
                 '-SS',
                 'SPATIAL',
                 SP)
  
  # paste the command together as a single char obj with spaces
  ss_command = paste(ss_command, collapse = ' ')
  # invoke the ss command
  system(ss_command)
  
  
  ss_command = c(CloudCompare,
                 '-SILENT',
                 '-C_EXPORT_FMT',
                 'LAS',
                 '-o',
                 '-GLOBAL_SHIFT',
                 GS_X,
                 GS_Y,
                 GS_Z,
                 ALS_TILE.f,
                 '-SS',
                 'SPATIAL',
                 SP)
  
  # paste the command together as a single char obj with spaces
  ss_command = paste(ss_command, collapse = ' ')
  # invoke the ss command
  system(ss_command)
  
  
  # ====================================================================
  # gather the newly subsampled point cloud
  # cloudcomapre will always append '_SPATIAL_SUBSAMPLE' to the output filename so search for it.
  DAP_TILE_SS = list.files(ICP_OUTPUT_DIR, pattern = glob2rx(paste('DAP_', i, '_SPATIAL_SUBSAMPLED*.las', sep='')))
  DAP_TILE_SS = paste(ICP_OUTPUT_DIR, DAP_TILE_SS, sep = '/')
  
  ALS_TILE_SS = list.files(ICP_OUTPUT_DIR, pattern = glob2rx(paste('ALS_', i, '_SPATIAL_SUBSAMPLED*.las', sep='')))
  ALS_TILE_SS = paste(ICP_OUTPUT_DIR, ALS_TILE_SS, sep = '/')
  
  ALS_TILE = ALS_TILE.f
  
  # make empty string for cluster directory
  DAP_cdir = ''
  # if not using only canopy points, the DAP_cdir is the ICP_OUTPUT_DIR
  if (CANOPY_ONLY == F){
    DAP_cdir = ICP_OUTPUT_DIR
  }
  
  if (CANOPY_ONLY == T){
    # create just canopy top points
    # DAP
    # make folders where tiled clouds reside
    DAP_cdir = paste(ICP_OUTPUT_DIR, "/cluster_DAP_", i, sep = '')
    ALS_cdir = paste(ICP_OUTPUT_DIR, "/cluster_ALS_", i, sep = '')
    
    # make into catalog object
    DAP_TILE.canopy = catalog(DAP_TILE_SS)
    
    # small tiles width
    opt_chunk_size(DAP_TILE.canopy) = SMALL_TILE_WINDOW
    # small tiles buffer width
    opt_chunk_buffer(DAP_TILE.canopy) = SMALL_TILE_BUFFER
    # cores to use for tiling procedure
    opt_cores(DAP_TILE.canopy) = SMALL_TILE_CORES
    # set output location for small tiles as inside the folders created above
    opt_output_files(DAP_TILE.canopy) = paste(DAP_cdir, "/part_{ID}", sep ='')
    
    if (!dir.exists(DAP_cdir)){
      dir.create(DAP_cdir)
    }
    
    # create the small tiles and store as object
    newctg_DAP = catalog_retile(DAP_TILE.canopy)
    
    # gather the files taht are the small tiles
    small_cells_DAP = list.files(DAP_cdir, full.names = T)
    for (f in small_cells_DAP){
      las = readLAS(f)
      # get the highest point in each tile
      las = lasfilter(las, Z >= quantile(Z, 0.95))
      # there will be new files with each only one point
      writeLAS(las, f)
    }
    
    # merge the all the single point las files into one las file
    DAP_merge_command = c(CloudCompare,
                          '-SILENT',
                          '-C_EXPORT_FMT',
                          'LAS',
                          '-o',
                          paste(small_cells_DAP, collapse = ' -o '),
                          '-MERGE_CLOUDS'
    )
    
    # run the merge command
    DAP_merge_command = paste(DAP_merge_command, collapse = ' ')
    system(DAP_merge_command)
    
    # ALS
    # do the exact same process as above for DAP on ALS
    ALS_TILE.canopy = catalog(ALS_TILE.f)
    
    opt_chunk_size(ALS_TILE.canopy) = SMALL_TILE_WINDOW
    opt_chunk_buffer(ALS_TILE.canopy) = SMALL_TILE_BUFFER
    opt_cores(ALS_TILE.canopy) = SMALL_TILE_CORES
    opt_output_files(ALS_TILE.canopy) = paste(ALS_cdir, "/part_{ID}", sep ='')
    
    if (!dir.exists(ALS_cdir)){
      dir.create(ALS_cdir)
    }
    
    newctg_ALS = catalog_retile(ALS_TILE.canopy)
    
    small_cells_ALS = list.files(ALS_cdir, full.names = T)
    # small_cells_ALS = small_cells_ALS[1:100]
    for (f in small_cells_ALS){
      las = readLAS(f)
      las = lasfilter(las, Z >= quantile(Z, 0.95))
      writeLAS(las, f)
    }
    
    ALS_merge_command = c(CloudCompare,
                          '-SILENT',
                          '-C_EXPORT_FMT',
                          'LAS',
                          '-o',
                          paste(small_cells_ALS, collapse = ' -o '),
                          '-MERGE_CLOUDS'
    )
    
    ALS_merge_command = paste(ALS_merge_command, collapse = ' ')
    system(ALS_merge_command)
    
    # gather the merged sparse canopy top point files 
    DAP_TILE_SS = list.files(DAP_cdir, pattern = 'MERGE', full.names = T)[[1]]
    ALS_TILE = list.files(ALS_cdir, pattern = 'MERGE', full.names = T)[[1]]
  }
  
  # ====================================================================
  # Starting the ICP process
  # first define a log file which the result of ICP are written.  Each ICP process will generate a logfile and 
  # a registration matrix file
  cc_logfile = paste(ICP_OUTPUT_DIR, '/', i, '_ICP.log', sep='')
  
  # build the ICP command as a list
  icp_command = c(CloudCompare,
                  '-SILENT',
                  '-LOG_FILE', 
                  cc_logfile,
                  '-C_EXPORT_FMT',
                  'LAS',
                  '-o',
                  '-GLOBAL_SHIFT',
                  'AUTO',
                  # GS_X,
                  # GS_Y,
                  # GS_Z,
                  DAP_TILE_SS,
                  '-o',
                  '-GLOBAL_SHIFT',
                  'AUTO',
                  # GS_X,
                  # GS_Y,
                  # GS_Z,
                  ALS_TILE_SS,
                  '-ICP',
                  '-OVERLAP',
                  overlap,
                  '-RANDOM_SAMPLING_LIMIT',
                  RSL)
  
  # paste the command into a sinlge line
  icp_command = paste(icp_command, collapse = ' ')
  # invoke the ICP command
  system(icp_command)
  
  # we can choose to remove the intermediate subsampled DAP
  # file.remove(DAP_TILE_SS)
  
  
  # =======================================
  # Gather the resulting ICP metrics from 
  # the MATRIX file and the log file
  # note the the RMS value is found in the log file not the matrix file
  RMS = grep("RMS", readLines(cc_logfile), value = T)
  RMS = as.numeric(strsplit(RMS, 'RMS:')[[1]][2])
  
  # get the matrix values
  #======================================================
  
  # build matrix
  if (CANOPY_ONLY == T){
    matrix.f = list.files(DAP_cdir, pattern = glob2rx(paste('*REGISTRATION_MATRIX*.txt', sep = '')), full.names = T)
  }
  if (CANOPY_ONLY == F){
    id = strsplit(basename(DAP_TILE.f), '.las')[[1]][1]
    matrix.f = list.files(DAP_cdir, pattern = glob2rx(paste(id, '*REGISTRATION_MATRIX*.txt', sep = '')), full.names = T)
  }
  
  matrix.df = read.table(matrix.f)
  matrix.m = matrix.df[-4,]
  m = as.matrix(matrix.m)
  
  # CREATE LIST OF RETURNED METRICS
  f.row = as.vector(m)
  f.row = c(RMS, f.row)
  f.row = as.character(f.row)
  f.row.char = paste(f.row, collapse = ',')
  
  return(f.row.char)
}

# ==================================================
# MAIN

# make cluster of set number of workers
cl = makeCluster(N_CLUSTERS)
# register cluster - see documentation for foreach looping
registerDoParallel(cl)

# this is for measuring the time taken for the whole process
library(tictoc)
tic('ICP foreach')

# for all the points in the initially defined grid...
# make sure to pass the necessary packages to the clusters
# if 'start' = 1 then we start at the beginning of the points and go sequentially
# this depends on how the points are ordered in Proj_gridpts
# note: must specify the packages which need to be loaded to each cluster using .packages argument
foreach (i = START_ROW:nrow(Proj_gridpts),
         .packages = c("lidR",
                       "sp",
                       "raster",
                       "spatstat")) %dopar% {

                         Moving_Window_ICP(i)

                       }
# gives us the end time of whole process
toc()
# stop any rogue clusters
stopCluster(cl)

# for testing in non-parallel
# for (i in 1:nrow(Proj_gridpts)){
#   Moving_Window_ICP(i)
# }





