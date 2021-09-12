library(terra)

# specify the nc file to use
fn <- paste0("NETCDF:", "/media/markus/Stash/data/nc/Bathymetry/GLO-MFC_001_024_mask_bathy.nc")
r <- rast(fn)

# inspect metadata
names(r)
varnames(r)

# inspect the depth layer
plot(r[["deptho"]])

# use the depth layer and resample to tdc
myr <- r[["deptho"]]

n <- 10
z <- rast(nrows = 180 * n, ncols = 360 * n)
ext(z) <- ext(-180, 180, -90, 90)

# 7 s
tictoc::tic()
r_tdc <- resample(myr, z, "bilinear")
tictoc::toc()

# output results to a TIFF file
myfn <- file.path(cnr_dir("tdc"), "deptho_0p100.tiff")

writeRaster(r_tdc, myfn, overwrite=TRUE, NAflag = -9999,
            gdal = c("COMPRESS=LZW", "TFW=YES", "of=COG","PROFILE=GeoTIFF"),
            datatype = "FLT8S")

# now compare tdc results to original results
plot(rast(myfn))
plot(myr)


# Use the same approach used when working with the GEBCO_2020 nc file

fn <- paste0("NETCDF:", "/media/markus/Stash/data/nc/Bathymetry/GEBCO_2020.nc")
r <- rast(fn)

# inspect metadata
names(r)
varnames(r)

# just one layer here...
myr <- r

n <- 10
z <- rast(nrows = 180 * n, ncols = 360 * n)
ext(z) <- ext(-180, 180, -90, 90)

# 500 s
tictoc::tic()
r_tdc <- resample(myr, z, "bilinear")
tictoc::toc()

# output results to a TIFF file
myfn <- file.path(cnr_dir("tdc"), "elevation_0p100.tiff")

writeRaster(r_tdc, myfn, overwrite=TRUE, NAflag = -9999,
            gdal = c("COMPRESS=LZW", "TFW=YES", "of=COG","PROFILE=GeoTIFF"),
            datatype = "FLT8S")

# compare tdc results to original results
plot(rast(myfn))
plot(myr)
