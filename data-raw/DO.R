# these are the groups of non-default and not-specifed crs data
sc <- paste0("NETCDF:", cnr_dir_meta("~/data/nc/BIO_001_029") %>% filter(var == "DO") %>% .$fp)
sc2 <- paste0("NETCDF:", cnr_dir_meta("~/data/nc/BIO_001_028") %>% filter(var == "DO") %>% .$fp)

rsc <- rast(sc)
crs(rsc) <- "EPSG:4258"

z <- rast(ncol = ncol(rsc), nrow = nrow(rsc), crs = cnr_crs("default"))
ext(z) <- ext(rsc)
res(z) <- res(rsc)

# # NB takes approx 4 hours
# tic()
# rscwgs <- project(rsc, z)
# time(rscwgs) <- time(rsc)
# toc()
#
# # could be done on flattened results instead, and happen in 5 minutes?
# rscwgs

# # approx 4 hours
# tic()
# rscwgs <- project(rsc, z)
# time(rscwgs) <- time(rsc)
# toc()

# does it need transform?
#range is 168:455 while earlier was 0.009:437, so appears to be needed...?

# writeCDF(rscwgs, filename = "/media/markus/Stash/data/nc/rscwgs.nc",
#          overwrite = TRUE,
#          varname = unique(varnames(rsc)),
#          longname = unique(longnames(rsc)),
#          unit = unique(units(rsc))
# )

# approx 10 mins, restrict sc to not contain the 2019 files?
tic()
cnr_flatten_save(sc, reproject_from = "EPSG:4258", prefix = "DO2", outdir = "/media/markus/Stash/data/nc")
toc()

tic()
cnr_flatten_save(sc2, prefix = "DO22", outdir = "/media/markus/Stash/data/nc")
toc()

# NOTE: by setting the CRS to the default, no error are spewed when plotting!

# do these two layers (one w 75 depth layers, one with 50 depth layers) have the same flattened info?
# it doesn't really seem so:

m1  <- rast("NETCDF:/home/markus/data/nc/BIO_001_029/global-reanalysis-bio-001-029-monthly_DO_2019_01.nc")
plot(m1[[1]])
m2 <- rast("NETCDF:/home/markus/data/nc/BIO_001_028/global-analysis-forecast-bio-001-028-monthly_DO_2019_01.nc")
plot(m2[[1]])

m1 <- rast("NETCDF:/media/markus/Stash/data/nc/DO2_2019_01.nc")
m2 <- rast("NETCDF:/media/markus/Stash/data/nc/DO22_2019_01.nc")
plot(m2)
plot(m1[[1]])
plot(m2[[1]])

#crs(m2)<- cnr_crs("default")
#rm2 <-cnr_reproject(m2, cnr_crs("default"))
#range(crds(rm2)[2,])
#crs(rm2)
#plot()
#cnr_crs()
#plot(m1)
#plot(m2)
#plot(rast("NETCDF:/media/markus/Stash/data/nc/DO_2017_02.nc"))

# unify all reprojected flattened layers into a single nc file

# according to the instruction only 2017, 2018 data should be used from BIO_001_029 (which is in EPSG:4258)
# but which still has data for 2019

fz1 <- dir("/media/markus/Stash/data/nc", pattern = "DO2_+", full.names = TRUE) %>%
  grep("_2019_", ., invert = TRUE, value = TRUE)

md1 <- cnr_dir_meta("/media/markus/Stash/data/nc", pat = "DO2_+") %>%
  filter(!(yr %in% "2019"))

fz2 <- dir("/media/markus/Stash/data/nc", pattern = "DO22_+", full.names = TRUE)
md2 <- cnr_dir_meta("/media/markus/Stash/data/nc", pat = "DO22_+")

fz <- c(fz1, fz2)
md <- bind_rows(md1, md2)

#View(md)
myr <- rast(fz)

prefix <- "DO"
names(myr) <- sprintf("%s_%s_%s", prefix, md$yr, md$mo)
time(myr) <- lubridate::as_datetime(lubridate::ymd(paste0(md$yr, "-", md$mo, "-", "16")))
crs(myr) <- cnr_crs("default")
varnames(myr) <- rep(prefix, nlyr(myr))
longnames(myr) <- "Bottom dissolved oxygen"

# NB: this is not needed as flattening earlier already did reprojection,
# so no need to reproject only the flattened layers, approx 500 sec

# z <- rast(ncol = ncol(myr), nrow = nrow(myr), crs = cnr_crs("default"))
# ext(z) <- ext(myr)
# res(z) <- res(myr)

# tic()
# myrwgs <- project(myr, z)
# time(myrwgs) <- time(myr)
# toc()
#
# writeCDF(myrwgs, filename = fn, overwrite = TRUE,
#          varname = unique(varnames(myr)),
#          longname = unique(longnames(myr)),
#          unit = unique(units(myr))
# )
# plot(myrwgs[[1]])
# plot(myr[[1]])
# na.omit(values(myrwgs[[1]])) == na.omit(values(myr[[1]]))

prefix <- "flattened_DO"
fn <- file.path("/media/markus/Stash/data/nc", sprintf("%s.nc", prefix))

writeCDF(myr, filename = fn, overwrite = TRUE,
         varname = unique(varnames(myr)),
         longname = unique(longnames(myr)),
         unit = unique(units(myr))
)

time(myr)


# resample to tdc

n <- 10
z <- rast(nrows = 180 * n, ncols = 360 * n, crs = cnr_crs("default"))

# this is pretty fast, only takes 33 s and we have results in tdc
tictoc::tic()
r <- resample(myr, z, "bilinear")
names(r) <- names(myr)
time(r) <- time(myr)
units(r) <- units(myr)
tictoc::toc()

prefix <- "flattened_DO"
fn <- file.path("/media/markus/Stash/data/nc", sprintf("%s_%s.nc", prefix, cnr_res(r)))
writeCDF(r, filename = fn, overwrite = TRUE,
         varname = unique(varnames(myr)),
         longname = unique(longnames(myr)),
         unit = unique(units(myr)))

# do the annual means etc and export to csv and tiffs

library(dplyr)
library(terra)
library(readr)
library(lubridate)
library(purrr)

longprefix <- "flattened_DO"
prefix <- "DO"
myfiles <- file.path("/media/markus/Stash/data/nc", sprintf("%s_%s.nc", longprefix, "0p100"))
r <- rast(paste0("NETCDF:", myfiles))
names(r) <- paste0(prefix, "_", year(time(r)), "_", sprintf("%02i", month(time(r))))

# write individual layers as monthly tiffs
fn <- sprintf(cnr_dir("tdc/months/%s_%s.tiff"), names(r), cnr_res(r))

writeRaster(r, fn, overwrite=TRUE, NAflag = -9999,
            gdal = c("COMPRESS=LZW", "TFW=YES", "of=COG","PROFILE=GeoTIFF"),
            datatype = "FLT8S")


# compute mean across years

years <- unique(year(time(r)))

# compute mean across years
tictoc::tic()

means <-
  years %>%
  purrr::map(function(x) cnr_annual_mean(r, x)) %>%
  do.call("c", .)

names(means) <- paste(sep = "_", prefix, years)

time(means) <- as_datetime(ymd(paste0(years, "0101")))
varnames(means) <- unique(varnames(r))
longnames(means) <- unique(longnames(r))
#units(tdc) <- unique(units(r))

tictoc::toc()

# plot all annual means
plot(means)

# output to CSV
means %>% to_df3() %>%
  write_csv("/media/markus/Stash/data/out/tdc/DOx_annuals_0p100.csv.gz", na = "")

# output to TIFF
means %>% writeRaster(
  filename = "/media/markus/Stash/data/out/tdc/DO_annuals_0p100.tiff",
  overwrite=TRUE, NAflag = -9999,
  gdal = c("COMPRESS=LZW", "TFW=YES", "of=COG","PROFILE=GeoTIFF"),
  datatype = "FLT8S"
)

tbs <- cnr_outer("DO", years = years, nm = c("var", "year", "fn"))

tbs %>% pmap(.f = function(var, year, ...)
  cnr_save_ym(r, y = year, varname = var, outdir = cnr_dir("tdc"))
)

gzfn <- dir(cnr_dir("tdc"), pattern = paste0(prefix, "_\\d{4}_0p100\\.csv"))
pfn <- gsub("\\.csv", "\\.parquet", gzfn)

cmd <- sprintf("duckdb :memory: \"copy '%s' to '%s' (format 'parquet');\"",
               gzfn, pfn)
fn <- cnr_dir(paste0("tdc/", prefix, "_monthlies.sh"))
writeLines(cmd, fn)

#' We can now export the yearly tables (with columns for the months) in .csv.gz-format to parquet
#' (this takes about 5 minutes using the shell script below)
readLines(fn)

"DO_2017_0p100.parquet" %>%
  cnr_view_parquet()


#' We already have .csv.gz files with annual means.
#' We can convert all .csv.gz files to parquet files

gzfn <- dir(cnr_dir("tdc"), pattern = paste0(prefix, "_annuals.*?\\.csv.gz$"))
pfn <- gsub("\\.csv\\.gz", "\\.parquet", gzfn)

cmd <- sprintf("duckdb :memory: \"copy '%s' to '%s' (format 'parquet');\"",
               gzfn, pfn)
cat(cmd)
fn <- cnr_dir(paste0("tdc/", prefix, "_annuals.sh"))
writeLines(cmd, fn)

"DO_annuals_0p100.parquet" %>%
  cnr_view_parquet()
