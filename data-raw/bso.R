
sc <- paste0("NETCDF:", cnr_dir_meta("~/data/nc/Phy_001-024") %>% filter(var == "salinity") %>% .$fp)

# idx_deepest <- function(x)
#   max(which(!is.nan(global(rast(sc[x]), "max", na.rm = TRUE)$max)))
#
# picks <- 1:length(sc) %>% map(idx_deepest)
#
# lyr_depth <- function(ifile, ilayer) {
#   readr::parse_double(gsub(".*?=(.*?)$", "\\1", names(rast(sc[ifile]))[[ilayer]]))
# }
#
# 1:length(picks) %>% map_dbl(function(x) lyr_depth(x, picks[[x]]))

# it appears the deepest depth layer across all of the files are layer 49 in each file
# it is measured at 5274.784
n <- 10
z <- rast(nrows = 180 * n, ncols = 360 * n)
#values(z) <- 1:ncell(z)

myrr <- app(rast(sc[1]), fun = function(x) max(which(!is.na(x))))
plot(myrr)

myrr2 <- app(rast(sc[1]), fun = function(x) length(which(!is.na(x))))

myrr3 <- app(rast(sc[1]), fun = function(x) rast(sc[1])[[max(which(!is.na(x)))]])


# for one of the bottom salinity files
s <- rast(sc[1])






cnr_flatten_save(sc, "bottom_salinity", outdir = "/media/markus/Stash/data/nc")

fz <- dir("/media/markus/Stash/data/nc", pattern = "bottom_salinity*", full.names = TRUE)
md <- cnr_dir_meta("/media/markus/Stash/data/nc", pat = "bottom_salinity*")

myr <- rast(fz)

names(myr) <- sprintf("%s_%s", md$yr, md$mo)
time(myr) <- lubridate::as_datetime(lubridate::ymd(paste0(md$yr, "-", md$mo, "-", "01")))

prefix <- "bottom_salinity"
fn <- file.path("/media/markus/Stash/data/nc", sprintf("%s.nc", prefix))
writeCDF(myr, filename = fn, overwrite = TRUE, varname = varnames(myr), longname = longnames(myr), unit = units(myr))

myr

# resample to tdc


n <- 10
z <- rast(nrows = 180 * n, ncols = 360 * n)

# 82 s
tictoc::tic()
r <- resample(myr, z, "bilinear")
names(r) <- names(myr)
time(r) <- time(myr)
units(r) <- units(myr)
tictoc::toc()


prefix <- "bottom_salinity"
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

longprefix <- "bottom_salinity"

myfiles <- file.path("/media/markus/Stash/data/nc",
  sprintf("%s_%s.nc", longprefix, "0p100"))

r <- rast(paste0("NETCDF:", myfiles))

prefix <- "bso"

names(r) <- paste0(prefix, year(time(r)), "_", sprintf("%02i", month(time(r))))

# write individual layers as monthly tiffs
fn <- sprintf(cnr_dir("tdc/months/%s_%s.tiff"), names(r), cnr_res(r))

writeRaster(r, fn, overwrite=TRUE, NAflag = -9999,
            gdal = c("COMPRESS=LZW", "TFW=YES", "of=COG","PROFILE=GeoTIFF"),
            datatype = "FLT8S")

# compute mean across years

years <- unique(year(time(r)))

tictoc::tic()

means <-
  years %>%
  purrr::map(function(x) cnr_annual_mean(r, x)) %>%
  do.call("c", .)

names(means) <- paste0(sep = "_", prefix, years)

time(means) <- as_datetime(ymd(paste0(years, "0101")))
varnames(means) <- unique(varnames(r))
longnames(means) <- unique(longnames(r))

tictoc::toc()

plot(means)

# output to CSV
means %>% to_df3() %>%
  write_csv(paste0("/media/markus/Stash/data/out/tdc/", prefix,
                   "_annuals_0p100.csv.gz", na = ""))

# output to TIFF
bso_means %>% writeRaster(filename =
  paste0("/media/markus/Stash/data/out/tdc/", prefix, "_annuals_0p100.tiff"),
  overwrite=TRUE, NAflag = -9999,
  gdal = c("COMPRESS=LZW", "TFW=YES", "of=COG","PROFILE=GeoTIFF"),
  datatype = "FLT8S"
)

# output yearly annuals as csv

tbs <- cnr_outer(prefix, years = years, nm = c("var", "year", "fn"))

tbs %>% pmap(.f = function(var, year, ...)
  cnr_save_ym(r, y = year, varname = var, outdir = cnr_dir("tdc"))
)

# convert those csv to parquet

gzfn <- dir(cnr_dir("tdc"), pattern = paste0(prefix, "_\\d{4}_0p100\\.csv"))
pfn <- gsub("\\.csv", "\\.parquet", gzfn)

cmd <- sprintf("duckdb :memory: \"copy '%s' to '%s' (format 'parquet');\"",
               gzfn, pfn)
writeLines(cmd, cnr_dir("tdc/bso_monthlies.sh"))
Sys.chmod(cnr_dir("tdc/bso_monthlies.sh"))

#' We can now export the yearly tables (with columns for the months) in .csv.gz-format to parquet
#' (this takes about 5 minutes using the shell script below)
readLines(cnr_dir("tdc/bso_monthlies.sh"))

#' and we can inspect the results
"bso_2017_0p100.parquet" %>%
  cnr_view_parquet()

#' We already have .csv.gz files with annual means.
#' We can convert all .csv.gz files to parquet files

gzfn <- dir(cnr_dir("tdc"), pattern = paste0(prefix, "_annuals.*?\\.csv.gz$"))
pfn <- gsub("\\.csv\\.gz", "\\.parquet", gzfn)

cmd <- sprintf("duckdb :memory: \"copy '%s' to '%s' (format 'parquet');\"",
               gzfn, pfn)
cat(cmd)
writeLines(cmd, cnr_dir("tdc/bso_annuals.sh"))
Sys.chmod(cnr_dir("tdc/bso_monthlies.sh"))

"bso_annuals_0p100.parquet" %>%
  cnr_view_parquet()






