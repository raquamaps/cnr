library(purrr)

# harmonize nppv layers
# all nc files for nppv, which have "non-default" CRS
nc <-
  filez %>% filter(var == "nppv") %>%
  pull(fp) %>%
  #grep("extra", invert = TRUE, value = TRUE, .) %>%
  grep("global", value = TRUE, .)


crsx <- function(x) cnr_crs("nppv", x)
ncfz <- nc[map_chr(nc, crsx) != cnr_crs()["default"]]
nppv <- rast(paste0("NETCDF:", ncfz))

# set all these layers to have the same (non-default) CRS
crs(nppv) <- cnr_crs("nppv", ncfz[1])

z <- rast(ncol = ncol(nppv), nrow = nrow(nppv), crs = cnr_crs("default"))
ext(z) <- ext(-180, 180, -90, 90)
res(z) <- 0.25

tic()
edc <- project(nppv, z)
time(edc) <- time(nppv)
toc()

plot(edc)

writeCDF(edc, filename = "/home/markus/data/nc/Phy_001-024/extra-monthly_nppv2.nc",
         overwrite = TRUE,
         varname = unique(varnames(nppv)),
         longname = unique(longnames(nppv)),
         unit = unique(units(nppv))
)

# for the other nppv layers, project also (to the same extents)
tic()
morez <- rast(paste0("NETCDF:", nc[-c(1,2)]))
crs(morez) <- cnr_crs("default")
edc2 <- project(morez, z)
time(edc2) <- time(morez)
varnames(edc2) <- unique(varnames(morez))
longnames(edc2) <- unique(longnames(morez))
units(edc2) <- unique(units(morez))
toc()

# combine all the layers into one single nppv layer

myr <- c(
  edc,
  #  rast("NETCDF:/home/markus/data/nc/Phy_001-024/extra-monthly_nppv2.nc"),
  edc2
)

tic()
tdc <- cnr_resample(myr, type = "rast")
time(tdc) <- time(myr)
varnames(tdc) <- unique(varnames(myr))
longnames(tdc) <- unique(longnames(myr))
units(tdc) <- unique(units(myr))
toc()

# this layer can now be saved in tenth degree cell resolution

writeCDF(tdc, filename = "/home/markus/data/nc/Phy_001-024/monthly_nppv_0p100.nc",
         overwrite = TRUE,
         varname = unique(varnames(tdc)),
         longname = unique(longnames(tdc)),
         unit = unique(units(tdc))
)


# write tdc as monthly tiffs

myvars <- "nppv"
r_all_tdc <- rast(paste0("NETCDF:", "/home/markus/data/nc/Phy_001-024/monthly_nppv_0p100.nc"))
names(r_all_tdc) <- paste0("nppv_", year(time(r_all_tdc)), "_", sprintf("%02i", month(time(r_all_tdc))))

# write individual layers as monthly tiffs
fn <- sprintf(cnr_dir("tdc/months/%s_%s.tiff"), names(r_all_tdc), cnr_res(r_all_tdc))

writeRaster(r_all_tdc, fn, overwrite=TRUE, NAflag = -9999,
            gdal = c("COMPRESS=LZW", "TFW=YES", "of=COG","PROFILE=GeoTIFF"),
            datatype = "FLT8S")

# we could write yearly tiffs, too... but would expect those to have annuals...
#idx <- intersect(which(year(time(r_all_tdc)) == 2017), grep("^so", names(r_all_tdc)))
#fn <- sprintf(cnr_dir("tdc/%s_%s.tiff"), "so_2017", cnr_res(r_all_tdc[[idx]]))
#names(r_all_tdc[[idx]])
#writeRaster(r_all_tdc[[idx]], fn, overwrite=TRUE, NAflag = -9999,
#            gdal = c("COMPRESS=LZW", "TFW=YES", "of=COG","PROFILE=GeoTIFF"),
#            datatype = "FLT8S")

# needs to be done for each variable
r <- r_all_tdc

for (myvar in myvars) {
  r2 <- r[[grep(paste0("^", myvar), names(r))]]
  years <- unique(year(time(r)))
  years <- years[which(years != 2016)]
  means <-
    years %>%
    purrr::map(function(x) cnr_annual_mean(r2, x)) %>%
    do.call("c", .)

  names(means) <- paste(sep = "_", myvar, years)

  fn <- sprintf(cnr_dir("tdc/%s_annuals_%s.tiff"),
                myvar, cnr_res(means))
  writeRaster(means, fn, overwrite=TRUE, NAflag = -9999,
              gdal = c("COMPRESS=LZW", "TFW=YES", "of=COG","PROFILE=GeoTIFF"),
              datatype = "FLT8S")

  means %>% cnr_save_tiff(outfile = fn)

  fn <-
    sprintf(cnr_dir("tdc/%s_annuals_%s.csv.gz"),
            myvar, cnr_res(means))

  means %>% to_df3() %>% write_csv(fn, na = "")

}

tbs <- cnr_outer(vars[which(vars == "nppv")], 2017:2021, c("var", "year", "fn"))

tbs %>% pmap(.f = function(var, year, ...)
  cnr_save_ym(r_all_tdc, y = year, varname = var, outdir = cnr_dir("tdc"))
)

gzfn <- dir(cnr_dir("tdc"), pattern = "nppv_\\d{4}_0p100\\.csv")
pfn <- gsub("\\.csv", "\\.parquet", gzfn)

cmd <- sprintf("duckdb :memory: \"copy '%s' to '%s' (format 'parquet');\"",
               gzfn, pfn)
writeLines(cmd, cnr_dir("tdc/nppv_monthlies.sh"))

#' We can now export the yearly tables (with columns for the months) in .csv.gz-format to parquet
#' (this takes about 5 minutes using the shell script below)
readLines(cnr_dir("tdc/nppv_monthlies.sh"))

#' ## Tables for annual means
#'
#' For each of the variables we already have .csv.gz files for annual means.
#'
#' We can convert all .csv.gz files to parquet files

gzfn <- dir(cnr_dir("tdc"), pattern = "nppv_annuals.*?\\.csv.gz$")
pfn <- gsub("\\.csv\\.gz", "\\.parquet", gzfn)

cmd <- sprintf("duckdb :memory: \"copy '%s' to '%s' (format 'parquet');\"",
               gzfn, pfn)
cat(cmd)
writeLines(cmd, cnr_dir("tdc/nppv_annuals.sh"))


