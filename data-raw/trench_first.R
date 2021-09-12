#' ---
#' title: "Processing of CNR ncdf files using 'terra'"
#' date: "`r Sys.Date()`"
#' output:
#'   rmdformats::robobook:
#'     self_contained: true
#'     highlight: kate
#' ---

#' This script has functions and steps for processing CNR files
#' using ['terra'](https://rspatial.org/terra).
#'
#' The script can be sourced or run or "knitted" into a HTML report
#'
#'#' Files to be processed, restricted to
#' surface-SST, sea ice, SalS; bottom-SBT
#'
library(cnr)
library(dplyr)

vars <- readLines(textConnection(
  "nppv
thetao
siconc
so
bottomT"
))

#' File listing to process, for a specific set of variables
filez <-
  cnr_dir_meta("/home/markus/data/nc/Phy_001-024") %>%
  filter(var %in% vars)

#' Display the list of ncdf files to process
filez %>%
  select(fn, yr, mo, var, id, sz) %>%
  knitr::kable()

#' This is the metadata records by file and variable
#' restricted to one of the variables
meta <-
  filez %>%
  select(fp, var, fn) %>%
  filter(var == vars[1]) %>%
  purrr::pmap_dfr(cnr_ncmeta)  %>%
  arrange(name)

#' We inspect to see if this variable has harmonizing metadata
#' across files
meta %>%
  select(-file) %>%
  unique() %>%
  knitr::kable()


# larger nc files with several full years of monthly data (60 months)
filez %>%
  filter(is.na(yr) & is.na(mo)) %>%
  select(fn, var, sz) %>%
  knitr::kable()

# not as large nc files from early 2021 (not full year of monthly data)
filez %>%
  filter(between(yr, 2021, 2021) & is.na(mo)) %>%
  select(fn, var, sz) %>%
  knitr::kable()

# smaller single nc files for recent 2021 months
filez %>%
  filter(!is.na(mo)) %>%
  select(fn, var, sz, mo) %>%
  arrange(desc(mo)) %>%
  knitr::kable()

# list of terra (multi-layer) rasters by variable

myfiles <- filez %>% filter(var != "nppv")

r_all <-
  myvars %>%
  purrr::map(function(x) cnr_unified(myfiles, x)) %>%
  setNames(., nm = myvars)

r_all %>% purrr::map(nlyr)

#' ## Step 2 - Resampling to tenth degree cells
#'

# terra raster with layers
# f_nppv %>% cnr_resample(transform = FALSE)

# output resamples data as HTML and CSV
# f_nppv %>% cnr_save_htmls(path = "~/data/out/tdc")  # 76 sec
# f_nppv %>% cnr_save_csvs(path = "~/data/out/tdc")  # 109 sec

#' ## Step 3 - Locations with missing data
#'
#' [This SO answer](https://stackoverflow.com/questions/67592103/how-to-unify-the-number-of-non-na-cells-in-raster-stack)
#' by the author of `terra` makes use of the mask function to filter out missing values

# about half of the cells (land) have no data
# na_ratio(cnr_resample(f_nppv, transform = FALSE))

#r <- cnr_raster(f_nppv)
#m <- any(is.na(r))
#x <- mask(r, m)

# gather all nppv data in one rast (with several layers)
myfiles <- filez %>% filter(grepl("^monthly_nppv.*?\\.nc$", fn)) %>% head(1)
r <- cnr_unified(myfiles, "nppv")

years <- unique(year(time(r)))

cnr_annual_mean <- function(r, y) {
  idx <- which(year(time(r)) == y)
  mean(r[[idx]], na.rm = TRUE)
}

# compute mean across years
nppv_means <-
  years %>%
  purrr::map(function(x) cnr_annual_mean(r, x)) %>%
  do.call("c", .)

nppv_means <- setNames(nppv_means, years)

# plot all annual means
plot(nppv_means)

#+ outputz, eval=FALSE

nppv_means %>%
  cnr_save_tiff(outfile = "~/data/out/nppv_annuals_orig.tiff")

# output to CSV
nppv_means %>%
  to_df3() %>%
  write_csv("~/data/out/nppv_annuals_orig.csv.gz", na = "")

# resample and output annual means
nppv10 <- cnr_resample(nppv_means, type = "rast")

nppv10 %>%
  cnr_save_tiff("~/data/out/tdc/nppv_annuals_orig.tif")


# compute global mean for full time period and all cells
global(nppv_means, na.rm = TRUE)

#+ notes


#' # Notes and questions

#' - variable naming for trench 1, what "vernacular" is preferred?

# "surface-SST, sea ice, SalS; bottom-SBT"
# "sea temperature (surface and bottom), surface salinity, primary production and sea ice concentration"

# nppv = approx primary production (total primary production of phyto, net_primary_production_of_biomass_expressed_as_carbon_per_unit_volume_in_sea_water)
# bottomT = bottom-SBT/sea temperature at bottom (sea_water_potential_temperature_at_sea_floor, sea floor potential temperature, bottom-SBT)
# so = SalS/surface salinity (sea water salinity)
# thetao = surface-SST/sea temperature surface (sea_water_potential_temperature)
# siconc = sea ice / sea ice concentration (ice concentration, sea_ice_area_fraction)

# |name          |variable |value                                                                              |
# |:-------------|:--------|:----------------------------------------------------------------------------------|
# |_FillValue    |bottomT  |-32767                                                                             |
# |add_offset    |bottomT  |21                                                                                 |
# |cell_methods  |bottomT  |area: mean                                                                         |
# |long_name     |bottomT  |Sea floor potential temperature                                                    |
# |scale_factor  |bottomT  |0.0007324442                                                                       |
# |standard_name |bottomT  |sea_water_potential_temperature_at_sea_floor                                       |
# |unit_long     |bottomT  |Degrees Celsius                                                                    |
# |units         |bottomT  |degrees_C                                                                          |
# |_ChunkSizes   |nppv     |1, 15, 137, 288                                                                    |
# |_ChunkSizes   |nppv     |1, 13, 171, 360                                                                    |
# |_FillValue    |nppv     |9.96921e+36                                                                        |
# |add_offset    |nppv     |0                                                                                  |
# |long_name     |nppv     |Total Primary Production of Phyto                                                  |
# |scale_factor  |nppv     |1                                                                                  |
# |standard_name |nppv     |net_primary_production_of_biomass_expressed_as_carbon_per_unit_volume_in_sea_water |
# |unit_long     |nppv     |milligrams of Carbon per cubic meter per day                                       |
# |units         |nppv     |mg m-3 day-1                                                                       |
# |_FillValue    |siconc   |-32767                                                                             |
# |add_offset    |siconc   |-3.814814e-05                                                                      |
# |cell_methods  |siconc   |area: mean where sea_ice                                                           |
# |long_name     |siconc   |Ice concentration                                                                  |
# |scale_factor  |siconc   |3.814814e-05                                                                       |
# |standard_name |siconc   |sea_ice_area_fraction                                                              |
# |unit_long     |siconc   |Fraction                                                                           |
# |units         |siconc   |1                                                                                  |
# |_FillValue    |so       |-32767                                                                             |
# |add_offset    |so       |-0.001525925                                                                       |
# |cell_methods  |so       |area: mean                                                                         |
# |long_name     |so       |Salinity                                                                           |
# |scale_factor  |so       |0.001525925                                                                        |
# |standard_name |so       |sea_water_salinity                                                                 |
# |unit_long     |so       |Practical Salinity Unit                                                            |
# |units         |so       |1e-3                                                                               |
# |_FillValue    |thetao   |-32767                                                                             |
# |add_offset    |thetao   |21                                                                                 |
# |cell_methods  |thetao   |area: mean                                                                         |
# |long_name     |thetao   |Temperature                                                                        |
# |scale_factor  |thetao   |0.0007324442                                                                       |
# |standard_name |thetao   |sea_water_potential_temperature                                                    |
# |unit_long     |thetao   |Degrees Celsius                                                                    |
# |units         |thetao   |degrees_C                                                                          |

#' Output from processing the variables
#'

#+ eval=FALSE
r_all <- vars %>% purrr::map(cnr_unified) %>% setNames(., nm = vars)
#system("sudo mkdir -p /media/markus/Reddish/data/out/orig/months /media/markus/Reddish/data/out/tdc/months /media/markus/Reddish/data/out/tdc/html")
#system("sudo chown -R markus:markus /media/markus/Reddish/data/out")
myvars <- vars[which(vars != "nppv")]

#for (myvar in myvars) {
#  cnr_process(r_all, myvar)
#}

r_all_tdc <-
  myvars %>%
  purrr::map(function(x) cnr_process2(r_all, x)) %>%
  do.call("c", .)

names(r_all_tdc)

# write individual layers as monthly tiffs
fn <- sprintf(cnr_dir("tdc/months/%s_%s.tiff"), names(r_all_tdc), cnr_res(r_all_tdc))

writeRaster(r_all_tdc, fn, overwrite=TRUE, NAflag = -9999,
            gdal = c("COMPRESS=LZW", "TFW=YES", "of=COG","PROFILE=GeoTIFF"),
            datatype = "FLT8S")

# we could write yearly tiffs, too... but would expect those to have annuals...
idx <- intersect(which(year(time(r_all_tdc)) == 2017), grep("^so", names(r_all_tdc)))
fn <- sprintf(cnr_dir("tdc/%s_%s.tiff"), "so_2017", cnr_res(r_all_tdc[[idx]]))
names(r_all_tdc[[idx]])
writeRaster(r_all_tdc[[idx]], fn, overwrite=TRUE, NAflag = -9999,
            gdal = c("COMPRESS=LZW", "TFW=YES", "of=COG","PROFILE=GeoTIFF"),
            datatype = "FLT8S")

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

# Usage:
# cnr_crs("default")
# cnr_crs("DO", fn = "001-029")
# purrr::map2_chr(filez$var, filez$fp, cnr_crs)
# filez %>% mutate(crs = map2_chr(.$var, .$fp, cnr_crs))


tbs <- cnr_outer(vars[which(vars != "nppv")], 2017:2021, c("var", "year", "fn"))

knitr::kable(tbs)

#' A couple of functions to output tables in .csv.gz format


#' Now, for each table representing such a variable and year combination,
#' we call the function to output .csv.gz files into a specific folder
#+ eval=FALSE
rz <-
  dir(cnr_dir("tdc"), pattern = "[bottomT|siconc|thetao|so]_0p100\\.tiff$", full.names = TRUE) %>%
  grep(pattern = "annuals", invert = TRUE, value = TRUE, x = .)

r_all_tdc <- rz %>% map(function(x) rast(x))

# NB: the order of myvars is important!
# previously the r_all_tdc$so slot was misnamed to r_all_tdc$thetao and vv

# myvars <- vars[vars != "nppv"]
# r_all_tdc <-
#   myvars %>%
#   map(function(x) cnr_transform(cnr_resample(r_all[[x]], type = "rast"), var = x))
#
# names(r_all_tdc) <- myvars
# myvars %>% walk(function(x) time(r_all_tdc[[x]]) <- time(r_all[[x]]))

tbs %>% pmap(.f = function(var, year, ...)
  cnr_save_ym(r_all_tdc, y = year, varname = var, outdir = cnr_dir("tdc"))
)

gzfn <- dir(cnr_dir("tdc"), pattern = "\\d{4}_0p100\\.csv")
pfn <- gsub("\\.csv", "\\.parquet", gzfn)

cmd <- sprintf("duckdb :memory: \"copy '%s' to '%s' (format 'parquet');\"",
               gzfn, pfn)
writeLines(cmd, cnr_dir("tdc/monthlies.sh"))

#' We can now export the yearly tables (with columns for the months) in .csv.gz-format to parquet
#' (this takes about 5 minutes using the shell script below)
readLines(cnr_dir("tdc/monthlies.sh"))

#' ## Tables for annual means
#'
#' For each of the variables we already have .csv.gz files for annual means.
#'
#' We can convert all .csv.gz files to parquet files

gzfn <- dir(cnr_dir("tdc"), pattern = "annuals.*?\\.csv.gz$")
pfn <- gsub("\\.csv\\.gz", "\\.parquet", gzfn)

cmd <- sprintf("duckdb :memory: \"copy '%s' to '%s' (format 'parquet');\"",
               gzfn, pfn)
cat(cmd)
writeLines(cmd, cnr_dir("tdc/annuals.sh"))

