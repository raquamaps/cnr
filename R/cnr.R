if (getRversion() >= "2.15.1")
  utils::globalVariables(c(
    "x", "y", "CtrLong", "CtrLat",
    "name", "variable", "value", "filez",
    "fp", "var", "fn",
    "yr", "mo",
    "is_match", "fn_re", ".", "id"
  ))

# regular expressions used to extract metadata from filenames
re_id <- "_(\\d{13})_"
re_var <- "_([a-zA-Z]+)(_|\\.)"
re_yr <- "_(\\d{4})(_|\\.)"
re_mo <- "_(\\d{2})\\."

#' Utility function to extract matches from regexes
#' @param x text
#' @param re regex
rex <- function(x, re) {
  res <- gsub(re, "\\1", regmatches(x, gregexpr(re, x)))
  res[which(res == "character(0)")] <- NA_character_
  return (res)
}

#' File path for scratch dir with nc files
#' @param subpath subpath
#' @param scratchdir default location of nc files
#' @export
cnr_dir <- function(subpath, scratchdir = "/media/markus/Stash/data/out") {

  if (!dir.exists(scratchdir)) {
    message("Did not find ", scratchdir, ", creating it...")
    dir.create(scratchdir, recursive = TRUE)
  }

  if (missing(subpath))
    return (scratchdir)

  fp <- file.path(scratchdir, subpath)

  if (!dir.exists(dirname(fp))) {
    message("Did not find ", fp, ", creating it...")
    dir.create(dirname(fp))
  }

  fp
}

#' Data frame with nc files listing with some metadata added,
#' extracted from the filenames
#' @param path path
#' @param pat pattern
#' @export
#' @importFrom dplyr tibble
#' @importFrom purrr map_chr
cnr_dir_meta <- function(path, pat = "*.nc$") {

  fn <- dir(path, pattern = pat)
  if (length(fn) < 0)
    stop("No files matching pattern found in that dir")

  rp <- file.path(path, fn)
  yr <- rex(fn, re_yr)
  mo <- rex(fn, re_mo)
  var <- rex(fn, re_var)
  id <- rex(fn, re_id)
  fp <- normalizePath(rp)

  sz <- fp %>% purrr::map_chr(
    function(x) format(structure(file.size(x), class = "object_size"),
      units = "auto", standard = "SI"))

  dplyr::tibble(
    fn, rp, yr, mo, var, id, fp, sz
    ) %>%
  arrange(var, -desc(yr), -desc(mo))
}

#' Metadata from nc file
#'
#' Inside a netcdf file there is metadata embedded, which can
#' be parsed with the `ncmeta` R package.
#' @param fp filepath
#' @param var variable
#' @param fn filename
#' @export
#' @import ncmeta
cnr_ncmeta <- function(fp, var, fn) {
  ncmeta::nc_atts(fp) %>%
    filter(variable == var) %>%
    mutate(file = fn) %>%
    select(-id)
}

#' Docs for terra::rast() says to prepend format
#' so netCDF files "do not open with the HDF5 driver".
#' @param file file
#' @export
nc_file <- function(file) {
  sprintf("NETCDF:\"%s\"", normalizePath(file))
}

#' Retrieve raster from a nc file (which may have several layers) where
#' translation of values is optional because
#'
#' 1. it takes time and could be delayed until after resampling
#' 1. results from terra seems to lose some original metadata in the translation
#' @param f file
#' @param translate boolean
#' @param rename_layers boolean
#' @export
#' @import terra tictoc
cnr_raster <- function(f, translate = FALSE, rename_layers = TRUE) {

  r <- rast(nc_file(f))

  if (isTRUE(rename_layers))
    names(r) <- cnr_filenames(r)

  if (!isTRUE(translate))
    return(r)

  tr <-
    ncmeta::nc_atts(f) %>%
    filter(name %in% c("add_offset", "scale_factor")) %>%
    pull(value) %>% unlist()

  tictoc::tic()
  r * tr["scale_factor"] + tr["add_offset"]
  tictoc::toc()
}

#' In addition to the metadata embedded in nc file names and to what `ncmeta`
#' can provide, there is also metadata available through inspection of
#' `terra` objects.
#'
#' Metadata from `terra` given a file
#' @param f filepath to nc file
#' @export
#' @import terra
#' @importFrom lubridate month year
cnr_terrameta <- function(f) {

  r <- rast(nc_file(f))

  list(
    time = time(r),  # time slices
    n_ts = length(r),
    mo = lubridate::month(time(r)),
    yr = lubridate::year(time(r)),
    varname = varnames(r),
    longname = longnames(r),
    nlyr = nlyr(r)  # layers (often the same as time slices)
  )
}

#' Resample raster to tenth degree resolution cells
#' with optional transformation of values
#' (using offset and scale from metadata)
#' @param f nc file path or rast
#' @param nth_degrees nth_degrees
#' @param transform boolean
#' @param type one of "file" or "rast"
#' @export
#' @import terra tictoc ncmeta
cnr_resample <- function(f, nth_degrees = 10, transform = TRUE,
                         type = c("file", "rast")) {

  # zones to resample data to
  n <- nth_degrees
  z <- rast(nrows = 180 * n, ncols = 360 * n)

  tictoc::tic()

  #orig <- rast(nc_file(f))
  if (missing(type)) type <- "file"

  orig <- switch(type,
                 file = rast(nc_file(f)),
                 rast = f)

  if (crs(orig) == "")
    crs(orig) <- "+proj=longlat +datum=WGS84 +no_defs000"

  r <- resample(orig, z, "bilinear")

  tictoc::toc()

  if (!isTRUE(transform) | type == "rast") return (r)

  tr <-
    ncmeta::nc_atts(f) %>%
    filter(name %in% c("add_offset", "scale_factor")) %>%
    pull(value) %>% unlist()

  if (is.null(tr)) return (r)

  r * tr["scale_factor"] + tr["add_offset"]

  # terra does not distinguish between NA and NaN, mostly uses NaN

  # TODO: determine deepest layer and return only that one?
  #lyr <-
  #  tibble(lyr = r@ptr$names, min = r@ptr$range_min, max = r@ptr$range_max, idx = 1:nlyr(r)) %>%
  #  filter(!is.nan(min) & !is.nan(max)) %>% tail(1) %>% pull(idx)

  #plot(mean(r[[lyr]]) * tr["scale_factor"] + tr["add_offset"])
  #plot(mean(max(range(r, na.rm = TRUE))))
  #mean(r, na.rm = TRUE) * tr["scale_factor"] + tr["add_offset"]

}

#' Filename(s) on the form variable_year_month
#' (useful for naming several output layers in a netcdf file)
#' @param md metadata from cnr_ncmeta
#' @export
cnr_filename <- function(md) {
  paste(sep = "_", md$varname, md$yr, sprintf("%02i", md$mo))
}

#' Filenames from rast layers
#' @param r rast
#' @param extension file extension
#' @export
#' @import terra
#' @importFrom lubridate month year
cnr_filenames <- function(r, extension) {

  md <- list(
    time = time(r),  # time slices
    mo = lubridate::month(time(r)),
    yr = lubridate::year(time(r)),
    varname = varnames(r),
    longname = longnames(r),
    nlyr = nlyr(r)  # layers (often the same as time slices)
  )

  fns <- paste(sep = "_", md$varname, md$yr, sprintf("%02i", md$mo))

  if (!missing(extension)) fns <- paste0(fns, extension)

  stopifnot(nlyr(r) == length(fns))

  fns

}

#' Output a terra layer to TIFF (using terra)
#' @param r rast
#' @param outfile outfile
#' @param ... args to send to terra::writeRaster
#' @export
#' @import terra
cnr_save_tiff <- function(r, outfile, ...) {

  #  writeRaster(r, outfile, overwrite=TRUE,
  #    #, "ESRI_XML_PAM=TRUE", "GDAL_GEOREF_SOURCES=PAM,WORLDFILE"),
  #    gdal=c("COMPRESS=LZW", "TFW=YES","of=COG"), datatype='INT1U')

  writeRaster(r, outfile, overwrite=TRUE, NAflag = -9999,
              gdal = c("COMPRESS=LZW", "TFW=YES", "of=COG","PROFILE=GeoTIFF"),
              datatype = "FLT8S", ...
  )

}

#' Output all layers in a nc file to TIFF (using terra)
#' @param f file
#' @param path output directory
#' @export
#' @importFrom utils flush.console
cnr_save_tiffs <- function(f, path = cnr_dir()) {
  stopifnot(dir.exists(path))
  rz <- cnr_raster(f)
  fz <- file.path(path, sprintf("%s.tif", cnr_filename(cnr_terrameta(f))))
  message("Saving ", nlyr(rz), " layers in ", f, " to ", path, "\n")
  tictoc::tic()
  for (x in 1:nlyr(rz)) {
    cat(".")
    flush.console()
    cnr_save_tiff(rz[[x]], fz[[x]])
  }
  cat("\n")
  tictoc::toc()

}

#' Output all layers in a nc file to TIFF (using terra)
#' @param r rast
#' @param path output directory
#' @export
cnr_save_tiffz <- function(r, path = cnr_dir()) {
  stopifnot(dir.exists(path))
  rz <- r
  fz <- file.path(path, sprintf("%s_%s.tiff", cnr_filenames(r), cnr_res(r)))
  message("Saving ", nlyr(rz), " layers to ", path, "\n")
  tictoc::tic()
  for (x in 1:nlyr(rz)) {
    cat(".")
    flush.console()
    cnr_save_tiff(rz[[x]], fz[[x]])
  }
  cat("\n")
  tictoc::toc()

}

# Usage:
# cnr_save_tiffs(filez[1, ]$fp)


#' Save outputs from raster to html file
#' @param r rast
#' @param outfile outfile
#' @export
#' @importFrom leaflet colorNumeric leaflet addLegend
cnr_save_html <- function(r, outfile) {

  #rt <- rast(ncol = ncol(r), nrow = nrow(r), crs="epsg:3857")
  #r2 <- terra::project(r, rt)
  #rl <- raster::raster(r2)#, "+proj=longlat +ellps=WGS84 +datum=WGS84")

  rl <- raster::raster(r)  # do not use terra raster
  raster::crs(rl) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  # NB: translation should be already done for input raster!
  val <- as.numeric(c(min(values(rl), na.rm=TRUE), max(values(rl), na.rm = TRUE)))
  pal <- colorNumeric(c("white", "darkred"), val, na.color = "transparent")
  cols <- c("#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#B10026")

  message("Processing conversion to leaflet projection for output to HTML")
  tictoc::tic()
  ras <- leaflet::projectRasterForLeaflet(rl, method = "bilinear")
  tictoc::toc()

  #pal <- leaflet::colorNumeric(cols,
  #  raster::values(ras),
  #  na.color = "transparent")
  #pal <- colorFactor(colors, levels = breaks, ordered = TRUE, na.color = "transparent")

  #pal <- leaflet::colorBin(cols, na.omit(unique(values(ras))),
  #                         bins = length(cols), pretty = TRUE, na.color = "#00000000")

  e <- raster::extent(rl)
  title <- basename(outfile)

  map <-
    leaflet() %>%
    leaflet::addProviderTiles(provider = "Esri.OceanBasemap") %>%
    leaflet::addRasterImage(ras, project = FALSE,
                            colors = pal, opacity = 0.9) %>%
    addLegend(values = raster::values(ras),
              title = title, pal = pal) %>%
    leaflet::fitBounds(lng1 = e@xmin, lat1 = e@ymin, lng2 = e@xmax, lat2 = e@ymax)

  htmlwidgets::saveWidget(map, file = outfile, selfcontained = TRUE)

}

#path = "~/data/out"
#r <- r1z[[1]]
#cnr_save_html(r1z[[1]], file.path(path, sprintf("%s.html", cnr_filename(cnr_terrameta(f1))[[1]])))

#' Save raster layers in nc file to a directory as html files
#' @param f file
#' @param path output dir
#' @param resample boolean
#' @export
cnr_save_htmls <- function(f, path = cnr_dir(), resample = TRUE) {

  stopifnot(dir.exists(path))

  if (isTRUE(resample)) {
    tictoc::tic()
    message("Resampling to tdc...")
    rz <- cnr_resample(f, transform = TRUE)
    tictoc::toc()
  } else {
    rz <- cnr_raster(f)
  }

  fz <- file.path(path, sprintf("%s.html", cnr_filename(cnr_terrameta(f))))
  n_layers <- nlyr(rz)
  message("Saving html: ", n_layers, " layers in ", f, " at ", path, "\n")

  tictoc::tic()
  for (x in 1:n_layers) {
    cat(".")
    flush.console()
    cnr_save_html(rz[[x]], fz[[x]])
  }
  cat("\n")
  tictoc::toc()
}

#myncf <- filez %>% filter(var == "nppv") %>% pull(fp) %>% tail(1)
#myncf <- filez %>% filter(var == "bottomT") %>% pull(fp) %>% tail(1)
#cnr_save_htmls(myncf)
#781.083 sec for 60 months of bottomT
#34 sec for 5 months of nppv 2021

#' Save mapview file for a raster to a specified file
#' @param r rast
#' @param outfile outfile
#' @export
cnr_save_mapview <- function(r, outfile) {

  ln <- names(r)
  rr <- raster::raster(r)
  crs(rr) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

  mv <- mapview::mapview(rr, layer.name = ln,
                         maxpixels =  ncol(r) * nrow(r))

  mapview::mapshot(mv, url = outfile, selfcontained = TRUE,
                   remove_controls = c("homeButton", "layersControl"))
}

# Usage:
# rl <- r_all$nppv[[1]]
# cnr_save_mapview(rl, outfile = sprintf("~/data/out/orig/%s.html", names(rl)))

#' Save raster layers as mapviews in a specified dir
#' @param r rast
#' @param outdir output directory
#' @export
cnr_save_mapviews <- function(r, outdir = "~/data/out/orig") {

  n <- nlyr(r)
  message("Processing ", n, " layers...")

  tic()

  fn <- paste0(file.path(outdir, names(r)), ".html")

  for (i in 1:n) {
    cnr_save_mapview(r[[i]], fn[i])
  }

  toc()

}

# Usage:
# cnr_save_mapviews(r_all$nppv)

#' Convert single raster to a data frame
#' @param r rast
#' @export
to_df3 <- function(r) {
  as.data.frame(r, xy = TRUE, na.rm = FALSE)
}


#' Output a terra layer to a compressed CSV file
#' (using .gz compression which is twice as fast compared to .zip)
#' @param r rast
#' @param outfile outfile
#' @param compress boolean
#' @export
#' @importFrom readr write_csv
cnr_save_csv <- function(r, outfile, compress = TRUE) {
  #readr::write_csv(to_df(r), outfile)
  #data <- to_df3(r, colnames = cnr_filenames(r))
  data <- to_df3(r)
  fn <- normalizePath(outfile, mustWork = FALSE)
  if (isTRUE(compress)) {
    fn <- paste0(fn, ".gz")
  }
  readr::write_csv(data, file = fn, na = "")
}

#' Write .csv files for variables across time
#' @param f file
#' @param path directory for output files
#' @export
cnr_save_csvs <- function(f, path = cnr_dir()) {
  stopifnot(dir.exists(path))
  rz <- cnr_resample(f, transform = TRUE)
  fz <- file.path(path, sprintf("%s.csv", cnr_filename(cnr_terrameta(f))))
  message("Saving csv: ", nlyr(rz), " layers in ", f, " at ", path, "\n")
  tictoc::tic()
  for (x in 1:nlyr(rz)) {
    cat(".")
    flush.console()
    cnr_save_csv(rz[[x]], fz[[x]], compress = FALSE)
  }
  cat("\n")
  tictoc::toc()
}

# Usage:
# cnr_save_csvs(myncf)
# 304.496 sec elapsed

#' Write .csv files for variables across time
#' @param r rast
#' @param path output directory
#' @export
cnr_save_csvz <- function(r, path = cnr_dir()) {
  stopifnot(dir.exists(path))
  rz <- r
  fz <- file.path(path, sprintf("%s_%s.csv", cnr_filenames(r), cnr_res(r)))
  message("Saving csvs: ", nlyr(rz), " layers to ", path, "\n")
  tictoc::tic()
  for (x in 1:nlyr(rz)) {
    cat(".")
    flush.console()
    cnr_save_csv(rz[[x]], fz[[x]], compress = FALSE)
  }
  cat("\n")
  tictoc::toc()
}

#' Unify several nc files into a raster with data
#' for a specific variable (multi-layer/"brick")
#' @param filez vector of file paths
#' @param varname varname
#' @export
cnr_unified <- function(filez, varname) {

  fnz <-
    filez %>%
    filter(var == varname) %>%
    pull(fp)

  # use terra metadata for yr and mo
  fp <-
    fnz %>%
    purrr::map_df(cnr_terrameta, .id = "id") %>%
    arrange(-desc(yr), -desc(mo)) %>%
    pull(id) %>%
    as.integer() %>%
    unique()

  fnz[fp] %>%
    purrr::map(cnr_raster) %>%
    do.call("c", .)

}

#' Ratio of NAs to total cell count in rast
#' @param r rast
#' @export
na_ratio <- function(r)
  global((is.na(r)), sum) / global((!is.na(r)), sum)

#' Mask template from marine tenth cell degrees dataset provided
#' @export
cnr_mask_tdc <- function() {

  tdc <-
    read_csv("~/data/Tdc_marine_cellsXY_rounded.txt", col_types = "ddd") %>%
    select(x = CtrLong, y = CtrLat)

  mtdc <- rast(nrows = 1800, ncols = 3600)
  values(mtdc) <- 0

  idx <- cellFromXY(mtdc, as.matrix(tdc))
  values(mtdc)[idx] <- NA
  values(mtdc)[-idx] <- 1
  mtdc
}

#' Find locations which do not match with the provided template
#' (note that there may possibly be mismatches due to approximations/rounding)
#' @param r rast
#' @param type one of rast or tibble
#' @export
cnr_mismatch_tdc <- function(r, type = c("rast", "tibble")) {

  mtdc <- cnr_mask_tdc()
  mr <- mask(r, mtdc)
  val <- values(mr)[!is.na(values(mr))]

  if (missing(type)) type <- "rast"

  if (type == "tibble") {
    data <-
      xyFromCell(mr, cells(mr)) %>%
      as_tibble() %>%
      mutate(var = val)
    return (data)
  }

  if (type == "rast") {
    message("Mismatches (n): ", length(val))
    return (mr)
  }

}

#' The nc metadata for scale factor and offset for a given variable
#' @param myvar variable name
#' @export
cnr_transform_params <- function(myvar) {

  meta <-
    filez %>%
    select(fp, var, fn) %>%
    purrr::pmap_dfr(cnr_ncmeta)  %>%
    arrange(name)

  res <-
    meta %>%
    select(-file) %>%
    unique() %>%
    filter(variable == myvar, name %in% c("scale_factor", "add_offset"))

  res %>% pull(value)

}

# Usage:
# vars %>% map(cnr_transform_params)

#' Transform given scale factor and offset
#' @param r rast
#' @param var variable
#' @export
cnr_transform <- function(r, var) {
  params <- cnr_transform_params(var)
  r * params$scale_factor + params$add_offset
}

# Usage:
# cnr_transform(cnr_resample(r_all$bottomT, type = "rast"), "bottomT")


#' Resolution string given a raster
#' @param r rast
#' @export
cnr_res <- function(r)
  paste0(collapse = "_", gsub("\\.", "p", unique(sprintf("%.03f", res(r)))))

#' Processing for a variable
#' @param r_all rast
#' @param myvar variable name
#' @export
#' @importFrom readr read_csv
cnr_process <- function(r_all, myvar) {

  tic()
  message("Processing ", myvar, "...")

  message("Writing tiff (original resolution, all layers)...")
  fn <- sprintf(cnr_dir("/orig/%s_%s.tiff"), myvar, cnr_res(r_all[[myvar]]))
  #ns <- sprintf(cnr_dir("/orig/%s_%s.tiff"), names(r_all[[myvar]]), cnr_res(r_all[[myvar]]))
  r_all[[myvar]] %>% cnr_save_tiff(fn)
  r_all[[myvar]] %>% cnr_save_tiffz(path = cnr_dir("/orig/months"))
  r_all[[myvar]] %>% cnr_save_csvz(path = cnr_dir("/orig/months"))

  message("Resampling and transforming variable(s)...")
  r <- cnr_transform(cnr_resample(r_all[[myvar]], type = "rast"), var = myvar)
  # metadata lost after resample -> reuse earlier object!
  #lyrnames <- cnr_filenames(r_all[[myvar]])
  #names(r) <- lyrnames
  time(r) <- time(r_all[[myvar]])

  message("Writing tiff (tenth degree cell resolution, all layers)...")
  fn <- sprintf(cnr_dir("/tdc/%s_%s.tiff"), myvar, cnr_res(r))
  r %>% cnr_save_tiff(outfile = fn)
  # TODO save_tiffz does not provide the proper varname! need to "reset" those
  varnames(r) <- myvar
  r %>% cnr_save_tiffz(path = cnr_dir("/tdc/months"))

  message("Computing annual means...")
  years <- unique(year(time(r)))
  means <-
    years %>%
    purrr::map(function(x) cnr_annual_mean(r, x)) %>%
    do.call("c", .)

  names(means) <- paste(sep = "_", myvar, years)

  message("Writing annual means as tiff....")
  fn <-
    sprintf(cnr_dir("/tdc/%s_annuals_%s.tiff"),
            myvar, cnr_res(means))

  means %>% cnr_save_tiff(outfile = fn)

  # TODO: fix to_df2 to not fiddle with the names
  # OOPS since resets the names for (data passed by reference, not by value!)

  message("Writing annual means as CSV...")

  fn <-
    sprintf(cnr_dir("/tdc/%s_annuals_%s.csv.gz"),
            myvar, cnr_res(means))

  means %>% to_df3() %>% write_csv(fn, na = "")

  #message("Writing annual means as HTML(mapview)...")
  #means %>% cnr_save_mapviews(outdir = cnr_dir("/tdc/html"))

  toc()
}


# get tdc variants (does this change by ref r_all?)

#' Process a list of raster with slots for variables
#' @param r_all rast
#' @param myvar variable to subset on
#' @export
cnr_process2 <- function(r_all, myvar) {
  r <- cnr_transform(cnr_resample(r_all[[myvar]], type = "rast"), var = myvar)
  time(r) <- time(r_all[[myvar]])
  r
}

#' Default CRS for specific nc files
#' @param var the variable
#' @param fn filename pattern
#' @importFrom purrr map_lgl
#' @importFrom stats setNames var
#' @importFrom dplyr `%>%` mutate filter arrange desc between first select as_tibble mutate pull
#' @export
cnr_crs <- function(var, fn) {

  crs <- readr::read_csv(trim_ws = TRUE, show_col_types = FALSE,
    "variable,crs,fn_re
    default,4326,
    thetao,4326,001-024
    siconc,4326,001-024
    so,4326,001-024
    bottomT,4326,001-024
    nppv,4326,001-028
    DO,4326,001-028
    nppv,4258,001-029
    DO,4258,001-029
    salinity column,4326,
    PP,32662,
    ") %>%
    dplyr::mutate(
      crs = paste0("EPSG:", crs)
    )

  if (missing(var))
    var <- unique(crs$variable)

  if (!missing(fn))
    crs <-
      crs %>%
      mutate(is_match = map_lgl(fn_re, function(x) grepl(x, fn))) %>%
      filter(is_match == TRUE)

  crs %>%
    filter(variable %in% var) %>%
    pull(crs) %>%
    setNames(., var)
}

#' A table of variable / year combinations
#' @param vars variable names
#' @param years vector of years
#' @param nm header names
#' @export
cnr_outer <- function(vars, years, nm = c("var", "year", "fn")) {
  xoy <- as.vector(outer(vars, years, function(x, y) paste0(x, "_", y)))
  myx <- as.vector(outer(vars, years, function(x, y) x))
  myy <- as.vector(outer(vars, years, function(x, y) y))
  setNames(data.frame(myx, myy, xoy), nm)
}


#' Index a set of raster layers by year
#' @param r rast
#' @param y year
#' @export
cnr_monthlies <- function(r, y) {
  idx <- which(year(time(r)) == y)
  res <- r[[idx]]
  #cols <- paste0(unique(year(time(res))), "_", tolower(month.abb[month(time(res))]))
  cols <- tolower(month.abb[month(time(res))])
  data <- to_df3(res)
  tb <- data
  names(tb) <- c("x", "y", cols)
  list(r = res, t = tb, cols = cols)
}

#' Write a set of raster layers by year with monthly columns as CSV
#' @param r rast
#' @param y year
#' @param varname variable name
#' @param outdir output directory
#' @param compress boolean
#' @export
cnr_save_ym <- function(r, y, varname, outdir, compress = FALSE) {
  #idx <- intersect(which(year(time(r)) == y), grep(paste0("^", var), names(r)))
  idx <- which(year(time(r)) == y & grepl(paste0("^", varname), names(r)))
  res <- r[[idx]]
  cols <- tolower(month.abb[month(time(res))])
  data <- to_df3(res)
  names(data) <- c("x", "y", cols)
  fn <- paste0(varname, "_", y, "_", cnr_res(r), ".csv")
  afn <- file.path(normalizePath(outdir, mustWork = FALSE), fn)
  if (isTRUE(compress)) {
    afn <- paste0(afn, ".gz")
  }
  readr::write_csv(data, file = afn, na = "")
}

#' Save yearly tables with monthly columns in an output directory as .csv.gz
#' @param r_unified rast
#' @param var variable
#' @param y year
#' @param outdir outdir
#' @param compress boolean
#' @export
cnr_save_monthlies <- function(r_unified, var, y, outdir, compress = TRUE) {
  ms <- cnr_monthlies(r_unified[[var]], year)
  data <- ms$t
  fn <- paste0(var, "_", y, "_", cnr_res(r_unified[[var]]), ".csv")
  afn <- file.path(normalizePath(outdir, mustWork = FALSE), fn)
  if (isTRUE(compress)) {
    afn <- paste0(afn, ".gz")
  }
  readr::write_csv(data, file = afn, na = "")
}

#' View first few rows in a parquet file
#' @param fn filename
#' @param n_row nr of rows to preview
#' @export
#' @importFrom readr read_csv
cnr_view_parquet <- function(fn, n_row = 10) {
  cnr_outdir <- cnr_dir("tdc")
  cmd <- sprintf(
    "duckdb :memory: -csv \"select * from '%s' limit %s;\"",
    file.path(cnr_outdir, fn), n_row)
  read_csv(paste0(system(cmd, intern = TRUE), "\n"), show_col_types = FALSE)
}

#' Query a parquet file for data within (inclusive) a certain extent
#' @param pf parquet file
#' @param ymin ymin
#' @param ymax ymax
#' @param xmin xmin
#' @param xmax xmax
#' @param pd parquet file directory
#' @export
cnr_extract_range <- function(
  pf, ymin, ymax, xmin, xmax, pd = cnr_dir("/tdc")) {

  cmd <- sprintf(
    paste0(
      "duckdb -csv :memory: \"select * from ",
      "'%s/%s' where round(y, 2) >= %s and round(y, 2) <= %s and ",
      "round(x, 2) >= %s and round(x, 2) <= %s \""
    ), pd, pf, xmin, xmax, ymin, ymax)

  message("Running sql query:\n ", cmd)

  readr::read_csv(paste0(system(cmd, intern = TRUE), "\n"))

}

#' Wrapper function to extract data at a specific coordinate
#' @param pf parquet file
#' @param y y
#' @param x x
#' @export
cnr_extract_xy <- function(pf, y, x) {
  cnr_extract_range(pf, ymin = y, ymax = y, xmin = x, xmax = x)
}

#' Wrapper function to extract data within specific bounds
#' @param data data
#' @param xmin xmin
#' @param xmax xmax
#' @param ymin ymin
#' @param ymax ymax
#' @export
cnr_ext_range_df <- function(data, xmin, xmax, ymin, ymax) {
  data %>% filter(between(x, xmin, xmax), between(y, ymin, ymax))
}

#' Extract coordinate data from data frame
#' @param data data
#' @param x x coord
#' @param y y coord
cnr_ext_df <- function(data, x, y) {
  cnr_ext_range_df(data, x, x, y, y)
}

#' Compute means for layers in rast in a given year
#' @param r rast
#' @param y year
#' @return rast with mean
#' @export
cnr_annual_mean <- function(r, y) {
  idx <- which(year(time(r)) == y)
  mean(r[[idx]], na.rm = TRUE)
}

#' Flatten layers in a raster across cells into a single layer
#' where each cell takes the first value which is not NaN (ie is numerical)
#' @param r rast
#' @return flattened raster
cnr_flatten_old <- function(r) {

  s <- r
  # reverse layers
  s <- s[[nlyr(s):1]]

  x <- app(s, fun = function(x) first(which(!is.nan(x))))

  idx <- as.vector(values(x))
  idx[which(is.nan(as.vector(values(x))))] <- nlyr(s) + 1

  # low indices indicate deep layers (due to reversing layers)

  # use a temporary fourth layer
  myl <- rast(nrows = nrow(s), ncols = ncol(s), nlyrs = 1)
  values(myl) <- NaN
  ext(myl) <- ext(s)
  time(myl) <- unique(time(s))

  values(myl) <- values(c(s, myl))[cbind(seq_along(idx), idx)]
  myl

}

cnr_depth_index <- function(r) {
  s <- r
  s <- s[[nlyr(s):1]]
  app(s, fun = function(x) first(which(!is.nan(x))))
}

#' Flatten layers in a raster across cells into a single layer
#' where each cell takes the first value which is not NaN (ie is numerical)
#' @param r rast
#' @return flattened raster
cnr_flatten <- function(r) {
  x <- cnr_depth_index(r)
  selectRange(r[[nlyr(r):1]], x)
}

#' Reproject a layer
#' where each cell takes the first value which is not NaN (ie is numerical)
#' @param r rast
#' @param crs_from valid EPSG crs string
#' @return flattened raster
cnr_reproject <- function(r, crs_from) {

  crs(r) <- crs_from
  z <- rast(ncol = ncol(r), nrow = nrow(r), crs = cnr_crs("default"))
  ext(z) <- ext(r)
  res(z) <- res(r)

  tic()
  rr <- project(r, z)
  toc()

  rr
}

#' Given ncfiles, flatten layers and save as individual nc files
#' @param ncfiles vector of ncdf files
#' @param reproject_from crs to project from
#' @param prefix variable name for prefixing file output
#' @param outdir location of resulting files
#' @return flattened raster
cnr_flatten_save <- function(ncfiles, reproject_from, prefix, outdir = cnr_dir("orig")) {

  for (f in ncfiles) {
    s <- rast(f)
    ts <- unique(time(s))
    fn <- file.path(outdir, sprintf("%s_%s_%s.nc",
      prefix, lubridate::year(ts), sprintf("%02i", lubridate::month(ts))))

    message("Flattening and saving ", fn)
    tictoc::tic()
    myl <- cnr_flatten(s)
    if (!missing(reproject_from)) {
      r <- cnr_reproject(myl, crs_from = reproject_from)
    } else {
      r <- myl
      crs(r) <- cnr_crs("default")
    }
    time(r) <- unique(time(s))
    varnames(r) <- varnames(s)
    longnames(r) <- longnames(s)
    units(r) <- unique(units(s))
    tictoc::toc()

    writeCDF(r, filename = fn,
             overwrite = TRUE,
             varname = unique(varnames(myl)),
             longname = unique(longnames(myl)),
             unit = unique(units(myl))
    )

  }

}
