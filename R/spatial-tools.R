# -----------------------------------------------------------------------------------------------------------------------
# Copyright (c) 2023, microgeo/Chaonan Li (licn@mtc.edu.cn).                                                           #
# The microgeo is distributed under the terms of the Modified BSD License.                                              #
# Full license is avaliable in the file LICENSE, distributed with this package.                                         #
# -----------------------------------------------------------------------------------------------------------------------

#' @title Show those numeric MODIS(Moderate-resolution Imaging Spectroradiometer, `https://modis.gsfc.nasa.gov/`) metrics
#' avaliable in microgeo R package.
#' @author Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
#' @description This function is used to show all numeric MODIS metrics avaliable in microgeo R package.
#' @return A `data.frame` of numeric MODIS metrics.
#' @examples
#' show_modis_num_metrics()
#' @export
show_modis_num_metrics = function(){
    prod.aval <- system.file("modis", "modis.products.Rds", package = "microgeo") %>% base::readRDS()
    ava.metrics <- c("NDVI", "EVI", "GPP", "PsnNet", "NPP", "LST", "ET", "PET", "LE","PLE")
    colname <- c("measure", "resolution", "type", "name", "sds", "unit", "scale.factor")
    prod.aval <- prod.aval[prod.aval$measure.name %in% ava.metrics %>% which,]
    rownames(prod.aval) <- prod.aval %>% nrow %>% seq
    colnames(prod.aval) <- colname
    return(prod.aval)
}

#' @title Show those classification MODIS (Moderate-resolution Imaging Spectroradiometer, `https://modis.gsfc.nasa.gov/`)
#' metrics avaliable in microgeo R package.
#' @author Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
#' @description This function is used to show all classification MODIS metrics avaliable in microgeo R package.
#' @return A `data.frame` of classification MODIS metrics.
#' @examples
#' show_modis_cla_metrics()
#' @export
show_modis_cla_metrics = function(){
    prod.aval <- system.file("modis", "modis.products.Rds", package = "microgeo") %>% base::readRDS()
    ava.metrics <- c("LC_Type1", "LC_Type2", "LC_Type3", "LC_Type4", "LC_Type5", "LC_Prop1",
                     "LC_Prop2", "LC_Prop3", "LC_Prop1_Assessment", "LC_Prop2_Assessment",
                     "LC_Prop3_Assessment")
    colname <- c("measure", "resolution", "type", "name", "sds", "unit", "scale.factor")
    prod.aval <- prod.aval[prod.aval$measure.name %in% ava.metrics %>% which,]
    rownames(prod.aval) <- prod.aval %>% nrow %>% seq
    colnames(prod.aval) <- colname
    return(prod.aval)
}

#' @title Get the variable names of all historically spatial data in a microgeo dataset
#' @author Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
#' @description This function is used to get the variable names of all historically spatial data in a microgeo dataset.
#' @param dataset A microgeo dataset with the class of `MicrogeoDataset`.
#' @return A `MicrogeoDataset` class with the following components:
#' \describe{
#'   \item{\code{object$mat}}{A `data.frame` of rarefied ASV/gene abundance.}
#'   \item{\code{object$ant}}{A `data.frame` of ASV/gene anootations.}
#'   \item{\code{object$met}}{A `data.frame` of sample information.}
#'   \item{\code{object$map}}{A `SpatialPolygonsDataFrame` of map.}
#'   \item{\code{object$phy}}{A phylogenetic tree with `newick` format (`phylo` class) if applicable.}
#'   \item{\code{object$env}}{A `data.frame` of measured environmental properties if applicable.}
#'   \item{\code{object$spa$rast}}{`SpatRaster` of all spatial data.}
#'   \item{\code{object$spa$unit}}{A `data.frame` of unit for historically spatial variables.}
#'   \item{\code{object$*}}{Other spatial and biogeographic traits if applicable.}
#' }
#' @seealso
#' \code{\link[microgeo:read_aliyun_map]{microgeo::read_aliyun_map()}}
#' \code{\link[microgeo:create_dataset]{microgeo::create_dataset()}}
#' \code{\link[microgeo:show_dataset]{microgeo::show_dataset()}}
#' \code{\link[microgeo:get_his_bioc]{microgeo::get_his_bioc()}}
#' @examples
#' # Create a microgeo dataset
#' data(qtp)
#' map <- read_aliyun_map(adcode = c(540000, 630000, 510000))
#' dataset.dts <- create_dataset(mat = qtp$asv, ant = qtp$tax, met = qtp$met, map = map,
#'                               phy = qtp$tre, env = qtp$env, lon = "longitude", lat = "latitude")
#' dataset.dts %>% show_dataset()
#' dataset.dts %<>% get_his_bioc(res = 2.5, out.dir = "test/microgeo_data")
#' dataset.dts %>% show_dataset()
#'
#' # Get the variable names of all historically spatial data in a microgeo dataset
#' dataset.dts %<>% get_spa_vars()
#' head(dataset.dts$spa$unit)
#' @export
get_spa_vars = function(dataset){
    dataset %>% check_dataset(); var.names <- c()
    unit.df <- system.file("unit", "spa.unit.Rds", package = "microgeo") %>% base::readRDS()
    if (!dataset$spa$rast$his %>% is.null) var.names <- c(var.names, dataset$spa$rast$his %>% names)
    if (!dataset$spa$rast$cla %>% is.null) var.names <- c(var.names, dataset$spa$rast$cla %>% names)
    if (var.names %>% length == 0) stop("No spatial data in your `microgeo` dataset!")
    res <- unit.df[which(unit.df$measure %in% var.names),]
    rownames(res) <- res %>% nrow %>% seq
    dataset$spa$unit <- res
    return(dataset)
}

#' @title Extract spatial data for each sample
#' @author Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
#' @description This function is used to extract spatial data for each sample based on longitudes and latitudes.
#' @param dataset A microgeo dataset with the class of `MicrogeoDataset`.
#' @param method Method for extracting values with points (`simple` or `bilinear`). With the `simple` values for the cell
#' a point falls in are returned. With a `bilinear` the returned values are interpolated from the values of the 4 nearest
#' raster cells. Default is `simple`.
#' @param type Which type of data would be extracted? Select one from `both`, `his` and `cla`. Default is `both`.
#' @param remove.na Remove rows with `NA` values. Default is `TRUE`.
#' @param ... Parameters parsed by \code{terra::extract()}.
#' \describe{
#'   \item{\code{object$mat}}{A data.frame of ASV/gene abundance.}
#'   \item{\code{object$ant}}{A data.frame of ASV/gene anootations.}
#'   \item{\code{object$met}}{A data.frame of sample information.}
#'   \item{\code{object$map}}{A `SpatialPolygonsDataFrame` of map.}
#'   \item{\code{object$phy}}{A phylogenetic tree with newick format if applicable.}
#'   \item{\code{object$env}}{A data.frame of measured environmental properties if applicable.}
#'   \item{\code{object$spa$rast}}{`SpatRaster` of all spatial data.}
#'   \item{\code{object$spa$unit}}{A `data.frame` of unit for historically spatial variables if applicable.}
#'   \item{\code{object$spa$tabs}}{A `data.frame` of historically spatial variables for each sample.}
#'   \item{\code{object$*}}{Other spatial and biogeographic traits if applicable.}
#' }
#' @seealso
#' \code{\link[terra:extract]{terra::extract()}}
#' \code{\link[microgeo:read_aliyun_map]{microgeo::read_aliyun_map()}}
#' \code{\link[microgeo:create_dataset]{microgeo::create_dataset()}}
#' \code{\link[microgeo:show_dataset]{microgeo::show_dataset()}}
#' \code{\link[microgeo:get_his_bioc]{microgeo::get_his_bioc()}}
#' \code{\link[microgeo:extract_data_from_spatraster]{microgeo::extract_data_from_spatraster()}}
#' \code{\link[microgeo:tidy_dataset]{microgeo::tidy_dataset()}}
#' @examples
#' # Create a microgeo dataset
#' data(qtp)
#' map <- read_aliyun_map(adcode = c(540000, 630000, 510000))
#' dataset.dts <- create_dataset(mat = qtp$asv, ant = qtp$tax, met = qtp$met, map = map,
#'                               phy = qtp$tre, env = qtp$env, lon = "longitude", lat = "latitude")
#' dataset.dts %>% show_dataset()
#' dataset.dts %<>% get_his_bioc(res = 2.5, out.dir = "test/microgeo_data")
#' dataset.dts %>% show_dataset()
#'
#' # Extract spatial data for each sample
#' dataset.dts %<>% extract_data_from_spatraster(type = 'his')
#' dataset.dts %<>% tidy_dataset()
#' dataset.dts %>% show_dataset()
#' head(dataset.dts$spa$tabs)
#' @export
extract_data_from_spatraster = function(dataset, method = c('simple', 'bilinear'),
                                        type = c('both', 'his', 'cla'), remove.na = TRUE, ...){
    dataset %>% check_dataset()
    method <- ifelse(method %>% length > 1, method[1], method)
    type <- ifelse(type %>% length > 1, type[1], type)
    if (type == 'both') type <- c('cla', 'his')
    extract.rst <- lapply(type, function(vname){
        data.object <- dataset$spa$rast[[vname]]
        if (data.object %>% is.null)
            paste0("No `", vname, "` data in your dataset! Try set <type> as `his` or `cla`!") %>% stop()
        ext.rs <- terra::extract(data.object, dataset$met[,c("longitude", "latitude")], method = method, ...)
        if (ext.rs %>% ncol == 2){ # only one variable
            data.colname <- colnames(ext.rs)[2]
            ext.rs <- data.frame(row.names = ext.rs %>% rownames, val = ext.rs[,2])
            rownames(ext.rs) <- dataset$met %>% rownames; colnames(ext.rs) <- data.colname
        }else{
            ext.rs <- ext.rs[,-1]; rownames(ext.rs) <- dataset$met %>% rownames
        }
        ext.rs
    }) %>% do.call("cbind", .)
    if (remove.na){
        extract.rst <- extract.rst %>% na.omit
        msg1 <- "Some samples were failed to be applied for extraction. use `remove.na = FALSE` to check them!"
        msg2 <- "All samples were failed to be applied for extraction. use `remove.na = FALSE` to check them!"
        if (extract.rst %>% nrow < dataset$met %>% nrow & extract.rst %>% nrow != 0) msg1 %>% warning()
        if (extract.rst %>% nrow == 0) msg2 %>% stop()
    }
    dataset$spa$tabs <- extract.rst
    show_stat_msg('spa$tabs')
    return(dataset)
}

#' @title Subset the classification SpatRaster.
#' @author Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
#' @description This function is designed to subset the SpatRaster returned by \code{microgeo::get_modis_cla_metrics()}.
#' @param spat.raster A classification SpatRaster returned by \code{microgeo::get_modis_cla_metrics()}.
#' @param use.class Class number. See \href{https://lpdaac.usgs.gov/documents/1409/MCD12_User_Guide_V61.pdf}{here}.
#' @return A `SpatRaster`.
#' @seealso
#' \code{\link[microgeo:read_aliyun_map]{microgeo::read_aliyun_map()}}
#' \code{\link[microgeo:create_dataset]{microgeo::create_dataset()}}
#' \code{\link[microgeo:show_dataset]{microgeo::show_dataset()}}
#' \code{\link[microgeo:get_modis_cla_metrics]{microgeo::get_modis_cla_metrics()}}
#' \code{\link[microgeo:plot_bmap]{microgeo::plot_bmap()}}
#' \code{\link[microgeo:add_spatraster]{microgeo::add_spatraster()}}
#' \code{\link[microgeo:add_north_arrow]{microgeo::add_north_arrow()}}
#' \code{\link[microgeo:add_scale_bar]{microgeo::add_scale_bar()}}
#' \code{\link[microgeo:get_modis_cla_metrics]{microgeo::get_modis_cla_metrics()}}
#' @examples
#' # Create a microgeo dataset
#' data(qtp)
#' map <- read_aliyun_map(adcode = c(540000, 630000, 510000))
#' dataset.dts <- create_dataset(mat = qtp$asv, ant = qtp$tax, met = qtp$met, map = map,
#'                               phy = qtp$tre, env = qtp$env, lon = "longitude", lat = "latitude")
#' dataset.dts %>% show_dataset()
#' dataset.dts %<>% get_modis_cla_metrics(username = "username", password = "passwd", out.dir = "test/microgeo_data")
#' dataset.dts %>% show_dataset()
#'
#' # Visualize all values of LC_Type1
#' dataset.dts$map %>%
#'     plot_bmap() %>%
#'     add_spatraster(spat.raster = dataset.dts$spa$rast$cla$LC_Type1,
#'                    color = c(RColorBrewer::brewer.pal(12, "Set3"), RColorBrewer::brewer.pal(9, "Set1"))) %>%
#'     add_north_arrow() %>% add_scale_bar()
#'
#' # Visualize Grasslands(10) and Barren(16)
#' gb.spat.raster <- subset_cla_spatraster(spat.raster = dataset.dts$spa$rast$cla$LC_Type1, use.class = c(10, 16))
#' dataset.dts$map %>%
#'     plot_bmap() %>%
#'     add_spatraster(spat.raster = gb.spat.raster,
#'                    color = RColorBrewer::brewer.pal(12, "Set3")[1:2]) %>%
#'     add_north_arrow() %>% add_scale_bar()
#'
#' # Visualize Grasslands(10)
#' g.spat.raster <- subset_cla_spatraster(spat.raster = dataset.dts$spa$rast$cla$LC_Type1, use.class = 10)
#' dataset.dts$map %>%
#'     plot_bmap() %>%
#'     add_spatraster(spat.raster = g.spat.raster, color = "green") %>%
#'     add_north_arrow() %>% add_scale_bar()
#' @export
subset_cla_spatraster = function(spat.raster, use.class){
    spat.raster %>% check_spatraster()
    if (use.class %>% length == 0) stop("At least one value is required in <use.class>!")
    index <- which(!terra::values(spat.raster) %in% use.class)
    terra::values(spat.raster)[index] <- NA
    terra::as.factor(spat.raster) %>% return()
}

#' @title Mask a SpatRaster through using another SpatRaster as a reference
#' @author Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
#' @description This function is used to mask a SpatRaster using another SpatRaster with classification values as a ref.
#' @param tar.spat Target `SpatRaster`.
#' @param ref.spat Reference `SpatRaster`.
#' @param use.class Class number. See \href{https://lpdaac.usgs.gov/documents/1409/MCD12_User_Guide_V61.pdf}{here}.
#' @return A `SpatRaster`.
#' @seealso
#' \code{\link[microgeo:read_aliyun_map]{microgeo::read_aliyun_map()}}
#' \code{\link[microgeo:create_dataset]{microgeo::create_dataset()}}
#' \code{\link[microgeo:show_dataset]{microgeo::show_dataset()}}
#' \code{\link[microgeo:get_modis_cla_metrics]{microgeo::get_modis_cla_metrics()}}
#' \code{\link[microgeo:rarefy_count_table]{microgeo::rarefy_count_table()}}
#' \code{\link[microgeo:tidy_dataset]{microgeo::tidy_dataset()}}
#' \code{\link[microgeo:calc_alpha_div]{microgeo::calc_alpha_div()}}
#' \code{\link[microgeo:interp_kri]{microgeo::interp_kri()}}
#' \code{\link[microgeo:plot_bmap]{microgeo::plot_bmap()}}
#' \code{\link[microgeo:add_spatraster]{microgeo::add_spatraster()}}
#' \code{\link[microgeo:add_north_arrow]{microgeo::add_north_arrow()}}
#' \code{\link[microgeo:add_scale_bar]{microgeo::add_scale_bar()}}
#' \code{\link[microgeo:add_crs]{microgeo::add_crs()}}
#' @examples
#' # Create a microgeo dataset
#' data(qtp)
#' map <- read_aliyun_map(adcode = c(540000, 630000, 510000))
#' dataset.dts <- create_dataset(mat = qtp$asv, ant = qtp$tax, met = qtp$met, map = map,
#'                               phy = qtp$tre, env = qtp$env, lon = "longitude", lat = "latitude")
#' dataset.dts %>% show_dataset()
#' dataset.dts %<>% get_modis_cla_metrics(username = "username", password = "passwd", out.dir = "test/microgeo_data")
#' dataset.dts %>% show_dataset()
#' dataset.dts %<>% rarefy_count_table()
#' dataset.dts %<>% tidy_dataset()
#' dataset.dts %<>% calc_alpha_div(measures = c("observed", "shannon"))
#' dataset.dts %>% show_dataset()
#'
#' # Perform kriging interpolation for alpha diversity indices
#' kri.rst.shannon <- interp_kri(map = dataset.dts$map,
#'                               met = dataset.dts$met,
#'                               dat = dataset.dts$div$alpha,
#'                               var = 'shannon', model = 'Mat', trim.dup = TRUE)
#'
#' # Mask the results by using Grasslands(10) and Barren(16)
#' kri.rst.shannon.masked <- mask_spatraster_by_cla(tar.spat = kri.rst.shannon,
#'                                                  ref.spat = dataset.dts$spa$rast$cla$LC_Type1,
#'                                                  use.class = c(10, 16))
#' dataset.dts$map %>% plot_bmap() %>%
#'     add_spatraster(spat.raster = kri.rst.shannon.masked) %>%
#'     add_scale_bar() %>% add_north_arrow() %>% add_crs()
#'
#' # Mask the results by using Grasslands(10)
#' kri.rst.shannon.masked <- mask_spatraster_by_cla(tar.spat = kri.rst.shannon,
#'                                                  ref.spat = dataset.dts$spa$rast$cla$LC_Type1,
#'                                                  use.class = 10)
#' dataset.dts$map %>% plot_bmap() %>%
#'     add_spatraster(spat.raster = kri.rst.shannon.masked) %>%
#'     add_scale_bar() %>% add_north_arrow() %>% add_crs()
#' @export
mask_spatraster_by_cla = function(tar.spat, ref.spat, use.class){
    tar.spat %>% check_spatraster(arg.name = "tar.spat")
    ref.spat %>% check_spatraster(arg.name = "ref.spat")
    if (length(use.class) == 0)
        stop("At least one value is required in <use.class>!")
    resampled.tar.spat <- terra::resample(x = tar.spat,
                                          y = ref.spat,
                                          method = 'bilinear') # resample the `tar.spat` by using `ref.spat` as a ref.
    index <- which(!ref.spat %>% terra::values() %in% use.class)
    if (index %>% length == 0) stop("Invalid <use.class>!")
    terra::values(resampled.tar.spat)[index] <- NA
    return(resampled.tar.spat)
}
