# -----------------------------------------------------------------------------------------------------------------------
# Copyright (c) 2023, microgeo/Chaonan Li (licn@mtc.edu.cn).                                                            #
# The microgeo is distributed under the terms of the Modified BSD License.                                              #
# Full license is avaliable in the file LICENSE, distributed with this package.                                         #
# -----------------------------------------------------------------------------------------------------------------------

#' @title Download aridity index from figshare
#' @author Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
#' @description This function is designed to download aridity index (AI) from the global aridity index database (figshare:
#' `10.6084/m9.figshare.7504448.v3`). Resolution is `30'`.
#' @param dataset A microgeo dataset with the class of `MicrogeoDataset`.
#' @param out.dir Directory path to save downloaded files. Default is `microgeo_data`
#' @return A `MicrogeoDataset` class with the following components:
#' \describe{
#'   \item{\code{object$mat}}{A `data.frame` of ASV/gene abundance.}
#'   \item{\code{object$ant}}{A `data.frame` of ASV/gene anootations.}
#'   \item{\code{object$met}}{A `data.frame` of sample information.}
#'   \item{\code{object$map}}{A `SpatialPolygonsDataFrame` of map.}
#'   \item{\code{object$phy}}{A phylogenetic tree with `newick` format if applicable.}
#'   \item{\code{object$env}}{A `data.frame` of measured environmental properties if applicable.}
#'   \item{\code{object$spa$rast$his}}{Historical `SpatRaster` of aridity index and other variables if applicable.}
#'   \item{\code{object$*}}{Spatial and biogeographic traits if applicable.}
#' }
#' @seealso
#' \code{\link[microgeo:read_aliyun_map]{microgeo::read_aliyun_map()}}
#' \code{\link[microgeo:create_dataset]{microgeo::create_dataset()}}
#' \code{\link[microgeo:show_dataset]{microgeo::show_dataset()}}
#' \code{\link[microgeo:plot_bmap]{microgeo::plot_bmap()}}
#' \code{\link[microgeo:add_spatraster]{microgeo::add_spatraster()}}
#' \code{\link[microgeo:add_label]{microgeo::add_label()}}
#' \code{\link[microgeo:add_north_arrow]{microgeo::add_north_arrow()}}
#' \code{\link[microgeo:add_scale_bar]{microgeo::add_scale_bar()}}
#' \code{\link[microgeo:add_crs]{microgeo::add_crs()}}
#' @examples
#' # Create a microgeo dataset
#' data(qtp)
#' showtext::showtext_auto(enable = TRUE)
#' map <- read_aliyun_map(adcode = c(540000, 630000, 510000))
#' dataset.dts <- create_dataset(mat = qtp$asv, ant = qtp$tax, met = qtp$met, map = map,
#'                               phy = qtp$tre, env = qtp$env, lon = "longitude", lat = "latitude")
#' dataset.dts %>% show_dataset()
#'
#' # Download aridity index for the research region
#' dataset.dts %<>% get_ai(out.dir = "test/microgeo_data")
#' dataset.dts %>% show_dataset()
#'
#' # Visualize the aridity index
#' dataset.dts$map %>% plot_bmap() %>%
#'    add_spatraster(spat.raster = dataset.dts$spa$rast$his$AI,
#'                   color = colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlGn"))(100)) %>%
#'    add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = 'X.CENTER', lat.var = 'Y.CENTER') %>%
#'    add_north_arrow() %>% add_scale_bar() %>% add_crs()
#' dataset.dts$map %>% plot_bmap() %>%
#'    add_spatraster(spat.raster = dataset.dts$spa$rast$his$AI,
#'                   color = RColorBrewer::brewer.pal(11, "RdYlGn")[c(1,3,5,9,11)],
#'                   breaks = c(0.03, 0.2, 0.5, 0.65), labels = c("HAR", "AR", "SER", "SHR", "HR")) %>%
#'    add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = 'X.CENTER', lat.var = 'Y.CENTER') %>%
#'    add_north_arrow() %>% add_scale_bar() %>% add_crs()
#' @export
get_ai = function(dataset, out.dir = "microgeo_data"){
    check_dataset(dataset)
    pat <- create_dir(dirpath = file.path(out.dir, "aridity_index"))
    ais <- download_aridity_index(outpath = pat)
    ais.crs <- terra::crs(ais, proj = TRUE, describe = TRUE, parse = TRUE)[1, 6]
    map.crs <- terra::crs(terra::vect(dataset$map), proj = TRUE, describe = TRUE, parse = TRUE)[1, 6]
    if (is.na(ais.crs) | ais.crs != map.crs){
        paste0("reprojecting the CRS of SpatRaster to epsg:", map.crs, ", it takes a while...") %>% show_comm_msg()
        ais <- terra::project(ais, as.character(terra::crs(dataset$map)))
    }
    ais <- terra::crop(ais, terra::ext(dataset$map)) %>% terra::mask(., terra::vect(dataset$map))
    names(ais) <- "AI"; ais <- ais * 1e-04
    dataset %<>% merge_spat_raster(spat.rast = ais)
    return(dataset)
}

#' @title Download elevation from WorldClim database
#' @author Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
#' @description This function is implemented to download the elevation data from the WorldClim database version 2.1 (URL:
#' `https://www.worldclim.org/`).
#' @param dataset A microgeo dataset with the class of `MicrogeoDataset`.
#' @param res Which resolution would be used? Valid values are `10, 5, 2.5, 0.5` (minutes of a degree). Default is `10`.
#' @param out.dir Directory path to save downloaded files. Default is `microgeo_data`
#' @return A `MicrogeoDataset` class with the following components:
#' \describe{
#'   \item{\code{object$mat}}{A `data.frame` of ASV/gene abundance.}
#'   \item{\code{object$ant}}{A `data.frame` of ASV/gene anootations.}
#'   \item{\code{object$met}}{A `data.frame` of sample information.}
#'   \item{\code{object$map}}{A `SpatialPolygonsDataFrame` of map.}
#'   \item{\code{object$phy}}{A phylogenetic tree with `newick` format if applicable.}
#'   \item{\code{object$env}}{A `data.frame` of measured environmental properties if applicable.}
#'   \item{\code{object$spa$rast$his}}{Historical `SpatRaster` of elevation and other variables if applicable.}
#'   \item{\code{object$*}}{Spatial and biogeographic traits if applicable.}
#' }
#' @seealso
#' \code{\link[microgeo:read_aliyun_map]{microgeo::read_aliyun_map()}}
#' \code{\link[microgeo:create_dataset]{microgeo::create_dataset()}}
#' \code{\link[microgeo:show_dataset]{microgeo::show_dataset()}}
#' \code{\link[microgeo:plot_bmap]{microgeo::plot_bmap()}}
#' \code{\link[microgeo:add_spatraster]{microgeo::add_spatraster()}}
#' \code{\link[microgeo:add_label]{microgeo::add_label()}}
#' \code{\link[microgeo:add_north_arrow]{microgeo::add_north_arrow()}}
#' \code{\link[microgeo:add_scale_bar]{microgeo::add_scale_bar()}}
#' \code{\link[microgeo:add_crs]{microgeo::add_crs()}}
#' @examples
#' # Create a microgeo dataset
#' data(qtp)
#' showtext::showtext_auto(enable = TRUE)
#' map <- read_aliyun_map(adcode = c(540000, 630000, 510000))
#' dataset.dts <- create_dataset(mat = qtp$asv, ant = qtp$tax, met = qtp$met, map = map,
#'                               phy = qtp$tre, env = qtp$env, lon = "longitude", lat = "latitude")
#' dataset.dts %>% show_dataset()
#'
#' # Download elevation for the research region
#' dataset.dts %<>% get_elev(res = 2.5, out.dir = "test/microgeo_data")
#' dataset.dts %>% show_dataset()
#'
#' # Visualize the elevation
#' dataset.dts$map %>% plot_bmap() %>%
#'     add_spatraster(spat.raster = dataset.dts$spa$rast$his$ELEV) %>%
#'     add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = 'X.CENTER', lat.var = 'Y.CENTER') %>%
#'     add_north_arrow() %>% add_scale_bar() %>% add_crs()
#' dataset.dts$map %>% plot_bmap() %>%
#'     add_spatraster(spat.raster = dataset.dts$spa$rast$his$ELEV, breaks = c(3000, 4000, 5000, 6000),
#'                    labels = c("<3000", "3000-4000", "4000-5000", "5000-6000", ">6000")) %>%
#'     add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = 'X.CENTER', lat.var = 'Y.CENTER') %>%
#'     add_north_arrow() %>% add_scale_bar() %>% add_crs()
#' @export
get_elev = function(dataset, res = c(10, 5, 2.5, 0.5), out.dir = "microgeo_data"){
    check_dataset(dataset)
    res <- ifelse(length(res) > 1, res[1], res)
    pat <- create_dir(dirpath = file.path(out.dir, "elevation"))
    ele <- download_elev(res = res, outpath = pat)
    ele.crs <- terra::crs(ele, proj = TRUE, describe = TRUE, parse = TRUE)[1, 6]
    map.crs <- terra::crs(terra::vect(dataset$map), proj = TRUE, describe = TRUE, parse = TRUE)[1, 6]
    if (is.na(ele.crs) | ele.crs != map.crs){
        paste0("reprojecting the CRS of SpatRaster to epsg:", map.crs, ", it takes a while...") %>% show_comm_msg()
        ele <- terra::project(ele, as.character(terra::crs(dataset$map)))
    }
    ele <- terra::crop(ele, terra::ext(dataset$map)) %>% terra::mask(., terra::vect(dataset$map))
    names(ele) <- "ELEV"; dataset %<>% merge_spat_raster(spat.rast = ele)
    return(dataset)
}

#' @title Download historical bioclimatic variables from WorldClim database
#' @author Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
#' @description This function is implemented to download historical bioclimatic variables from WorldClim database version
#' 2.1 (`https://www.worldclim.org/`).
#' @param dataset A microgeo dataset with the class of `MicrogeoDataset`.
#' @param res Which resolution would be used? Valid values are `10, 5, 2.5, 0.5` (minutes of a degree). Default is `10`.
#' @param out.dir Directory path to save downloaded files. Default is `microgeo_data`
#' @return A `MicrogeoDataset` class with the following components:
#' \describe{
#'   \item{\code{object$mat}}{A `data.frame` of ASV/gene abundance.}
#'   \item{\code{object$ant}}{A `data.frame` of ASV/gene anootations.}
#'   \item{\code{object$met}}{A `data.frame` of sample information.}
#'   \item{\code{object$map}}{A `SpatialPolygonsDataFrame` of map.}
#'   \item{\code{object$phy}}{A phylogenetic tree with `newick` format if applicable.}
#'   \item{\code{object$env}}{A `data.frame` of measured environmental properties if applicable.}
#'   \item{\code{object$spa$rast$his}}{Historical `SpatRaster` of bioclimatic variables and other vars if applicable.}
#'   \item{\code{object$*}}{Spatial and biogeographic traits if applicable.}
#' }
#' @seealso
#' \code{\link[microgeo:read_aliyun_map]{microgeo::read_aliyun_map()}}
#' \code{\link[microgeo:create_dataset]{microgeo::create_dataset()}}
#' \code{\link[microgeo:show_dataset]{microgeo::show_dataset()}}
#' \code{\link[microgeo:plot_bmap]{microgeo::plot_bmap()}}
#' \code{\link[microgeo:add_spatraster]{microgeo::add_spatraster()}}
#' \code{\link[microgeo:add_label]{microgeo::add_label()}}
#' \code{\link[microgeo:add_north_arrow]{microgeo::add_north_arrow()}}
#' \code{\link[microgeo:add_scale_bar]{microgeo::add_scale_bar()}}
#' \code{\link[microgeo:add_crs]{microgeo::add_crs()}}
#' @examples
#' # Create a microgeo dataset
#' data(qtp)
#' showtext::showtext_auto(enable = TRUE)
#' map <- read_aliyun_map(adcode = c(540000, 630000, 510000))
#' dataset.dts <- create_dataset(mat = qtp$asv, ant = qtp$tax, met = qtp$met, map = map,
#'                               phy = qtp$tre, env = qtp$env, lon = "longitude", lat = "latitude")
#' dataset.dts %>% show_dataset()
#'
#' # Download 19 historical bioclimatic variables for the research region
#' dataset.dts %<>% get_his_bioc(res = 2.5, out.dir = "test/microgeo_data")
#' dataset.dts %>% show_dataset()
#'
#' # Visualize the Bio12 (Mean annual precipitation). Just an example.
#' dataset.dts$map %>% plot_bmap() %>%
#'    add_spatraster(spat.raster = dataset.dts$spa$rast$his$Bio12,
#'                   color = colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlGn"))(100)) %>%
#'    add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = 'X.CENTER', lat.var = 'Y.CENTER') %>%
#'    add_north_arrow() %>% add_scale_bar() %>% add_crs()
#' dataset.dts$map %>% plot_bmap() %>%
#'    add_spatraster(spat.raster = dataset.dts$spa$rast$his$Bio12, breaks = c(100, 300, 500, 700),
#'                   color = RColorBrewer::brewer.pal(11, "RdYlGn")[c(1,3,5,9,11)],
#'                   labels = c("A", "B", "C", "D", "E")) %>%
#'    add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = 'X.CENTER', lat.var = 'Y.CENTER') %>%
#'    add_north_arrow() %>% add_scale_bar() %>% add_crs()
#'
#' # Visualize the Bio1 (Mean annual temperature). Just an example.
#' dataset.dts$map %>% plot_bmap() %>%
#'    add_spatraster(spat.raster = dataset.dts$spa$rast$his$Bio1,
#'                   color = colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlGn"))(100)) %>%
#'    add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = 'X.CENTER', lat.var = 'Y.CENTER') %>%
#'    add_north_arrow() %>% add_scale_bar() %>% add_crs()
#' dataset.dts$map %>% plot_bmap() %>%
#'    add_spatraster(spat.raster = dataset.dts$spa$rast$his$Bio1, breaks = c(-2, -1, 1, 2),
#'                   color = RColorBrewer::brewer.pal(11, "RdYlGn")[c(1,3,5,9,11)],
#'                   labels = c("A", "B", "C", "D", "E")) %>%
#'    add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = 'X.CENTER', lat.var = 'Y.CENTER') %>%
#'    add_north_arrow() %>% add_scale_bar() %>% add_crs()
#' @export
get_his_bioc = function(dataset, res = c(10, 5, 2.5, 0.5), out.dir = "microgeo_data"){
    check_dataset(dataset)
    res <- ifelse(length(res) > 1, res[1], res)
    pat <- create_dir(dirpath = file.path(out.dir, "his_bioclimatic_vars"))
    bio <- download_his_bioc(res = res, outpath = pat)
    bio.crs <- terra::crs(bio, proj = TRUE, describe = TRUE, parse = TRUE)[1, 6]
    map.crs <- terra::crs(terra::vect(dataset$map), proj = TRUE, describe = TRUE, parse = TRUE)[1, 6]
    if (is.na(bio.crs) | bio.crs != map.crs){
        paste0("reprojecting the CRS of SpatRaster to epsg:", map.crs, ", it takes a while...") %>% show_comm_msg()
        bio <- terra::project(bio, as.character(terra::crs(dataset$map)))
    }
    bio <- terra::crop(bio, terra::ext(dataset$map)) %>% terra::mask(., terra::vect(dataset$map))
    names(bio) <- paste0("Bio", seq(19))
    dataset %<>% merge_spat_raster(spat.rast = bio)
    return(dataset)
}

#' @title Download future bioclimatic variables from WorldClim database
#' @author Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
#' @description This function is implemented to download future bioclimatic variables from WorldClim database version 2.1
#' (`https://www.worldclim.org/`).
#' @param dataset A microgeo dataset with the class of `MicrogeoDataset`.
#' @param gcm Which model abbreviation of future bioclimatic data would be used? Default is `BCC-CSM2-MR`. There are many
#' climate model are avaliable, Visit `https://www.worldclim.org/data/cmip6/cmip6climate.html` for more details.
#' @param sce Which Shared Socio-economic Pathway code of future bioclimatic data would be applied? Select one or mutiple
#' from `ssp126, ssp245, ssp370, ssp585`.
#' @param yea Which time period would be applied? Select from `2021-2040, 2041-2060, 2061-2080, 2081-2100`.
#' @param res Which resolution would be used? Valid values are `10, 5, 2.5, 0.5` (minutes of a degree). Default is `10`.
#' @param out.dir Directory path to save downloaded files. Default is `microgeo_data`
#' @return A `MicrogeoDataset` class with the following components:
#' \describe{
#'   \item{\code{object$mat}}{A `data.frame` of ASV/gene abundance.}
#'   \item{\code{object$ant}}{A `data.frame` of ASV/gene anootations.}
#'   \item{\code{object$met}}{A `data.frame` of sample information.}
#'   \item{\code{object$map}}{A `SpatialPolygonsDataFrame` of map.}
#'   \item{\code{object$phy}}{A phylogenetic tree with `newick` format if applicable.}
#'   \item{\code{object$env}}{A `data.frame` of measured environmental properties if applicable.}
#'   \item{\code{object$spa$rast$fut}}{Future `SpatRaster` of bioclimatic variables.}
#'   \item{\code{object$*}}{Spatial and biogeographic traits if applicable.}
#' }
#' @seealso
#' \code{\link[microgeo:read_aliyun_map]{microgeo::read_aliyun_map()}}
#' \code{\link[microgeo:create_dataset]{microgeo::create_dataset()}}
#' \code{\link[microgeo:show_dataset]{microgeo::show_dataset()}}
#' \code{\link[microgeo:get_his_bioc]{microgeo::get_his_bioc()}}
#' \code{\link[microgeo:plot_bmap]{microgeo::plot_bmap()}}
#' \code{\link[microgeo:add_spatraster]{microgeo::add_spatraster()}}
#' \code{\link[microgeo:add_label]{microgeo::add_label()}}
#' \code{\link[microgeo:add_north_arrow]{microgeo::add_north_arrow()}}
#' \code{\link[microgeo:add_scale_bar]{microgeo::add_scale_bar()}}
#' \code{\link[microgeo:add_crs]{microgeo::add_crs()}}
#' @examples
#' # Create a microgeo dataset
#' data(qtp)
#' showtext::showtext_auto(enable = TRUE)
#' map <- read_aliyun_map(adcode = c(540000, 630000, 510000))
#' dataset.dts <- create_dataset(mat = qtp$asv, ant = qtp$tax, met = qtp$met, map = map,
#'                               phy = qtp$tre, env = qtp$env, lon = "longitude", lat = "latitude")
#' dataset.dts %>% show_dataset()
#'
#' # Download 19 historical and future bioclimatic variables for the research region
#' dataset.dts %<>% get_his_bioc(res = 2.5, out.dir = "test/microgeo_data")
#' dataset.dts %<>% get_fut_bioc(res = 2.5, out.dir = "test/microgeo_data")
#' dataset.dts %>% show_dataset()
#'
#' # Visualize the Bio12 (Mean annual precipitation) of `BCC-CSM2-MR|ssp585|2061-2080`. Just an example.
#' dataset.dts$map %>% plot_bmap() %>%
#'    add_spatraster(spat.raster = dataset.dts$spa$rast$fut$`BCC-CSM2-MR|ssp585|2061-2080`$Bio12,
#'                   color = colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlGn"))(100)) %>%
#'    add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = 'X.CENTER', lat.var = 'Y.CENTER') %>%
#'    add_north_arrow() %>% add_scale_bar() %>% add_crs()
#' dataset.dts$map %>% plot_bmap() %>%
#'    add_spatraster(spat.raster = dataset.dts$spa$rast$fut$`BCC-CSM2-MR|ssp585|2061-2080`$Bio12,
#'                   breaks = c(100, 300, 500, 700),
#'                   color = RColorBrewer::brewer.pal(11, "RdYlGn")[c(1,3,5,9,11)],
#'                   labels = c("A", "B", "C", "D", "E")) %>%
#'    add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = 'X.CENTER', lat.var = 'Y.CENTER') %>%
#'    add_north_arrow() %>% add_scale_bar() %>% add_crs()
#'
#' # Visualize the Bio1 (Mean annual temperature) of `BCC-CSM2-MR|ssp585|2061-2080`. Just an example.
#' dataset.dts$map %>% plot_bmap() %>%
#'    add_spatraster(spat.raster = dataset.dts$spa$rast$fut$`BCC-CSM2-MR|ssp585|2061-2080`$Bio1,
#'                   color = colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlGn"))(100)) %>%
#'    add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = 'X.CENTER', lat.var = 'Y.CENTER') %>%
#'    add_north_arrow() %>% add_scale_bar() %>% add_crs()
#' dataset.dts$map %>% plot_bmap() %>%
#'    add_spatraster(spat.raster = dataset.dts$spa$rast$fut$`BCC-CSM2-MR|ssp585|2061-2080`$Bio1,
#'                   breaks = c(-2, -1, 1, 2),
#'                   color = RColorBrewer::brewer.pal(11, "RdYlGn")[c(1,3,5,9, 11)],
#'                   labels = c("A", "B", "C", "D", "E")) %>%
#'    add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = 'X.CENTER', lat.var = 'Y.CENTER') %>%
#'    add_north_arrow() %>% add_scale_bar() %>% add_crs()
#' @export
get_fut_bioc = function(dataset, gcm = c("BCC-CSM2-MR", "ACCESS-CM2", "CNRM-CM6-1", "CNRM-ESM2-1", "CanESM5"),
                        sce = c('ssp126', 'ssp585'), yea = c("2021-2040", "2061-2080"),
                        res = c(10, 5, 2.5), out.dir = "microgeo_data"){
    check_dataset(dataset)
    check.his.bio.rst <- lapply(X = seq(19), function(x) is.null(dataset$spa$rast$his[[paste0("Bio", x)]])) %>% unlist()
    if (TRUE %in% check.his.bio.rst) stop("historical bioclimatic variables are required, please run `get_his_bioc()`!")
    gcm <- ifelse(length(gcm) > 1, gcm[1], gcm)
    res <- ifelse(length(res) > 1, res[1], res)
    pat <- file.path(out.dir, "fut_bioclimatic_vars") %>% create_dir(., recursive = T)
    obj <- lapply(sce, function(s) {rst = c();
    for (y in yea) rst <- c(rst, paste(gcm, s, y, sep = '|')); rst}) %>% unlist()
    fut <- lapply(obj, function(j) {
        dnrst <- download_fut_bioc(arg = j, res = res, outpath = pat)
        dnrst.crs <- terra::crs(dnrst, proj = TRUE, describe = TRUE, parse = TRUE)[1, 6]
        mapss.crs <- terra::crs(terra::vect(dataset$map), proj = TRUE, describe = TRUE, parse = TRUE)[1, 6]
        if (is.na(dnrst.crs) | dnrst.crs != mapss.crs){
            paste0("reprojecting the CRS of SpatRaster to epsg:", mapss.crs, ", it takes a while...") %>%
                show_comm_msg()
            dnrst <- terra::project(dnrst, as.character(terra::crs(dataset$map)))
        }
        dnrst <- terra::crop(dnrst, terra::ext(dataset$map)) %>%
            terra::mask(., terra::vect(dataset$map)); names(dnrst) <- paste0("Bio", seq(19))
        dnrst <- terra::resample(x = dnrst, y = dataset$spa$rast$his, method = 'bilinear')
    }); names(fut) <- obj
    dataset$spa$rast$fut <- fut
    paste0("spa$rast$fut [", "19 variables; ", length(obj), " groups]") %>% show_stat_msg()
    return(dataset)
}

#' @title Download numeric MODIS metrics from EOSDIS
#' @author Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
#' @description This function is implemented to download numeric MODIS metrics from `https://cmr.earthdata.nasa.gov`.
#' @param dataset A microgeo dataset with the class of `MicrogeoDataset`.
#' @param username Username of your (free) EOSDIS account. Sign up \href{https://urs.earthdata.nasa.gov/users/new}{here}.
#' @param password Password of your (free) EOSDIS account. Sign up \href{https://urs.earthdata.nasa.gov/users/new}{here}.
#' @param measures Measures. Use \code{\link[microgeo:show_modis_num_metrics]{microgeo::show_modis_num_metrics()}} to see
#' all avaliable measures. More details are avaliable at \href{https://modis.gsfc.nasa.gov/data/}{MODIS Web}. The default
#' is all avaliable measures.
#' @param prod.res Spatial resolution of MODIS image product. Valid resolution are `1000, 500, 250`(meter). It only works
#' when the measures are `NDVI` and/or `ENV`. Other MODIS image product harbor a fixed resolution. Default is `1000` m.
#' @param prod.typ Type of MODIS image product. Select one from `"Terra", "Aqua"`. Default is `Terra`.
#' @param date.ran Which date range would be applied? Default is `c("2019-08-01|2019-09-01", "2020-08-01|2020-09-01")`.
#' @param nums.job How many threads would be used for the remote-sense image merging? Default is `30`. This would consume
#' a large amounts of computational resources!
#' @param out.dir Directory path to save downloaded files. Default is `microgeo_data`
#' @return A `MicrogeoDataset` class with the following components:
#' \describe{
#'   \item{\code{object$mat}}{A `data.frame` of ASV/gene abundance.}
#'   \item{\code{object$ant}}{A `data.frame` of ASV/gene anootations.}
#'   \item{\code{object$met}}{A `data.frame` of sample information.}
#'   \item{\code{object$map}}{A `SpatialPolygonsDataFrame` of map.}
#'   \item{\code{object$phy}}{A phylogenetic tree with `newick` format if applicable.}
#'   \item{\code{object$env}}{A `data.frame` of measured environmental properties if applicable.}
#'   \item{\code{object$spa$rast$his}}{Historical `SpatRaster` of numeric MODIS metrics and other vars if applicable.}
#'   \item{\code{object$*}}{Spatial and biogeographic traits if applicable.}
#' }
#' @seealso
#' \code{\link[microgeo:show_modis_num_metrics]{microgeo::show_modis_num_metrics()}}
#' \code{\link[microgeo:read_aliyun_map]{microgeo::read_aliyun_map()}}
#' \code{\link[microgeo:create_dataset]{microgeo::create_dataset()}}
#' \code{\link[microgeo:show_dataset]{microgeo::show_dataset()}}
#' \code{\link[microgeo:plot_bmap]{microgeo::plot_bmap()}}
#' \code{\link[microgeo:add_spatraster]{microgeo::add_spatraster()}}
#' \code{\link[microgeo:add_label]{microgeo::add_label()}}
#' \code{\link[microgeo:add_north_arrow]{microgeo::add_north_arrow()}}
#' \code{\link[microgeo:add_scale_bar]{microgeo::add_scale_bar()}}
#' \code{\link[microgeo:add_crs]{microgeo::add_crs()}}
#' @examples
#' # Create a microgeo dataset
#' data(qtp)
#' showtext::showtext_auto(enable = TRUE)
#' map <- read_aliyun_map(adcode = c(540000, 630000, 510000))
#' dataset.dts <- create_dataset(mat = qtp$asv, ant = qtp$tax, met = qtp$met, map = map,
#'                               phy = qtp$tre, env = qtp$env, lon = "longitude", lat = "latitude")
#' dataset.dts %>% show_dataset()
#'
#' # Download numeric MODIS metrics of research region
#' dataset.dts %<>% get_modis_num_metrics(username = "username", password = "password",
#'                                        prod.typ = "Terra", out.dir = "test/microgeo_data")
#' dataset.dts %>% show_dataset()
#'
#' # Visualize NDVI
#' dataset.dts$map %>% plot_bmap() %>%
#'    add_spatraster(spat.raster = dataset.dts$spa$rast$his$NDVI,
#'                   color = colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlGn"))(100)) %>%
#'    add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = 'X.CENTER', lat.var = 'Y.CENTER') %>%
#'    add_north_arrow() %>% add_scale_bar() %>% add_crs()
#'
#' # Visualize EVI
#' dataset.dts$map %>% plot_bmap() %>%
#'    add_spatraster(spat.raster = dataset.dts$spa$rast$his$EVI,
#'                   color = colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlGn"))(100)) %>%
#'    add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = 'X.CENTER', lat.var = 'Y.CENTER') %>%
#'    add_north_arrow() %>% add_scale_bar() %>% add_crs()
#' @export
get_modis_num_metrics = function(dataset, username, password, measures = c("NDVI", "EVI"),
                                 prod.res = c(1000, 500, 250), prod.typ = c("Terra", "Aqua"),
                                 date.ran = c("2019-08-01|2019-09-01", "2020-08-01|2020-09-01"),
                                 nums.job = 30, out.dir = "microgeo_data"){
    check_dataset(dataset)
    if (username %>% is.null) stop("The <username> of EOSDIS is required!")
    if (password %>% is.null) stop("The <password> of EOSDIS is required!")
    prod.res <- ifelse(length(prod.res) > 1, prod.res[1], prod.res)
    prod.typ <- ifelse(length(prod.typ) > 1, prod.typ[1], prod.typ)
    modis.prod <- get_modis_products(dataset = dataset, measures = measures,
                                     prod.res = prod.res, prod.typ = prod.typ, date.ran = date.ran)
    query.list <- que_modis_products(modis.prod = modis.prod)
    if (nrow(query.list) == 0) stop("Failed to search any image data!")
    downs.path <- dow_modis_products(username = username, password = password, search.res = query.list,
                                     outpath = out.dir)
    merge.ptvs <- ptv_modis_products(prod.list = modis.prod$prod.list, hdfs.path = downs.path)
    modis.data <- meg_modis_products(prod.list = modis.prod$prod.list, ptvdata = merge.ptvs,
                                     hdfpath = downs.path, threads = nums.job, outpath = out.dir)
    modis.dbss <- avg_modis_products(modis.data = modis.data, date.ran = date.ran, threads = nums.job,
                                     outpath = out.dir)
    for (mea in modis.dbss$measure){
        tif.path <- modis.dbss[which(modis.dbss$measure == mea),]$dbpath
        tmp.plnt <- terra::rast(tif.path)
        plnt.crs <- terra::crs(tmp.plnt, proj = TRUE, describe = TRUE, parse = TRUE)[1, 6]
        maps.crs <- terra::crs(terra::vect(dataset$map), proj = TRUE, describe = TRUE, parse = TRUE)[1, 6]
        if (is.na(plnt.crs) | plnt.crs != maps.crs){
            paste0("reprojecting the CRS of SpatRaster to epsg:", maps.crs, ", it takes a while...") %>% show_comm_msg()
            tmp.plnt <- terra::project(tmp.plnt, as.character(terra::crs(dataset$map)))
        }
        tmp.plnt <- terra::crop(tmp.plnt, terra::ext(dataset$map)) %>% terra::mask(., terra::vect(dataset$map))
        names(tmp.plnt) <- mea; dataset %<>% merge_spat_raster(spat.rast = tmp.plnt)
    }
    return(dataset)
}

#' @title Download classification MODIS metrics from EOSDIS
#' @author Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
#' @description This function is implemented to download some classification MODIS metrics from the EOSDIS website ( URL:
#' `https://cmr.earthdata.nasa.gov` ). When visualizing the `SpatRaster` obtained from this function, we do not recommend
#' using `add_crs()`. This is because the data has a high resolution, which would significantly slow down the plotting.
#' @param dataset A microgeo dataset with the class of `MicrogeoDataset`.
#' @param username Username of your (free) EOSDIS account. Sign up \href{https://urs.earthdata.nasa.gov/users/new}{here}.
#' @param password Password of your (free) EOSDIS account. Sign up \href{https://urs.earthdata.nasa.gov/users/new}{here}.
#' @param measures Measures. Use \code{\link[microgeo:show_modis_cla_metrics]{microgeo::show_modis_cla_metrics()}} to see
#' all avaliable  measures. More details are avaliable at \href{https://lpdaac.usgs.gov/products/mcd12q1v006/}{MODIS} and
#' \href{https://lpdaac.usgs.gov/documents/1409/MCD12_User_Guide_V61.pdf}{MCD12_User_Guide_V61}. Default is `"LC_Type1"`.
#' @param year Download data for which year? If it is `least`, we will search data by using value obtained by subtracting
#' two years from the current year. Default is `least`.
#' @param nums.job How many threads would be used for the remote-sense image merging? Default is `30`. This would consume
#' a large amounts of computational resources!
#' @param out.dir Directory path to save downloaded files. Default is `microgeo_data`
#' @return A `MicrogeoDataset` class with the following components:
#' \describe{
#'   \item{\code{object$mat}}{A `data.frame` of ASV/gene abundance.}
#'   \item{\code{object$ant}}{A `data.frame` of ASV/gene anootations.}
#'   \item{\code{object$met}}{A `data.frame` of sample information.}
#'   \item{\code{object$map}}{A `SpatialPolygonsDataFrame` of map.}
#'   \item{\code{object$phy}}{A phylogenetic tree with `newick` format if applicable.}
#'   \item{\code{object$env}}{A `data.frame` of measured environmental properties if applicable.}
#'   \item{\code{object$spa$rast$cla}}{Historical `SpatRaster` of classification MODIS metrics.}
#'   \item{\code{object$*}}{Spatial and biogeographic traits if applicable.}
#' }
#' @seealso
#' \code{\link[microgeo:show_modis_cla_metrics]{microgeo::show_modis_cla_metrics()}}
#' \code{\link[microgeo:read_aliyun_map]{microgeo::read_aliyun_map()}}
#' \code{\link[microgeo:create_dataset]{microgeo::create_dataset()}}
#' \code{\link[microgeo:show_dataset]{microgeo::show_dataset()}}
#' \code{\link[microgeo:plot_bmap]{microgeo::plot_bmap()}}
#' \code{\link[microgeo:add_spatraster]{microgeo::add_spatraster()}}
#' \code{\link[microgeo:add_label]{microgeo::add_label()}}
#' \code{\link[microgeo:add_north_arrow]{microgeo::add_north_arrow()}}
#' \code{\link[microgeo:add_scale_bar]{microgeo::add_scale_bar()}}
#' @examples
#' # Create a microgeo dataset
#' data(qtp)
#' showtext::showtext_auto(enable = TRUE)
#' map <- read_aliyun_map(adcode = c(540000, 630000, 510000))
#' dataset.dts <- create_dataset(mat = qtp$asv, ant = qtp$tax, met = qtp$met, map = map,
#'                               phy = qtp$tre, env = qtp$env, lon = "longitude", lat = "latitude")
#' dataset.dts %>% show_dataset()
#'
#' # Download classification MODIS metrics of research region
#' dataset.dts %<>% get_modis_cla_metrics(username = "username",
#'                                        password = "password", out.dir = "test/microgeo_data")
#' dataset.dts %>% show_dataset()
#'
#' # Visualize LC_Type1
#' dataset.dts$map %>% plot_bmap() %>%
#'    add_spatraster(spat.raster = dataset.dts$spa$rast$cla$LC_Type1,
#'                   color = c(RColorBrewer::brewer.pal(12, "Set3"), RColorBrewer::brewer.pal(9, "Set1"))) %>%
#'    add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = 'X.CENTER', lat.var = 'Y.CENTER') %>%
#'    add_north_arrow() %>% add_scale_bar()
#' @export
get_modis_cla_metrics = function(dataset, username, password, measures = "LC_Type1", year = "least",
                                 nums.job = 30, out.dir = "microgeo_data"){
    check_dataset(dataset)
    if (username %>% is.null) stop("The <username> of EOSDIS is required!")
    if (password %>% is.null) stop("The <password> of EOSDIS is required!")
    if (!is.numeric(year) & year != "least") stop("The <year> must be a number if it is not 'least'!")
    if (year == 'least'){
        year <- date() %>% strsplit(" ") %>% unlist()
        year <- year[length(year)] %>% as.numeric(); year <- year - 2
    }
    date.ran <- paste0(year, "-01-01|", year, "-12-31")
    modis.prod <- get_modis_products(dataset = dataset, measures = measures, prod.res = prod.res,
                                     prod.typ = "Combined", date.ran = date.ran)
    query.list <- que_modis_products(modis.prod = modis.prod)
    if (nrow(query.list) == 0) stop("Failed to search image data!")
    downs.path <- dow_modis_products(username = username, password = password,
                                     search.res = query.list, outpath = out.dir)
    merge.ptvs <- ptv_modis_products(prod.list = modis.prod$prod.list, hdfs.path = downs.path)
    modis.data <- meg_modis_products(prod.list = modis.prod$prod.list, ptvdata = merge.ptvs,
                                     hdfpath = downs.path, threads = nums.job, outpath = out.dir)
    modis.dbss <- col_modis_products(modis.data = modis.data, date.ran = date.ran, outpath = out.dir)
    for (mea in modis.dbss$measure){
        tif.path <- modis.dbss[which(modis.dbss$measure == mea),]$dbpath
        tmp.plnt <- terra::rast(tif.path)
        plnt.crs <- terra::crs(tmp.plnt, proj = TRUE, describe = TRUE, parse = TRUE)[1, 6]
        maps.crs <- terra::crs(terra::vect(dataset$map), proj = TRUE, describe = TRUE, parse = TRUE)[1, 6]
        if (is.na(plnt.crs) | plnt.crs != maps.crs){
            paste0("reprojecting the CRS of SpatRaster to epsg:", maps.crs, ", it takes a while...") %>% show_comm_msg()
            tmp.plnt <- terra::project(tmp.plnt, as.character(terra::crs(dataset$map)))
        }
        tmp.plnt <- terra::crop(tmp.plnt, terra::ext(dataset$map)) %>% terra::mask(., terra::vect(dataset$map))
        names(tmp.plnt) <- mea; dataset %<>% merge_spat_raster(spat.rast = terra::as.factor(tmp.plnt), type = 'cla')
    }
    return(dataset)
}

#' @title Download soil metrics from SoilGRIDS
#' @author Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
#' @description This function is used to download soil metrics from SoilGRIDS (`https://www.isric.org/explore/soilgrids`)
#' by using \code{geodata::soil_world()}.
#' @param dataset A microgeo dataset with the class of `MicrogeoDataset`.
#' @param measures Variables to be downloaded. Default is all avaliable variables. See \code{geodata::soil_world()}.
#' @param depth Soil depth. One of `c(5, 15, 30, 60, 100, 200)`, which means the depth ranges of 0-5, 5-15, 15-30, 30-60,
#' 60-100, 100-200 cm. Default is `5`.
#' @param out.dir Directory path to save downloaded files. Default is `microgeo_data`
#' @return A `MicrogeoDataset` class with the following components:
#' \describe{
#'   \item{\code{object$mat}}{A `data.frame` of ASV/gene abundance.}
#'   \item{\code{object$ant}}{A `data.frame` of ASV/gene anootations.}
#'   \item{\code{object$met}}{A `data.frame` of sample information.}
#'   \item{\code{object$map}}{A `SpatialPolygonsDataFrame` of map.}
#'   \item{\code{object$phy}}{A phylogenetic tree with `newick` format if applicable.}
#'   \item{\code{object$env}}{A `data.frame` of measured environmental properties if applicable.}
#'   \item{\code{object$spa$rast$his}}{Historical `SpatRaster` of soil metrics (SoilGRIDS) and other vars if applicable.}
#'   \item{\code{object$*}}{Spatial and biogeographic traits if applicable.}
#' }
#' @seealso
#' \code{\link[geodata:soil_world]{geodata::soil_world()}}
#' \code{\link[microgeo:read_aliyun_map]{microgeo::read_aliyun_map()}}
#' \code{\link[microgeo:create_dataset]{microgeo::create_dataset()}}
#' \code{\link[microgeo:show_dataset]{microgeo::show_dataset()}}
#' \code{\link[microgeo:plot_bmap]{microgeo::plot_bmap()}}
#' \code{\link[microgeo:add_spatraster]{microgeo::add_spatraster()}}
#' \code{\link[microgeo:add_label]{microgeo::add_label()}}
#' \code{\link[microgeo:add_north_arrow]{microgeo::add_north_arrow()}}
#' \code{\link[microgeo:add_scale_bar]{microgeo::add_scale_bar()}}
#' \code{\link[microgeo:add_crs]{microgeo::add_crs()}}
#' @examples
#' # Create a microgeo dataset
#' data(qtp)
#' showtext::showtext_auto(enable = TRUE)
#' map <- read_aliyun_map(adcode = c(540000, 630000, 510000))
#' dataset.dts <- create_dataset(mat = qtp$asv, ant = qtp$tax, met = qtp$met, map = map,
#'                               phy = qtp$tre, env = qtp$env, lon = "longitude", lat = "latitude")
#' dataset.dts %>% show_dataset()
#'
#' # Download soil metrics from the SoilGRIDS for the research region
#' dataset.dts %<>% get_soilgrid(depth = 5, out.dir = "test/microgeo_data")
#' dataset.dts %>% show_dataset()
#'
#' # Visualize phh2o
#' dataset.dts$map %>% plot_bmap() %>%
#'    add_spatraster(spat.raster = dataset.dts$spa$rast$his$phh2o,
#'                   color = colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlGn"))(100)) %>%
#'    add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = 'X.CENTER', lat.var = 'Y.CENTER') %>%
#'    add_north_arrow() %>% add_scale_bar() %>% add_crs()
#'
#' # Visualize soc
#' dataset.dts$map %>% plot_bmap() %>%
#'    add_spatraster(spat.raster = dataset.dts$spa$rast$his$soc,
#'                   color = colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlGn"))(100)) %>%
#'    add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = 'X.CENTER', lat.var = 'Y.CENTER') %>%
#'    add_north_arrow() %>% add_scale_bar() %>% add_crs()
#'
#' # Visualize nitrogen
#' dataset.dts$map %>% plot_bmap() %>%
#'    add_spatraster(spat.raster = dataset.dts$spa$rast$his$nitrogen,
#'                   color = colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlGn"))(100)) %>%
#'    add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = 'X.CENTER', lat.var = 'Y.CENTER') %>%
#'    add_north_arrow() %>% add_scale_bar() %>% add_crs()
#' @export
get_soilgrid = function(dataset, measures = c("bdod", "cfvo", "nitrogen", "phh2o", "sand", "silt", "soc", "ocd"),
                        depth = c(5, 15, 30, 60, 100, 200), out.dir = "microgeo_data"){ # no "clay" data
    check_dataset(dataset)
    pat.soil <- create_dir(dirpath = file.path(out.dir, "soilgrid_products"))
    depth <- ifelse(depth %>% length > 1, depth[1], depth)
    for (mea in measures) {
        soilrst <- geodata::soil_world(var = mea, depth = depth, stat = "mean", path = pat.soil)
        spa.crs <- terra::crs(soilrst, proj = TRUE, describe = TRUE, parse = TRUE)[1, 6]
        map.crs <- terra::crs(terra::vect(dataset$map), proj = TRUE, describe = TRUE, parse = TRUE)[1, 6]
        if (spa.crs %>% is.na | spa.crs != map.crs){
            paste0("reprojecting the CRS of SpatRaster to epsg:", map.crs, ", it takes a while...") %>% show_comm_msg()
            soilrst <- terra::project(soilrst, as.character(terra::crs(dataset$map)))
        }
        soilrst <- terra::crop(soilrst, terra::ext(dataset$map)) %>% terra::mask(., terra::vect(dataset$map))
        names(soilrst) <- mea; dataset %<>% merge_spat_raster(spat.rast = soilrst)
    }
    return(dataset)
}

#' @title Process soil metrics of CHINA
#' @author Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
#' @description This function could process soil metrics downloaded from `http://globalchange.bnu.edu.cn/research/soil2`.
#' Because of the limitations in copyrights, the spatial data of China soil properties should be manually downloaded from
#' \href{http://globalchange.bnu.edu.cn/research/soil2}{here}.
#' @param dataset A microgeo dataset with the class of `MicrogeoDataset`.
#' @param measures Variables to be processed. Default is all avaliable variables.
#' @param depth Soil depth. It should be one of `0.045, 0.091, 0.166, 0.289, 0.493, 0.829, 1.383, 2.296`, which means the
#' depth ranges of 0-0.045, 0.045-0.091, 0.091-0.166, 0.166-0.289, 0.289-0.493, 0.493-0.829, 0.829-1.383 & 1.383-2.296 m.
#' Default is `0.045`.
#' @param out.dir Directory path to save downloaded files. Default is `microgeo_data`
#' @return A `MicrogeoDataset` class with the following components:
#' \describe{
#'   \item{\code{object$mat}}{A `data.frame` of ASV/gene abundance.}
#'   \item{\code{object$ant}}{A `data.frame` of ASV/gene anootations.}
#'   \item{\code{object$met}}{A `data.frame` of sample information.}
#'   \item{\code{object$map}}{A `SpatialPolygonsDataFrame` of map.}
#'   \item{\code{object$phy}}{A phylogenetic tree with `newick` format if applicable.}
#'   \item{\code{object$env}}{A `data.frame` of measured environmental properties if applicable.}
#'   \item{\code{object$spa$rast$his}}{Historical `SpatRaster` of soil metrics (CHINA) and other vars if applicable.}
#'   \item{\code{object$*}}{Spatial and biogeographic traits if applicable.}
#' }
#' @seealso
#' \code{\link[microgeo:read_aliyun_map]{microgeo::read_aliyun_map()}}
#' \code{\link[microgeo:create_dataset]{microgeo::create_dataset()}}
#' \code{\link[microgeo:show_dataset]{microgeo::show_dataset()}}
#' \code{\link[microgeo:plot_bmap]{microgeo::plot_bmap()}}
#' \code{\link[microgeo:add_spatraster]{microgeo::add_spatraster()}}
#' \code{\link[microgeo:add_label]{microgeo::add_label()}}
#' \code{\link[microgeo:add_north_arrow]{microgeo::add_north_arrow()}}
#' \code{\link[microgeo:add_scale_bar]{microgeo::add_scale_bar()}}
#' \code{\link[microgeo:add_crs]{microgeo::add_crs()}}
#' @examples
#' # Create a microgeo dataset
#' data(qtp)
#' showtext::showtext_auto(enable = TRUE)
#' map <- read_aliyun_map(adcode = c(540000, 630000, 510000))
#' dataset.dts <- create_dataset(mat = qtp$asv, ant = qtp$tax, met = qtp$met, map = map,
#'                               phy = qtp$tre, env = qtp$env, lon = "longitude", lat = "latitude")
#' dataset.dts %>% show_dataset()
#'
#' # Process soil metrics of CHINA for a research region
#' dataset.dts %<>% get_soilcn(depth = 0.045, out.dir = "test/microgeo_data")
#' dataset.dts %>% show_dataset()
#'
#' # Visualize PH
#' dataset.dts$map %>% plot_bmap() %>%
#'    add_spatraster(spat.raster = dataset.dts$spa$rast$his$PH,
#'                   color = colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlGn"))(100)) %>%
#'    add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = 'X.CENTER', lat.var = 'Y.CENTER') %>%
#'    add_north_arrow() %>% add_scale_bar() %>% add_crs()
#'
#' # Visualize SOM
#' dataset.dts$map %>% plot_bmap() %>%
#'    add_spatraster(spat.raster = dataset.dts$spa$rast$his$SOM,
#'                   color = colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlGn"))(100)) %>%
#'    add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = 'X.CENTER', lat.var = 'Y.CENTER') %>%
#'    add_north_arrow() %>% add_scale_bar() %>% add_crs()
#'
#' # Visualize TN
#' dataset.dts$map %>% plot_bmap() %>%
#'    add_spatraster(spat.raster = dataset.dts$spa$rast$his$TN,
#'                   color = colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlGn"))(100)) %>%
#'    add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = 'X.CENTER', lat.var = 'Y.CENTER') %>%
#'    add_north_arrow() %>% add_scale_bar() %>% add_crs()
#'
#' # Visualize AP
#' dataset.dts$map %>% plot_bmap() %>%
#'    add_spatraster(spat.raster = dataset.dts$spa$rast$his$AP,
#'                   color = colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlGn"))(100)) %>%
#'    add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = 'X.CENTER', lat.var = 'Y.CENTER') %>%
#'    add_north_arrow() %>% add_scale_bar() %>% add_crs()
#' @export
get_soilcn = function(dataset, measures = c("PH", "SOM", "TN", "TP", "TK", "AN", "AP", "AK", "CEC", "H",
                                            "AL", "CA", "MG", "K", "NA", "SA", "SI", "CL", "GRAV", "BD", "POR"),
                      depth = c(0.045, 0.091, 0.166, 0.289, 0.493, 0.829, 1.383, 2.296), out.dir = "microgeo_data"){
    check_dataset(dataset)
    pat.soil <- create_dir(dirpath = file.path(out.dir, "soilchina_products"))
    depth <- ifelse(length(depth) > 1, depth[1] * 100, depth * 100)
    if (depth == 9.1)   depth = 9.1000004
    if (depth == 49.3)  depth = 49.299999
    if (depth == 82.9)  depth = 82.900002
    if (depth == 229.6) depth = 229.60001
    if (!pat.soil %>% file.exists) {
        msg1 <- paste0("Please download all database files from",
                       " <http://globalchange.bnu.edu.cn/research/soil2> and unzip them into: ")
        paste0(msg1, pat.soil, ". The file suffix may be `.nc`.") %>% stop()
    }
    lostdb.files <- lapply(measures, function(measure){
        dbfile <- file.path(pat.soil, paste0(measure, '.nc'))
        res <- ifelse(!dbfile %>% file.exists, dbfile, NA)
    }) %>% unlist() %>% unique() %>% na.omit() %>% as.vector()
    if (lostdb.files %>% length > 0)
        paste(lostdb.files, collapse = ", ") %>% paste0("Can not find database files: ", .) %>% stop()
    for (i in seq(length(measures))){
        measure <- measures[i]
        soilrst <- terra::rast(file.path(pat.soil, paste0(measure, ".nc"))) %>% suppressWarnings()
        title   <- paste0(measure, "_depth=", depth)
        soilrst <- soilrst[[title]]
        spa.crs <- terra::crs(soilrst, proj = TRUE, describe = TRUE, parse = TRUE)[1, 6]
        map.crs <- terra::crs(terra::vect(dataset$map), proj = TRUE, describe = TRUE, parse = TRUE)[1, 6]
        if (spa.crs %>% is.na | spa.crs != map.crs){
            paste0("reprojecting the CRS of SpatRaster to epsg:", map.crs, ", it takes a while...") %>% show_comm_msg()
            soilrst <- terra::project(soilrst, as.character(terra::crs(dataset$map)))
        }
        soilrst <- terra::crop(soilrst, terra::ext(dataset$map)) %>% terra::mask(., terra::vect(dataset$map)) %>%
            suppressWarnings()
        names(soilrst) <- measure; dataset %<>% merge_spat_raster(spat.rast = soilrst)
    }
    return(dataset)
}
