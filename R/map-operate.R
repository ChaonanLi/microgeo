# -----------------------------------------------------------------------------------------------------------------------
# Copyright (c) 2023, microgeo/Chaonan Li (licn@mtc.edu.cn).                                                            #
# The microgeo is distributed under the terms of the GPL-3 License.                                                     #
# Full license is avaliable in the file LICENSE, distributed with this package.                                         #
# -----------------------------------------------------------------------------------------------------------------------

#' @title Read one or multiple areas of China map from the DataV.GeoAtlas (Aliyun, China)
#' @author Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
#' @description This function is implemented to download one or multiple area(s) of China map from Aliyun DataV.GeoAtlas.
#' Undoubtedly, it is user friendly for Chinese since the boundary data is relatively accurate and is very likely to have
#' undergone official review. It is a better choice to use \code{microgeo::read_aliyun_map()} when you publishing several
#' graphics in official publications if you are unsure whether your map dataset is accurate.
#' @param adcode Adcode of the map area(s). For example, the `650000` means Xinjiang, and the `c(540000, 630000, 510000)`
#' represents Xizang, Qinghai, and Sichuan. You can visit the `http://datav.aliyun.com/portal/school/atlas/area_selector`
#' to select one or mutiple adcodes for your interested areas. Default is `c(540000, 630000, 510000)`.
#' @param crs.string Coordinate Reference System (CRS). Default is `"EPSG:4326"`.
#' @return A `SpatialPolygonsDataFrame`.
#' @seealso \code{\link[microgeo:plot_bmap]{microgeo::plot_bmap()}}
#' @examples
#' ## Example 1: using one adcode
#' map <- read_aliyun_map(adcode = 650000)
#' head(map@data)
#' map %>% plot_bmap()
#'
#' ## Example 1: using mutiple adcodes
#' map <- read_aliyun_map(adcode = c(540000, 630000, 510000))
#' head(map@data)
#' map %>% plot_bmap()
#' @export
read_aliyun_map = function(adcode = c(540000, 630000, 510000), crs.string = "EPSG:4326"){
    base.url <- "https://geo.datav.aliyun.com/areas_v3/bound/geojson?code="
    if (adcode %>% length == 1){
        map <- paste0(base.url, adcode) %>% sf::read_sf() %>% sf::as_Spatial()
        # address bugs in windows
        if (length(sf::st_is_valid(sf::st_as_sf(map))) > 1 || !sf::st_is_valid(sf::st_as_sf(map))){
            unio.sp <- rgeos::gUnaryUnion(map) %>% suppressWarnings() %>%
                suppressMessages() %>% sf::st_as_sf() %>% as(., "Spatial")
            map <- sp::SpatialPolygonsDataFrame(unio.sp, data = data.frame(map@data, row.names = rownames(map@data)))
        }
    }else{
        map <- lapply(adcode, function(code){
            map.tmp <- paste0(base.url, code) %>% sf::read_sf() %>% sf::as_Spatial()
            # address bugs in windows
            if (length(sf::st_is_valid(sf::st_as_sf(map.tmp))) > 1 || !sf::st_is_valid(sf::st_as_sf(map.tmp))){
                unio.sp <- rgeos::gUnaryUnion(map.tmp) %>% suppressWarnings() %>%
                    suppressMessages() %>% sf::st_as_sf() %>% as(., "Spatial")
                map.tmp <- sp::SpatialPolygonsDataFrame(unio.sp, data = data.frame(map.tmp@data, row.names = rownames(map.tmp@data)))
            }
            map.tmp <- map.tmp %>% sf::st_as_sf() %>% dplyr::as_tibble()
        }) %>% do.call('rbind', .) %>% sf::st_as_sf() %>% sf::as_Spatial()
    }
    map@data <- lapply(X = map@data %>% nrow %>% seq, function(x){
        data.frame(TYPE = 'DataV.GeoAtlas', FMTS = 'microgeo', NAME = map@data$name[x],
                   X.CENTER = map@data$centroid[[x]][1], Y.CENTER = map@data$centroid[[x]][2])
    }) %>% do.call("rbind", .)
    if (terra::crs(map %>% terra::vect(), proj = TRUE, describe = TRUE, parse = TRUE)[1, 6] %>% is.na){
        sp::proj4string(map) <- crs.string %>% sp::CRS()
    }else{
        map %<>% sp::spTransform(., crs.string %>% sp::CRS())
    }
    return(map)
}

#' @title Read a geographic map from the GADM (Global Administrative Areas)
#' @author Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
#' @description This function is designed to read a map from the GADM (`https://gadm.org/`) using \code{geodata::gadm()}.
#' Particularly, the Chinese territorial claims are not reflected in the boundary data provided by the GADM, and the data
#' for provincial/city/district boundaries may not necessarily be the most up-to-date versions. Thus, additional cautions
#' should be exercised when publishing graphics that utilize GADM data in official publications. The author of microgeo R
#' package is not responsible for any political disputes and instability caused by the use of this R package by anyone.
#' @param iso Three-letter ISO code or a full country name. Multiple names would be all downloaded and rbinded.
#' @param res An integer indicating the level of detail. Only for the version 4.1 of `GADM`. It must be either `1` (high)
#' or `2` (low). Default is `1`.
#' @param level The level of administrative subdivision requested. (starting with `0` for country, then `1` for the first
#' level of subdivision). Default is `1`.
#' @param version Either `latest` or GADM version number (can be `3.6`, `4.0` or `4.1`). Default is `latest`.
#' @param out.dir Path for storing the downloaded data. Default is `GADM`.
#' @param crs.string Coordinate Reference System (CRS). Default is `"EPSG:4326"`.
#' @param ... Parameters parsed by \code{geodata::gadm()}.
#' @return A `SpatialPolygonsDataFrame`.
#' @seealso
#' \code{\link[geodata:gadm]{geodata::gadm()}}
#' \code{\link[microgeo:plot_bmap]{microgeo::plot_bmap()}}
#' \code{\link[microgeo:trans_map_fmt]{microgeo::trans_map_fmt()}}
#' @examples
#' ## Example 1: the map of Australia
#' map <- read_gadm_map(iso = 'Australia', out.dir = 'test')
#' map %<>% trans_map_fmt(var = 'NAME_1')
#' head(map@data)
#' map %>% plot_bmap()
#'
#' ## Example 2: the map of Washington State, USA
#' map <- read_gadm_map(iso = 'USA', out.dir = 'test')
#' map %<>% terra::subset(map$NAME_1 == 'Washington')
#' map %<>% trans_map_fmt(var = 'NAME_1')
#' head(map@data)
#' map %>% plot_bmap()
#' @export
read_gadm_map = function(iso, res = 1, level = 1, version = 3.6, out.dir = "GADM", crs.string = "EPSG:4326", ...){
    map <- geodata::gadm(country = iso, resolution = res, level = level, version = version, path = out.dir, ...) %>%
        suppressWarnings()
    if (terra::crs(map) == '')
        stop("The crs is not avaliable in downloaded SpatVector! Try to change version as 3.6.")
    map <- map %>% sf::st_as_sf() %>% sf::as_Spatial()
    if (terra::crs(map %>% terra::vect(), proj = TRUE, describe = TRUE, parse = TRUE)[1, 6] %>% is.na){
        sp::proj4string(map) <- crs.string %>% sp::CRS()
    }else{
        map %<>% sp::spTransform(., crs.string %>% sp::CRS())
    }
    return(map)
}

#' @title Read a geographic map from the ESRI Shapefile
#' @author Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
#' @description This function is used to read a geographic map from the ESRI Shapefile via using \code{terra::vect()}.
#' @param shpfile A file path of ESRI Shapefile.
#' @param crs.string Coordinate Reference System (CRS). Default is `"EPSG:4326"`.
#' @param ... Parameters parsed by \code{terra::vect()}.
#' @return A `SpatialPolygonsDataFrame`.
#' @seealso
#' \code{\link[terra::vect]{terra::vect()}}
#' \code{\link[microgeo:plot_bmap]{microgeo::plot_bmap()}}
#' \code{\link[microgeo:trans_map_fmt]{microgeo::trans_map_fmt()}}
#' @examples
#' ##### Example 1: subset the research areas from a ESRI Shapefile of China map #####
#' map <- system.file("shapefiles/china-map", "china.shp", package = "microgeo") %>% read_shp_map()
#' map %<>% terra::subset(map$FENAME %in% c("Xizang Zizhiqu", "Qinghai Sheng", "Sichuan Sheng"))
#' map %<>% trans_map_fmt(var = 'FCNAME')
#' head(map@data)
#' map %>% plot_bmap()
#'
#' ##### Example 2: directly load the research areas from a ESRI Shapefile #####
#' map <- system.file("shapefiles/qtp-map", "DBATP_Polygon.shp", package = "microgeo") %>% read_shp_map()
#' map@data$NAME <- "Qinghai-Tibet Plateau" # Add a column to meet the demands of `microgeo::trans_map_fmt()`
#' map %<>% trans_map_fmt(var = 'NAME')
#' head(map@data)
#' map %>% plot_bmap()
#' @export
read_shp_map = function(shpfile, crs.string = "EPSG:4326", ...){
    map <- terra::vect(shpfile, ...) %>% sf::st_as_sf(y) %>% as(., "Spatial")
    if (terra::crs(map %>% terra::vect(), proj = TRUE, describe = TRUE, parse = TRUE)[1, 6] %>% is.na){
        sp::proj4string(map) <- crs.string %>% sp::CRS()
    }else{
        map %<>% sp::spTransform(., crs.string %>% sp::CRS())
    }
}

#' @title Convert the non-standardized `SpatialPolygonsDataFrame` to a microgeo-compatible one
#' @author Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
#' @description This function is designed to transform the `SpatialPolygonsDataFrame` to a microgeo-compatible one. It is
#' used for a `SpatialPolygonsDataFrame` returned by \code{microgeo::read_gadm_map()} or \code{microgeo::read_shp_map()}.
#' @param map A geographic map with the class of `SpatialPolygonsDataFrame`. Only accepts the `SpatialPolygonsDataFrame`
#' returned by \code{microgeo::read_gadm_map()} and \code{microgeo::read_shp_map()}.
#' @param var Which variable in the <map> would be used as the name(s) of research area(s)?
#' @return A `SpatialPolygonsDataFrame`.
#' @seealso
#' \code{\link[microgeo:read_gadm_map]{microgeo::read_gadm_map()}}
#' \code{\link[microgeo:read_shp_map]{microgeo::read_shp_map()}}
#' @examples
#' ## Example 1: the map of Australia
#' map <- read_gadm_map(iso = 'Australia', out.dir = 'test')
#' head(map@data)
#' map %<>% trans_map_fmt(var = 'NAME_1')
#' head(map@data)
#'
#' ## Example 2: the map of Washington State, USA
#' map <- read_gadm_map(iso = 'USA', out.dir = 'test')
#' head(map@data)
#' map %<>% terra::subset(map$NAME_1 == 'Washington')
#' head(map@data)
#' map %<>% trans_map_fmt(var = 'NAME_1')
#' head(map@data)
#' @export
trans_map_fmt = function(map, var){
    is.aliyun.map <- map %>% check_aliyun_map()
    if (is.aliyun.map) stop('Can not accept a `SpatialPolygonsDataFrame` returned by `read_aliyun_map()`!')
    map %>% check_map_var(., var)
    centroids <- map %>% sp::coordinates()
    map@data <- lapply(X = map@data %>% nrow %>% seq, function(x){
        data.frame(TYPE  = 'Non.DataV.GeoAtlas', FMTS  = 'microgeo', NAME = map@data[,var][x],
                   X.CENTER = centroids[x, 1], Y.CENTER = centroids[x, 2])
    }) %>% do.call("rbind", .)
    return(map)
}

#' @title Add grids to a map
#' @author Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
#' @description This function is used to add grids to a map based on a specified spatial resolution.
#' @param map A geographic map with the class of `SpatialPolygonsDataFrame`.
#' @param res Resolution of grids (minutes of a degree). A larger value would generate larger cells. Default is `1.5`.
#' @return A `SpatialPolygonsDataFrame`.
#' @seealso
#' \code{\link[microgeo:read_aliyun_map]{microgeo::read_aliyun_map()}}
#' \code{\link[microgeo:read_gadm_map]{microgeo::read_gadm_map()}}
#' \code{\link[microgeo:read_shp_map]{microgeo::read_shp_map()}}
#' \code{\link[microgeo:trans_map_fmt]{microgeo::trans_map_fmt()}}
#' \code{\link[microgeo:plot_bmap]{microgeo::plot_bmap()}}
#' @examples
#' ##### Example 1: add grids on the research areas subsetted from China map (`DataV.GeoAtlas`) #####
#' map <- read_aliyun_map(adcode = c("540000", "630000", "510000"))
#' head(map@data)
#' map %>% plot_bmap()
#' map %<>% grid_map(res = 1.5)
#' head(map@data)
#' map %>% plot_bmap()
#'
#' ##### Example 2: add grids on the map of Australia #####
#' map <- read_gadm_map(iso = 'Australia', out.dir = 'test')
#' map %<>% trans_map_fmt(var = 'NAME_1')
#' head(map@data)
#' map %>% plot_bmap()
#' map %<>% grid_map(res = 1.5)
#' head(map@data)
#' map %>% plot_bmap()
#'
#' ##### Example 3: add grids on the research areas subsetted from China map (`ESRI Shapefile`) #####
#' map <- system.file("shapefiles/china-map", "china.shp", package = "microgeo") %>% read_shp_map()
#' map %<>% terra::subset(map$FENAME %in% c("Xizang Zizhiqu", "Qinghai Sheng", "Sichuan Sheng"))
#' map %<>% trans_map_fmt(var = 'FCNAME')
#' head(map@data)
#' map %>% plot_bmap()
#' map %<>% grid_map(res = 1.5)
#' head(map@data)
#' map %>% plot_bmap()
#'
#' ##### Example 4: add grids on the research areas that are directly loaded from a `ESRI Shapefile` #####
#' map <- system.file("shapefiles/qtp-map", "DBATP_Polygon.shp", package = "microgeo") %>% read_shp_map()
#' map@data$NAME <- "Qinghai-Tibet Plateau" # Add a column to meet the demands of `microgeo::trans_map_fmt()`
#' map %<>% trans_map_fmt(var = 'NAME')
#' head(map@data)
#' map %>% plot_bmap()
#' map %<>% grid_map(res = 1.5)
#' head(map@data)
#' map %>% plot_bmap()
#' @export
grid_map = function(map, res = 1.5){
    map %>% check_mapdata()
    map.grid <- tryCatch({
        grid <- map %>% raster::extent() %>% raster::raster()
        raster::res(grid) <- res
        sp::proj4string(grid) <- map %>% sp::proj4string()
        gridpolygon <- grid %>% raster::rasterToPolygons()
        map.grid <- terra::intersect(map, gridpolygon)
        map.grid@data$grids <- map.grid@data %>% rownames()
        map.grid
    }, error = function(e){
        map <- sp::SpatialPolygonsDataFrame(map %>% rgeos::gUnaryUnion(), data = data.frame(TYPE = "Repaired")) %>%
            suppressWarnings() %>% suppressMessages()
        grid <- map %>% raster::extent() %>% raster::raster()
        raster::res(grid) <- res
        sp::proj4string(grid) <- map %>% sp::proj4string()
        gridpolygon <- grid %>% raster::rasterToPolygons()
        map.grid <- terra::intersect(map, gridpolygon)
        map.grid@data$grids <- map.grid@data %>% rownames()
        map.grid
    })
    map.grid@data <- map.grid@data[, colSums(map.grid@data %>% is.na) != map.grid@data %>% nrow] %>% na.omit()
    centroids <- sp::coordinates(map.grid)
    map.grid@data <- lapply(X = map.grid@data %>% nrow %>% seq, function(x){
        bbox.dat <- map.grid@polygons[[x]] %>% sp::bbox() %>% as.data.frame()
        data.frame(TYPE = "Gridded.Map", FMTS = 'microgeo', NAME = map.grid@data[x,]$grids,
                   X.CENTER = centroids[x, 1], Y.CENTER = centroids[x, 2])
    }) %>% do.call("rbind", .)
    return(map.grid)
}

#' @title Merge a `data.frame` with a map
#' @author Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
#' @description This function is used to merge a `data.frame` with a map. Values in this `data.frame` must be numeric.
#' @param map Geographic map with a class of `SpatialPolygonsDataFrame`.
#' @param dat A `data.frame` containing the variables to be merged with the <map>. Such a `data.frame` should not include
#' sample names (please set sample names as the row names). You can include any number of variables in this `data.frame`.
#' @param met A `data.frame` of sample information. Row names should be sample ids and the column names must be variables
#' (e.g., `longitude`, `latitude` and `group`). Longitude and latitude must be included in this `data.frame` as we mainly
#' focus on the spatial patterns of microbial traits.
#' @param med Method to aggregate numeric variables. Select one from `mean` and `median`. Default is `mean`.
#' @return A `SpatialPolygonsDataFrame`.
#' @seealso
#' \code{\link[microgeo:read_aliyun_map]{microgeo::read_aliyun_map()}}
#' \code{\link[microgeo:create_dataset]{microgeo::create_dataset()}}
#' \code{\link[microgeo:rarefy_count_table]{microgeo::rarefy_count_table()}}
#' \code{\link[microgeo:tidy_dataset]{microgeo::tidy_dataset()}}
#' \code{\link[microgeo:calc_alpha_div]{microgeo::calc_alpha_div()}}
#' \code{\link[microgeo:show_dataset]{microgeo::show_dataset()}}
#' \code{\link[microgeo:grid_map]{microgeo::grid_map()}}
#' @examples
#' # Create a standard microgeo dataset
#' data(qtp)
#' map <- read_aliyun_map(adcode = c(540000, 630000, 510000))
#' dataset.dts <- create_dataset(mat = qtp$asv, ant = qtp$tax, met = qtp$met, map = map,
#'                               phy = qtp$tre, env = qtp$env, lon = "longitude", lat = "latitude")
#' dataset.dts %<>% rarefy_count_table()
#' dataset.dts %<>% tidy_dataset()
#' dataset.dts %<>% calc_alpha_div(measures = c("observed", "shannon"))
#' dataset.dts %>% show_dataset()
#'
#' # Merge alpha diversity indices to a common map
#' common.map.mean <- merge_dfs_to_map(map = dataset.dts$map, dat = dataset.dts$div$alpha,
#'                                     met = dataset.dts$met, med = 'mean')
#' common.map.median <- merge_dfs_to_map(map = dataset.dts$map, dat = dataset.dts$div$alpha,
#'                                     met = dataset.dts$met, med = 'median')
#' View(common.map.mean@data)
#' View(common.map.median@data)
#'
#' # Merge alpha diversity indices to a gridded map
#' gridded.map <- dataset.dts$map %>% grid_map(res = 1.5)
#' gridded.map.mean <- merge_dfs_to_map(map = gridded.map, dat = dataset.dts$div$alpha,
#'                                      met = dataset.dts$met, med = 'mean')
#' gridded.map.median <- merge_dfs_to_map(map = gridded.map, dat = dataset.dts$div$alpha,
#'                                      met = dataset.dts$met, med = 'median')
#' View(gridded.map.mean@data)
#' View(gridded.map.median@data)
#' @export
merge_dfs_to_map = function(map, dat, met, med = c('mean', 'median')){

    # Check arguments for data merging
    map %>% check_mapdata()
    dat %>% check_merging_data(met = met, map = map)
    dat.use <- met[,c('longitude', 'latitude')] %>% cbind(., dat)
    med <- ifelse(med %>% length > 1, med[1], med)
    if (!med %in% c('mean', 'median')) 'The <med> must be one of `mean` and `median`!' %>% stop()
    var.names <- names(dat.use)[which(!dat.use %>% names %in% c('longitude', 'latitude'))]
    num.chk <- lapply(var.names, function(var) dat[,var] %>% is.numeric) %>% unlist()
    if (FALSE %in% num.chk) "All variables should be numeric!" %>% stop()

    # Add a field of `NAME` (in <map>) to <dat.use>
    dat.use$NAME <- rep(NA, dat.use %>% nrow)
    mapgeo <- map %>% sf::st_as_sf() %>% sf::st_geometry()
    for (i in map@data %>% nrow %>% seq){
        mapgeo.make.valid <- mapgeo[[i]] %>% sf::st_make_valid()
        test.rst <- lapply(X = dat.use %>% nrow %>% seq, function(x){
            point <- c(dat.use$longitude[x], dat.use$latitude[x]) %>% sf::st_point()
            is_within <- sf::st_within(point, mapgeo.make.valid)
            ifelse(is_within[[1]] %>% length == 0, 0, 1)
        }) %>% unlist()
        dat.use$NAME[which(test.rst == 1)] <- map@data[i,]$NAME
    }

    # Merge <dat.use> with the map
    grid.data.rs <- lapply(map@data$NAME, function(name){
        tmp <- dat.use[which(dat.use$NAME == name),]
        sample.num <- ifelse(tmp %>% nrow == 0, NA, tmp %>% nrow)
        sample.names <- ifelse(tmp %>% nrow == 0, NA, paste0(tmp %>% rownames, collapse = '|'))
        if (med == 'mean'){
            dd <- matrix(nrow = 1, ncol = var.names %>% length * 3) %>% as.data.frame()
            colnames(dd) <- c(paste0(var.names, "_mean"), paste0(var.names, "_sd"), paste0(var.names, "_se"))
        }else{
            dd <- matrix(ncol = var.names %>% length, nrow = 1) %>% as.data.frame()
            colnames(dd) <- var.names %>% paste0(., '_median')
        }
        dd <- data.frame(dd, sample.num = sample.num, sample.names = sample.names)
        if (tmp %>% nrow > 0){
            for (var in var.names){
                if (med == 'mean'){
                    val <- ifelse(sample.num > 1,
                                  paste0(mean(tmp[, var]), '_', sd(tmp[, var]), '_', sd(tmp[, var])/sqrt(nrow(tmp))),
                                  paste0(tmp[, var], "_", NA, "_", NA))
                    avg.rst <- strsplit(val, split = '_', fixed = T) %>% unlist() %>% as.numeric() %>% suppressWarnings()
                    dd[, paste0(var, '_mean')] <- avg.rst[1]
                    dd[, paste0(var, '_sd')] <- avg.rst[2]
                    dd[, paste0(var, '_se')] <- avg.rst[3]
                }else{
                    dd[, paste0(var, '_median')] <- median(tmp[, var])
                }
            }
        }
        dd
    }) %>% do.call('rbind', .)
    map@data <- cbind(map@data, grid.data.rs)
    if (map@data$sample.num %>% na.omit %>% sum != met %>% nrow) {
        paste0("only ", map@data$sample.num %>% na.omit %>% sum,
               " out of ", met %>% nrow, " samples fall with the map areas!") %>% warning()
    }
    return(map)
}

#' @title Merge a distance `matrix` with a map
#' @author Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
#' @description This function is used to merge a distance matrix with a map. Values in this `matrix` must be numeric.
#' @param map Geographic map with a class of `SpatialPolygonsDataFrame`.
#' @param dat A symmetric distance matrix.
#' @param met A `data.frame` of sample information. Row names should be sample ids and the column names must be variables
#' (e.g., `longitude`, `latitude` and `group`). Longitude and latitude must be included in this `data.frame` as we mainly
#' focus on the spatial patterns of microbial traits.
#' @param var A string to name the variable merged with a map. If it is `NULL`, the default, we will use the `dist_*`.
#' @param med Method to aggregate the variable. Select one from `mean`, `median`, `dtc.mean` (mean distances to centroid)
#' and `dtc.median` (median distances to centroid). Default is `mean`.
#' @return A `SpatialPolygonsDataFrame`.
#' @seealso
#' \code{\link[vegan:betadisper]{vegan:betadisper()}}
#' \code{\link[microgeo:read_aliyun_map]{microgeo::read_aliyun_map()}}
#' \code{\link[microgeo:create_dataset]{microgeo::create_dataset()}}
#' \code{\link[microgeo:rarefy_count_table]{microgeo::rarefy_count_table()}}
#' \code{\link[microgeo:tidy_dataset]{microgeo::tidy_dataset()}}
#' \code{\link[microgeo:calc_beta_div]{microgeo::calc_beta_div()}}
#' \code{\link[microgeo:show_dataset]{microgeo::show_dataset()}}
#' \code{\link[microgeo:grid_map]{microgeo::grid_map()}}
#' @examples
#' # Create a standard microgeo dataset
#' data(qtp)
#' map <- read_aliyun_map(adcode = c(540000, 630000, 510000))
#' dataset.dts <- create_dataset(mat = qtp$asv, ant = qtp$tax, met = qtp$met, map = map,
#'                               phy = qtp$tre, env = qtp$env, lon = "longitude", lat = "latitude")
#' dataset.dts %<>% rarefy_count_table()
#' dataset.dts %<>% tidy_dataset()
#' dataset.dts %<>% calc_beta_div(measures = c("bray", "jaccard"))
#' dataset.dts %>% show_dataset()
#'
#' # Merge distance matrix to a common map
#' map.mean <- merge_mtx_to_map(map = dataset.dts$map, dat = dataset.dts$div$beta$bray,
#'                              met = dataset.dts$met, var = 'bray', med = 'mean')
#' map.median <- merge_mtx_to_map(map = dataset.dts$map, dat = dataset.dts$div$beta$bray,
#'                                met = dataset.dts$met, var = 'bray', med = 'median')
#' map.dtc.mean <- merge_mtx_to_map(map = dataset.dts$map, dat = dataset.dts$div$beta$bray,
#'                                  met = dataset.dts$met, var = 'bray', med = 'dtc.mean')
#' map.dtc.median <- merge_mtx_to_map(map = dataset.dts$map, dat = dataset.dts$div$beta$bray,
#'                                  met = dataset.dts$met, var = 'bray', med = 'dtc.median')
#' View(map.mean@data)
#' View(map.median@data)
#' View(map.dtc.mean@data)
#' View(map.dtc.median@data)
#'
#' # Merge distance matrix to a gridded map
#' gridded.map <- dataset.dts$map %>% grid_map(res = 1.5)
#' gridded.map.mean <- merge_mtx_to_map(map = gridded.map, dat = dataset.dts$div$beta$bray,
#'                                      met = dataset.dts$met, var = 'bray', med = 'mean')
#' gridded.map.median <- merge_mtx_to_map(map = gridded.map, dat = dataset.dts$div$beta$bray,
#'                                        met = dataset.dts$met, var = 'bray', med = 'median')
#' gridded.map.dtc.mean <- merge_mtx_to_map(map = gridded.map, dat = dataset.dts$div$beta$bray,
#'                                          met = dataset.dts$met, var = 'bray', med = 'dtc.mean')
#' gridded.map.dtc.median <- merge_mtx_to_map(map = gridded.map, dat = dataset.dts$div$beta$bray,
#'                                            met = dataset.dts$met, var = 'bray', med = 'dtc.median')
#' View(gridded.map.mean@data)
#' View(gridded.map.median@data)
#' View(gridded.map.dtc.mean@data)
#' View(gridded.map.dtc.median@data)
#' @export
merge_mtx_to_map = function(map, dat, met, var = NULL, med = c('mean', 'median', 'dtc.mean', 'dtc.median')){

    # Check the arguments
    map %>% check_mapdata()
    dat %>% check_merging_data(met = met, map = map)
    med <- ifelse(med %>% length > 1, med[1], med)
    if (!med %in% c('mean', 'median', 'dtc.mean', 'dtc.median')){
        'The <med> must be one of `mean`, `median`, `dtc.mean` and `dtc.median`!' %>% stop()
    }
    if ((dat %>% is.matrix || dat %>% is.data.frame) && dat %>% as.matrix %>% unname %>% isSymmetric) {
        dat %<>% as.matrix()
    }else{
        'The <dat> must be a symmetric distance matrix!' %>% stop()
    }
    num.chk <- lapply(dat %>% colnames, function(var) dat[,var] %>% is.numeric) %>% unlist()
    if (FALSE %in% num.chk) "All values in <dat> should be numeric!" %>% stop()

    # Add a field of `NAME` (in <map>) to <met>
    met$NAME <- rep(NA, met %>% nrow)
    mapgeo <- map %>% sf::st_as_sf() %>% sf::st_geometry()
    for (i in map@data %>% nrow %>% seq){
        mapgeo.make.valid <- mapgeo[[i]] %>% sf::st_make_valid()
        test.rst <- lapply(X = met %>% nrow %>% seq, function(x){
            point <- c(met$longitude[x], met$latitude[x]) %>% sf::st_point()
            is_within <- sf::st_within(point, mapgeo.make.valid)
            ifelse(is_within[[1]] %>% length == 0, 0, 1)
        }) %>% unlist()
        met$NAME[which(test.rst == 1)] <- map@data[i,]$NAME
    }

    # Merge distance matrix with the map
    grid.data.rs <- lapply(map@data$NAME, function(name){
        tmp <- met[which(met$NAME == name),]
        sample.num   <- ifelse(tmp %>% nrow == 0, NA, tmp %>% nrow)
        sample.names <- ifelse(tmp %>% nrow == 0, NA, tmp %>% rownames %>% paste0(., collapse = '|'))
        if (med == 'median'){
            dd <- data.frame(dist_median = NA, sample.num = sample.num, sample.names = sample.names)
        }else if (med == 'dtc.median'){
            dd <- data.frame(dist_dtc.median = NA, sample.num = sample.num, sample.names = sample.names)
        }else if (med == 'mean'){
            dd <- data.frame(dist_mean = NA, dist_sd = NA, dist_se = NA, sample.num = sample.num,
                             sample.names = sample.names)
        }else if (med == 'dtc.mean'){
            dd <- data.frame(dist_dtc.mean = NA, dist_dtc.sd = NA, dist_dtc.se = NA, sample.num = sample.num,
                             sample.names = sample.names)
        }
        if (tmp %>% nrow > 0){
            samples <- met[which(met$NAME == name),] %>% rownames()
            dist.M <- dat[samples, samples]
            if (med == 'median'){
                dd$dist_median <- dist.M %>% as.dist() %>% median()
            }else if (med == 'dtc.median'){
                if (samples %>% length < 3){
                    dd$dist_dtc.median <- NA
                }else{
                    dd$dist_dtc.median <- vegan::betadisper(d = dist.M %>% as.dist,
                                                            group = rep("A", dist.M %>% nrow),
                                                            type = 'centroid')$distances %>% median()
                }
            }else if (med == 'mean'){
                dd$dist_mean <- dist.M %>% as.dist() %>% mean()
                dd$dist_sd <- dist.M %>% as.dist() %>% sd()
                dd$dist_se <- dd$dist_sd/sqrt(dist.M %>% as.dist %>% length)
            }else if (med == 'dtc.mean'){
                if (samples %>% length < 3){
                    dd$dist_dtc.mean <- dd$dist_dtc.sd <- dd$dist_dtc.se <- NA
                }else{
                    dd$dist_dtc.mean <- vegan::betadisper(d = dist.M %>% as.dist, group = rep("A", dist.M %>% nrow),
                                                          type = 'centroid')$distances %>% mean()
                    dd$dist_dtc.sd   <- vegan::betadisper(d = dist.M %>% as.dist, group = rep("A", dist.M %>% nrow),
                                                          type = 'centroid')$distances %>% sd()
                    dd$dist_dtc.se   <- dd$dist_dtc.sd/sqrt(
                        vegan::betadisper(d = dist.M %>% as.dist, group = rep("A", dist.M %>% nrow),
                                          type = 'centroid')$distances %>% length())
                }
            }
        }
        if (med == 'median' & !var %>% is.null){
            colnames(dd)[which(dd %>% colnames == 'dist_median')] <- paste0(var, '_median')
        }
        if (med == 'dtc.median' & !var %>% is.null){
            colnames(dd)[which(dd %>% colnames == 'dist_dtc.median')] <- paste0(var, '_dtc.median')
        }
        if (med == 'mean' & !var %>% is.null){
            colnames(dd)[which(dd %>% colnames == 'dist_mean')] <- paste0(var, '_mean')
            colnames(dd)[which(dd %>% colnames == 'dist_sd')] <- paste0(var, '_sd')
            colnames(dd)[which(dd %>% colnames == 'dist_se')] <- paste0(var, '_se')
        }
        if (med == 'dtc.mean' & !var %>% is.null){
            colnames(dd)[which(dd %>% colnames == 'dist_dtc.mean')] <- paste0(var, '_dtc.mean')
            colnames(dd)[which(dd %>% colnames == 'dist_dtc.sd')] <- paste0(var, '_dtc.sd')
            colnames(dd)[which(dd %>% colnames == 'dist_dtc.se')] <- paste0(var, '_dtc.se')
        }
        dd
    }) %>% do.call('rbind', .)
    map@data <- cbind(map@data, grid.data.rs)
    if (map@data$sample.num %>% na.omit %>% sum != met %>% nrow) {
        paste0("only ", map@data$sample.num %>% na.omit %>% sum,
               " out of ", met %>% nrow, " samples fall with the map areas!") %>% warning()
    }
    return(map)
}

#' @title Extract the metadata table from a geographic map
#' @author Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
#' @description This function is designed to extract metadata table from a geographic map. The results returned by such a
#' function would be useful for the statistical analysis based on administrative regions or grid names.
#' @param map A geographic map with the class of `SpatialPolygonsDataFrame`.
#' @param met Original metadata table. See the examples for more details about such a table.
#' @return A `data.frame` of new metadata table extracted from the geographic map.
#' @seealso
#' \code{\link[microgeo:read_aliyun_map]{microgeo::read_aliyun_map()}}
#' \code{\link[microgeo:grid_map]{microgeo::grid_map()}}
#' \code{\link[microgeo:create_dataset]{microgeo::create_dataset()}}
#' \code{\link[microgeo:rarefy_count_table]{microgeo::rarefy_count_table()}}
#' \code{\link[microgeo:tidy_dataset]{microgeo::tidy_dataset()}}
#' \code{\link[microgeo:calc_alpha_div]{microgeo::calc_alpha_div()}}
#' \code{\link[microgeo:show_dataset]{microgeo::show_dataset()}}
#' \code{\link[microgeo:merge_dfs_to_map]{microgeo::merge_dfs_to_map()}}
#' @examples
#' ##### Example 1: extract the metadata from a map retrieved from the DataV.GeoAtlas #####
#' data(qtp)
#' names(qtp)
#' head(qtp$met)
#' map <- read_aliyun_map(adcode = c(540000, 630000, 510000))
#' metadata <- map %>% extract_metadata_from_map(met = qtp$met)
#' head(metadata) # Now we can compare any metrics among administrative regions based on such a new metadata table
#'
#' ##### Example 2: extract the metadata from a gridded map #####
#' map %<>% grid_map(res = 1.5)
#' metadata <- map %>% extract_metadata_from_map(met = qtp$met)
#' head(metadata) # Now we can compare any metrics among grids based on such a new metadata table
#'
#' ##### Example 3: extract the metadata from a map with additional data #####
#' dataset.dts <- create_dataset(mat = qtp$asv, ant = qtp$tax, met = qtp$met, map = map,
#'                               phy = qtp$tre, env = qtp$env, lon = "longitude", lat = "latitude")
#' dataset.dts %<>% rarefy_count_table()
#' dataset.dts %<>% tidy_dataset()
#' dataset.dts %<>% calc_alpha_div(measures = c("observed", "shannon"))
#' dataset.dts %>%  show_dataset()
#' map <- merge_dfs_to_map(map = dataset.dts$map, dat = dataset.dts$div$alpha,
#'                         met = dataset.dts$met, agg.method = 'mean')
#' metadata <- map %>% extract_metadata_from_map(met = dataset.dts$met)
#' head(metadata) # Now we can compare the alpha diversity metrics among administrative regions
#' @export
extract_metadata_from_map = function(map, met){
    map %>% check_mapdata(); met %>% check_metdata()
    met$NAME <- met %>% nrow() %>% rep(NA, .)
    mapgeo <- map %>% sf::st_as_sf() %>% sf::st_geometry()
    for (i in map@data %>% nrow %>% seq){
        mapgeo.make.valid <- mapgeo[[i]] %>% sf::st_make_valid()
        test.rst <- lapply(X = met %>% nrow %>% seq, function(x){
            is_within <- c(met$longitude[x], met$latitude[x]) %>% sf::st_point() %>% sf::st_within(., mapgeo.make.valid)
            ifelse(is_within[[1]] %>% length == 0, 0, 1)
        }) %>% unlist()
        met$NAME[which(test.rst == 1)] <- map@data[i,]$NAME
    }
    supp.dat <- lapply(met$NAME, function(name){
        tmpdat <- map@data[which(map@data$NAME == name),]
        if ('sample.names' %in% colnames(tmpdat)) tmpdat <- tmpdat[,-which(tmpdat %>% colnames == 'sample.names')]
        tmpdat[,-which(tmpdat %>% colnames == 'NAME')]
    }) %>% do.call('rbind', .)
    met.new <- cbind(met, supp.dat)
    return(met.new)
}
