# -----------------------------------------------------------------------------------------------------------------------
# Copyright (c) 2023, microgeo/Chaonan Li (licn@mtc.edu.cn).                                                            #
# The microgeo is distributed under the terms of the Modified BSD License.                                              #
# Full license is avaliable in the file LICENSE, distributed with this package.                                         #
# -----------------------------------------------------------------------------------------------------------------------

#' @title Nearest neighbour interpolation
#' @author Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
#' @description This function is used to perform nearest neighbour interpolation (Voronoi diagram).
#' @param map A `SpatialPolygonsDataFrame` of geographic map.
#' @param met A `data.frame` of sample information. Row names should be sample ids and the column names must be variables
#' (e.g., `longitude`, `latitude` and `group`). Longitude and latitude must be included in this `data.frame` as we mainly
#' focus on the spatial patterns of microbial traits.
#' @param dat A `data.frame` containing the variables to be interpolated.
#' @param var Name of variable to be interpolated.
#' @param trim.dup Whether to randomly remove sampling sites with duplicated coordinates? Default is `FALSE`.
#' @return An R `list` with the following components:
#' \describe{
#'   \item{\code{object$res}}{Data for ggplot2 visualization.}
#'   \item{\code{object$var}}{Name of interpolated variable.}
#'   \item{\code{object$type}}{Type of interpolation, which is `'nen'`.}
#' }
#' @seealso
#' \code{\link[sf:st_voronoi]{sf::st_voronoi()}}
#' \code{\link[microgeo:read_aliyun_map]{microgeo::read_aliyun_map()}}
#' \code{\link[microgeo:create_dataset]{microgeo::create_dataset()}}
#' \code{\link[microgeo:show_dataset]{microgeo::show_dataset()}}
#' \code{\link[microgeo:rarefy_count_table]{microgeo::rarefy_count_table()}}
#' \code{\link[microgeo:tidy_dataset]{microgeo::tidy_dataset()}}
#' \code{\link[microgeo:calc_alpha_div]{microgeo::calc_alpha_div()}}
#' \code{\link[microgeo:plot_nmap]{microgeo::plot_nmap()}}
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
#'
#' # Tidy microgeo dataset, and calculate microbial alpha diversity indices
#' dataset.dts %<>% rarefy_count_table()
#' dataset.dts %<>% tidy_dataset()
#' dataset.dts %<>% calc_alpha_div(measures = c("observed", "shannon"))
#' dataset.dts %>% show_dataset()
#'
#' # Nearest neighbour interpolation for alpha diversity indices
#' observed.rst <- interp_nen(map = dataset.dts$map, met = dataset.dts$met,
#'                            dat = dataset.dts$div$alpha, var = 'observed', trim.dup = TRUE)
#' shannon.rst <- interp_nen(map = dataset.dts$map, met = dataset.dts$met,
#'                           dat = dataset.dts$div$alpha, var = 'shannon', trim.dup = TRUE)
#'
#' # Visualize the interpolated results
#' observed.rst %>% plot_nmap() %>% add_north_arrow() %>% add_scale_bar() %>% add_crs()
#' shannon.rst %>% plot_nmap() %>% add_north_arrow() %>% add_scale_bar() %>% add_crs()
#' @export
interp_nen = function(map, met, dat, var, trim.dup = FALSE){
    map %>% check_mapdata(); map.sfs <- sf::st_as_sf(map)
    map.crs <- terra::crs(terra::vect(map), proj = TRUE, describe = TRUE, parse = TRUE)[1, 6]
    if (length(sf::st_is_valid(sf::st_as_sf(map))) > 1 || !sf::st_is_valid(sf::st_as_sf(map))){
        map.sfs <- rgeos::gUnaryUnion(map) %>% suppressMessages() %>% suppressWarnings() %>% sf::st_as_sf()
    }
    map.sfs %<>% sf::st_simplify(preserveTopology = TRUE, dTolerance = 200)
    use.dat <- get_interpolation_data(met = met, dat = dat, var = var, trim.dup = trim.dup)
    use.dat.sf <- sf::st_as_sf(use.dat, coords = c("longitude", "latitude"), crs = map.crs)
    use.dat.sf.gm <- sf::st_union(use.dat.sf)
    use.dat_sf_voronoi <- sf::st_voronoi(use.dat.sf.gm) %>%
        sf::st_sf() %>%
        sf::st_cast() %>%
        sf::st_intersection(., sf::st_union(sf::st_buffer(map.sfs, dist = 500))) %>%
        sf::st_sf() %>%
        sf::st_join(use.dat.sf) %>% suppressWarnings()
    res <- list(res = use.dat_sf_voronoi, var = var, type = 'nen')
    return(res)
}

#' @title Polynomial fit (2nd order) interpolation
#' @author Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
#' @description This function is used to perform polynomial fit (2nd order) interpolation.
#' @param map A `SpatialPolygonsDataFrame` of geographic map.
#' @param met A `data.frame` of sample information. Row names should be sample ids and the column names must be variables
#' (e.g., `longitude`, `latitude` and `group`). Longitude and latitude must be included in this `data.frame` as we mainly
#' focus on the spatial patterns of microbial traits.
#' @param dat A `data.frame` containing the variables to be interpolated.
#' @param var Name of variable to be interpolated.
#' @param n Approximate sample size. Default is `4000`. See \code{sp::spsample()} for more details.
#' @param trim.dup Whether to randomly remove sampling sites with duplicated coordinates? Default is `FALSE`.
#' @return A `SpatRaster`.
#' @seealso
#' \code{\link[stats:lm]{stats::lm()}}
#' \code{\link[sp:spsample]{sp::spsample()}}
#' \code{\link[microgeo:read_aliyun_map]{microgeo::read_aliyun_map()}}
#' \code{\link[microgeo:create_dataset]{microgeo::create_dataset()}}
#' \code{\link[microgeo:show_dataset]{microgeo::show_dataset()}}
#' \code{\link[microgeo:rarefy_count_table]{microgeo::rarefy_count_table()}}
#' \code{\link[microgeo:tidy_dataset]{microgeo::tidy_dataset()}}
#' \code{\link[microgeo:calc_alpha_div]{microgeo::calc_alpha_div()}}
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
#' # Tidy microgeo dataset, and calculate microbial alpha diversity indices
#' dataset.dts %<>% rarefy_count_table()
#' dataset.dts %<>% tidy_dataset()
#' dataset.dts %<>% calc_alpha_div(measures = c("observed", "shannon"))
#' dataset.dts %>% show_dataset()
#'
#' # 2nd polynomial fit interpolation for alpha diversity indices
#' pol.rst.observed <- interp_pol(map = dataset.dts$map, met = dataset.dts$met,
#'                                dat = dataset.dts$div$alpha, var = 'observed', trim.dup = TRUE)
#' pol.rst.shannon  <- interp_pol(map = dataset.dts$map, met = dataset.dts$met,
#'                                dat = dataset.dts$div$alpha, var = 'shannon', trim.dup = TRUE)
#'
#' # Visualize the interpolated results
#' dataset.dts$map %>% plot_bmap() %>%
#'     add_spatraster(spat.raster = pol.rst.observed) %>%
#'     add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = "X.CENTER", lat.var = "Y.CENTER") %>%
#'     add_scale_bar() %>% add_north_arrow() %>% add_crs()
#' dataset.dts$map %>% plot_bmap() %>%
#'     add_spatraster(spat.raster = pol.rst.shannon) %>%
#'     add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = "X.CENTER", lat.var = "Y.CENTER") %>%
#'     add_scale_bar() %>% add_north_arrow() %>% add_crs()
#' @export
interp_pol = function(map, met, dat, var, n = 4000, trim.dup = FALSE){
    map %>% check_mapdata(); map.sfs <- sf::st_as_sf(map)
    map.crs <- terra::crs(terra::vect(map), proj = TRUE, describe = TRUE, parse = TRUE)[1, 6]
    if (length(sf::st_is_valid(sf::st_as_sf(map))) > 1 || !sf::st_is_valid(sf::st_as_sf(map))){
        map.sfs <- rgeos::gUnaryUnion(map) %>% suppressMessages() %>% suppressWarnings() %>% sf::st_as_sf()
    }
    map.sfs %<>% sf::st_simplify(preserveTopology = TRUE, dTolerance = 200)
    map.sp <- sf::as_Spatial(map.sfs)
    grd <- as.data.frame(sp::spsample(map.sp, "regular", n = n))
    names(grd) <- c("X", "Y"); sp::coordinates(grd) <- c("X", "Y")
    sp::gridded(grd) <- TRUE; sp::fullgrid(grd) <- TRUE
    use.dat <- get_interpolation_data(met = met, dat = dat, var = var, trim.dup = trim.dup)
    use.dat.sf <- sf::st_as_sf(use.dat, coords = c("longitude", "latitude"), crs = map.crs)
    use.dat.sp <- sf::as_Spatial(use.dat.sf); use.dat.sf.xy <- as(use.dat.sf, "Spatial")
    use.dat.sf.xy <- sp::spTransform(use.dat.sf.xy, terra::crs(use.dat.sp))
    sp::proj4string(grd) <- as.character(terra::crs(use.dat.sp))
    use.dat.sf.xy@bbox <- sf::as_Spatial(map.sfs)@bbox
    use.dat.sf.xy$X <- sp::coordinates(use.dat.sf.xy)[,1]
    use.dat.sf.xy$Y <- sp::coordinates(use.dat.sf.xy)[,2]
    f <- as.formula(target ~ X + Y + I(X*X)+I(Y*Y) + I(X*Y))
    lm.rst <- lm(f, data = use.dat.sf.xy)
    dat.res <- sp::SpatialGridDataFrame(grd, data.frame(var1.pred = predict(lm.rst, newdata = grd))) # regression model
    r <- raster::raster(dat.res); r.m <- raster::mask(r, sf::as_Spatial(map.sfs)) %>% terra::rast()
    terra::crs(r.m) <- as.character(terra::crs(map)); names(r.m) <- var
    return(r.m)
}

#' @title Inverse distance weighting (IDW) interpolation
#' @author Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
#' @description This function is implemented to perform the inverse distance weighting (IDW) interpolation.
#' @param map A `SpatialPolygonsDataFrame` of geographic map.
#' @param met A `data.frame` of sample information. Row names should be sample ids and the column names must be variables
#' (e.g., `longitude`, `latitude` and `group`). Longitude and latitude must be included in this `data.frame` as we mainly
#' focus on the spatial patterns of microbial traits.
#' @param dat A `data.frame` containing the variables to be interpolated.
#' @param var Name of variable to be interpolated.
#' @param grid.step Step for grid expanding (minutes of a degree). Default is `0.1`.
#' @param type Please select one from `regular` and `hexagonal`. `regular` for regular (systematically aligned) sampling,
#' while `hexagonal` for sampling on a hexagonal lattice. Default is `regular`. See \code{sp::spsample()} for details.
#' @param n Approximate sample size. A large value represents a high resolution. Only work when `type='regular'`. Default
#' is `50000`. See \code{sp::spsample()} for more details.
#' @param cellsize Cell size. Only work when `type='hexagonal'` Default is `1`. See \code{sp::spsample()} for details.
#' @param idp Specify the inverse distance weighting power. Default is `2`. See \code{gstat::idw()} for details.
#' @param trim.dup Whether to randomly remove sampling sites with duplicated coordinates? Default is `FALSE`.
#' @param ... Parameters parsed by \code{gstat::idw()}.
#' @return
#' A `SpatRaster` if the <type> is `regular`; An R `list` with the following components if the <type> is `hexagonal`.
#' \describe{
#'   \item{\code{object$res}}{Data for ggplot2 visualization.}
#'   \item{\code{object$var}}{Name of interpolated variable.}
#'   \item{\code{object$type}}{Type of interpolation, which is `'idw_hex'`.}
#' }
#' @seealso
#' \code{\link[sp:spsample]{sp::spsample()}}
#' \code{\link[gstat:idw]{gstat::idw()}}
#' \code{\link[microgeo:read_aliyun_map]{microgeo::read_aliyun_map()}}
#' \code{\link[microgeo:create_dataset]{microgeo::create_dataset()}}
#' \code{\link[microgeo:show_dataset]{microgeo::show_dataset()}}
#' \code{\link[microgeo:rarefy_count_table]{microgeo::rarefy_count_table()}}
#' \code{\link[microgeo:tidy_dataset]{microgeo::tidy_dataset()}}
#' \code{\link[microgeo:calc_alpha_div]{microgeo::calc_alpha_div()}}
#' \code{\link[microgeo:plot_bmap]{microgeo::plot_bmap()}}
#' \code{\link[microgeo:plot_imap]{microgeo::plot_bmap()}}
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
#' # Tidy microgeo dataset, and calculate microbial alpha diversity indices
#' dataset.dts %<>% rarefy_count_table()
#' dataset.dts %<>% tidy_dataset()
#' dataset.dts %<>% calc_alpha_div(measures = c("observed", "shannon"))
#' dataset.dts %>% show_dataset()
#'
#' # Inverse distance weighting (IDW) interpolation for alpha diversity indices(type = 'regular')
#' idw.rst.observed <- interp_idw(map = dataset.dts$map, met = dataset.dts$met,
#'                                dat = dataset.dts$div$alpha, var = 'observed',
#'                                type = 'regular', trim.dup = TRUE)
#' idw.rst.shannon  <- interp_idw(map = dataset.dts$map, met = dataset.dts$met,
#'                                dat = dataset.dts$div$alpha, var = 'shannon',
#'                                type = 'regular', trim.dup = TRUE)
#'
#' # Visualize the interpolated results
#' dataset.dts$map %>% plot_bmap() %>%
#'     add_spatraster(spat.raster = idw.rst.observed) %>%
#'     add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = "X.CENTER", lat.var = "Y.CENTER") %>%
#'     add_scale_bar() %>% add_north_arrow() %>% add_crs()
#' dataset.dts$map %>% plot_bmap() %>%
#'     add_spatraster(spat.raster = idw.rst.observed) %>%
#'     add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = "X.CENTER", lat.var = "Y.CENTER") %>%
#'     add_scale_bar() %>% add_north_arrow() %>% add_crs()
#'
#' # Inverse distance weighting (IDW) interpolation for alpha diversity indices (type = 'hexagonal')
#' observed.hex <- interp_idw(map = dataset.dts$map, met = dataset.dts$met,
#'                            dat = dataset.dts$div$alpha, var = 'observed',
#'                            type = 'hexagonal', trim.dup = TRUE)
#' shannon.hex <- interp_idw(map = dataset.dts$map, met = dataset.dts$met,
#'                           dat = dataset.dts$div$alpha, var = 'shannon',
#'                           type = 'hexagonal', trim.dup = TRUE)
#'
#' # Visualize the interpolated results
#' observed.hex %>% plot_imap() %>% add_scale_bar() %>% add_north_arrow() %>% add_crs()
#' shannon.hex %>% plot_imap() %>% add_scale_bar() %>% add_north_arrow() %>% add_crs()
#' @export
interp_idw = function(map, met, dat, var, grid.step = 0.1, type = c('regular', 'hexagonal'), n = 50000,
                      cellsize = 1, idp = 2, trim.dup = FALSE, ...){
    map %>% check_mapdata(); type <- ifelse(length(type) > 1, type[1], type)
    if (!type %in% c('regular', 'hexagonal')) stop('<type> must be one of `regular` and `hexagonal`!')
    map.crs <- terra::crs(terra::vect(map), proj = TRUE, describe = TRUE, parse = TRUE)[1, 6]
    map.sfs <- sf::st_as_sf(map)
    if (length(sf::st_is_valid(sf::st_as_sf(map))) > 1 || !sf::st_is_valid(sf::st_as_sf(map))){
        map.sfs <- rgeos::gUnaryUnion(map) %>% suppressMessages() %>% suppressWarnings() %>% sf::st_as_sf()
    }
    map.sfs %<>% sf::st_simplify(preserveTopology = TRUE, dTolerance = 200)
    grd.extent <- sf::st_bbox(map.sfs)
    x.range <- as.numeric(c(grd.extent[1], grd.extent[3])) # min/max longitude of the interpolation area
    y.range <- as.numeric(c(grd.extent[2], grd.extent[4])) # min/max latitude of the interpolation area
    grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = grid.step),
                       y = seq(from = y.range[1], to = y.range[2], by = grid.step))
    sp::coordinates(grd) <- ~x + y; sp::gridded(grd) <- TRUE
    use.dat <- get_interpolation_data(met = met, dat = dat, var = var, trim.dup = trim.dup)
    use.dat.sf <- sf::st_as_sf(use.dat, coords = c("longitude", "latitude"), crs = map.crs)
    use.dat.sp <- sf::as_Spatial(use.dat.sf); map.sp <- sf::as_Spatial(map.sfs)
    if (type == 'regular'){
        idw.grid <- sp::spsample(map.sp, type = "regular", n = n)
        idw.grid <- sp::spTransform(idw.grid, terra::crs(use.dat.sp))
        rslt.pt <- gstat::idw(target ~ 1, use.dat.sp, newdata = idw.grid, idp = idp, ...)
        rslt.pt@data <- data.frame(row.names = seq(nrow(rslt.pt@data)),
                                   var1.pred = rslt.pt@data[,which(colnames(rslt.pt@data) != 'var1.var')])
        rslt.pt %<>% raster::rasterFromXYZ() %>% terra::rast()
        terra::crs(rslt.pt) <- terra::crs(map) %>% as.character()
        names(rslt.pt) <- var; return(rslt.pt)
    }
    if (type == 'hexagonal'){
        grid.hex <- sp::spsample(map.sp, type = "hexagonal", cellsize = cellsize)
        grid.hex <- sp::HexPoints2SpatialPolygons(grid.hex)
        grid.hex <- sp::spTransform(grid.hex, terra::crs(use.dat.sp))
        p.idw.hex <- gstat::idw(target ~ 1, use.dat.sp, newdata = grid.hex, idp = idp, ...)
        rslt.hex <- sf::st_as_sf(p.idw.hex)
        colnames(rslt.hex)[which(colnames(rslt.hex) == 'var1.pred')] <- 'target'
        list(res = rslt.hex, var = var, type = 'idw_hex') %>% return()
    }
}

#' @title Kriging interpolation
#' @author Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
#' @description This function is implemented to perform the kriging interpolation.
#' @param map A `SpatialPolygonsDataFrame` of geographic map.
#' @param met A `data.frame` of sample information. Row names should be sample ids and the column names must be variables
#' (e.g., `longitude`, `latitude` and `group`). Longitude and latitude must be included in this `data.frame` as we mainly
#' focus on the spatial patterns of microbial traits.
#' @param dat A `data.frame` containing the variables to be interpolated.
#' @param var Name of variable to be interpolated.
#' @param model Model for kriging interpolation. e.g., `Mat`, `Exp`, `Sph` and `Gau`. Calling \code{gstat::vgm()} without
#' a argument returns a `data.frame` with available models. See \code{gstat::vgm()} for more details. Default is `Mat`.
#' @param n Approximate sample size. A large value means a high resolution. Default is `4000`. See \code{sp::spsample()}.
#' @param trim.dup Whether to randomly remove sampling sites with duplicated coordinates? Default is `FALSE`.
#' @param test.model Whether to test selected model (<model>) rather than directly perform kriging interpolation? Default
#' is `FALSE`.
#' @param ... Parameters parsed by \code{gstat::krige()}.
#' @return A `SpatRaster` if the <test.model> is `FALSE`.
#' @seealso
#' \code{\link[sp:spsample]{sp::spsample()}}
#' \code{\link[gstat:vgm]{gstat::vgm()}}
#' \code{\link[gstat:fit.variogram]{gstat::fit.variogram()}}
#' \code{\link[gstat:krige]{gstat::krige()}}
#' \code{\link[microgeo:read_aliyun_map]{microgeo::read_aliyun_map()}}
#' \code{\link[microgeo:create_dataset]{microgeo::create_dataset()}}
#' \code{\link[microgeo:show_dataset]{microgeo::show_dataset()}}
#' \code{\link[microgeo:rarefy_count_table]{microgeo::rarefy_count_table()}}
#' \code{\link[microgeo:tidy_dataset]{microgeo::tidy_dataset()}}
#' \code{\link[microgeo:calc_alpha_div]{microgeo::calc_alpha_div()}}
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
#' # Tidy microgeo dataset, and calculate microbial alpha diversity indices
#' dataset.dts %<>% rarefy_count_table()
#' dataset.dts %<>% tidy_dataset()
#' dataset.dts %<>% calc_alpha_div(measures = c("observed", "shannon"))
#' dataset.dts %>% show_dataset()
#'
#' # Test the `Sph` and `Mat` models
#' interp_kri(map = dataset.dts$map, met = dataset.dts$met, dat = dataset.dts$div$alpha,
#'            var = 'observed', model = 'Sph', test.model = TRUE) # error
#' interp_kri(map = dataset.dts$map, met = dataset.dts$met, dat = dataset.dts$div$alpha,
#'            var = 'observed', model = 'Sph', test.model = TRUE, trim.dup = TRUE)
#' interp_kri(map = dataset.dts$map, met = dataset.dts$met, dat = dataset.dts$div$alpha,
#'            var = 'shannon', model = 'Mat', test.model = TRUE) # error
#' interp_kri(map = dataset.dts$map, met = dataset.dts$met, dat = dataset.dts$div$alpha,
#'            var = 'shannon', model = 'Mat', test.model = TRUE, trim.dup = TRUE)
#'
#' # Perform kriging interpolation for alpha diversity indices
#' kri.rst.observed <- interp_kri(map = dataset.dts$map, met = dataset.dts$met,
#'                                dat = dataset.dts$div$alpha, var = 'observed',
#'                                model = 'Sph', trim.dup = TRUE)
#' kri.rst.shannon <- interp_kri(map = dataset.dts$map, met = dataset.dts$met,
#'                                dat = dataset.dts$div$alpha, var = 'shannon',
#'                                model = 'Mat', trim.dup = TRUE)
#'
#' # Visualize the interpolated results
#' dataset.dts$map %>% plot_bmap() %>%
#'     add_spatraster(spat.raster = kri.rst.observed) %>%
#'     add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = "X.CENTER", lat.var = "Y.CENTER") %>%
#'     add_scale_bar() %>% add_north_arrow() %>% add_crs()
#' dataset.dts$map %>% plot_bmap() %>%
#'     add_spatraster(spat.raster = kri.rst.shannon) %>%
#'     add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = "X.CENTER", lat.var = "Y.CENTER") %>%
#'     add_scale_bar() %>% add_north_arrow() %>% add_crs()
#' @export
interp_kri = function(map, met, dat, var, model = c('Mat', 'Exp', 'Sph', 'Gau'), n = 4000,
                      trim.dup = FALSE, test.model = FALSE, ...){

    # prepare data for kriging interpolation
    map %>% check_mapdata(); model <- ifelse(length(model) > 1, model[1], model)
    #if (!model %in% c('Mat', 'Exp', 'Sph', 'Gau')) stop('The <model> must be one of `Mat`, `Exp`, `Sph` and `Gau`!')
    #map.crs <- terra::crs(terra::vect(map), describe = TRUE)[1, 3] %>% as.numeric()
    map.crs <- terra::crs(terra::vect(map), proj = TRUE, describe = TRUE, parse = TRUE)[1, 6]
    map.sfs <- sf::st_as_sf(map)
    if (length(sf::st_is_valid(sf::st_as_sf(map))) > 1 || !sf::st_is_valid(sf::st_as_sf(map))){
        map.sfs <- rgeos::gUnaryUnion(map) %>% suppressMessages() %>% suppressWarnings() %>% sf::st_as_sf()
    }
    map.sfs %<>% sf::st_simplify(preserveTopology = TRUE, dTolerance = 200)
    map.sp <- sf::as_Spatial(map.sfs)
    grd <- as.data.frame(sp::spsample(map.sp, "regular", n = n))
    names(grd) <- c("X", "Y"); sp::coordinates(grd) <- c("X", "Y")
    sp::gridded(grd) <- TRUE; sp::fullgrid(grd) <- TRUE
    use.dat <- get_interpolation_data(met = met, dat = dat, var = var, trim.dup = trim.dup)
    use.dat.sf <- sf::st_as_sf(use.dat, coords = c("longitude", "latitude"), crs = map.crs)
    use.dat.sp <- sf::as_Spatial(use.dat.sf); use.dat.sf.xy <- as(use.dat.sf, "Spatial")
    use.dat.sf.xy <- sp::spTransform(use.dat.sf.xy, terra::crs(use.dat.sp))
    sp::proj4string(grd) <- as.character(terra::crs(use.dat.sp))
    use.dat.sf.xy@bbox <- sf::as_Spatial(map.sfs)@bbox
    use.dat.sf.xy$X <- sp::coordinates(use.dat.sf.xy)[,1]
    use.dat.sf.xy$Y <- sp::coordinates(use.dat.sf.xy)[,2]

    # define and fit model
    p.semivariog <- gstat::variogram(target ~ 1, locations = use.dat.sf.xy, data = use.dat.sf.xy) # plot(p.semivariog)
    psill  <- p.semivariog$gamma[nrow(p.semivariog)]
    nugget <- p.semivariog$gamma[1]
    range  <- ceiling(p.semivariog$dist[nrow(p.semivariog)])
    p.model.variog <- gstat::vgm(psill = psill, model = model, nugget = nugget, range = range)
    p.fit.variog <- gstat::fit.variogram(p.semivariog, p.model.variog) %>% suppressWarnings()
    if (test.model){
        plot(p.semivariog, p.fit.variog)
    }else{
        p.krig <- gstat::krige(formula = target ~ 1, locations = use.dat.sf.xy, newdata = grd,
                               model = p.model.variog, ...)
        r <- raster::raster(p.krig); r.m <- raster::mask(r, sf::as_Spatial(map.sfs)) %>% terra::rast()
        terra::crs(r.m) <- as.character(terra::crs(map)); names(r.m) <- var
        return(r.m)
    }
}
