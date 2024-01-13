# -----------------------------------------------------------------------------------------------------------------------
# Copyright (c) 2023, microgeo/Chaonan Li (licn@mtc.edu.cn).                                                            #
# The microgeo is distributed under the terms of the Modified BSD License.                                              #
# Full license is avaliable in the file LICENSE, distributed with this package.                                         #
# -----------------------------------------------------------------------------------------------------------------------

#' @description Show an error message on the screen.
#' @param msg An error message to be shown.
stop = function(msg) msg %>% cli::cli_abort(call. = TRUE)

#' @description Show a warning message on the screen.
#' @param msg A warning message to be shown.
warning = function(msg){
    list(color = "#FF0000") %>% list(span.emph = .) %>% cli::cli_div(theme = .)
    prefix <- "WARN"; paste0("[", Sys.time(), "] {.emph {prefix}} ==> ", msg) %>% cli::cli_alert_warning()
}

#' @description Show a common message on the screen.
#' @param msg A common message to be shown.
show_comm_msg = function(msg){
    list(color = "#0000FF") %>% list(span.emph = .) %>% cli::cli_div(theme = .)
    prefix <- "INFO"; paste0("[", Sys.time(), "] {.emph {prefix}} ==> ", msg) %>% cli::cli_alert_info()
}

#' @description Show a warning message on the screen to indicate the existence of a file path.
#' @param filepath A file path to be shown.
show_exis_war = function(filepath){
    list(color = "#FF0000") %>% list(span.emph = .) %>% cli::cli_div(theme = .)
    msg <- "] {.emph {prefix}} ==> file {.path {filepath}} exists, we directly use it!"
    prefix <- "WARN"; paste0("[", Sys.time(), msg) %>% cli::cli_alert_warning()
}

#' @description Show a common message on the screen to indicate the success of results saving.
#' @param str A string indicating the position of results in the microgeo dataset.
show_stat_msg = function(str){
    list(color = "#33FF00") %>% list(span.emph = .) %>% cli::cli_div(theme = .)
    msg <- "] {.emph {prefix}} ==> results have been saved to: object$"
    prefix <- "SAVE"; paste0("[", Sys.time(), msg, str) %>% cli::cli_alert_success()
}

#' @description Show a common message on the screen to indicate the success of file loading.
#' @param filepath A file path to be shown.
show_load_msg = function(filepath){
    list(color = "#33FF00") %>% list(span.emph = .) %>% cli::cli_div(theme = .)
    msg <- "] {.emph {prefix}} ==> results have been loaded from: {.path {filepath}}"
    prefix <- "LOAD"; paste0("[", Sys.time(), msg) %>% cli::cli_alert_success()
}

#' @description Create a directory according to specified path.
#' @param dirpath A directory path to be created.
#' @param recursive Should elements of the path other than the last be created? Default is `TRUE`.
#' @return A created directory path.
create_dir = function(dirpath, recursive = TRUE){
    if (!dirpath %>% dir.exists) dirpath %>% dir.create(recursive = recursive)
    return(dirpath)
}

#' @description Check the class of a geographic map dataset.
#' @param map A geographic map dataset. It is expected to be a `SpatialPolygonsDataFrame`.
check_mapdata = function(map){
    if (map %>% class != 'SpatialPolygonsDataFrame') 'The <map> must be a `SpatialPolygonsDataFrame`!' %>% stop()
    msg <- "The <map> must be a `SpatialPolygonsDataFrame` returned by `microgeo::read_aliyun_map()`,
    `microgeo::trans_map_fmt()`, `microgeo::grid_map()` or `merge_*_to_map()`!"
    if (!'FMTS' %in% names(map@data) || map@data$FMTS %>% unique != 'microgeo') msg %>% stop()
}

#' @description Check the esential elements of a metadata table.
#' @param met A `data.frame` of sample info., and row names are sample ids while the column names are observed variables.
#' @param lon Column name of longitude in <met>. Default is `longitude`.
#' @param lat Column name of latitude in <met>. Default is `latitude`.
check_metdata = function(met, lon = "longitude", lat = "latitude"){
    if (met %>% class != "data.frame") stop("The <met> must be a `data.frame`!")
    if (!lon %in% colnames(met)) paste0("The `", lon, "` not in <met> table!") %>% stop()
    if (!lat %in% colnames(met)) paste0("The `", lat, "` not in <met> table!") %>% stop()
}

#' @description Check the class of microgeo dataset.
#' @param dataset A microgeo dataset. It is expected to be a `MicrogeoDataset`.
check_dataset = function(dataset){
    msg <- "The <dataset> must be a `MicrogeoDataset`!"
    if (dataset %>% class %>% length != 2) msg %>% stop()
    if (class(dataset)[1] != "list" || class(dataset)[2] != "MicrogeoDataset") msg %>% stop()
}

#' @description Check the class of ggplot2 map layer.
#' @param map.layer A ggplot2 object. It is expected to be a `gg` class.
check_ggplot2_object = function(map.layer){
    msg <- "The <map.layer> must be a `ggplot2` object!"
    if (map.layer %>% class %>% length != 2) msg %>% stop()
    if (class(map.layer)[1] != "gg" || class(map.layer)[2] != "ggplot") msg %>% stop()
}

#' @description Check whether the map object is returned by `microgeo::read_aliyun_map()`.
#' @param map A geographic map with the class of `SpatialPolygonsDataFrame`.
#' @return `TRUE` if the map object is returned by `microgeo::read_aliyun_map()`;
check_aliyun_map = function(map){
    is.aliyun.map <- FALSE
    if ('FMTS' %in% names(map@data) && unique(map@data$FMTS) == 'microgeo' &&
        unique(map@data$TYPE) == 'DataV.GeoAtlas'){
        is.aliyun.map <- TRUE
    }
    return (is.aliyun.map)
}

#' @description Check whether a variable exists in the `SpatialPolygonsDataFrame`.
#' @param map A geographic map object with the class of `SpatialPolygonsDataFrame`.
#' @param var A variable name to be checked.
check_map_var = function(map, var){
    msg <- '` not in <map>! Use `head(<your map object>@data)` to check the avaliable variables!'
    if (!var %in% names(map@data)) paste0('The `', var, msg) %>% stop()
}

#' @description Check whether the ggplot2 theme is valid.
#' @param gg.theme A ggplot2 theme like `theme_bw()`.
check_ggplot2_theme = function(gg.theme) {
    if (!gg.theme %>% ggplot2::is.theme()) "Invalid ggplot2 theme!" %>% stop()
}

#' @description Get the position of legend for ggplot2 visualization.
#' @param legend.position A character or numeric vector indicating the position of legend.
#' @return The position of ggplot2 legend.
get_legend_position = function(legend.position){
    # legend position: c(0.2, 0.3) or c("right", "left", "top", "bottom")
    if (legend.position %>% length > 1 && !legend.position %>% is.numeric)
        legend.position <- legend.position[1]
    return(legend.position)
}

#' @description Add a rectangle for ggplot2 legend.
#' @param p.map A ggplot2 object with the class of `gg`.
add_legend_rect = function(p.map){
    legend.background.obj <- element_rect(fill = NA, size = 0.2, linetype = "solid", colour = "gray30")
    p.map <- p.map + theme(legend.background = legend.background.obj)
    return(p.map)
}

#' @description Check the data used for `microgeo::add_label()`.
#' @param dat A `data.frame` containing coordinates and the variable to be checked.
#' @param lon.var Variable of longitude. Default is `longitude`.
#' @param lat.var Variable of latitude. Default is `latitude`.
#' @param lab.var Variable of labels to be checked.
#' @param lab.var.must.be.num Should the <lab.var> be numeric? Default is `FALSE`.
check_label_data = function(dat, lon.var, lat.var, lab.var, lab.var.must.be.num = FALSE){
    if (dat %>% class != "data.frame") "The <dat> must be a `data.frame`!" %>% stop()
    if (!lon.var %in% colnames(dat)) paste0("The `", lon.var, "` not in <dat>!") %>% stop()
    if (!lat.var %in% colnames(dat)) paste0("The `", lat.var, "` not in <dat>!") %>% stop()
    if (!lab.var %in% colnames(dat)) paste0("The `", lab.var, "` not in <dat>!") %>% stop()
    if (lab.var.must.be.num && !dat[,lab.var] %>% is.numeric) 'The <lab.var> must be numeric!' %>% stop()
}

#' @description Check the `SpatRaster`.
#' @param spat.raster A `SpatRaster` to be checked.
#' @param arg.name Argument name in the error message. Default is `NULL`.
check_spatraster = function(spat.raster, arg.name = NULL){
    arg.name <- ifelse(arg.name %>% is.null, "spat.raster", arg.name)
    if (spat.raster %>% class != "SpatRaster")
        paste0('The <', arg.name,'> must be a `SpatRaster`!') %>% stop()
}

#' @description Check the data to be merged with a map object.
#' @param dat A `data.frame` or `matrix` containing the variables to be merged with a map.
#' @param met A `data.frame` of sample information. Row names should be sample ids and the column names must be variables
#' (e.g., `longitude`, `latitude` and `group`). Longitude and latitude must be included in this `data.frame` as we mainly
#' focus on the spatial patterns of microbial traits.
#' @param map Geographic map with a class of `SpatialPolygonsDataFrame`.
check_merging_data = function(dat, met, map){
    met %>% check_metdata(lon = 'longitude', lat = 'latitude') # 'longitude' and 'latitude' must be the colnames of <met>
    check.length <- unique(rownames(dat) == rownames(met)) %>% length
    check.logics <- !unique(rownames(dat) == rownames(met))
    if (check.length > 1 || check.logics) 'Sample ids in <dat> can not be matched to those in <met>!' %>% stop()
    is.matrix <- ifelse((dat %>% is.matrix || dat %>% is.data.frame) &&
                            isSymmetric(dat %>% as.matrix %>% unname), TRUE, FALSE)
    if (is.matrix){
        check.length <- unique(colnames(dat) == rownames(met)) %>% length()
        check.logics <- !unique(colnames(dat) == rownames(met))
        if (check.length > 1 || check.logics) 'Sample ids in <dat> can not be matched to those in <met>!' %>% stop()
    }
    map %>% check_mapdata()
    if (!'NAME' %in% names(map@data)) "Invalid `SpatialPolygonsDataFrame`!" %>% stop()
    if (map@polygons %>% length != map@data %>% nrow) "Invalid `SpatialPolygonsDataFrame`!" %>% stop()
}

#' @description Sort samples to ensure the order is same as that in metadata table.
#' @param dat A `data.frame` to be applied for sorting.
#' @param met A `data.frame` of metadata serving as a reference.
#' @param is.matrix Is the <dat> a matrix? Default is `FALSE`.
#' @return A `data.frame` (<is.matrix> is `FALSE`) or `matrix` (<is.matrix> is `TRUE`).
sort_samples = function(dat, met, is.matrix = FALSE){

    # if the <dat> is a data.frame with samples as the row names: is.matrix = FALSE
    if (!is.matrix){
        idx <- sapply(met %>% rownames, function(x){ which(rownames(dat) == x) })
        colname <- colnames(dat); dat <- dat[idx,] %>% as.data.frame()
        if (dat %>% ncol == 1){
            dat <- data.frame(row.names = rownames(met), val = dat)
            colnames(dat) <- colname
        }
        if (length(unique(rownames(met) == rownames(dat))) > 1 || !unique(rownames(met) == rownames(dat)))
            stop("Failed to sort sample ids, please check your data!")

    # if the <dat> is a matrix with samples as the row and column names: is.matrix = TRUE
    }else{
        idx <- sapply(rownames(met), function(x){ which(rownames(dat) == x) })
        dat <- dat[idx, idx] %>% as.data.frame()
        if (length(unique(rownames(met) == rownames(dat))) > 1 || !unique(rownames(met) == rownames(dat)))
            stop("Failed to sort samples, please check your data!")
        if (length(unique(rownames(met) == colnames(dat))) > 1 || !unique(rownames(met) == colnames(dat)))
            stop("Failed to sort samples, please check your data!")
        dat %<>% as.matrix()
    }
    return(dat)
}

#' @description Create a `data.frame` for spatial interpolation
#' @param dat A `data.frame` containing variables to be interpolated. Row names should be sample ids and the column names
#' are observed variables.
#' @param met A `data.frame` of sample information. Row names should be sample ids and the column names must be variables
#' (e.g., `longitude`, `latitude` and `group`). Longitude and latitude must be included in this `data.frame` as we mainly
#' focus on the spatial patterns of microbial traits.
#' @param var Which variable would be used for interpolation (in <dat>)?
#' @param trim.dup Whether to randomly remove sampling sites with duplicated coordinates? Default is `FALSE`.
#' @return A `data.frame` with the `longitude`, `latitude` and `target` (<var>) as the column names.
get_interpolation_data = function(dat, met, var, trim.dup = FALSE){
    met %>% check_metdata();
    if (!var %in% colnames(dat)) paste0('Can not find `', var, '` in <dat>!') %>% stop()
    check.length <- unique(rownames(met) == rownames(dat)) %>% length
    check.logics <- !unique(rownames(met) == rownames(dat))
    if (check.length > 1 || check.logics) stop('Sample ids in <dat> can not be matched to those in <met>!')
    use.dat <- data.frame(row.names = rownames(met), longitude = met$longitude,
                          latitude  = met$latitude, target = dat[,var])
    if (use.dat[, c('longitude', 'latitude')] %>% unique() %>% nrow() < nrow(met)){ # there are duplicated coordinates
        if (!trim.dup) # rise an error if trim.dup = FALSE
            stop('Duplicated coordinates were detected! you can use `trim.dup = TRUE` to ignore this error,
                 but it would randomly remove several sampling sites!')
        use.dat.0 <- use.dat[,c('longitude', 'latitude')] %>% unique()
        use.dat   <- use.dat[rownames(use.dat.0),]
        paste0('only use ', nrow(use.dat), ' out of ', nrow(met),' sampling sites for interpolation!') %>% warning()
    }
    return(use.dat)
}

#' @description Format file size
#' @param size_in_bytes File size in bytes
#' @return A formatted file size with correct unit
format_file_size = function(size_in_bytes) {
    units <- c("B", "KB", "MB", "GB", "TB")
    unit_index <- max(0, floor(log(size_in_bytes, base = 1024)))
    size_in_unit <- size_in_bytes / (1024 ^ unit_index)
    formatted_size <- sprintf("%.2f %s", size_in_unit, units[unit_index + 1])
    return(formatted_size)
}

#' @description Download file from public database
#' @param url A remote file URL.
#' @param destfile A local file path.
#' @param show.size Show downloading size? Default is `TRUE`.
download_remote_file = function(url, destfile, show.size = TRUE){
    if (show.size){
        head_response <- httr::HEAD(url)
        if (head_response$status_code == 200 && "content-length" %in% names(httr::headers(head_response))){
            size_in_bytes <- as.numeric(httr::headers(head_response)$`content-length`)
            formatted_size <- format_file_size(size_in_bytes)
            paste0("Trying URL (", formatted_size, "): ", url) %>% show_comm_msg()
        }else{
            paste0("Trying URL: ", url) %>% show_comm_msg()
        }
    }else{
        paste0("Trying URL: ", url) %>% show_comm_msg()
    }
    response <- tryCatch({
        httr::GET(url, httr::progress(),
                  httr::progress(),
                  httr::write_disk(destfile, overwrite = TRUE),
                  httr::timeout(600))
    }, error = function(e){
        paste0("Failed to download file: `", url,
               "`. Please manually download this file through using a browser, and copy it to ", dirname(destfile)) %>% stop()
    })
    if (response$status_code != 200)
        paste0("Failed to download file: `", url,
               "`. Please manually download this file through using a browser, and copy it to ", dirname(destfile)) %>% stop()
}
