# -----------------------------------------------------------------------------------------------------------------------
# Copyright (c) 2023, microgeo/Chaonan Li (licn@mtc.edu.cn).                                                            #
# The microgeo is distributed under the terms of the Modified BSD License.                                              #
# Full license is avaliable in the file LICENSE, distributed with this package.                                         #
# -----------------------------------------------------------------------------------------------------------------------

#' @description Merge a `SpatRaster` with a microgeo dataset.
#' @param dataset A microgeo dataset with the class of `MicrogeoDataset`.
#' @param spat.rast A `SpatRaster` to be merged with microgeo dataset.
#' @param type Which type of `SpatRaster` would be merged? Default is `his`. Also could be `cla`.
#' @return A `MicrogeoDataset` class containing all avaliable data for microgeo R package.
merge_spat_raster = function(dataset, spat.rast, type = "his"){
    dataset %>% check_dataset()
    spat.rast %>% check_spatraster(arg.name = "spat.rast")
    if (dataset$spa$rast[[type]] %>% is.null){
        dataset$spa$rast[[type]] <- spat.rast
    }else{
        vars.in.new.rast <- spat.rast %>% names()
        vars.in.old.rast <- dataset$spa$rast[[type]] %>% names()
        intersect.vars   <- intersect(vars.in.new.rast, vars.in.old.rast)
        subset.rast.vars <- vars.in.old.rast[which(!vars.in.old.rast %in% intersect.vars)]
        if (subset.rast.vars %>% length > 0){
            dataset$spa$rast[[type]] <- dataset$spa$rast[[type]][[subset.rast.vars]]
            re.rast <- terra::resample(x = spat.rast, y = dataset$spa$rast[[type]], method = 'bilinear')
            dataset$spa$rast[[type]] <- c(dataset$spa$rast[[type]], re.rast)
        }else{
            dataset$spa$rast[[type]] <- spat.rast
        }
    }
    msg <- ifelse(spat.rast %>% names %>% length > 1,
                  names(spat.rast) %>% length() %>% paste0('spa$rast$', type, '(', ., " variables)"),
                  names(spat.rast) %>% paste0('spa$rast$', type, "$", .))
    msg %>% show_stat_msg()
    return (dataset)
}

#' @description Download aridity index from global aridity index database.
#' @param outpath Output directory path.
#' @return A `SpatRaster` of aridity index.
download_aridity_index = function(outpath){
    ai.tiff <- outpath %>% file.path(., "ai_et0.tif")
    if (ai.tiff %>% file.exists) terra::rast(ai.tiff) %>% return ()
    downloadURL <- "https://figshare.com/ndownloader/files/14118800"
    downloadPAT <- outpath %>% file.path(., "global-ai_et0.zip")
    if (!downloadPAT %>% file.exists) download.file(url = downloadURL, destfile = downloadPAT)
    unziped.dirpath <- outpath %>% file.path(., "global-ai_et0")
    unzip(zipfile = downloadPAT, overwrite = T, exdir = unziped.dirpath)
    copy.status <- file.copy(from = unziped.dirpath %>% file.path(., 'ai_et0/ai_et0.tif'), to = ai.tiff, overwrite = T)
    unlink(unziped.dirpath, recursive = T)
    terra::rast(ai.tiff) %>% return()
}

#' @description Download elevation from WorldClim database version 2.1.
#' @param res Which resolution would be used (minutes of a degree)?
#' @param outpath Output directory path.
#' @return A `SpatRaster` of elevation.
download_elev = function(res, outpath){
    res.part <- ifelse(res == 0.5, '30s', paste0(res, 'm'))
    elev.path <- paste0("wc2.1_", res.part) %>% file.path(outpath, .) %>% create_dir(dirpath = .)
    elev.tiff <- paste0("wc2.1_", res.part, '_elev.tif') %>% file.path(elev.path, .)
    if (elev.tiff %>% file.exists) terra::rast(elev.tiff) %>% return ()
    zipfile <- paste0("wc2.1_", res.part, '_elev.zip')
    downloadURL <- file.path("https://biogeo.ucdavis.edu/data/worldclim/v2.1/base", zipfile)
    downloadPAT <- file.path(elev.path, zipfile)
    if (!downloadPAT %>% file.exists) download.file(url = downloadURL, destfile = downloadPAT)
    unzip(zipfile = downloadPAT, overwrite = T, exdir = elev.path)
    terra::rast(elev.tiff) %>% return()
}

#' @description Download bioclimatic variables from WorldClim database version 2.1.
#' @param res Which resolution would be used (minutes of a degree)?
#' @param outpath Output directory path.
#' @return A `SpatRaster` of historical bioclimatic variables.
download_his_bioc = function(res, outpath){
    res.part <- ifelse(res == 0.5, '30s', paste0(res, 'm'))
    bioc.path <- paste0("wc2.1_", res.part) %>% file.path(outpath, .) %>% create_dir(dirpath = .)
    zipfile <- paste0("wc2.1_", res.part, '_bio.zip')
    zippath <- file.path(bioc.path, zipfile)
    downloadURL <- file.path("https://biogeo.ucdavis.edu/data/worldclim/v2.1/base", zipfile)
    if (!zippath %>% file.exists) download.file(url = downloadURL, destfile = zippath)
    unzip(zipfile = zippath, overwrite = T, exdir = bioc.path)
    file.list <- paste0("wc2.1_", res.part, '_bio_') %>% paste0(., seq(19), '.tif') %>% file.path(bioc.path, .)
    terra::rast(file.list) %>% return()
}

#' @description Download future bioclimatic variables from WorldClim database version 2.1.
#' @param arg A string of additional arguments (`gcm|sec|yea`).
#' @param res Which resolution would be used (minutes of a degree)?
#' @param outpath Output directory path.
#' @return A `SpatRaster` of future bioclimatic variables.
download_fut_bioc = function(arg, res, outpath){
    dat.part  <- strsplit(arg, split = "|", fixed = T) %>% unlist()
    res.part  <- ifelse(res == 0.5, '30s', paste0(res, 'm'))
    bioc.path <- paste0("wc2.1_", res.part) %>% file.path(outpath, .) %>% create_dir(dirpath = .)
    tiffile   <- paste0("wc2.1_", res.part, '_bioc_', dat.part[1], '_', dat.part[2], '_', dat.part[3], '.tif')
    tifpath   <- file.path(bioc.path, tiffile)
    downloadURL <- paste0("https://geodata.ucdavis.edu/cmip6/",
                          res.part, '/', dat.part[1], '/', dat.part[2], '/', tiffile)
    if (!tifpath %>% file.exists) download.file(url = downloadURL, destfile = tifpath)
    terra::rast(tifpath) %>% return()
}

#' @description Create a `data.frame` for all avaliable MODIS product in microgeo R package.
#' @param dataset A microgeo dataset with the class of `MicrogeoDataset`.
#' @param measures Which measure would be downloaded?
#' @param prod.res The resolution of MODIS product. Only for `NDVI` and `EVI`.
#' @param prod.typ The type of MODIS product.
#' @param date.ran Which date range would be applied for MODIS product downloading?
#' @return An R `list` of MODIS data information.
get_modis_products = function(dataset, measures, prod.res, prod.typ, date.ran){
    show_comm_msg("preparing MODIS product list for searching..."); check_dataset(dataset)
    prod.aval <- system.file("modis", "modis.products.Rds", package = "microgeo") %>% base::readRDS()
    prod.list <- lapply(measures, function(measure){
        prod.res <- ifelse(measure %in% c("NDVI", "EVI"), prod.res,
                           prod.aval[which(prod.aval$measure.name == measure),]$prod.resolution %>% unique())
        prod.aval[which(prod.aval$measure.name == measure &
                        prod.aval$prod.resolution == prod.res &
                        prod.aval$prod.type == prod.typ),]
    }) %>% do.call("rbind", .) # a data.frame
    modis.list <- lapply(date.ran, function(date.r){
        bbox <- terra::ext(dataset$map) %>% as.vector() %>% matrix(ncol=2) %>% t() %>%
                as.vector() %>% paste(collapse = ",")
        se.date <- strsplit(date.r, '|', fixed = T) %>% unlist()
        unique(prod.list$product.name) %>% paste(., date.r, bbox, sep = "|")
    }) %>% unlist() # a vector
    list(modis.list = modis.list, prod.list = prod.list) %>% return()
}

#' @description Search avaliable MODIS products from `https://cmr.earthdata.nasa.gov`.
#' @param modis.prod An R list containing the searching items returned by \code{get_modis_products()}.
#' @return A data.frame of MODIS data information.
que_modis_products = function(modis.prod){
    show_comm_msg("searching avaliable MODIS products...")
    modis.list <- modis.prod$modis.list; prod.list  <- modis.prod$prod.list
    search.res <- lapply(X = modis.list %>% length %>% seq, function(x){
        d.vector <- strsplit(modis.list[x], '|', fixed = T) %>% unlist()
        date.range <- paste0(as.Date(d.vector[2]), "T00:00:00Z", ",", as.Date(d.vector[3]), "T00:00:00Z")
        kwargs.s <- list(short_name = d.vector[1], temporal = date.range, downloadable = "true",
                         bounding_box = d.vector[4])
        url <- "https://cmr.earthdata.nasa.gov/search/granules"; page.num <- 1; results <- NULL
        paste0('current product (', x, '/', length(modis.list), '): ',
               d.vector[1], " (", prod.list[which(prod.list$product.name == d.vector[1]),]$measure.name %>%
               unique() %>% paste(., collapse = "|"), "--> ", d.vector[2],
               " to ", d.vector[3], ")") %>% show_comm_msg()
        while (TRUE){
            rsts  <- httr::GET(url = url, httr::add_headers(Accept = "text/csv"),
                               query = c(kwargs.s, page_num = page.num)) %>% httr::stop_for_status()
            page  <- utils::read.csv(text = httr::content(rsts, as="text"),
                                     check.names = FALSE, stringsAsFactors = FALSE)
            if (page %>% nrow == 0) break(); page.num <- page.num + 1; results <- rbind(results, page)
        }
        if (results %>% is.null)
            paste0("Faild in searching: ", d.vector[1],
                   " (", prod.list[which(prod.list$product.name == d.vector[1]),]$measure.name %>%
                   unique() %>% paste(., collapse = "|"), "--> ", d.vector[2],
                   " to ", d.vector[3], "). In the specified time period, there may be no data.") %>% stop()
        p2v.df <- lapply(results$`Online Access URLs`, function(url){
            info.parts <- basename(url) %>% strsplit(., split = '.', fixed = T) %>% unlist()
            data.frame(prod.name = info.parts[1], version.name = info.parts[4], count = 1, url = url)
        }) %>% do.call("rbind", .)
        p2v.df.agg <- aggregate(p2v.df$count, by = list(p2v.df$prod.name, p2v.df$version.name), FUN = sum)
        colnames(p2v.df.agg) <- c("product.name", "product.version", "hdf.file.count")
        use.version <- ifelse("061" %in% p2v.df.agg$product.version, "061", "006")
        p2v.df <- p2v.df[which(p2v.df$version.name == use.version),]
        final.rst <- results[which(results$`Online Access URLs` %in% p2v.df$url),]
    }) %>% do.call("rbind", .)
    file.size <- round(sum(search.res$Size)/1024, 2) # total size of files
    paste0("find ", nrow(search.res), " files with ",
           file.size, " GB in total...") %>% show_comm_msg()
    return(search.res)
}

#' @description Download avaliable MODIS products from `https://cmr.earthdata.nasa.gov`.
#' @param username Username of your EOSDIS account. You can sign up \href{https://urs.earthdata.nasa.gov/users/new}{here}.
#' @param password Password of your EOSDIS account. You can sign up \href{https://urs.earthdata.nasa.gov/users/new}{here}.
#' @param search.res A `data.frame` returned by \code{que_modis_products()}.
#' @param outpath Output directory path.
#' @return A directory path of downloaded MODIS products (HDF files).
dow_modis_products = function(username, password, search.res, outpath){
    show_comm_msg("downloading all avaliable MODIS products[skip if the file exists]...")
    hdfpath <- file.path(outpath, "modis_products/hdf") %>% create_dir()
    res <- lapply(X = 1:nrow(search.res), function(x){
        info.part <- basename(search.res[x, 5]) %>% strsplit(., split = '.', fixed = TRUE) %>% unlist()
        tmp.path <- create_dir(dirpath = file.path(hdfpath, info.part[1]), recursive = T)
        destfile <- file.path(tmp.path, basename(search.res[x, 5]))
        if (!file.exists(destfile)){
            paste0("current file (", x, "/", search.res %>% nrow, "): ",
                   basename(search.res[x, 5]), " (", round(search.res[x, 9], 2), " MB)") %>% show_comm_msg()
            f <- httr::GET(search.res[x, 5], httr::authenticate(username, password),
                           httr::progress(), httr::write_disk(destfile, overwrite = TRUE))
        }
    })
    return(hdfpath)
}

#' @description  Create a PTV (Product, Time, Version) dataframe for remote-sensing image merging.
#' @param prod.list A `data.frame` created by \code{get_modis_products()}.
#' @param hdfs.path Directory path of downloaded MODIS HDF files. It is the value returned by \code{dow_modis_products()}.
#' @return A `data.frame` of PTV (Product, Time, Version).
ptv_modis_products = function(prod.list, hdfs.path){
    show_comm_msg("preparing the PTVs (Product, Time, Version) for merging remote-sensing images...")
    product.names <- prod.list$product.name %>% unique()
    ptv <- lapply(product.names, function(product.name){
        prod.path <- file.path(hdfs.path, product.name)
        res1 <- lapply(list.files(prod.path), function(hdf){
            name.list <- strsplit(hdf, split = ".", fixed = T) %>% unlist()
            pattern <- paste0("^", name.list[1], "\\.", name.list[2], ".*\\.", name.list[4], ".*\\hdf$")
            if (!name.list[1] == product.name) stop("Error in product name!")
            measures.v <- prod.list[which(prod.list$product.name == name.list[1]),]$measure.name %>% unique()
            res2 <- lapply(measures.v, function(measure){
                sds.name <- prod.list[which(prod.list$measure.name == measure),]$sds.name %>% unique()
                unit.name <- prod.list[which(prod.list$measure.name == measure),]$unit.name %>% unique()
                scale.factor <- prod.list[which(prod.list$measure.name == measure),]$scale.factor %>% unique()
                res3 <- data.frame(measure = measure, unit.name = unit.name, scale.factor = scale.factor,
                        product = name.list[1], sds.name = sds.name, time = name.list[2], version = name.list[4],
                        pattern = pattern)
            }) %>% do.call("rbind", .)
            res2
        }) %>% do.call("rbind", .) %>% na.omit()
        res1
    }) %>% do.call("rbind", .) %>% unique()
    return(ptv)
}

#' @description Merge MODIS HDF files downloaded from `https://cmr.earthdata.nasa.gov`.
#' @param prod.list A data.frame created by \code{get_modis_products()}.
#' @param ptvdata A data.frame created by \code{ptv_modis_products()}.
#' @param hdfpath Directory path of downloaded MODIS HDF files. It is the value returned by \code{dow_modis_products()}.
#' @param threads How many threads would be used for merging?
#' @param outpath Output directory path.
#' @return A `data.frame` of merged remote-sensing image information.
meg_modis_products = function(prod.list, ptvdata, hdfpath, threads, outpath){
    show_comm_msg("converting hdf files to tif files...")
    product.names <- prod.list$product.name %>% unique(); modis.dat.all.rst <- NULL
    for (i in seq(length(product.names))){
        product.name <- product.names[i];
        hdfpath.tmp <- file.path(hdfpath, product.name); hdf.files <- list.files(hdfpath.tmp)
        ptvdata.tmp <- ptvdata[which(ptvdata$product == product.name),]
        tifpath <- create_dir(file.path(outpath, "modis_products/tif", product.name), recursive = T)
        threads.use <- ifelse(threads > nrow(ptvdata.tmp), nrow(ptvdata.tmp), threads)
        paste0('current product (', i, '/', length(product.names), '): ',
               product.name, " (convert ", length(hdf.files), " hdf files into ", nrow(ptvdata.tmp),
               " tif files using ", threads.use, " threads)") %>% show_comm_msg()
        cl <- parallel::makeCluster(threads.use)
        parts <- split(x = 1:nrow(ptvdata.tmp), f = 1:threads.use) %>% suppressWarnings()
        parallel::clusterExport(cl = cl, varlist = c("ptvdata.tmp", "parts", "tifpath", "hdfpath.tmp", "%>%"),
                                envir = environment())
        parallelX <- parallel::parLapply(cl = cl, X = 1:threads.use, fun = function(x){
            dat.parts <- ptvdata.tmp[parts[[x]], ]
            res <- lapply(X = 1:nrow(dat.parts), function(y){
                dat.part <- dat.parts[y, ]
                files <- list.files(hdfpath.tmp, pattern = dat.part$pattern)
                files <- file.path(hdfpath.tmp, files)
                outfile.name <- file.path(tifpath, paste(dat.part$product, dat.part$time,
                                                         dat.part$version, gsub(" ", "_", dat.part$sds.name), sep = '.'))
                outfile.name <- paste0(outfile.name, ".tif")
                if (!file.exists(outfile.name)){
                    merg.list <- lapply(files, function(file){
                        sds.object <- MODIS::getSds(file); sds.names <- sds.object$SDSnames
                        hdf.layer <- terra::rast(sds.object$SDS4gdal[which(sds.names == dat.part$sds.name)])
                        if (dat.part$scale.factor != "NA") hdf.layer <- hdf.layer * as.numeric(dat.part$scale.factor)
                        hdf.layer
                    })
                    layer.text <- paste(paste0(rep("merg.list[[", length(merg.list)),
                                               seq(length(merg.list)), rep("]]", length(merg.list))), collapse = ", ")
                    layer.cmds <- paste0("terra::mosaic(", layer.text, ", fun = 'mean')")
                    merged.layer <- parse(text = layer.cmds) %>% eval()
                    terra::writeRaster(x = merged.layer, filename = outfile.name, overwrite = TRUE)
                }
                res <- data.frame(measure.name = dat.part$measure, product.name = dat.part$product,
                                  time.name = dat.part$time,
                                  band.name = dat.part$sds.name, veri.name = dat.part$version,
                                  unit.name = dat.part$unit.name,
                                  tif.path = outfile.name)
            }) %>% do.call("rbind", .) %>% as.data.frame()
        })
        parallel::stopCluster(cl)
        modis.dat <- do.call("rbind", parallelX)
        modis.dat.all.rst <- rbind(modis.dat.all.rst, modis.dat)
    }
    return(modis.dat.all.rst)
}

#' @description Copy MODIS remote-sensing images into one dir. It only works when the data in images are classification
#' @param modis.data A `data.frame` created by \code{meg_modis_products()}.
#' @param date.ran Which date range would be applied for MODIS product downloading?
#' @param outpath Output directory path.
#' @return A `data.frame` of remote-sensing image information.
col_modis_products = function(modis.data, date.ran, outpath){
    show_comm_msg("collecting all merged image files...")
    system_info <- Sys.info()
    is_windows <- tolower(system_info["sysname"]) == "windows"
    savepath <- file.path(outpath, "modis_products", "avg") %>% create_dir()
    measure2veris <- paste(modis.data$measure.name, modis.data$veri.name, sep = '_') %>% unique()
    res <- lapply(X = seq(length(measure2veris)), function(x){
        measure2veri <- measure2veris[x]
        paste0('current measure (', x, '/', length(measure2veris), '): ', measure2veri) %>% show_comm_msg()
        str.part <- strsplit(measure2veri, "_", fixed = T) %>% unlist()
        savefile <- file.path(savepath, paste0(measure2veri, "_", paste(date.ran, collapse = "_"), ".tif"))
        tmp.dfs  <- modis.data[which(modis.data$measure.name == paste(str.part[1:(length(str.part) - 1)],
            collapse = "_") & modis.data$veri.name == str.part[length(str.part)]),]
        if (is_windows){
            copy.res <- file.copy(from = normalizePath(tmp.dfs$tif.path), to = savefile, overwrite = T)
        }else{
            copy.res <- file.copy(from = tmp.dfs$tif.path, to = savefile, overwrite = T)
        }
        if (!copy.res) stop("Failed to copy file!")
        unit  <-  tmp.dfs$unit.name %>% unique(); unit <- ifelse(unit == 'no unit', "", unit)
        title <- paste(tmp.dfs$product.name, gsub(" ", "_", tmp.dfs$band.name),
                       tmp.dfs$veri.name, sep = ".") %>% unique()
        title <- ifelse(unit != "", paste0(title, " (", unit, ")"), title)
        res0  <- data.frame(measure = paste(str.part[1:(length(str.part) - 1)], collapse = "_"),
                            title = title, dbpath = savefile)
    }) %>% do.call("rbind", .)
    return (res)
}

#' @description Calculate the mean values of MODIS remote-sensing images based on specified date range. It works when the
#' data in images are numeric！
#' @param modis.data A `data.frame` created by \code{meg_modis_products()}.
#' @param date.ran Which date range would be applied for MODIS product downloading?
#' @param threads How many threads would be used?
#' @param outpath Output directory path.
#' @return A `data.frame` of remote-sensing image information.
avg_modis_products = function(modis.data, date.ran, threads, outpath){
    show_comm_msg("calculating average values based on date range...")
    savepath <- file.path(outpath, "modis_products/avg") %>% create_dir()
    measure2veris <- paste(modis.data$measure.name, modis.data$veri.name, sep = '_') %>% unique()
    res <- lapply(X = seq(length(measure2veris)), function(x){
        measure2veri <- measure2veris[x]
        paste0('current measure (', x, '/', length(measure2veris), '): ',
               measure2veri, "(", threads, " threads)") %>% show_comm_msg()
        str.part <- strsplit(measure2veri, "_", fixed = T) %>% unlist()
        savefile <- file.path(savepath, paste0(measure2veri, "_", paste(date.ran, collapse = "_"), ".tif"))
        tmp.dfs  <- modis.data[which(modis.data$measure.name == paste(str.part[1:(length(str.part) - 1)],
            collapse = "_") & modis.data$veri.name == str.part[length(str.part)]),]
        if (!file.exists(savefile)) {
            files <- tmp.dfs$tif.path
            stk <- raster::stack(); for (i in seq(length(files))) stk <- raster::addLayer(stk, files[i])
            raster::beginCluster(n = threads)
            stk.avg <- raster::clusterR(stk, fun = raster::calc, args = list(fun = mean, na.rm = TRUE))
            raster::endCluster()
            raster::writeRaster(stk.avg, savefile, overwrite = TRUE)
        }
        unit  <-  tmp.dfs$unit.name %>% unique(); unit <- ifelse(unit == 'no unit', "", unit)
        title <- paste(tmp.dfs$product.name, gsub(" ", "_", tmp.dfs$band.name),
                       tmp.dfs$veri.name, sep = ".") %>% unique()
        title <- ifelse(unit != "", paste0(title, " (", unit, ")"), title)
        res0  <- data.frame(measure = paste(str.part[1:(length(str.part) - 1)], collapse = "_"),
                            title = title, dbpath = savefile)
    }) %>% do.call("rbind", .)
    return (res)
}
