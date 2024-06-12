# -----------------------------------------------------------------------------------------------------------------------
# Copyright (c) 2023, microgeo/Chaonan Li (licn@mtc.edu.cn).                                                            #
# The microgeo is distributed under the terms of the GPL-3 License.                                                     #
# Full license is avaliable in the file LICENSE, distributed with this package.                                         #
# -----------------------------------------------------------------------------------------------------------------------

#' @title Create a machine learning (ML) model
#' @author Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
#' @description This function is implemented to create a machine learning model through using \code{caret::train()}.
#' @param y.data A `data.frame` containing variables to be used for ML modeling. For example, a `data.frame` of microbial
#' alpha-diversity indices (e.g., Shannon-Wiener diversity index).
#' @param x.data A `data.frame` containing the prediction variables, e.g., a `data.frame` of 19 bioclimatic variables.
#' @param var Variable (in <y.data>) to be used for ML modeling. For example, the `Shannon-Wiener` index if applicable.
#' @param method ML method. Please see `http://topepo.github.io/caret/train-models-by-tag.html` for more options. Default
#' is `rf`, which means a `Random Forest` model.
#' @param type Type of ML model. Select one from `regression` and `classification`. Default is `regression`.
#' @param testing.ratio How many samples would be applied for testing? Default is `0.2`.
#' @param threads Threads used for model training (logical CPUs). Default is `60`.
#' @param remove.outlier Whether to remove those outliers before creating a machine learning model? If it is `TRUE`, only
#' those values (<y.data>) ranging from `lower quartile - 1.5 × interquartile` to `upper quartile + 1.5 × interquartile`
#' would be applied for model training. It only works when `type = 'regression'` Default is `TRUE`.
#' @param ... Parameters parsed by \code{caret::train()}.
#' @return An R `list` with the following components:
#' \describe{
#'   \item{\code{object$model}}{ML model returned by \code{caret::train()}.}
#'   \item{\code{object$model.method}}{ML method.}
#'   \item{\code{object$model.type}}{ML type (`regression` or `classification`).}
#'   \item{\code{object$tests.dat}}{A `data.frame` of testing data for ML model.}
#'   \item{\code{object$train.dat}}{A `data.frame` of training data for ML model.}
#'   \item{\code{object$var.name}}{Name of variable used for ML modeling.}
#' }
#' @seealso
#' \code{\link[caret:train]{caret::train()}}
#' \code{\link[microgeo:read_aliyun_map]{microgeo::read_aliyun_map()}}
#' \code{\link[microgeo:create_dataset]{microgeo::create_dataset()}}
#' \code{\link[microgeo:show_dataset]{microgeo::show_dataset()}}
#' \code{\link[microgeo:get_ai]{microgeo::get_ai()}}
#' \code{\link[microgeo:get_his_bioc]{microgeo::get_his_bioc()}}
#' \code{\link[microgeo:extract_data_from_spatraster]{microgeo::extract_data_from_spatraster()}}
#' \code{\link[microgeo:rarefy_count_table]{microgeo::rarefy_count_table()}}
#' \code{\link[microgeo:tidy_dataset]{microgeo::tidy_dataset()}}
#' \code{\link[microgeo:calc_rel_abund]{microgeo::calc_rel_abund()}}
#' \code{\link[microgeo:calc_markers]{microgeo::calc_markers()}}
#' @examples
#' # Create a microgeo dataset
#' data(qtp)
#' map <- read_aliyun_map(adcode = c(540000, 630000, 510000))
#' dataset.dts <- create_dataset(mat = qtp$asv, ant = qtp$tax, met = qtp$met, map = map,
#'                               phy = qtp$tre, env = qtp$env, lon = "longitude", lat = "latitude")
#' dataset.dts %>% show_dataset()
#'
#' # Download aridity index and 19 historical bioclimatic variables
#' dataset.dts %<>% get_ai(out.dir = "test/microgeo_data")
#' dataset.dts %<>% get_his_bioc(res = 2.5, out.dir = "test/microgeo_data")
#' dataset.dts %<>% extract_data_from_spatraster(type = 'his')
#' dataset.dts %>% show_dataset()
#'
#' # Tidy up microgeo dataset
#' dataset.dts %<>% rarefy_count_table()
#' dataset.dts %<>% tidy_dataset()
#' dataset.dts %>% show_dataset()
#'
#' # Calculate the relative abundance and the ecological markers at `Family` level
#' dataset.dts %<>% calc_rel_abund()
#' dataset.dts %<>% calc_markers(use.dat = 'spa', use.var = 'AI', annotation.level = 'Family', r.thres = 0.3)
#' head(dataset.dts$abd$mar$correlation)
#' dataset.dts %>% show_dataset()
#'
#' # Create a regression random forest model for soil pH
#' rf.rst.reg <- create_ml_model(y.data = dataset.dts$env,
#'                               x.data = dataset.dts$spa$tabs[,paste0("Bio", seq(19))], var = 'pH', method = 'rf')
#' print(rf.rst.reg$model)
#'
#' # Create a binary classification random forest model for the family of f__A21b
#' family.bins <- data.frame(row.names = rownames(dataset.dts$abd$mar$abundance),
#'                           f__A21b = dataset.dts$abd$mar$abundance$f__A21b)
#' family.bins$f__A21b <- ifelse(family.bins$f__A21b > 0, "presence", "absence") # two classifications
#' family.bins$f__A21b <- as.factor(family.bins$f__A21b)
#' rf.rst.cla.bin <- create_ml_model(y.data = family.bins,
#'                                   x.data = dataset.dts$spa$tabs[,paste0("Bio", seq(19))],
#'                                   var = 'f__A21b', method = 'rf', type = 'classification')
#' print(rf.rst.cla.bin$model)
#'
#' # Create a mutiple class classification random forest model for the family of f__A21b
#' family.mutiple <- data.frame(row.names = rownames(dataset.dts$abd$mar$abundance),
#'                              f__A21b = dataset.dts$abd$mar$abundance$f__A21b)
#' family.mutiple$f__A21b <- cut(family.mutiple$f__A21b, breaks = c(-Inf, 0.05, 0.2, 0.9, Inf),
#'                               labels = c("H", "A", "S", "Y")) # four classifications
#' family.mutiple$f__A21b <- as.factor(family.mutiple$f__A21b)
#' rf.rst.cla.mutiple <- create_ml_model(y.data = family.mutiple,
#'                                       x.data = dataset.dts$spa$tabs[,paste0("Bio", seq(19))],
#'                                       var = 'f__A21b', method = 'rf', type = 'classification')
#' print(rf.rst.cla.mutiple$model)
#' @export
create_ml_model = function(y.data, x.data, var, method = 'rf', type = c('regression', 'classification'),
                           testing.ratio = 0.2, threads = 60, remove.outlier = TRUE, ...){

    # check the data and arguments for ML modeling
    if (class(y.data) != "data.frame") stop("The <y.data> must be a dataframe!")
    if (class(x.data) != "data.frame") stop("The <x.data> must be a dataframe!")
    if (!var %in% colnames(y.data)) paste0('Can not find `', var, '` in <y.data>, please check it!') %>% stop()
    type <- ifelse(length(type) == 1, type, type[1])
    if (type == 'regression' && !is.numeric(y.data[,var])) {
        paste0('The `', var, '` in <y.data> should be numeric if the <type> is `regression`!') %>% stop()
    }
    if (type == 'classification' && !is.factor(y.data[,var])) {
        paste0('The `', var, '` in <y.data> should be factor if the <type> is `classification`!') %>% stop()
    }
    pred.dat.chk <- lapply(X = seq(ncol(x.data)), function(x) is.numeric(x.data[,x])) %>% unlist()
    if (FALSE %in% pred.dat.chk) stop("All data in <x.data> must be numeric!")
    check.length <- length(unique(rownames(y.data) == rownames(x.data)))
    check.logics <- !unique(rownames(y.data) == rownames(x.data))
    if (check.length > 1 | check.logics) stop('Sample ids in <y.data> can not be matched to those in <x.data>!')
    type <- ifelse(length(type) > 1, type[1], type)
    if (!type %in% c('regression', 'classification')) stop('The <type> must be one of `regression` and `classification`!')
    dw.raw <- data.frame(target = y.data[,var], x.data)

    # remove outliers if <type> is `regression` and <remove.outlier> is `TRUE`.
    if (type == 'regression' & remove.outlier){
        dw <- lapply(X = colnames(dw.raw), function(variable){
            dat.tmps <- data.frame(row.names = rownames(dw.raw),
                                   vals = dw.raw[, which(colnames(dw.raw) == variable)])
            Q <- quantile(dat.tmps$vals, probs = c(.25, .75), na.rm = FALSE)
            iqr <- IQR(dat.tmps$vals)
            up <- Q[2]  + 1.5*iqr %>% as.vector() %>% as.numeric() # upper range
            low <- Q[1] - 1.5*iqr %>% as.vector() %>% as.numeric() # lower range
            dat.tmps$vals[which(dat.tmps$vals < low | dat.tmps$vals > up)] <- NA
            colnames(dat.tmps) <- variable; dat.tmps
        }) %>% do.call("cbind", .) %>% na.omit()
    }else{ dw <- dw.raw }

    # create a machine learning model
    set.seed(1234); i <- base::sample(nrow(dw), testing.ratio * nrow(dw))
    tests.dat <- as.data.frame(dw[i,]); train.dat <- as.data.frame(dw[-i,])
    if (type == 'classification'){ # requires at least two classifications
        if (length(unique(train.dat$target)) < 2)
            stop('At least two classifications are required for trainning dataset!')
        if (length(unique(tests.dat$target)) < 2)
            stop('At least two classifications are required for testing dataset!')
    }
    cl <- parallel::makePSOCKcluster(threads); doParallel::registerDoParallel(cl); set.seed(1234)
    if (type == 'regression') { # regression model
        model <- caret::train(target ~ ., data = train.dat, method = method, ...)
    }else{ # classification model
        model <- caret::train(x = train.dat[, 2:ncol(train.dat)], y = train.dat$target, method = method, ...)
    }
    parallel::stopCluster(cl)
    dat <- list(model = model, model.method = method, model.type = type,
                tests.dat = tests.dat, train.dat = train.dat, var.name = var)
    return(dat)
}

#' @title Evaluate machine learning model created by \code{microgeo::create_ml_model()}
#' @author Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
#' @description This function is desinged to evaluate a ML model created by \code{microgeo::create_ml_model()}.
#' @param model.dat A machine learning model object returned by \code{microgeo::create_ml_model()}.
#' @param only.show.class Only show observed classifications in ROC curves (delete micro-/macro- average)? Only works for
#' classification model. Default is `TRUE`.
#' @return A `ggplot2` object.
#' @seealso
#' \code{\link[psych:corr.test]{psych::corr.test()}}
#' \code{\link[multiROC:multi_roc]{multiROC::multi_roc()}}
#' \code{\link[multiROC:plot_roc_data]{multiROC::plot_roc_data()}}
#' \code{\link[microgeo:create_ml_model]{microgeo::create_ml_model()}}
#' \code{\link[microgeo:read_aliyun_map]{microgeo::read_aliyun_map()}}
#' \code{\link[microgeo:create_dataset]{microgeo::create_dataset()}}
#' \code{\link[microgeo:show_dataset]{microgeo::show_dataset()}}
#' \code{\link[microgeo:get_ai]{microgeo::get_ai()}}
#' \code{\link[microgeo:get_his_bioc]{microgeo::get_his_bioc()}}
#' \code{\link[microgeo:extract_data_from_spatraster]{microgeo::extract_data_from_spatraster()}}
#' \code{\link[microgeo:rarefy_count_table]{microgeo::rarefy_count_table()}}
#' \code{\link[microgeo:tidy_dataset]{microgeo::tidy_dataset()}}
#' \code{\link[microgeo:calc_rel_abund]{microgeo::calc_rel_abund()}}
#' \code{\link[microgeo:calc_markers]{microgeo::calc_markers()}}
#' @examples
#' # Create a microgeo dataset
#' data(qtp)
#' map <- read_aliyun_map(adcode = c(540000, 630000, 510000))
#' dataset.dts <- create_dataset(mat = qtp$asv, ant = qtp$tax, met = qtp$met, map = map,
#'                               phy = qtp$tre, env = qtp$env, lon = "longitude", lat = "latitude")
#' dataset.dts %>% show_dataset()
#'
#' # Download aridity index and 19 historical bioclimatic variables
#' dataset.dts %<>% get_ai(out.dir = "test/microgeo_data")
#' dataset.dts %<>% get_his_bioc(res = 2.5, out.dir = "test/microgeo_data")
#' dataset.dts %<>% extract_data_from_spatraster(type = 'his')
#'
#' # Tidy up microgeo dataset
#' dataset.dts %<>% rarefy_count_table()
#' dataset.dts %<>% tidy_dataset()
#' dataset.dts %>% show_dataset()
#'
#' # Calculate the relative abundance and the ecologycal markers at `Family` level
#' dataset.dts %<>% calc_rel_abund()
#' dataset.dts %<>% calc_markers(use.dat = 'spa', use.var = 'AI', annotation.level = 'Family', r.thres = 0.3)
#' head(dataset.dts$abd$mar$correlation)
#'
#' # Create a regression random forest model for soil pH
#' rf.rst.reg <- create_ml_model(y.data = dataset.dts$env,
#'                               x.data = dataset.dts$spa$tabs[,paste0("Bio", seq(19))], var = 'pH', method = 'rf')
#' rf.rst.reg %>% evaluate_ml_model()
#'
#' # Create a binary classification random forest model for the family of f__A21b
#' family.bins <- data.frame(row.names = rownames(dataset.dts$abd$mar$abundance),
#'                           f__A21b = dataset.dts$abd$mar$abundance$f__A21b)
#' family.bins$f__A21b <- ifelse(family.bins$f__A21b > 0, "presence", "absence") # two classifications
#' family.bins$f__A21b <- as.factor(family.bins$f__A21b)
#' rf.rst.cla.bin <- create_ml_model(y.data = family.bins,
#'                                   x.data = dataset.dts$spa$tabs[,paste0("Bio", seq(19))],
#'                                   var = 'f__A21b', method = 'rf', type = 'classification')
#' rf.rst.cla.bin %>% evaluate_ml_model()
#'
#' # Create a mutiple class classification random forest model for the family of f__A21b
#' family.mutiple <- data.frame(row.names = rownames(dataset.dts$abd$mar$abundance),
#'                              f__A21b = dataset.dts$abd$mar$abundance$f__A21b)
#' family.mutiple$f__A21b <- cut(family.mutiple$f__A21b, breaks = c(-Inf, 0.05, 0.2, 0.9, Inf),
#'                               labels = c("H", "A", "S", "Y")) # four classifications
#' family.mutiple$f__A21b <- as.factor(family.mutiple$f__A21b)
#' rf.rst.cla.mutiple <- create_ml_model(y.data = family.mutiple,
#'                                       x.data = dataset.dts$spa$tabs[,paste0("Bio", seq(19))],
#'                                       var = 'f__A21b', method = 'rf', type = 'classification')
#' rf.rst.cla.mutiple %>% evaluate_ml_model()
#' @export
evaluate_ml_model = function(model.dat, only.show.class = TRUE){
    if (model.dat$model.type == 'regression'){
        model <- model.dat$model; tests.dat <- model.dat$tests.dat; train.dat <- model.dat$train.dat
        obs.val <- tests.dat$target; prd.dat <- tests.dat[,which(colnames(tests.dat) != "target")]
        predictions <- stats::predict(model, newdata = prd.dat); rmse <- sqrt(mean((obs.val - predictions)^2))
        r <- psych::corr.test(obs.val, predictions, method = 'pearson', adjust = 'fdr')$r
        p <- psych::corr.test(obs.val, predictions, method = 'pearson', adjust = 'fdr')$p
        plotdata  <- data.frame(obs = obs.val, prd = predictions)
        labedata1 <- data.frame(x = -Inf, y = Inf,
                                label = paste0("Training size: ",
                                               nrow(train.dat), " | ",
                                               "Testing  size: ", nrow(tests.dat)))
        labedata2 <- data.frame(x = -Inf, y = Inf,
                                label = paste0("RMSE = ",
                                               round(rmse, 3),
                                               ", Pearson R = ", round(r, 3), ", P = ", round(p, 3)))
        ggplot(plotdata, aes(x = obs, y = prd)) + geom_point(shape = 21, size = 3, fill = 'gray60', alpha = 0.7) +
            geom_smooth(method = 'lm', formula = y ~ x) + xlab(paste0("#Observed ", model.dat$var.name)) +
            ylab(paste0("#Predicted ", model.dat$var.name)) +
            geom_text(data = labedata1, aes(x = x, y = y, label = label),
                      vjust= 2, hjust = -0.02, size = 5, color = 'black') +
            geom_text(data = labedata2, aes(x = x, y = y, label = label),
                      vjust= 4, hjust = -0.02, size = 5, color = 'black') +
            theme_bw() + labs(title = paste0(model.dat$model.method, " model for ", model.dat$var.name)) +
            theme(axis.text = element_text(size = 12),
                  axis.title = element_text(size = 14), plot.title = element_text(hjust = 0.5, face = 'bold'))
    }else if(model.dat$model.type == 'classification'){
        model <- model.dat$model; tests.dat <- model.dat$tests.dat; train.dat <- model.dat$train.dat
        labedata <- data.frame(x = -Inf, y = Inf,
                               label = paste0("Training size: ",
                                              nrow(train.dat), " | ",
                                              "Testing  size: ", nrow(tests.dat)))
        predicted <- predict(model, newdata = tests.dat[,which(colnames(tests.dat) != "target")], type = "prob")
        label.df <- lapply(colnames(predicted), function(label) as.character(tests.dat$target)) %>%
            do.call('cbind', .); colnames(label.df) <- colnames(predicted)
        for (i in seq(ncol(label.df)))
            label.df[,i] <- ifelse(label.df[,i] == colnames(label.df)[i], 1, 0) %>% as.numeric()
        label.df <- as.data.frame(label.df); rownames(label.df) <- rownames(tests.dat)
        colnames(predicted) <- paste0(colnames(predicted), "_pred_", model.dat$model.method)
        colnames(label.df) <- paste0(colnames(label.df), "_true"); final.df <- cbind(label.df, predicted)
        roc.res  <- multiROC::multi_roc(final.df, force_diag = F) %>% suppressWarnings()
        plot.roc.df <- multiROC::plot_roc_data(roc.res)
        if (only.show.class) plot.roc.df <- plot.roc.df[which(plot.roc.df$Group %in% unique(tests.dat$target)),]
        plot.roc.df$Group <- paste0(plot.roc.df$Group, ' (AUC = ', round(plot.roc.df$AUC, 3), ")")
        ggplot(plot.roc.df, aes(x = 1 - Specificity, y = Sensitivity)) +
            geom_path(aes(color = Group), linewidth = 1) +
            geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), colour = 'grey', linetype = 'dotdash') +
            geom_text(data = labedata,
                      aes(x = x, y = y, label = label), vjust= 2, hjust = -0.02, size = 5, color = 'black') +
            theme_bw() + scale_color_manual(values = RColorBrewer::brewer.pal(9, "Set1")) +
            labs(title = paste0(model.dat$model.method, " model for ", model.dat$var.name)) +
            theme(plot.title = element_text(hjust = 0.5, face = 'bold'), legend.justification = c(1, 0),
                legend.position = c(.999, .001), legend.title=element_blank(),
                legend.background = element_rect(fill = NULL, size = 0.5, linetype = "solid", colour = "black"),
                axis.text = element_text(size = 12), axis.title = element_text(size = 14),
                legend.text = element_text(size = 12))
    }
}

#' @title Predict values on a map
#' @author Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
#' @description This function is implemented to predict the value on a geographic map based on the machine learning model
#' created by \code{microgeo::create_ml_model()()}. The results generated by this function can be saved into a file using
#' \code{terra::writeRaster()}. See \code{terra::writeRaster()} for more details.
#' @param model.dat A machine learning model object returned by \code{microgeo::create_ml_model()}.
#' @param spat.raster A `SpatRaster` object used as prediction variables. Variable names in such a `SpatRaster` should be
#' same as those in ML modeling (<x.data>).
#' @param ... Parameters parsed by \code{terra::predict()}.
#' @return A `SpatRaster`.
#' @seealso
#' \code{\link[terra:predict]{terra::predict()}}
#' \code{\link[terra:writeRaster]{terra::writeRaster()}}
#' \code{\link[microgeo:create_ml_model]{microgeo::create_ml_model()}}
#' \code{\link[microgeo:evaluate_ml_model]{microgeo::evaluate_ml_model()}}
#' \code{\link[microgeo:read_aliyun_map]{microgeo::read_aliyun_map()}}
#' \code{\link[microgeo:create_dataset]{microgeo::create_dataset()}}
#' \code{\link[microgeo:show_dataset]{microgeo::show_dataset()}}
#' \code{\link[microgeo:get_ai]{microgeo::get_ai()}}
#' \code{\link[microgeo:get_his_bioc]{microgeo::get_his_bioc()}}
#' \code{\link[microgeo:get_fut_bioc]{microgeo::get_fut_bioc()}}
#' \code{\link[microgeo:get_modis_cla_metrics]{microgeo::get_modis_cla_metrics()}}
#' \code{\link[microgeo:extract_data_from_spatraster]{microgeo::extract_data_from_spatraster()}}
#' \code{\link[microgeo:rarefy_count_table]{microgeo::rarefy_count_table()}}
#' \code{\link[microgeo:tidy_dataset]{microgeo::tidy_dataset()}}
#' \code{\link[microgeo:calc_rel_abund]{microgeo::calc_rel_abund()}}
#' \code{\link[microgeo:calc_markers]{microgeo::calc_markers()}}
#' \code{\link[microgeo:mask_spatraster_by_cla]{microgeo::mask_spatraster_by_cla()}}
#' \code{\link[microgeo:plot_bmap]{microgeo::plot_nmap()}}
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
#' # Download aridity index and 19 historical bioclimatic variables
#' dataset.dts %<>% get_ai(out.dir = "test/microgeo_data") # Aridity index
#' dataset.dts %<>% get_his_bioc(res = 2.5, out.dir = "test/microgeo_data")
#' dataset.dts %<>% get_fut_bioc(res = 2.5, out.dir = "test/microgeo_data", gcm = "BCC-CSM2-MR")
#' dataset.dts %<>% get_modis_cla_metrics(username = "bioape.lichaonan", password = "Lichaonan@19910905", out.dir = "test/microgeo_data")
#' dataset.dts %<>% extract_data_from_spatraster(type = 'his')
#'
#' # Tidy up microgeo dataset
#' dataset.dts %<>% rarefy_count_table()
#' dataset.dts %<>% tidy_dataset()
#' dataset.dts %>% show_dataset()
#'
#' # Calculate the relative abundance and the ecologycal markers at `Family` level
#' dataset.dts %<>% calc_rel_abund()
#' dataset.dts %<>% calc_markers(use.dat = 'spa', use.var = 'AI', annotation.level = 'Family', r.thres = 0.3)
#' head(dataset.dts$abd$mar$correlation)
#'
#' # Create a regression random forest model for soil pH
#' rf.rst.reg <- create_ml_model(y.data = dataset.dts$env,
#'                               x.data = dataset.dts$spa$tabs[,paste0("Bio", seq(19))], var = 'pH', method = 'rf')
#' rf.rst.reg %>% evaluate_ml_model()
#'
#' # Predict soil pH using 19 historically bioclimatic variables
#' rf.rst.reg.pred.his <- rf.rst.reg %>%
#'     predict_ml_geomap(spat.raster = dataset.dts$spa$rast$his[[paste0("Bio", seq(19))]])
#' dataset.dts$map %>% plot_bmap() %>%
#'     add_spatraster(spat.raster = rf.rst.reg.pred.his) %>%
#'     add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = "X.CENTER", lat.var = "Y.CENTER") %>%
#'     add_scale_bar() %>% add_north_arrow() %>% add_crs()
#'
#' # Mask the results by using Grasslands(10) and Barren(16)
#' rf.rst.reg.pred.his.masked.gb <- mask_spatraster_by_cla(tar.spat = rf.rst.reg.pred.his,
#'                                                         ref.spat = dataset.dts$spa$rast$cla$LC_Type1,
#'                                                         use.class = c(10, 16))
#' dataset.dts$map %>% plot_bmap() %>%
#'     add_spatraster(spat.raster = rf.rst.reg.pred.his.masked.gb) %>%
#'     add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = "X.CENTER", lat.var = "Y.CENTER") %>%
#'     add_scale_bar() %>% add_north_arrow() %>% add_crs()
#'
#' # Mask the results by using Grasslands(10)
#' rf.rst.reg.pred.his.masked.g <- mask_spatraster_by_cla(tar.spat = rf.rst.reg.pred.his,
#'                                                        ref.spat = dataset.dts$spa$rast$cla$LC_Type1,
#'                                                        use.class = 10)
#' dataset.dts$map %>% plot_bmap() %>%
#'     add_spatraster(spat.raster = rf.rst.reg.pred.his.masked.g) %>%
#'     add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = "X.CENTER", lat.var = "Y.CENTER") %>%
#'     add_scale_bar() %>% add_north_arrow() %>% add_crs()
#'
#' # Mask the results by using Forest (1,2,3,4,5)
#' rf.rst.reg.pred.his.masked.f <- mask_spatraster_by_cla(tar.spat = rf.rst.reg.pred.his,
#'                                                        ref.spat = dataset.dts$spa$rast$cla$LC_Type1,
#'                                                        use.class = c(1,2,3,4,5))
#' dataset.dts$map %>% plot_bmap(bg.color = "gray40", gd.color = "white") %>%
#'     add_spatraster(spat.raster = rf.rst.reg.pred.his.masked.f, border.color = 'white', border.size = 1) %>%
#'     add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = "X.CENTER", lat.var = "Y.CENTER") %>%
#'     add_scale_bar() %>% add_north_arrow() %>% add_crs()
#'
#' # Predict soil pH using 19 future bioclimatic variables [`BCC-CSM2-MR|ssp585|2061-2080`]
#' rf.rst.reg.pred.fut <- rf.rst.reg %>%
#'     predict_ml_geomap(spat.raster = dataset.dts$spa$rast$fut$`BCC-CSM2-MR|ssp585|2061-2080`)
#' plot_bmap(map = dataset.dts$map) %>%
#'     add_spatraster(spat.raster = rf.rst.reg.pred.fut) %>%
#'     add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = "X.CENTER", lat.var = "Y.CENTER") %>%
#'     add_scale_bar() %>% add_north_arrow() %>% add_crs()
#'
#' # Create a binary classification random forest model for the family of f__A21b
#' family.bins <- data.frame(row.names = rownames(dataset.dts$abd$mar$abundance),
#'                           f__A21b = dataset.dts$abd$mar$abundance$f__A21b)
#' family.bins$f__A21b <- ifelse(family.bins$f__A21b > 0, "presence", "absence") # two classifications
#' family.bins$f__A21b <- as.factor(family.bins$f__A21b)
#' rf.rst.cla.bin <- create_ml_model(y.data = family.bins,
#'                                   x.data = dataset.dts$spa$tabs[,paste0("Bio", seq(19))],
#'                                   var = 'f__A21b', method = 'rf', type = 'classification')
#' rf.rst.cla.bin %>% evaluate_ml_model()
#'
#' # Predict the presence/absence probability of family f__A21b using 19 historically bioclimatic variables
#' rf.rst.cla.bin.pred.his <- rf.rst.cla.bin %>%
#'     predict_ml_geomap(spat.raster = dataset.dts$spa$rast$his[[paste0("Bio", seq(19))]])
#' dataset.dts$map %>% plot_bmap() %>%
#'     add_spatraster(spat.raster = rf.rst.cla.bin.pred.his$presence) %>%
#'     add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = "X.CENTER", lat.var = "Y.CENTER") %>%
#'     add_scale_bar() %>% add_north_arrow() %>% add_crs()
#' dataset.dts$map %>% plot_bmap() %>%
#'     add_spatraster(spat.raster = rf.rst.cla.bin.pred.his$absence) %>%
#'     add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = "X.CENTER", lat.var = "Y.CENTER") %>%
#'     add_scale_bar() %>% add_north_arrow() %>% add_crs()
#'
#' # Predict the presence/absence probability of family f__A21b
#' # using 19 future bioclimatic variables [`BCC-CSM2-MR|ssp585|2061-2080`]
#' rf.rst.cla.bin.pred.fut <- rf.rst.cla.bin %>%
#'     predict_ml_geomap(spat.raster = dataset.dts$spa$rast$fut$`BCC-CSM2-MR|ssp585|2061-2080`)
#' dataset.dts$map %>% plot_bmap() %>%
#'     add_spatraster(spat.raster = rf.rst.cla.bin.pred.fut$presence) %>%
#'     add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = "X.CENTER", lat.var = "Y.CENTER") %>%
#'     add_scale_bar() %>% add_north_arrow() %>% add_crs()
#' dataset.dts$map %>% plot_bmap() %>%
#'     add_spatraster(spat.raster = rf.rst.cla.bin.pred.fut$absence) %>%
#'     add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = "X.CENTER", lat.var = "Y.CENTER") %>%
#'     add_scale_bar() %>% add_north_arrow() %>% add_crs()
#'
#' # Create a mutiple class classification random forest model for the family of f__A21b
#' family.mutiple <- data.frame(row.names = rownames(dataset.dts$abd$mar$abundance),
#'                              f__A21b = dataset.dts$abd$mar$abundance$f__A21b)
#' family.mutiple$f__A21b <- cut(family.mutiple$f__A21b, breaks = c(-Inf, 0.05, 0.2, 0.9, Inf),
#'                               labels = c("H", "A", "S", "Y")) # four classifications
#' family.mutiple$f__A21b <- as.factor(family.mutiple$f__A21b) # must convert the classification to factor
#' rf.rst.cla.mutiple <- create_ml_model(y.data = family.mutiple,
#'                                       x.data = dataset.dts$spa$tabs[,paste0("Bio", seq(19))],
#'                                       var = 'f__A21b', method = 'rf', type = 'classification')
#' rf.rst.cla.mutiple %>% evaluate_ml_model()
#'
#' # Predict the probability of each classification for family f__A21b using 19 historically bioclimatic variables
#' rf.rst.cla.mutiple.pred.his <- rf.rst.cla.mutiple %>%
#'     predict_ml_geomap(spat.raster = dataset.dts$spa$rast$his[[paste0("Bio", seq(19))]])
#' dataset.dts$map %>% plot_bmap() %>%
#'     add_spatraster(spat.raster = rf.rst.cla.mutiple.pred.his$H) %>%
#'     add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = "X.CENTER", lat.var = "Y.CENTER") %>%
#'     add_scale_bar() %>% add_north_arrow() %>% add_crs()
#' dataset.dts$map %>% plot_bmap() %>%
#'     add_spatraster(spat.raster = rf.rst.cla.mutiple.pred.his$A) %>%
#'     add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = "X.CENTER", lat.var = "Y.CENTER") %>%
#'     add_scale_bar() %>% add_north_arrow() %>% add_crs()
#' dataset.dts$map %>% plot_bmap() %>%
#'     add_spatraster(spat.raster = rf.rst.cla.mutiple.pred.his$S) %>%
#'     add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = "X.CENTER", lat.var = "Y.CENTER") %>%
#'     add_scale_bar() %>% add_north_arrow() %>% add_crs()
#' dataset.dts$map %>% plot_bmap() %>%
#'     add_spatraster(spat.raster = rf.rst.cla.mutiple.pred.his$Y) %>%
#'     add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = "X.CENTER", lat.var = "Y.CENTER") %>%
#'     add_scale_bar() %>% add_north_arrow() %>% add_crs()
#'
#' # Predict the probability of each classification for family f__A21b
#' # using 19 future bioclimatic variables [`BCC-CSM2-MR|ssp585|2061-2080`]
#' rf.rst.cla.mutiple.pred.fut <- rf.rst.cla.mutiple %>%
#'     predict_ml_geomap(spat.raster = dataset.dts$spa$rast$fut$`BCC-CSM2-MR|ssp585|2061-2080`)
#' dataset.dts$map %>% plot_bmap() %>%
#'     add_spatraster(spat.raster = rf.rst.cla.mutiple.pred.fut$H) %>%
#'     add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = "X.CENTER", lat.var = "Y.CENTER") %>%
#'     add_scale_bar() %>% add_north_arrow() %>% add_crs()
#' dataset.dts$map %>% plot_bmap() %>%
#'     add_spatraster(spat.raster = rf.rst.cla.mutiple.pred.fut$A) %>%
#'     add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = "X.CENTER", lat.var = "Y.CENTER") %>%
#'     add_scale_bar() %>% add_north_arrow() %>% add_crs()
#' dataset.dts$map %>% plot_bmap() %>%
#'     add_spatraster(spat.raster = rf.rst.cla.mutiple.pred.fut$S) %>%
#'     add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = "X.CENTER", lat.var = "Y.CENTER") %>%
#'     add_scale_bar() %>% add_north_arrow() %>% add_crs()
#' dataset.dts$map %>% plot_bmap() %>%
#'     add_spatraster(spat.raster = rf.rst.cla.mutiple.pred.fut$Y) %>%
#'     add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = "X.CENTER", lat.var = "Y.CENTER") %>%
#'     add_scale_bar() %>% add_north_arrow() %>% add_crs()
#'
#' # Save the predicted result info file. Just an example
#' print(rf.rst.cla.mutiple.pred.fut)
#' terra::writeRaster(rf.rst.cla.mutiple.pred.fut, file = 'test/rf.rst.cla.mutiple.pred.fut.tif')
#' @export
predict_ml_geomap = function(model.dat, spat.raster, ...){
    spat.raster %>% check_spatraster(); model <- model.dat$model
    train.vars <- model.dat$train.dat[, which(colnames(model.dat$train.dat) != 'target')]
    check.length <- length(unique(colnames(train.vars) == names(spat.raster)))
    check.logics <- !unique(colnames(train.vars) == names(spat.raster))
    if (check.length > 1 || check.logics)
        stop('Invalid `SpatRaster` for prediction!')
    if (model.dat$model.type == 'regression'){
        rp <- terra::predict(spat.raster, model, na.rm = TRUE, ...)
        names(rp) <- model.dat$var.name
    }else{
        # use probability if the model type is classification
        rp <- terra::predict(spat.raster, model, na.rm = TRUE, type = 'prob', ...)
    }
    return(rp)
}
