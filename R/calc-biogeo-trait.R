# -----------------------------------------------------------------------------------------------------------------------
# Copyright (c) 2023, microgeo/Chaonan Li (licn@mtc.edu.cn).                                                            #
# The microgeo is distributed under the terms of the GPL-3 License.                                                     #
# Full license is avaliable in the file LICENSE, distributed with this package.                                         #
# -----------------------------------------------------------------------------------------------------------------------

#' @title Calculate relatively taxonomic/gene abundance
#' @author Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
#' @description This function is used to calculate the relatively taxonomic/gene abundance based on classification levels.
#' @param dataset A microgeo dataset with the class of `MicrogeoDataset`.
#' @param use.perct Whether to use the percentages for relatively taxonomic/gene abundance. Default is `TRUE`.
#' @return A `MicrogeoDataset` class with the following components:
#' \describe{
#'   \item{\code{object$mat}}{A `data.frame` of ASV/gene abundance.}
#'   \item{\code{object$ant}}{A `data.frame` of ASV/gene anootations.}
#'   \item{\code{object$met}}{A `data.frame` of sample information.}
#'   \item{\code{object$map}}{A `SpatialPolygonsDataFrame` of map.}
#'   \item{\code{object$phy}}{A phylogenetic tree with `newick` format if applicable.}
#'   \item{\code{object$env}}{A `data.frame` of measured environmental properties if applicable.}
#'   \item{\code{object$abd$raw}}{A `list` of taxonomic/gene abundance tables.}
#'   \item{\code{object$*}}{Spatial and biogeographic traits if applicable.}
#' }
#' @seealso
#' \code{\link[microgeo:read_aliyun_map]{microgeo::read_aliyun_map()}}
#' \code{\link[microgeo:create_dataset]{microgeo::create_dataset()}}
#' \code{\link[microgeo:show_dataset]{microgeo::show_dataset()}}
#' \code{\link[microgeo:rarefy_count_table]{microgeo::rarefy_count_table()}}
#' \code{\link[microgeo:tidy_dataset]{microgeo::tidy_dataset()}}
#' @examples
#' # Create a microgeo dataset
#' data(qtp)
#' map <- read_aliyun_map(adcode = c(540000, 630000, 510000))
#' dataset.dts <- create_dataset(mat = qtp$asv, ant = qtp$tax, met = qtp$met, map = map,
#'                               phy = qtp$tre, env = qtp$env, lon = "longitude", lat = "latitude")
#' dataset.dts %>% show_dataset()
#'
#' # Tidy up microgeo dataset
#' dataset.dts %<>% rarefy_count_table()
#' dataset.dts %<>% tidy_dataset()
#' dataset.dts %>% show_dataset()
#'
#' # Calculate relative abundance
#' dataset.dts %<>% calc_rel_abund()
#' dataset.dts %>% show_dataset()
#'
#' # Check the results
#' head(dataset.dts$abd$raw$Phylum[,1:5])
#' @export
calc_rel_abund = function(dataset, use.perct = TRUE){
    dataset %>% check_dataset(); mat <- dataset$mat; ant <- dataset$ant; met <- dataset$met
    check.length <- unique(rownames(mat) == rownames(ant)) %>% length
    check.logics <- !unique(rownames(mat) == rownames(ant))
    if (check.length > 1 || check.logics) stop('The ASV/gene ids in <mat> can not be matched to those in <ant>!')
    classifications <- colnames(ant)
    abund.res <- lapply(classifications, function(classification){
        abund.agg <- aggregate(mat, by = list(ant[,which(colnames(ant) == classification)]), FUN = sum)
        rownames(abund.agg) <- abund.agg[,1]; abund.agg <- abund.agg[,-1] %>% t()
        abund.agg.abund <- abund.agg/rowSums(abund.agg)
        if (use.perct){
            abund.agg.abund <- data.frame(row.names = rownames(abund.agg.abund), abund.agg.abund * 100)
        }else{
            abund.agg.abund <- data.frame(row.names = rownames(abund.agg.abund), abund.agg.abund)
        }
        abund.agg.abund %<>% sort_samples(met = met)
    })
    names(abund.res) <- classifications; dataset$abd$raw <- abund.res
    show_stat_msg('abd$raw'); return(dataset)
}

#' @title Identify ecological markers highly correlated with one or mutiple abiotic properties
#' @author Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
#' @description This function is implemented to identify ecological markers highly correlate with one or mutiple abiotic
#' properties. Particularly, it requires the relative abundance data calculated by \code{microgeo::calc_rel_abund()}.
#' @param dataset A microgeo dataset with the class of `MicrogeoDataset`.
#' @param use.var Variables used for correlations. If there are only one variable, then the correlation coefficient would
#' be calculated by using \code{psych::corr.test()}. If there are mutiple variables, the correlation coefficient would be
#' calculated by using \code{vegan::mantel()}.
#' @param use.dat Which type of data would be used for correlations? Select one from `env` (measured properties) and `spa`
#' (the data extracted from spatial dataset). Default is `env`.
#' @param annotation.level Which annotation level would be used? Default is `Phylum`.
#' @param cor.method Correlation method, as accepted by cor: `pearson`, `spearman` or `kendall`. Default is `pearson`.
#' @param dist.method Method for distance calculations if there are mutiple variables (<use.var>). Default is `bray`. See
#' \code{vegan::vegdist()} for more options.
#' @param p.adjust.method A method to adjust P value. Select from `holm`, `hochberg`, `hommel`, `bonferroni`, `BH`, `BY`,
#' `fdr`, `none`. Default is `none`.
#' @param r.thres Threshold of correlation coefficient to identify markers. Default is `0.1`.
#' @param p.thres Threshold of P value to identify markers. Default is `0.05`.
#' @param mantel.permutations Permutations for mantel test if there are mutiple variables (<use.var>). Default is `999`.
#' @param mantel.parallel Parallel number for mantel test if there are mutiple variables (<use.var>). Default is `2`.
#' @return A `MicrogeoDataset` class with the following components:
#' \describe{
#'   \item{\code{object$mat}}{A `data.frame` of ASV/gene abundance.}
#'   \item{\code{object$ant}}{A `data.frame` of ASV/gene anootations.}
#'   \item{\code{object$met}}{A `data.frame` of sample information.}
#'   \item{\code{object$map}}{A `SpatialPolygonsDataFrame` of map.}
#'   \item{\code{object$phy}}{A phylogenetic tree with `newick` format if applicable.}
#'   \item{\code{object$env}}{A `data.frame` of measured environmental properties if applicable.}
#'   \item{\code{object$abd$raw}}{A `list` of taxonomic/gene abundance tables.}
#'   \item{\code{object$abd$mar}}{A `list` of taxonomic/gene marker tables.}
#'   \item{\code{object$*}}{Spatial and biogeographic traits if applicable.}
#' }
#' @seealso
#' \code{\link[psych:corr.test]{psych::corr.test()}}
#' \code{\link[vegan:mantel]{vegan::mantel()}}
#' \code{\link[microgeo:read_aliyun_map]{microgeo::read_aliyun_map()}}
#' \code{\link[microgeo:create_dataset]{microgeo::create_dataset()}}
#' \code{\link[microgeo:show_dataset]{microgeo::show_dataset()}}
#' \code{\link[microgeo:get_his_bioc]{microgeo::get_his_bioc()}}
#' \code{\link[microgeo:extract_data_from_spatraster]{microgeo::extract_data_from_spatraster()}}
#' \code{\link[microgeo:rarefy_count_table]{microgeo::rarefy_count_table()}}
#' \code{\link[microgeo:tidy_dataset]{microgeo::tidy_dataset()}}
#' \code{\link[microgeo:calc_rel_abund]{microgeo::calc_rel_abund()}}
#' @examples
#' # Create a microgeo dataset
#' data(qtp)
#' map <- read_aliyun_map(adcode = c(540000, 630000, 510000))
#' dataset.dts <- create_dataset(mat = qtp$asv, ant = qtp$tax, met = qtp$met, map = map,
#'                               phy = qtp$tre, env = qtp$env, lon = "longitude", lat = "latitude")
#' dataset.dts %>% show_dataset()
#'
#' # Extract 19 historical bioclimatic variables
#' dataset.dts %<>% get_his_bioc(res = 2.5, out.dir = "test/microgeo_data")
#' dataset.dts %<>% extract_data_from_spatraster(type = 'his')
#' dataset.dts %>% show_dataset()
#'
#' # Tidy up microgeo dataset
#' dataset.dts %<>% rarefy_count_table()
#' dataset.dts %<>% tidy_dataset()
#' dataset.dts %>% show_dataset()
#'
#' # Calculate relative abundance
#' dataset.dts %<>% calc_rel_abund()
#' dataset.dts %>% show_dataset()
#'
#' # Identify ecological markers based on soil pH in <env>
#' # Correlation coefficients would be calculated by `psych::corr.test()`
#' dataset.dts %<>% calc_markers(use.var = 'pH', annotation.level = 'Phylum', r.thres = 0.1)
#' dataset.dts %>% show_dataset()
#' head(dataset.dts$abd$mar$correlation)
#'
#' # Identify ecological markers based on Bio12 in <spa>
#' # Correlation coefficients would be calculated by `psych::corr.test()`
#' dataset.dts %<>% calc_markers(use.var = 'Bio12', use.dat = 'spa', annotation.level = 'Phylum', r.thres = 0.1)
#' dataset.dts %>% show_dataset()
#' head(dataset.dts$abd$mar$correlation)
#'
#' # Identify ecological markers based on soil pH and TOC in <env>
#' # Correlation coefficients would be calculated by `vegan::mantel()`
#' dataset.dts %<>% calc_markers(use.var = c('pH', 'TOC'), annotation.level = 'Phylum', r.thres = 0.1)
#' dataset.dts %>% show_dataset()
#' head(dataset.dts$abd$mar$correlation)
#' @export
calc_markers = function(dataset, use.var, use.dat = c('env', 'spa'), annotation.level = 'Phylum',
                        cor.method = c('pearson', 'spearman', 'kendall'),
                        dist.method = 'bray', p.adjust.method = 'none',
                        r.thres = 0.1, p.thres = 0.05,
                        mantel.permutations = 999,
                        mantel.parallel = 2){

    # check dataset and arguments
    dataset %>% check_dataset()
    if (dataset$abd$raw %>% is.null)
        stop("Please calculate the relative abundance of ASVs/genes by using `calc_rel_abund()`!")
    if (dataset$abd$raw[[annotation.level]] %>% is.null)
        paste0('No `', annotation.level, '` in the annotation table, please check your dataset!') %>% stop()
    cor.method <- ifelse(cor.method %>% length > 1, cor.method[1], cor.method)
    if (!cor.method %in% c('pearson', 'spearman', 'kendall'))
        stop("The <cor.method> must be one of `pearson`, `spearman` and `kendall`!")
    use.dat <- ifelse(use.dat %>% length > 1, use.dat[1], use.dat)
    if (!use.dat %in% c('env', 'spa')) stop("The <use.dat> must be one of `env` and `spa`!")
    if (dataset[[use.dat]] %>% is.null) paste0("No <", use.dat, "> data in your dataset!") %>% stop()
    if (use.dat == 'env'){
        var.detect.res <- intersect(dataset[[use.dat]] %>% colnames, use.var)
    }else{
        if (dataset[[use.dat]]$tabs %>% is.null)
            paste0("No extracted <", use.dat,
                   "> data in your dataset, please run `extract_data_from_spatraster()`!") %>% stop()
        var.detect.res <- intersect(dataset[[use.dat]]$tabs %>% colnames, use.var)
    }
    if (var.detect.res %>% length != use.var %>% length)
        paste0("Some variables in <use.var> can not be found in your <", use.dat, "> data!") %>% stop()

    # prepare <env> data
    if (use.var %>% length > 1) {
        if (use.dat == 'env'){
            env.dat <- dataset[[use.dat]][,use.var]
        }else{
            env.dat <- dataset[[use.dat]]$tabs[,use.var]
        }
    }else{
        if (use.dat == 'env'){
            env.dat <- data.frame(row.names = rownames(dataset[[use.dat]]),
                                  dataset[[use.dat]][,use.var])
        }else{
            env.dat <- data.frame(row.names = rownames(dataset[[use.dat]]$tabs),
                                  dataset[[use.dat]]$tabs[,use.var])
        }
        colnames(env.dat) <- use.var
    }

    # match samples
    abd.dat <- dataset$abd$raw[[annotation.level]]
    check.length <- unique(env.dat %>% rownames == abd.dat %>% rownames) %>% length
    check.logics <- !unique(env.dat %>% rownames == abd.dat %>% rownames)
    if (check.length > 1 | check.logics)
        stop('Sample ids in relative abundance table can not be matched to those in env/spa table!')

    # calculate correlation coefficient
    if (use.var %>% length == 1){
        cor.res <- psych::corr.test(abd.dat, env.dat, method = cor.method, adjust = p.adjust.method)
        cor.r <- cor.res$r; cor.p <- cor.res$p.adj
        cor.df <- as.data.frame(cbind(cor.r, cor.p)); colnames(cor.df) <- c("r", "p")
        cor.df <- data.frame(var = cor.df %>% rownames, cor.df)
    }else{
        paste0("mantel test would take a while...") %>% show_comm_msg()
        annot.names <- colnames(abd.dat)
        cor.df <- lapply(annot.names %>% length %>% seq, function(ct){
            annot.n <- annot.names[ct]
            paste0('current annot. name(', ct, '/', length(annot.names), '): ',
                   annot.n) %>% show_comm_msg()
            tmp.abund <- data.frame(row.names = abd.dat %>% rownames, abd.dat[,annot.n]); colnames(tmp.abund) <- annot.n
            tmp.cor.res <- vegan::mantel(vegan::vegdist(tmp.abund + 0.000000001, method = dist.method),
                                         dist(env.dat), method = cor.method, permutations = mantel.permutations,
                                         parallel = mantel.parallel) # add 0.000000001 to avoid empty row error
            res <- data.frame(var = annot.n, r = tmp.cor.res$statistic, p = tmp.cor.res$signif)
        }) %>% do.call("rbind", .)
        if (p.adjust.method != 'none') cor.df$p <- stats::p.adjust(cor.df$p, method = p.adjust.method)
        rownames(cor.df) <- cor.df$var
    }

    # filter results
    cor.filter <- cor.df[which(cor.df$r %>% abs >= r.thres & cor.df$p <= p.thres),]
    if (cor.filter %>% nrow == 0)
        paste0("Can not find any markers when `|r.thres| >= ", r.thres,
               " & p.thres < ", p.thres, "`!") %>% stop()
    colnames(abd.dat) <- lapply(colnames(abd.dat),function(n) gsub(pattern = " ", replacement = ".", x = n)) %>% unlist()
    if (cor.filter$var %>% length > 1){
        cor.filter <- cor.filter[order(cor.filter$r %>% abs, decreasing = T),]
        cor.filter$var <- factor(cor.filter$var, levels = cor.filter$var %>% as.character)
        marker.abund <- abd.dat[,cor.filter$var %>% as.character]
    }else{
        marker.abund <- data.frame(row.names = abd.dat %>% rownames, abd.dat[,cor.filter$var])
        colnames(marker.abund) <- cor.filter$var
    }
    paste0('found ', ncol(marker.abund),
           ' ecological markers at `', annotation.level, '` level..') %>% show_comm_msg()
    dataset$abd$mar <- list(abundance = marker.abund,
                            correlation = cor.filter, annotation.level = annotation.level)
    show_stat_msg('abd$mar')
    return(dataset)
}

#' @title Calculate alpha diversity indices
#' @author Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
#' @description This function is used to calculate alpha diversity indices based on a ASV/gene count/abundance matrices.
#' @param dataset A microgeo dataset with the class of `MicrogeoDataset`.
#' @param measures Diversity measures. Default is: `"observed", "chao1", "ace", "shannon", "simpson", "invsimpson", "pd"`.
#' @return A `MicrogeoDataset` class with the following components:
#' \describe{
#'   \item{\code{object$mat}}{A `data.frame` of ASV/gene abundance.}
#'   \item{\code{object$ant}}{A `data.frame` of ASV/gene anootations.}
#'   \item{\code{object$met}}{A `data.frame` of sample information.}
#'   \item{\code{object$map}}{A `SpatialPolygonsDataFrame` of map.}
#'   \item{\code{object$phy}}{A phylogenetic tree with `newick` format if applicable.}
#'   \item{\code{object$env}}{A `data.frame` of measured environmental properties if applicable.}
#'   \item{\code{object$div$alpha}}{A `data.frame` of alpha diversity indices.}
#'   \item{\code{object$*}}{Spatial and biogeographic traits if applicable.}
#' }
#' @seealso
#' \code{\link[vegan:estimateR]{vegan::estimateR()}}
#' \code{\link[vegan:diversity]{vegan::diversity()}}
#' \code{\link[picante:pd]{picante::pd()}}
#' \code{\link[microgeo:read_aliyun_map]{microgeo::read_aliyun_map()}}
#' \code{\link[microgeo:create_dataset]{microgeo::create_dataset()}}
#' \code{\link[microgeo:show_dataset]{microgeo::show_dataset()}}
#' \code{\link[microgeo:rarefy_count_table]{microgeo::rarefy_count_table()}}
#' \code{\link[microgeo:tidy_dataset]{microgeo::tidy_dataset()}}
#' @examples
#' # Create a microgeo dataset
#' data(qtp)
#' map <- read_aliyun_map(adcode = c(540000, 630000, 510000))
#' dataset.dts <- create_dataset(mat = qtp$asv, ant = qtp$tax, met = qtp$met, map = map,
#'                               phy = qtp$tre, env = qtp$env, lon = "longitude", lat = "latitude")
#' dataset.dts %>% show_dataset()
#'
#' # Tidy up microgeo dataset
#' dataset.dts %<>% rarefy_count_table()
#' dataset.dts %<>% tidy_dataset()
#' dataset.dts %>% show_dataset()
#'
#' # Calculate alpha diversity indices
#' dataset.dts %<>% calc_alpha_div(measures = c("observed", "shannon"))
#' dataset.dts %>%  show_dataset()
#'
#' # Check the results
#' head(dataset.dts$div$alpha)
#' @export
calc_alpha_div = function(dataset, measures = c('observed', 'chao1', 'ace', 'shannon', 'simpson', 'invsimpson', 'pd')){

    # calculate alpha diversity indices by using `vegan::estimateR`.
    dataset %>% check_dataset(); estimateR.indices <- NULL
    estimateR.measures <- measures[which(measures %in% c("observed", "chao1", "ace"))]
    if (estimateR.measures %>% length > 0){
        estimateR.indices <- dataset$mat %>% t() %>% vegan::estimateR() %>% t()
        estimateR.indices <- estimateR.indices[,c(1,2,4)]
        colnames(estimateR.indices) <- c("observed", "chao1", "ace")
        if (estimateR.measures %>% length > 1) {
            estimateR.indices <- estimateR.indices[,which(colnames(estimateR.indices) %in% estimateR.measures)]
        }else{
            estimateR.indices <- data.frame(
                row.names = estimateR.indices %>% rownames,
                val = estimateR.indices[,which(colnames(estimateR.indices) %in% estimateR.measures)]
            )
            colnames(estimateR.indices) <- estimateR.measures
        }
    }

    # calculate alpha diversity indices by using `vegan::diversity`
    diversity.measures <- measures[which(measures %in% c("shannon", "simpson", "invsimpson"))]; diversity.indices <- NULL
    if (diversity.measures %>% length > 0){
        diversity.indices <- lapply(diversity.measures, function(diversity.measure){
            res <- vegan::diversity(dataset$mat %>% t, index = diversity.measure) %>% as.data.frame()
            colnames(res) <- diversity.measure; res
        })
        diversity.indices <- do.call("cbind", diversity.indices)
    }

    # calculate alpha diversity indices by using `picante::pd`
    pd.measures <- measures[which(measures == "pd")]; pd.indices <- NULL
    if (pd.measures %>% length > 0){
        if (dataset$phy %>% is.null) stop("A rooted phylogenetic tree is required for pd calculation!")
        if (dataset$phy$edge.length %>% is.null) stop("Phylogenetic tree has no branch lengths, cannot compute pd!")
        pd <- picante::pd(dataset$mat %>% t, dataset$phy); colnames(pd) <- c("pd", 'sr')
        pd.indices <- data.frame(row.names = pd %>% rownames, pd = pd$pd)
    }

    # combine all avaliable alpha diversity indices
    samples <- dataset$mat %>% t() %>% rownames()
    alpha.div <- data.frame(row.names = samples)
    if (!estimateR.indices %>% is.null) alpha.div <- cbind(alpha.div, estimateR.indices)
    if (!diversity.indices %>% is.null) alpha.div <- cbind(alpha.div, diversity.indices)
    if (!pd.indices %>% is.null) alpha.div <- cbind(alpha.div, pd.indices)
    alpha.div %<>% sort_samples(met = dataset$met)
    dataset$div$alpha <- alpha.div
    show_stat_msg('div$alpha')
    return(dataset)
}

#' @title Calculate beta diversity matrices
#' @author Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
#' @description This function is used to calculate beta diversity distances based on a ASV/gene count/abundance matrix.
#' @param dataset A microgeo dataset with the class of `MicrogeoDataset`.
#' @param measures Diversity distance measures. Default is: `'bray', 'jaccard', 'euclidean', 'unifrac'`.
#' @return A `MicrogeoDataset` class with the following components:
#' \describe{
#'   \item{\code{object$mat}}{A `data.frame` of ASV/gene abundance.}
#'   \item{\code{object$ant}}{A `data.frame` of ASV/gene anootations.}
#'   \item{\code{object$met}}{A `data.frame` of sample information.}
#'   \item{\code{object$map}}{A `SpatialPolygonsDataFrame` of map.}
#'   \item{\code{object$phy}}{A phylogenetic tree with `newick` format if applicable.}
#'   \item{\code{object$env}}{A `data.frame` of measured environmental properties if applicable.}
#'   \item{\code{object$div$beta}}{A `list` of beta diversity distance matrices.}
#'   \item{\code{object$*}}{Spatial and biogeographic traits if applicable.}
#' }
#' @seealso
#' \code{\link[vegan:vegdist]{vegan::vegdist()}}
#' \code{\link[GUniFrac:GUniFrac]{GUniFrac::GUniFrac()}}
#' \code{\link[microgeo:read_aliyun_map]{microgeo::read_aliyun_map()}}
#' \code{\link[microgeo:create_dataset]{microgeo::create_dataset()}}
#' \code{\link[microgeo:show_dataset]{microgeo::show_dataset()}}
#' \code{\link[microgeo:rarefy_count_table]{microgeo::rarefy_count_table()}}
#' \code{\link[microgeo:tidy_dataset]{microgeo::tidy_dataset()}}
#' @examples
#' # Create a microgeo dataset
#' data(qtp)
#' map <- read_aliyun_map(adcode = c(540000, 630000, 510000))
#' dataset.dts <- create_dataset(mat = qtp$asv, ant = qtp$tax, met = qtp$met, map = map,
#'                               phy = qtp$tre, env = qtp$env, lon = "longitude", lat = "latitude")
#' dataset.dts %>% show_dataset()
#'
#' # Tidy up microgeo dataset
#' dataset.dts %<>% rarefy_count_table()
#' dataset.dts %<>% tidy_dataset()
#' dataset.dts %>% show_dataset()
#'
#' # Calculate beta diversity indices
#' dataset.dts %<>% calc_beta_div(measures = c("bray", "jaccard"))
#' dataset.dts %>%  show_dataset()
#'
#' # Check the results
#' names(dataset.dts$div$bet)
#' dataset.dts$div$bet$bray[1:4, 1:4]
#' @export
calc_beta_div = function(dataset, measures = c('bray', 'jaccard', 'euclidean', 'unifrac')){

    # calculate taxonomic beta diversity indices
    dataset %>% check_dataset(); taxa.beta.res <- list()
    taxa.beta.idx <- measures[which(measures %in% c('bray', 'jaccard', 'euclidean'))]
    if (taxa.beta.idx %>% length > 0){
        taxa.beta.res <- lapply(taxa.beta.idx, function(idx){
            binary <- ifelse(idx == 'jaccard', TRUE, FALSE)
            dist.matrix <- as.matrix(vegan::vegdist(t(dataset$mat), method = idx, binary = binary))
        })
        names(taxa.beta.res) <- taxa.beta.idx
    }

    # calculate phylogenetic beta diversity indices
    phyl.beta.res <- list()
    if ('unifrac' %in% measures){
        if (dataset$phy %>% is.null) stop("A rooted phylogenetic tree is required for unifrac calculation!")
        if (dataset$phy$edge.length %>% is.null) stop("Phylogenetic tree has no branch lengths, cannot compute unifrac!")
        unifrac.obj <- GUniFrac::GUniFrac(dataset$mat %>% t, phy, alpha = c(0, 0.5, 1))
        weighted.unifrac <- as.matrix(unifrac.obj$unifracs[,, "d_1"])
        unweighted.unifrac <- as.matrix(unifrac.obj$unifracs[,, "d_UW"])
        phyl.beta.res <- list(weighted.unifrac = weighted.unifrac, unweighted.unifrac = unweighted.unifrac)
    }

    # combine all distance matrices into one list and save to RDS
    beta.div.matrix <- append(taxa.beta.res, phyl.beta.res)
    beta.div.matrix.sorted <- lapply(beta.div.matrix, function(div.matrix){
        div.matrix %<>% sort_samples(met = dataset$met, is.matrix = TRUE)
    })
    names(beta.div.matrix.sorted) <- beta.div.matrix %>% names()
    dataset$div$beta <- beta.div.matrix.sorted
    show_stat_msg('div$beta'); return(dataset)
}

#' @title Calculate the matrices related to microbial community assembly
#' @author Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
#' @description This function is used to calculate matrices related to microbial community assembly based on null model.
#' @param dataset A microgeo dataset with the class of `MicrogeoDataset`.
#' @param type Which type of matrices would be used? Select one from `'alpha.phylo', 'beta.phylo'`. A `alpha.phylo` means
#' alpha-type phylogenetic null model (e.g., ses.MNTD), and a `beta.phylo` means beta-type phylogenetic null model (e.g.,
#' betaNTI). Default is `alpha.phylo`.
#' @param model Null model used to calculate `alpha.phylo` indices. Default is `taxa.labels`.
#' @param runs How many runs would be used for calculation? Defaults is `999`.
#' @param nworker How many threads would be used for calculation? Defaults is `2`.
#' @param memory.G How many memory would be used for calculation? Defaults is `10` (GB).
#' @param sig.bNTI Cutoffs for betaNTI. Required if the <type> is `beta.phylo`. Default is `2`.
#' @param sig.rc Cutoffs for RC. Required if the <type> is `beta.phylo`. Default is `0.95`.
#' @param abundance.weighted Use abundance weighted? Defaults is `TRUE`.
#' @param out.dir A directory saving output files. Default is `calc_comm_asmb.rst`.
#' @return A `MicrogeoDataset` class with the following components:
#' \describe{
#'   \item{\code{object$mat}}{A `data.frame` of ASV/gene abundance.}
#'   \item{\code{object$ant}}{A `data.frame` of ASV/gene anootations.}
#'   \item{\code{object$met}}{A `data.frame` of sample information.}
#'   \item{\code{object$map}}{A `SpatialPolygonsDataFrame` of map.}
#'   \item{\code{object$phy}}{A phylogenetic tree with `newick` format if applicable.}
#'   \item{\code{object$env}}{A `data.frame` of measured environmental properties if applicable.}
#'   \item{\code{object$asb}}{Matrices related to microbial community assembly.}
#'   \item{\code{object$*}}{Spatial and biogeographic traits if applicable.}
#' }
#' @seealso
#' \code{\link[iCAMP:pdist.big]{iCAMP::pdist.big()}}
#' \code{\link[iCAMP:qpen]{iCAMP::qpen()}}
#' \code{\link[picante:ses.mntd]{picante::ses.mntd()}}
#' \code{\link[microgeo:read_aliyun_map]{microgeo::read_aliyun_map()}}
#' \code{\link[microgeo:create_dataset]{microgeo::create_dataset()}}
#' \code{\link[microgeo:show_dataset]{microgeo::show_dataset()}}
#' \code{\link[microgeo:rarefy_count_table]{microgeo::rarefy_count_table()}}
#' \code{\link[microgeo:tidy_dataset]{microgeo::tidy_dataset()}}
#' @examples
#' # Create a microgeo dataset
#' data(qtp)
#' map <- read_aliyun_map(adcode = c(540000, 630000, 510000))
#' dataset.dts <- create_dataset(mat = qtp$asv, ant = qtp$tax, met = qtp$met, map = map,
#'                               phy = qtp$tre, env = qtp$env, lon = "longitude", lat = "latitude")
#' dataset.dts %>% show_dataset()
#'
#' # Tidy up microgeo dataset
#' dataset.dts %<>% rarefy_count_table()
#' dataset.dts %<>% tidy_dataset()
#' dataset.dts %>% show_dataset()
#'
#' # Calculate `alpha.phylo` null model
#' # runs = 9 just for an example
#' # runs = 999 would be better
#' dataset.dts %<>% calc_phylo_asmb(type = 'alpha.phylo', runs = 9, out.dir = 'test/calc_comm_asmb.rst')
#' dataset.dts %>%  show_dataset()
#' head(dataset.dts$asb$alpha.phylo)
#'
#' # Calculate `beta.phylo` null model
#' # runs = 9 just for an example
#' # runs = 999 would be better
#' dataset.dts %<>% calc_phylo_asmb(type = 'beta.phylo', runs = 9, out.dir = 'test/calc_comm_asmb.rst')
#' dataset.dts %>%  show_dataset()
#' names(dataset.dts$asb$beta.phylo$dis) # distance matrices
#' head(dataset.dts$asb$beta.phylo$raw$result) # raw result of `iCAMP::qpen()`
#' @export
calc_phylo_asmb = function(
    dataset, type = c('alpha.phylo', 'beta.phylo'),
    model = c('taxa.labels', 'richness', 'frequency', 'sample.pool', 'phylogeny.pool', 'independentswap', 'trialswap'),
    runs = 999, nworker = 2, memory.G = 10, sig.bNTI = 2, sig.rc = 0.95, abundance.weighted = TRUE,
    out.dir = 'calc_comm_asmb.rst'){

    # check dataset and arguments
    dataset %>% check_dataset(); mat <- dataset$mat; met <- dataset$met; phy = dataset$phy
    if (dataset$phy %>% is.null) stop("A rooted phylogenetic tree is required for null model!")
    if (dataset$phy$edge.length %>% is.null) stop("Phylogenetic tree has no branch lengths, apply it for null model!")
    type <- ifelse(type %>% length > 1, type[1], type); model <- ifelse(model %>% length > 1, model[1], model)
    if (!type %in% c('alpha.phylo', 'beta.phylo')) stop('The <type> must be one of `alpha.phylo` and `beta.phylo`!')

    # set the output file path
    if (type == 'alpha.phylo'){
        outfile <- file.path(out.dir, 'alpha_phylo_asmb.Rds')
        tmppath <- file.path(out.dir, 'alpha.phylo.tmp') %>% create_dir(recursive = T)
    }else if (type == 'beta.phylo'){
        outfile <- file.path(out.dir, 'beta_phylo_asmb.Rds')
        tmppath <- file.path(out.dir, 'beta.phylo.tmp') %>% create_dir(recursive = T)
    }

    # `alpha.phylo` null model
    if (type == 'alpha.phylo'){ #
        if (!outfile %>% file.exists){
            show_comm_msg("calculating phylogenetic distances, please wait...")
            phylo.dist <- iCAMP::pdist.big(
                tree = phy,
                wd = tmppath,
                tree.asbig = FALSE,
                output = TRUE,
                nworker = nworker,
                nworker.pd = nworker,
                memory.G = memory.G
            ) %>% suppressWarnings() %>% suppressMessages()
            asv.idx <- sapply(mat %>% rownames, function(x){which(phylo.dist %>% rownames == x)})
            phylo.dist.sig <- as.matrix(phylo.dist[asv.idx, asv.idx])
            if (length(unique(rownames(mat) == rownames(phylo.dist.sig))) > 1 ||
                !unique(rownames(mat) == rownames(phylo.dist.sig))){
                stop("The ASV/gene ids in <mat> can not be matched to those in phylogenetic distance matrix!")
            }
            if (length(unique(rownames(mat) == colnames(phylo.dist.sig))) > 1 ||
                !unique(rownames(mat) == colnames(phylo.dist.sig))){
                stop("The ASV/gene ids in <mat> can not be matched to those in phylogenetic distance matrix!")
            }
            show_comm_msg("calculating community assembly indices[alpha.phylo], please wait....")
            ses.mntd.obj <- picante::ses.mntd(
                samp = t(mat),
                dis = phylo.dist.sig,
                null.model = model,
                abundance.weighted = abundance.weighted,
                runs = runs
            ) %>% suppressWarnings() %>% suppressMessages()
            saveRDS(ses.mntd.obj, file = outfile)
            unlink(tmppath, recursive = T)
        }else{
            show_exis_war(outfile)
            ses.mntd.obj <- readRDS(file = outfile)
            show_load_msg(outfile)
        }
        ses.mntd.obj %<>% sort_samples(met = met)
        dataset$asb$alpha.phylo <- ses.mntd.obj
        show_stat_msg('asb$alpha.phylo')
        return(dataset)
    }

    # `beta.phylo` null model
    if (type == 'beta.phylo'){
        if (!outfile %>% file.exists){
            show_comm_msg("calculating community assembly indices[beta.phylo], please wait....")
            betaNTI.obj <- iCAMP::qpen(
                comm = t(mat),
                tree = phy,
                ab.weight = abundance.weighted,
                rand.time = runs,
                sig.bNTI = sig.bNTI,
                sig.rc = sig.rc,
                nworker = nworker,
                memory.G = memory.G,
                project = "beta.phylo",
                wd = tmppath,
                output.detail = TRUE
            ) %>% suppressWarnings() %>% suppressMessages()
            show_comm_msg("saving results to files...")
            saveRDS(betaNTI.obj, file = outfile)
            unlink(tmppath, recursive = T)
        }else{
            show_exis_war(outfile)
            betaNTI.obj <- readRDS(file = outfile)
            show_load_msg(outfile)
        }
        show_comm_msg("extracting distance matrix...")
        mntd.matirx.obs <- betaNTI.obj$bMNTD
        bc.matirx.obs <- betaNTI.obj$BC
        rc.matirx <- nti.matirx <- matrix(nrow = nrow(mntd.matirx.obs),
                                          ncol = ncol(mntd.matirx.obs), data = 0)
        rownames(rc.matirx) <- rownames(nti.matirx) <- rownames(mntd.matirx.obs)
        colnames(rc.matirx) <- colnames(nti.matirx) <- colnames(mntd.matirx.obs)
        null.rst <- betaNTI.obj$result
        for (i in null.rst %>% nrow %>% seq){
            sample1 <- null.rst[i,2] %>% as.character()
            sample2 <- null.rst[i,1] %>% as.character()
            nti <- null.rst[i,5] %>% as.numeric()
            rc <-  null.rst[i,6] %>% as.numeric()
            nti.matirx[which(rownames(nti.matirx) == sample1), which(colnames(nti.matirx) == sample2)] <- nti
            nti.matirx[which(rownames(nti.matirx) == sample2), which(colnames(nti.matirx) == sample1)] <- nti
            rc.matirx [which(rownames(rc.matirx)  == sample1), which(colnames(rc.matirx) == sample2)]  <- rc
            rc.matirx [which(rownames(rc.matirx)  == sample2), which(colnames(rc.matirx) == sample1)]  <- rc
        }

        # integrate results
        res.dat <- list(
            raw = betaNTI.obj,
            dis = list(b.mntd = sort_samples(dat = mntd.matirx.obs, met = met, is.matrix = T),
                       bc = sort_samples(dat = bc.matirx.obs, met = met, is.matrix = T),
                       rc = sort_samples(dat = rc.matirx, met = met, is.matrix = T),
                       b.nti = sort_samples(dat = nti.matirx, met = met, is.matrix = T)))
        dataset$asb$beta.phylo <- res.dat
        show_stat_msg('asb$beta.phylo')
        return(dataset)
    }
}
