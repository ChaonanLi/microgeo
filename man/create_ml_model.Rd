% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trait-ml-model.R
\name{create_ml_model}
\alias{create_ml_model}
\title{Create a machine learning (ML) model}
\usage{
create_ml_model(
  y.data,
  x.data,
  var,
  method = "rf",
  type = c("regression", "classification"),
  testing.ratio = 0.2,
  threads = 60,
  remove.outlier = TRUE,
  ...
)
}
\arguments{
\item{y.data}{A `data.frame` containing variables to be used for ML modeling. For example, a `data.frame` of microbial
alpha-diversity indices (e.g., Shannon-Wiener diversity index).}

\item{x.data}{A `data.frame` containing the prediction variables, e.g., a `data.frame` of 19 bioclimatic variables.}

\item{var}{Variable (in <y.data>) to be used for ML modeling. For example, the `Shannon-Wiener` index if applicable.}

\item{method}{ML method. Please see `http://topepo.github.io/caret/train-models-by-tag.html` for more options. Default
is `rf`, which means a `Random Forest` model.}

\item{type}{Type of ML model. Select one from `regression` and `classification`. Default is `regression`.}

\item{testing.ratio}{How many samples would be applied for testing? Default is `0.2`.}

\item{threads}{Threads used for model training (logical CPUs). Default is `60`.}

\item{remove.outlier}{Whether to remove those outliers before creating a machine learning model? If it is `TRUE`, only
those values (<y.data>) ranging from `lower quartile - 1.5 × interquartile` to `upper quartile + 1.5 × interquartile`
would be applied for model training. It only works when `type = 'regression'` Default is `TRUE`.}

\item{...}{Parameters parsed by \code{caret::train()}.}
}
\value{
An R `list` with the following components:
\describe{
  \item{\code{object$model}}{ML model returned by \code{caret::train()}.}
  \item{\code{object$model.method}}{ML method.}
  \item{\code{object$model.type}}{ML type (`regression` or `classification`).}
  \item{\code{object$tests.dat}}{A `data.frame` of testing data for ML model.}
  \item{\code{object$train.dat}}{A `data.frame` of training data for ML model.}
  \item{\code{object$var.name}}{Name of variable used for ML modeling.}
}
}
\description{
This function is implemented to create a machine learning model through using \code{caret::train()}.
}
\examples{
# Create a microgeo dataset
data(qtp)
map <- read_aliyun_map(adcode = c(540000, 630000, 510000))
dataset.dts <- create_dataset(mat = qtp$asv, ant = qtp$tax, met = qtp$met, map = map,
                              phy = qtp$tre, env = qtp$env, lon = "longitude", lat = "latitude")
dataset.dts \%>\% show_dataset()

# Download aridity index and 19 historical bioclimatic variables
dataset.dts \%<>\% get_ai(out.dir = "test/microgeo_data")
dataset.dts \%<>\% get_his_bioc(res = 2.5, out.dir = "test/microgeo_data")
dataset.dts \%<>\% extract_data_from_spatraster(type = 'his')
dataset.dts \%>\% show_dataset()

# Tidy up microgeo dataset
dataset.dts \%<>\% rarefy_count_table()
dataset.dts \%<>\% tidy_dataset()
dataset.dts \%>\% show_dataset()

# Calculate the relative abundance and the ecological markers at `Family` level
dataset.dts \%<>\% calc_rel_abund()
dataset.dts \%<>\% calc_markers(use.dat = 'spa', use.var = 'AI', annotation.level = 'Family', r.thres = 0.3)
head(dataset.dts$abd$mar$correlation)
dataset.dts \%>\% show_dataset()

# Create a regression random forest model for soil pH
rf.rst.reg <- create_ml_model(y.data = dataset.dts$env,
                              x.data = dataset.dts$spa$tabs[,paste0("Bio", seq(19))], var = 'pH', method = 'rf')
print(rf.rst.reg$model)

# Create a binary classification random forest model for the family of f__A21b
family.bins <- data.frame(row.names = rownames(dataset.dts$abd$mar$abundance),
                          f__A21b = dataset.dts$abd$mar$abundance$f__A21b)
family.bins$f__A21b <- ifelse(family.bins$f__A21b > 0, "presence", "absence") # two classifications
family.bins$f__A21b <- as.factor(family.bins$f__A21b)
rf.rst.cla.bin <- create_ml_model(y.data = family.bins,
                                  x.data = dataset.dts$spa$tabs[,paste0("Bio", seq(19))],
                                  var = 'f__A21b', method = 'rf', type = 'classification')
print(rf.rst.cla.bin$model)

# Create a mutiple class classification random forest model for the family of f__A21b
family.mutiple <- data.frame(row.names = rownames(dataset.dts$abd$mar$abundance),
                             f__A21b = dataset.dts$abd$mar$abundance$f__A21b)
family.mutiple$f__A21b <- cut(family.mutiple$f__A21b, breaks = c(-Inf, 0.05, 0.2, 0.9, Inf),
                              labels = c("H", "A", "S", "Y")) # four classifications
family.mutiple$f__A21b <- as.factor(family.mutiple$f__A21b)
rf.rst.cla.mutiple <- create_ml_model(y.data = family.mutiple,
                                      x.data = dataset.dts$spa$tabs[,paste0("Bio", seq(19))],
                                      var = 'f__A21b', method = 'rf', type = 'classification')
print(rf.rst.cla.mutiple$model)
}
\seealso{
\code{\link[caret:train]{caret::train()}}
\code{\link[microgeo:read_aliyun_map]{microgeo::read_aliyun_map()}}
\code{\link[microgeo:create_dataset]{microgeo::create_dataset()}}
\code{\link[microgeo:show_dataset]{microgeo::show_dataset()}}
\code{\link[microgeo:get_ai]{microgeo::get_ai()}}
\code{\link[microgeo:get_his_bioc]{microgeo::get_his_bioc()}}
\code{\link[microgeo:extract_data_from_spatraster]{microgeo::extract_data_from_spatraster()}}
\code{\link[microgeo:rarefy_count_table]{microgeo::rarefy_count_table()}}
\code{\link[microgeo:tidy_dataset]{microgeo::tidy_dataset()}}
\code{\link[microgeo:calc_rel_abund]{microgeo::calc_rel_abund()}}
\code{\link[microgeo:calc_markers]{microgeo::calc_markers()}}
}
\author{
Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
}
