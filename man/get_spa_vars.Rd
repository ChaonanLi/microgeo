% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial-tools.R
\name{get_spa_vars}
\alias{get_spa_vars}
\title{Get the variable names of all historically spatial data in a microgeo dataset}
\usage{
get_spa_vars(dataset)
}
\arguments{
\item{dataset}{A microgeo dataset with the class of `MicrogeoDataset`.}
}
\value{
A `MicrogeoDataset` class with the following components:
\describe{
  \item{\code{object$mat}}{A `data.frame` of rarefied ASV/gene abundance.}
  \item{\code{object$ant}}{A `data.frame` of ASV/gene anootations.}
  \item{\code{object$met}}{A `data.frame` of sample information.}
  \item{\code{object$map}}{A `SpatialPolygonsDataFrame` of map.}
  \item{\code{object$phy}}{A phylogenetic tree with `newick` format (`phylo` class) if applicable.}
  \item{\code{object$env}}{A `data.frame` of measured environmental properties if applicable.}
  \item{\code{object$spa$rast}}{`SpatRaster` of all spatial data.}
  \item{\code{object$spa$unit}}{A `data.frame` of unit for historically spatial variables.}
  \item{\code{object$*}}{Other spatial and biogeographic traits if applicable.}
}
}
\description{
This function is used to get the variable names of all historically spatial data in a microgeo dataset.
}
\examples{
# Create a microgeo dataset
data(qtp)
map <- read_aliyun_map(adcode = c(540000, 630000, 510000))
dataset.dts <- create_dataset(mat = qtp$asv, ant = qtp$tax, met = qtp$met, map = map,
                              phy = qtp$tre, env = qtp$env, lon = "longitude", lat = "latitude")
dataset.dts \%>\% show_dataset()
dataset.dts \%<>\% get_his_bioc(res = 2.5, out.dir = "test/microgeo_data")
dataset.dts \%>\% show_dataset()

# Get the variable names of all historically spatial data in a microgeo dataset
dataset.dts \%<>\% get_spa_vars()
head(dataset.dts$spa$unit)
}
\seealso{
\code{\link[microgeo:read_aliyun_map]{microgeo::read_aliyun_map()}}
\code{\link[microgeo:create_dataset]{microgeo::create_dataset()}}
\code{\link[microgeo:show_dataset]{microgeo::show_dataset()}}
\code{\link[microgeo:get_his_bioc]{microgeo::get_his_bioc()}}
}
\author{
Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
}
