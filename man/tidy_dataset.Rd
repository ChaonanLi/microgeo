% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/microgeo-dataset.R
\name{tidy_dataset}
\alias{tidy_dataset}
\title{Tidy up a microgeo dataset to ensure the valid order of sample and ASV/gene ids}
\usage{
tidy_dataset(dataset)
}
\arguments{
\item{dataset}{A microgeo dataset with the class of `MicrogeoDataset`.}
}
\value{
A `MicrogeoDataset` class with the following components:
\describe{
  \item{\code{object$mat}}{A `data.frame` of ASV/gene abundance.}
  \item{\code{object$ant}}{A `data.frame` of ASV/gene anootations.}
  \item{\code{object$met}}{A `data.frame` of sample information.}
  \item{\code{object$map}}{A `SpatialPolygonsDataFrame` of map.}
  \item{\code{object$phy}}{A phylogenetic tree with `newick` format (`phylo` class) if applicable.}
  \item{\code{object$env}}{A `data.frame` of measured environmental properties if applicable.}
  \item{\code{object$*}}{Spatial and biogeographic traits if applicable.}
}
}
\description{
This function is designed to tidy up a dataset, which ensures the valid order of sample and ASV/gene ids.
When unsure of the dataset's ID order, always use this function.
}
\examples{
data(qtp)
map <- read_aliyun_map(adcode = c(540000, 630000, 510000))
map \%>\% plot_bmap() \%>\%
    add_label(dat = map@data, lab.var = 'NAME', lon.var = 'X.CENTER', lat.var = 'Y.CENTER') \%>\%
    add_scale_bar() \%>\% add_north_arrow() \%>\% add_crs()
dataset.dts <- create_dataset(mat = qtp$asv, ant = qtp$tax, met = qtp$met, map = map,
                              phy = qtp$tre, env = qtp$env, lon = "longitude", lat = "latitude")
dataset.dts \%<>\% tidy_dataset()
dataset.dts \%>\% show_dataset()
unique(rownames(dataset.dts$mat) == rownames(dataset.dts$ant)) # return `TRUE`
unique(colnames(dataset.dts$mat) == rownames(dataset.dts$met)) # return `TRUE`
unique(colnames(dataset.dts$mat) == rownames(dataset.dts$env)) # return `TRUE`
}
\seealso{
\code{\link[microgeo:read_aliyun_map]{microgeo::read_aliyun_map()}}
\code{\link[microgeo:plot_bmap]{microgeo::plot_bmap()}}
\code{\link[microgeo:add_label]{microgeo::add_label()}}
\code{\link[microgeo:add_scale_bar]{microgeo::add_scale_bar()}}
\code{\link[microgeo:add_north_arrow]{microgeo::add_north_arrow()}}
\code{\link[microgeo:add_crs]{microgeo::add_crs()}}
\code{\link[microgeo:create_dataset]{microgeo::create_dataset()}}
\code{\link[microgeo:show_dataset]{microgeo::show_dataset()}}
\code{\link[ape:is.rooted]{ape::is.rooted()}}
\code{\link[ape:multi2di]{ape::multi2di()}}
\code{\link[ape:drop.tip]{ape::drop.tip()}}
}
\author{
Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
}
