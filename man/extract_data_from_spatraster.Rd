% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial-tools.R
\name{extract_data_from_spatraster}
\alias{extract_data_from_spatraster}
\title{Extract spatial data for each sample}
\usage{
extract_data_from_spatraster(
  dataset,
  method = c("simple", "bilinear"),
  type = c("both", "his", "cla"),
  remove.na = TRUE,
  ...
)
}
\arguments{
\item{dataset}{A microgeo dataset with the class of `MicrogeoDataset`.}

\item{method}{Method for extracting values with points (`simple` or `bilinear`). With the `simple` values for the cell
a point falls in are returned. With a `bilinear` the returned values are interpolated from the values of the 4 nearest
raster cells. Default is `simple`.}

\item{type}{Which type of data would be extracted? Select one from `both`, `his` and `cla`. Default is `both`.}

\item{remove.na}{Remove rows with `NA` values. Default is `TRUE`.}

\item{...}{Parameters parsed by \code{terra::extract()}.
\describe{
  \item{\code{object$mat}}{A data.frame of ASV/gene abundance.}
  \item{\code{object$ant}}{A data.frame of ASV/gene anootations.}
  \item{\code{object$met}}{A data.frame of sample information.}
  \item{\code{object$map}}{A `SpatialPolygonsDataFrame` of map.}
  \item{\code{object$phy}}{A phylogenetic tree with newick format if applicable.}
  \item{\code{object$env}}{A data.frame of measured environmental properties if applicable.}
  \item{\code{object$spa$rast}}{`SpatRaster` of all spatial data.}
  \item{\code{object$spa$unit}}{A `data.frame` of unit for historically spatial variables if applicable.}
  \item{\code{object$spa$tabs}}{A `data.frame` of historically spatial variables for each sample.}
  \item{\code{object$*}}{Other spatial and biogeographic traits if applicable.}
}}
}
\description{
This function is used to extract spatial data for each sample based on longitudes and latitudes.
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

# Extract spatial data for each sample
dataset.dts \%<>\% extract_data_from_spatraster(type = 'his')
dataset.dts \%<>\% tidy_dataset()
dataset.dts \%>\% show_dataset()
head(dataset.dts$spa$tabs)
}
\seealso{
\code{\link[terra:extract]{terra::extract()}}
\code{\link[microgeo:read_aliyun_map]{microgeo::read_aliyun_map()}}
\code{\link[microgeo:create_dataset]{microgeo::create_dataset()}}
\code{\link[microgeo:show_dataset]{microgeo::show_dataset()}}
\code{\link[microgeo:get_his_bioc]{microgeo::get_his_bioc()}}
\code{\link[microgeo:extract_data_from_spatraster]{microgeo::extract_data_from_spatraster()}}
\code{\link[microgeo:tidy_dataset]{microgeo::tidy_dataset()}}
}
\author{
Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
}
