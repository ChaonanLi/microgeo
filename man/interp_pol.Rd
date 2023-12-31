% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trait-interp-model.R
\name{interp_pol}
\alias{interp_pol}
\title{Polynomial fit (2nd order) interpolation}
\usage{
interp_pol(map, met, dat, var, n = 4000, trim.dup = FALSE)
}
\arguments{
\item{map}{A `SpatialPolygonsDataFrame` of geographic map.}

\item{met}{A `data.frame` of sample information. Row names should be sample ids and the column names must be variables
(e.g., `longitude`, `latitude` and `group`). Longitude and latitude must be included in this `data.frame` as we mainly
focus on the spatial patterns of microbial traits.}

\item{dat}{A `data.frame` containing the variables to be interpolated.}

\item{var}{Name of variable to be interpolated.}

\item{n}{Approximate sample size. Default is `4000`. See \code{sp::spsample()} for more details.}

\item{trim.dup}{Whether to randomly remove sampling sites with duplicated coordinates? Default is `FALSE`.}
}
\value{
A `SpatRaster`.
}
\description{
This function is used to perform polynomial fit (2nd order) interpolation.
}
\examples{
# Create a microgeo dataset
data(qtp)
showtext::showtext_auto(enable = TRUE)
map <- read_aliyun_map(adcode = c(540000, 630000, 510000))
dataset.dts <- create_dataset(mat = qtp$asv, ant = qtp$tax, met = qtp$met, map = map,
                              phy = qtp$tre, env = qtp$env, lon = "longitude", lat = "latitude")
dataset.dts \%>\% show_dataset()

# Tidy microgeo dataset, and calculate microbial alpha diversity indices
dataset.dts \%<>\% rarefy_count_table()
dataset.dts \%<>\% tidy_dataset()
dataset.dts \%<>\% calc_alpha_div(measures = c("observed", "shannon"))
dataset.dts \%>\% show_dataset()

# 2nd polynomial fit interpolation for alpha diversity indices
pol.rst.observed <- interp_pol(map = dataset.dts$map, met = dataset.dts$met,
                               dat = dataset.dts$div$alpha, var = 'observed', trim.dup = TRUE)
pol.rst.shannon  <- interp_pol(map = dataset.dts$map, met = dataset.dts$met,
                               dat = dataset.dts$div$alpha, var = 'shannon', trim.dup = TRUE)

# Visualize the interpolated results
dataset.dts$map \%>\% plot_bmap() \%>\%
    add_spatraster(spat.raster = pol.rst.observed) \%>\%
    add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = "X.CENTER", lat.var = "Y.CENTER") \%>\%
    add_scale_bar() \%>\% add_north_arrow() \%>\% add_crs()
dataset.dts$map \%>\% plot_bmap() \%>\%
    add_spatraster(spat.raster = pol.rst.shannon) \%>\%
    add_label(dat = dataset.dts$map@data, lab.var = 'NAME', lon.var = "X.CENTER", lat.var = "Y.CENTER") \%>\%
    add_scale_bar() \%>\% add_north_arrow() \%>\% add_crs()
}
\seealso{
\code{\link[stats:lm]{stats::lm()}}
\code{\link[sp:spsample]{sp::spsample()}}
\code{\link[microgeo:read_aliyun_map]{microgeo::read_aliyun_map()}}
\code{\link[microgeo:create_dataset]{microgeo::create_dataset()}}
\code{\link[microgeo:show_dataset]{microgeo::show_dataset()}}
\code{\link[microgeo:rarefy_count_table]{microgeo::rarefy_count_table()}}
\code{\link[microgeo:tidy_dataset]{microgeo::tidy_dataset()}}
\code{\link[microgeo:calc_alpha_div]{microgeo::calc_alpha_div()}}
\code{\link[microgeo:plot_bmap]{microgeo::plot_bmap()}}
\code{\link[microgeo:add_spatraster]{microgeo::add_spatraster()}}
\code{\link[microgeo:add_label]{microgeo::add_label()}}
\code{\link[microgeo:add_north_arrow]{microgeo::add_north_arrow()}}
\code{\link[microgeo:add_scale_bar]{microgeo::add_scale_bar()}}
\code{\link[microgeo:add_crs]{microgeo::add_crs()}}
}
\author{
Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
}
