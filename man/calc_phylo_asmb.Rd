% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc-biogeo-trait.R
\name{calc_phylo_asmb}
\alias{calc_phylo_asmb}
\title{Calculate the matrices related to microbial community assembly}
\usage{
calc_phylo_asmb(
  dataset,
  type = c("alpha.phylo", "beta.phylo"),
  model = c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool",
    "independentswap", "trialswap"),
  runs = 999,
  nworker = 2,
  memory.G = 10,
  sig.bNTI = 2,
  sig.rc = 0.95,
  abundance.weighted = TRUE,
  out.dir = "calc_comm_asmb.rst"
)
}
\arguments{
\item{dataset}{A microgeo dataset with the class of `MicrogeoDataset`.}

\item{type}{Which type of matrices would be used? Select one from `'alpha.phylo', 'beta.phylo'`. A `alpha.phylo` means
alpha-type phylogenetic null model (e.g., ses.MNTD), and a `beta.phylo` means beta-type phylogenetic null model (e.g.,
betaNTI). Default is `alpha.phylo`.}

\item{model}{Null model used to calculate `alpha.phylo` indices. Default is `taxa.labels`.}

\item{runs}{How many runs would be used for calculation? Defaults is `999`.}

\item{nworker}{How many threads would be used for calculation? Defaults is `2`.}

\item{memory.G}{How many memory would be used for calculation? Defaults is `10` (GB).}

\item{sig.bNTI}{Cutoffs for betaNTI. Required if the <type> is `beta.phylo`. Default is `2`.}

\item{sig.rc}{Cutoffs for RC. Required if the <type> is `beta.phylo`. Default is `0.95`.}

\item{abundance.weighted}{Use abundance weighted? Defaults is `TRUE`.}

\item{out.dir}{A directory saving output files. Default is `calc_comm_asmb.rst`.}
}
\value{
A `MicrogeoDataset` class with the following components:
\describe{
  \item{\code{object$mat}}{A `data.frame` of ASV/gene abundance.}
  \item{\code{object$ant}}{A `data.frame` of ASV/gene anootations.}
  \item{\code{object$met}}{A `data.frame` of sample information.}
  \item{\code{object$map}}{A `SpatialPolygonsDataFrame` of map.}
  \item{\code{object$phy}}{A phylogenetic tree with `newick` format if applicable.}
  \item{\code{object$env}}{A `data.frame` of measured environmental properties if applicable.}
  \item{\code{object$asb}}{Matrices related to microbial community assembly.}
  \item{\code{object$*}}{Spatial and biogeographic traits if applicable.}
}
}
\description{
This function is used to calculate matrices related to microbial community assembly based on null model.
}
\examples{
# Create a microgeo dataset
data(qtp)
map <- read_aliyun_map(adcode = c(540000, 630000, 510000))
dataset.dts <- create_dataset(mat = qtp$asv, ant = qtp$tax, met = qtp$met, map = map,
                              phy = qtp$tre, env = qtp$env, lon = "longitude", lat = "latitude")
dataset.dts \%>\% show_dataset()

# Tidy up microgeo dataset
dataset.dts \%<>\% rarefy_count_table()
dataset.dts \%<>\% tidy_dataset()
dataset.dts \%>\% show_dataset()

# Calculate `alpha.phylo` null model
# runs = 9 just for an example
# runs = 999 would be better
dataset.dts \%<>\% calc_phylo_asmb(type = 'alpha.phylo', runs = 9, out.dir = 'test/calc_comm_asmb.rst')
dataset.dts \%>\%  show_dataset()
head(dataset.dts$asb$alpha.phylo)

# Calculate `beta.phylo` null model
# runs = 9 just for an example
# runs = 999 would be better
dataset.dts \%<>\% calc_phylo_asmb(type = 'beta.phylo', runs = 9, out.dir = 'test/calc_comm_asmb.rst')
dataset.dts \%>\%  show_dataset()
names(dataset.dts$asb$beta.phylo$dis) # distance matrices
head(dataset.dts$asb$beta.phylo$raw$result) # raw result of `iCAMP::qpen()`
}
\seealso{
\code{\link[iCAMP:pdist.big]{iCAMP::pdist.big()}}
\code{\link[iCAMP:qpen]{iCAMP::qpen()}}
\code{\link[picante:ses.mntd]{picante::ses.mntd()}}
\code{\link[microgeo:read_aliyun_map]{microgeo::read_aliyun_map()}}
\code{\link[microgeo:create_dataset]{microgeo::create_dataset()}}
\code{\link[microgeo:show_dataset]{microgeo::show_dataset()}}
\code{\link[microgeo:rarefy_count_table]{microgeo::rarefy_count_table()}}
\code{\link[microgeo:tidy_dataset]{microgeo::tidy_dataset()}}
}
\author{
Li Chaonan (Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University)
}
