# -----------------------------------------------------------------------------------------------------------------------
# Copyright (c) 2023, microgeo/Chaonan Li (licn@mtc.edu.cn).                                                            #
# The microgeo is distributed under the terms of the Modified BSD License.                                              #
# Full license is avaliable in the file LICENSE, distributed with this package.                                         #
# -----------------------------------------------------------------------------------------------------------------------
# For Ubuntu OS: we need to install three libraries
# Just run `apt update && apt install -y libgdal-dev libudunits2-dev uuid-dev libharfbuzz-dev libfribidi-dev libv8-dev` in terminal
for (pkg in c('IRkernel', 'Cairo', 'remotes')){
    if (!suppressMessages(require(pkg, character.only = TRUE)))
        install.packages(pkg, dependencies = TRUE, repos = "http://cran.rstudio.com/")
}
IRkernel::installspec(user = FALSE)
if (!suppressMessages(require('microgeo', character.only = TRUE)))
    remotes::install_github('ChaonanLi/microgeo')
    source(system.file("scripts", "install-extra-pkgs.R", package = "microgeo"))
