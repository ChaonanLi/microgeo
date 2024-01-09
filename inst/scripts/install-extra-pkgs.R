# -----------------------------------------------------------------------------------------------------------------------
# Copyright (c) 2023, microgeo/Chaonan Li (licn@mtc.edu.cn).                                                            #
# The microgeo is distributed under the terms of the Modified BSD License.                                              #
# Full license is avaliable in the file LICENSE, distributed with this package.                                         #
# -----------------------------------------------------------------------------------------------------------------------
# Install several required R packages after the `microgeo` has been successfully installed
# These packages are not avaliable at CRAN.
if (!suppressMessages(require('maptools', character.only = TRUE)))
    install.packages(system.file("pkgs", "maptools_1.1-8.tar.gz", package = "microgeo"),
                     repos = NULL, type = "source", dependencies = TRUE)
if (!suppressMessages(require('MODIS', character.only = TRUE)))
    install.packages(system.file("pkgs", "MODIS-1.2.11.tar.gz", package = "microgeo"),
                     repos = NULL, type = "source", dependencies = TRUE)
if (!suppressMessages(require('rgeos', character.only = TRUE)))
    install.packages(system.file("pkgs", "rgeos_0.6-4.tar.gz", package = "microgeo"),
                     repos = NULL, type = "source", dependencies = TRUE)
