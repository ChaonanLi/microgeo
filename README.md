## microgeo: An R package rapidly visualizing the biogeographic traits of soil microbes onto maps

Created by 

* [Li, Chaonan (李超男)](https://www.researchgate.net/profile/Chaonan-Li-5) / licn@mtc.edu.cn / [Ecological Security and Protection Key Laboratory of Sichuan Province, Mianyang Normal University](https://zdsys.mtc.edu.cn/)

* Li, Chi (刘驰) / liuchi0426@126.com /
[College of Resources and Environment, Fujian Agriculture and Forestry University](https://zhxy.fafu.edu.cn/main.htm)

* Liao, Haijun (廖海君) / liaohj@mtc.edu.cn /
[Engineering Research Center of Chuanxibei RHS Construction at Mianyang Normal University of Sichuan Province](https://rhs.mtc.edu.cn/)

Reviewed by 

* [Li, Xiangzhen (李香真)](https://www.researchgate.net/profile/Xiangzhen-Li-2) / lixz@fafu.edu.cn /
[College of Resources and Environment, Fujian Agriculture and Forestry University](https://zhxy.fafu.edu.cn/main.htm)


### 1. Introduction

With the accumulation of massive microbial sequence data, there are growing interests in the analysis of datasets deposited in public databases (e.g., [NCBI Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra/?term=) and [EMBL-EBI European Nucleotide Archive](https://www.ebi.ac.uk/ena/browser/home)) to explain the biodiversity patterns and species distributions at a large spatial scale. Yet, analyzing these public datasets faces tremendous challenges not only in processing a large amount of sequence data but also in effectively exploring and visualizing microbial traits summaried from sequence data in-depth. A large number of open-source R packages can provide diverse methods for complex statistical analysis and visualization, for example, [phyloseq](https://github.com/joey711/phyloseq), [microbiome](https://github.com/microbiome/microbiome), [microbiomeSeq](http://www.github.com/umerijaz/microbiomeSeq), [ampvis2](https://madsalbertsen.github.io/ampvis2/reference/index.html), [MicrobiomeR](https://github.com/vallenderlab/MicrobiomeR) and [Rhea](https://pubmed.ncbi.nlm.nih.gov/28097056/). Our recently developed [microeco](https://github.com/ChiLiubio/microeco) R package even provides almost all common analysis methods of microbiome, thereby effectively addressed many challenges of statistical analysis. However, these packages tend to provide abstract and statistically-oriented results, and can not intuitively render microbial traits or the results generated by themselves onto maps. Although the [ArcGIS](https://www.arcgis.com/index.html) and the R packages like [ggplot2](https://github.com/tidyverse/ggplot2), [ggspatial](https://github.com/paleolimbot/ggspatial), [terra](https://github.com/rspatial/terra), [sf](https://github.com/r-spatial/sf/) and [sp](https://github.com/edzer/sp) are capable of flexibly rendering microbial traits on maps, they are relatively large and come with a steep learning curve, and the time cost of learning many geographic information system (GIS) tools or R packages only for the purpose of visualizing microbial traits onto a map is not acceptable. Hence, there is a pressing need for a lightweight, flexible, and user-friendly R package to visualize the outcomes derived from aforementioned packages that are focused on statistical analysis, thereby enabling the analysis of microbiome dataset in conjunction with the metadata obtained from geographic maps.

An common approach in analyzing the spatial pattern of microbial traits at a large spatial scale is to compare the measures summaried from sequence data through classifying samples into different groups based on biomes or ecosystems ([Luke R. Thompson et al., 2017](https://www.nature.com/articles/nature24621) and [Salvador Ramírez-Flandes et al., 2018](https://www.pnas.org/doi/10.1073/pnas.1817554116)) or fit these measures to specific environmental/spatial gradients (e.g., soil pH, temperature, longitude and latitude) ([Luke R. Thompson et al., 2017](https://www.nature.com/articles/nature24621) and [Mohammad Bahram et al., 2018](https://www.nature.com/articles/s41586-018-0386-6)). However, such an approach can only yield statistically-oriented results, which are far less intuitive compared to directly displaying these outcomes on a map. Moreover, the sampling sites of the analyzed data may be insufficient to cover the whole study area, even with a grid-based sampling strategy. In this light, the results obtained from the aforementioned method are difficult to accurately describe the spatial patterns of microbial traits, especially at a large spatial scale. Several studies attempt to visualize microbial traits onto maps by using spatial interpolation ([Bin Ma et al., 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5029158/) and [Shuo Jiao et al., 2022](https://onlinelibrary.wiley.com/doi/10.1002/imt2.2)), machine learning ([Shang Wang et al., 2021](https://www.sciencedirect.com/science/article/pii/S2095927321000578) and [Tarciso C. C. Leão et al., 2021](https://onlinelibrary.wiley.com/doi/10.1111/geb.13365)) or other ecological models ([Carlos A. Guerra et al., 2021](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7610617/)). These approaches can estimate microbial traits at the unknown sites by using a model or the known values at nearby sites, thereby describing the regional patterns of microbial traits in the whole study areas. Yet, although these approaches are attractive, applying them to the microbiome data analysis still remains a lot of significant challenges. Besides, analyzing the spatial patterns of microbial traits summarized from public sequence data also poses a significant challenge due to the limited availability of environmental datasets that correspond to microbiome data, and this dearth of paired datasets is particularly prevalent in the majority of currently accessible public repositories. The [geodata](https://github.com/rspatial/geodata) R package provides some method to download spatial, climate and soil properties from public repositories, and the [MODIS](https://modis.gsfc.nasa.gov/) hosts many remote-sense images of plant metrics and land cover. To a certain extent, this can largely mitigate the difficult of limited environmental data, though the accuracy of their absolute values is relatively lower compared to those measured values. However, it is extremely challenging for most R language beginners to effectively integrate these datasets for microbial biogeography analysis, as the whole process involves manipulating diverse datasets.

Here we present the [microgeo](https://github.com/ChaonanLi/microgeo) R package which warps the [ggplot2](https://github.com/tidyverse/ggplot2) and [ggspatial](https://github.com/paleolimbot/ggspatial), and other R packages related to GIS and machine learning such as [gstat](https://github.com/r-spatial/gstat/), [raster](https://github.com/rspatial/raster), [terra](https://github.com/rspatial/terra), [sf](https://github.com/r-spatial/sf/), [sp](https://github.com/edzer/sp) and [caret](https://github.com/topepo/caret/). It permits the microbial biogeographical trait calculation, spatial data collection as well as spatial interpolation, machine learning modeling/prediction and visualization for microbial biogeographical traits and spatial data. Specifically, the [microgeo](https://github.com/ChaonanLi/microgeo) provides the flexible visualization methods for microbial traits onto maps, e.g.,  gridded-based visualization, spatial interpolation and machine learning, and it is not limited to those traits calculated by [microgeo](https://github.com/ChaonanLi/microgeo) package itself. Users can rapidly visualize any microbial traits calculated by other tools onto maps, and also can analyze microbiome dataset in conjunction with the data obtained from [microgeo](https://github.com/ChaonanLi/microgeo) R package.

The [microgeo](https://github.com/ChaonanLi/microgeo) is currently undergoing continuous developments and updates. We welcome any ideas and suggestions. If you find any errors during use, please submit them to [GitHub Issues](https://github.com/ChaonanLi/microgeo/issues).

<img src="images/microgeo-workflow.png" width = "100%" alt="Fig-1" align=center />

### 2. Installation and dependencies

The [microgeo](https://github.com/ChaonanLi/microgeo) requires ![](https://img.shields.io/badge/R->=4.1.0-orange), ![](https://img.shields.io/badge/GDAL-green) and following dependent R packages. For installing R and GDAL, you can refer to the methods described on their official websites.

To install [microgeo](https://github.com/ChaonanLi/microgeo), just type following codes in your R console. The installation may take a few minutes due to the dependencies. Please be patient and wait for the process to complete.

```{.r}
if (!suppressMessages(require('remotes', character.only = TRUE))) 
    install.packages('remotes', dependencies = TRUE, repos = "http://cran.rstudio.com/")
if (!suppressMessages(require('microgeo', character.only = TRUE)))
    # remotes::install_github('ChaonanLi/microgeo') 
    remotes::install_git("https://gitee.com/bioape/microgeo") # Specifically for Mainland Chinese users; 
source(system.file("scripts", "install-extra-pkgs.R", package = "microgeo")) # Install additional R packages required by `microgeo`
```

For the convenience of users to quickly use the microgeo R package and to facilitate the deployment of the microgeo R package on servers, we have built a Docker image based on Ubuntu 22.04 and R 4.3.2. You can quickly set up a microgeo running environment by using the following commands:

```{.shell}
docker pull microgeo-jupyterlab:ubuntu-22.04-R4.3.2

docker run -itd --rm -p 10000:8888 -v $PWD:$PWD --name microgeo-jupyterlab-ubuntu-22.04-R4.32 microgeo-jupyterlab:ubuntu-22.04-R4.3.2
```

Then, you can visit the http://[IP]:10000 or http://127.0.0.1:10000 in your browser to use microgeo R package.

### 3. Citation

If you use [microgeo](https://github.com/ChaonanLi/microgeo) for data processing and publication of a research paper, please cite: https://github.com/ChaonanLi/microgeo.

### 4. Usages and examples

To make it more convenient for users to use [microgeo](https://github.com/ChaonanLi/microgeo), we not only provide detailed examples in the help document section of each function (which can be viewed by using `?<function_name()>`), but also provide a more detailed usage [tutorial](https://github.com/ChaonanLi/microgeo/tutorial/*.ipynb) by using JupyterLab Notebook. If the [microgeo](https://github.com/ChaonanLi/microgeo) package and its dependencies have been successfully installed, users can directly run these examples in this Notebook. If the JupyterLab is not avaliable on your PC or server, users also can copy the codes from the [HTML version of tutorial](https://github.com/ChaonanLi/microgeo/tutorial/*.html) to the R console to run the examples.

It is recommended to use [JupyterLab](https://jupyter.org/install) or [RStudio](https://posit.co/downloads/) as your IDE (Integrated Development Environment). However, if you prefer the black screen with white text, using R's command interface directly is also fine.
