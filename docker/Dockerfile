FROM ubuntu-22.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y tzdata python3 python3-pip libffi-dev libssl-dev libzmq3-dev build-essential \
    gfortran wget libreadline-dev libcurl4-openssl-dev libpcre2-dev libbz2-dev libv8-dev liblzma-dev libcairo2-dev \
    libxt-dev libgdal-dev libudunits2-dev uuid-dev libharfbuzz-dev libfribidi-dev cmake libgit2-dev && apt-get clean && \
    wget https://mirrors.pku.edu.cn/CRAN/src/base/R-4/R-4.3.2.tar.gz && \
    tar -xf R-4.3.2.tar.gz && cd R-4.3.2 && ./configure --prefix=/usr/local --enable-R-shlib && \
    make -j && make install && pip3 install --upgrade setuptools_scm -i https://mirrors.aliyun.com/pypi/simple/ && \
    pip3 install jupyterlab -i https://mirrors.aliyun.com/pypi/simple/ && \
    Rscript -e "install.packages(c('IRkernel', 'Cairo'), repos='https://cloud.r-project.org/', dependencies = TRUE); \
    IRkernel::installspec(user = FALSE)" && echo "options(bitmapType = 'cairo')" > ~/.Rprofile && Rscript \
    -e "install.packages('fpc', version = '2.2.10', repos = 'http://cran.rstudio.com/'); \
    install.packages('ggdendro', version = '0.1.23', repos = 'http://cran.rstudio.com/'); \
    install.packages('factoextra', version = '1.0.7', repos = 'http://cran.rstudio.com/'); \
    install.packages(c('remotes', 'cowplot'), dependencies = TRUE, repos = 'http://cran.rstudio.com/'); \
    remotes::install_git('https://gitee.com/bioape/microgeo'); \
    install.packages(system.file('pkgs', 'rgdal_1.6-7.tar.gz', package = 'microgeo'), repos = NULL, type = 'source', dependencies = TRUE); \
    install.packages(system.file('pkgs', 'maptools_1.1-8.tar.gz', package = 'microgeo'), repos = NULL, type = 'source', dependencies = TRUE); \
    install.packages(system.file('pkgs', 'MODIS-1.2.11.tar.gz', package = 'microgeo'), repos = NULL, type = 'source', dependencies = TRUE); \
    install.packages(system.file('pkgs', 'rgeos_0.6-4.tar.gz', package = 'microgeo'), repos = NULL, type = 'source', dependencies = TRUE)"
WORKDIR /workspace
EXPOSE 8888
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--no-browser", "--allow-root", "--NotebookApp.token=''"]
