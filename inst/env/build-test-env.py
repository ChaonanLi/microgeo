#!/usr/bin/python3
# -----------------------------------------------------------------------------------------------------------------------
# Copyright (c) 2023, microgeo/Chaonan Li (licn@mtc.edu.cn).                                                            #
# The microgeo is distributed under the terms of the Modified BSD License.                                              #
# Full license is avaliable in the file LICENSE, distributed with this package.                                         #
# -----------------------------------------------------------------------------------------------------------------------
# Build testing environment for microgeo R package 

import os

def create_directories(base_path):
    ubuntu_versions = ['18.04', '20.04', '22.04', '23.04']
    r_versions = ['4.1.0', '4.1.1', '4.1.2', '4.1.3', '4.2.0', '4.2.1', '4.2.2', '4.2.3', '4.3.0', '4.3.1', '4.3.2']
    for ubuntu_version in ubuntu_versions:
        for r_version in r_versions:
            directory_path = os.path.join(base_path, f'ubuntu-{ubuntu_version}', f'R{r_version}')
            os.makedirs(directory_path, exist_ok=True)
            write_dockerfile(directory_path, ubuntu_version, r_version)
            build_docker_images(base_path)

def write_dockerfile(directory, ubuntu_version, r_version):
    r_version_prefix = r_version.split(".")[0]
    if ubuntu_version != '23.04':
        dockerfile_content = f"""FROM ubuntu:{ubuntu_version}
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y python3 python3-pip \\
    libffi-dev libssl-dev libzmq3-dev build-essential gfortran \\
    wget libreadline-dev libcurl4-openssl-dev libpcre2-dev libbz2-dev libv8-dev \\
    liblzma-dev libcairo2-dev libxt-dev libgdal-dev libudunits2-dev uuid-dev libharfbuzz-dev libfribidi-dev && apt-get clean && \\
    wget https://mirrors.pku.edu.cn/CRAN/src/base/R-{r_version_prefix}/R-{r_version}.tar.gz && tar -xf R-{r_version}.tar.gz && \\
    cd R-{r_version} && ./configure --prefix=/usr/local --enable-R-shlib && make -j && make install && \\
    pip3 install --upgrade setuptools_scm -i https://mirrors.aliyun.com/pypi/simple/ && \\
    pip3 install jupyterlab -i https://mirrors.aliyun.com/pypi/simple/ 
COPY ./install-packages.r /install-packages.r
RUN Rscript /install-packages.r
RUN echo "options(bitmapType='cairo')" > ~/.Rprofile
WORKDIR /workspace
EXPOSE 8888
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--no-browser", "--allow-root", "--NotebookApp.token=''"]
"""
    else:
        dockerfile_content = f"""FROM ubuntu:{ubuntu_version}
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y python3 python3-pip pipx \\
    libffi-dev libssl-dev libzmq3-dev build-essential gfortran \\
    wget libreadline-dev libcurl4-openssl-dev libpcre2-dev libbz2-dev libv8-dev \\
    liblzma-dev libcairo2-dev libxt-dev libgdal-dev libudunits2-dev uuid-dev libharfbuzz-dev libfribidi-dev && apt-get clean && \\
    wget https://mirrors.pku.edu.cn/CRAN/src/base/R-{r_version_prefix}/R-{r_version}.tar.gz && tar -xf R-{r_version}.tar.gz && \\
    cd R-{r_version} && ./configure --prefix=/usr/local --enable-R-shlib && make -j && make install && \\
    pipx ensurepath && pipx install jupyterlab -i https://mirrors.aliyun.com/pypi/simple/ 
COPY ./install-packages.r /install-packages.r
RUN Rscript /install-packages.r
RUN echo "options(bitmapType='cairo')" > ~/.Rprofile
WORKDIR /workspace
EXPOSE 8888
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--no-browser", "--allow-root", "--NotebookApp.token=''"]
"""
    with open(os.path.join(directory, 'Dockerfile'), 'w') as f:
        f.write(dockerfile_content)
    os.system(f"cp ./install-packages.r {directory}")

def build_docker_images(base_path):
    for root, dirs, files in os.walk(base_path):
        if 'Dockerfile' in files:
            os.chdir(root)
            imagename = root.split("/")[-2] + "-" + root.split("/")[-1]
            try:
                os.system(f'docker build -t jupyterlab:{imagename} .')
            except Exception as e:
                pass
            os.chdir("../../")
            
if __name__ == "__main__":
    base_path = input("Enter the base path: ")
    create_directories(base_path)
