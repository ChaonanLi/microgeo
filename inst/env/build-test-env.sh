# -----------------------------------------------------------------------------------------------------------------------
# Copyright (c) 2023, microgeo/Chaonan Li (licn@mtc.edu.cn).                                                            #
# The microgeo is distributed under the terms of the Modified BSD License.                                              #
# Full license is avaliable in the file LICENSE, distributed with this package.                                         #
# -----------------------------------------------------------------------------------------------------------------------
# Build testing environments
# run `build-test-env.py` in `./env`
docker run -itd --rm -p 10000:8888 -v $PWD/tutorial:/workspace/tutorial \
    -v $PWD/images:/workspace/images \
    -v $PWD/test:/workspace/test \
    --name ubuntu-18.04-R4.1.0 jupyterlab:ubuntu-18.04-R4.1.0














# Build testing environments
# We use four ubuntu versions as the testing environments
docker build -t jupyter-lab:ubuntu-18.04-R3.6.3-1bionic .




docker build -t jupyter-lab:ubuntu-20.04 .
docker build -t jupyter-lab:ubuntu-22.04 .
docker build -t jupyter-lab:ubuntu-23.04 .

# Run:
docker run -itd --rm -p 18040:8888 -v $PWD/tutorial:/workspace/tutorial -v $PWD/images:/workspace/images --name ubuntu-18.04-R3.6.3-1bionic jupyter-lab:ubuntu-18.04-R3.6.3-1bionic
docker run -itd --rm -p 20040:8888 -v $PWD/tutorial:/workspace/tutorial -v $PWD/images:/workspace/images --name jupyter-lab-ubuntu-20.04 jupyter-lab:ubuntu-20.04
docker run -itd --rm -p 22040:8888 -v $PWD/tutorial:/workspace/tutorial -v $PWD/images:/workspace/images --name jupyter-lab-ubuntu-22.04 jupyter-lab:ubuntu-22.04
docker run -itd --rm -p 23040:8888 -v $PWD/tutorial:/workspace/tutorial -v $PWD/images:/workspace/images --name jupyter-lab-ubuntu-23.04 jupyter-lab:ubuntu-23.04








apt update && apt install -y build-essential gfortran wget libreadline-dev libcurl4-openssl-dev libpcre2-dev libbz2-dev liblzma-dev libpcre2-dev libcairo2-dev libxt-dev
wget https://cran.r-project.org/src/base/R-4/R-4.0.0.tar.gz && tar -xf R-4.0.0.tar.gz && cd R-4.0.0
./configure --prefix=/usr/local --enable-R-shlib && make -j && make install



R --version



