FROM bicore/10x_multiome_qc:4.2.0

MAINTAINER Anish Chakka anish.chakka@alleninstitute.org

ENV DEBIAN_FRONTEND=noninteractive 

RUN apt-get update && apt-get install -y \
 	python3 \
 	python3-pip \
	bedtools

# install MACS2
RUN python3 -m pip install MACS2

# copy the package
COPY peakCallingPipeline_0.15.tar.gz /tools/peakCallingPipeline_0.15.tar.gz
COPY PeakRankR_0.0.0.9000.tar.gz /tools/PeakRankR_0.0.0.9000.tar.gz

# install dependency packages
RUN install2.r -e -s DescTools

## Install R depends 
RUN R -e 'devtools::install_github("PhanstielLab/bedtoolsr")'
RUN R -e 'install.packages("tidyr")'

RUN R -e 'install.packages("BiocManager", update=FALSE)' 
RUN R -e 'BiocManager::install(c( "BiocParallel" )'

# install the package
RUN R --quiet -e 'install.packages("/tools/PeakRankR_0.0.0.9000.tar.gz", repos=NULL, type="source")'
RUN R --quiet -e 'install.packages("/tools/peakCallingPipeline_0.15.tar.gz", repos=NULL, type="source")'

## Clean up
RUN rm -rf /var/lib/apt/lists/*
RUN rm -rf /tmp/downloaded_packages

## Strip binary installed lybraries from RSPM
## https://github.com/rocker-org/rocker-versioned2/issues/340
RUN strip /usr/local/lib/R/site-library/*/libs/*.so