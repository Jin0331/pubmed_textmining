FROM rocker/tidyverse:4.1.3
MAINTAINER Jinoo <jinoo@wmbio.co>

USER root
WORKDIR /home/rstudio

# package install
RUN R -e 'install.packages(c("httr", "xml2", "rentrez", "RMariaDB", "parallel", "arules", "BiocManager"))'
RUN R -e 'BiocManager::install(c("org.Hs.eg.db", "AnnotationDbi"), ask = FALSE, force = TRUE)'

# DIR copy
ADD RAW_DATA RAW_DATA
ADD dict dict
ADD R R

CMD [ "echo", "Run Textmining" ]