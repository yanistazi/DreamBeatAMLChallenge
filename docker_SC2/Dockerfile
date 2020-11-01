## Start from this Docker image
FROM ubuntu:18.04

## Install R in Docker image

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && apt-get -y install --no-install-recommends --no-install-suggests ca-certificates software-properties-common gnupg2 gnupg1 \
      && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
      && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' \
      && apt-get install -y r-base 


## Install R packages in Docker image

#RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
#RUN Rscript -e "install.packages('reshape2'); install.packages('caret'); install.packages('xgboost'); install.packages('glmnet'); install.packages('e1071'); install.packages('randomForest'); install.packages('randomForestSRC')"
RUN Rscript -e "install.packages('reshape2'); install.packages('glmnet'); install.packages('randomForestSRC')"

## Copy your files into Docker image
ADD docker /docker
RUN mkdir features input output  
RUN chmod a+x /docker/submit.R

## Make Docker container executable
ENTRYPOINT ["Rscript", "/docker/submit.R"]
