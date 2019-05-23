# Base Image
FROM r-base:3.6.0

# Metadata
LABEL base.image="who-multi-omics-results:v.20190430"
LABEL version="1"
LABEL software="WHO multi-omics Results"
LABEL software.version="1.0.0"
LABEL description="Generate results combining all different analyses results sucha s mutation, expression, translocation, CNV, etc."
LABEL tags="CNV Survival Mutation Expression Translocation Subgroups Cell of Origin"

# Maintainer
MAINTAINER DaveLab <lab.dave@gmail.com>

# update the OS related packages
RUN apt-get update -y &&\
    apt-get install git -y

# install required dependencies for QCParser
RUN R --vanilla -e 'install.packages(c("optparse", "ggplot2", "jsonlite", "survival", "here"), repos="http://cran.us.r-project.org")'

# clone who-multi-omics-results repo
RUN git clone https://github.com/labdave/who-multi-omics-results.git

# add who-multi-omics-results repo to SYSPATH
ENV PATH /who-multi-omics-results:$PATH

# change the permission of the repo
RUN chmod 777 -R /who-multi-omics-results

# make the meta_analysis_CLI as default script to execute
CMD ["meta_analysis_CLI.R"]
