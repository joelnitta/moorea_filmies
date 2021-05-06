FROM rocker/verse:4.0.3

ARG DEBIAN_FRONTEND=noninteractive

###############################################
### Install other dependencies with apt-get ###
###############################################

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
  libpoppler-cpp-dev

####################################
### Install R packages with renv ###
####################################

# Create directory for renv project library
RUN mkdir renv

# Modify Rprofile.site so renv uses /renv for project library, and doesn't use the cache
RUN echo 'Sys.setenv(RENV_PATHS_LIBRARY = "/renv")' >> /usr/local/lib/R/etc/Rprofile.site

# Initialize a 'dummy' project and restore the renv library.
# Since the library path is specified as above, the library will be restored to /renv
RUN mkdir tmp/project

COPY ./renv.lock tmp/project

WORKDIR tmp/project

# Don't use cache (the symlinks won't work from Rstudio server)
RUN Rscript -e 'install.packages("renv"); renv::consent(provided = TRUE); renv::settings$use.cache(FALSE); renv::init(bare = TRUE); renv::restore()'

###########################################
### Install latex packages with tinytex ###
###########################################

COPY install_latex.R .

RUN Rscript install_latex.R

WORKDIR /home/rstudio/
