FROM openanalytics/r-base:4.0.3

# system libraries of general use
RUN apt-get update && apt-get install -y \
    sudo \
    build-essential \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libssl1.1 \
    libxml2-dev \
    libproj-dev \ 
    libgdal-dev \
    libharfbuzz-dev \
    libfribidi-dev

# basic shiny functionality
RUN R -e "install.packages(c('shiny'), repos='https://cloud.r-project.org/')"
RUN R -e "install.packages(c('shinyBS', 'htmlwidgets', 'shinycssloaders', 'shinybusy', 'DT', 'stringr' , 'stats', 'plyr', 'MESS', 'dplyr', 'RColorBrewer', 'ggsci', 'ggplot2', 'multcompView', 'gridExtra', 'grid', 'igraph', 'visNetwork', 'pheatmap', 'tidyr', 'pastecs', 'zip', 'eulerr', 'VennDiagram'), repos='https://cloud.r-project.org/')"


# copy the app to the image
RUN mkdir /root/PROMISed
COPY PROMISed /root/PROMISed

COPY Rprofile.site /usr/lib/R/etc/

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/root/PROMISed')"]
