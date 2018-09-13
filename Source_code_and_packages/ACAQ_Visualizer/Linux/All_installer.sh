#!/bin/bash

sudo apt-get install r-base

wget -q https://cran.rstudio.com/src/contrib/shiny_1.0.5.tar.gz
wget -q https://cran.rstudio.com/src/contrib/shinyBS_0.61.tar.gz
wget -q https://cran.rstudio.com/src/contrib/shinyFiles_0.6.2.tar.gz
wget -q https://cran.rstudio.com/src/contrib/ggplot2_2.2.1.tar.gz
wget -q https://cran.rstudio.com/src/contrib/gridExtra_2.3.tar.gz
wget -q https://cran.rstudio.com/src/contrib/corrplot_0.84.tar.gz
wget -q https://cran.rstudio.com/src/contrib/DT_0.4.tar.gz

sudo R CMD INSTALL shiny_1.0.5.tar.gz
sudo R CMD INSTALL shinyBS_0.61.tar.gz
sudo R CMD INSTALL shinyFiles_0.6.2.tar.gz
sudo R CMD INSTALL ggplot2_2.2.1.tar.gz
sudo R CMD INSTALL gridExtra_2.3.tar.gz
sudo R CMD INSTALL corrplot_0.84.tar.gz
sudo R CMD INSTALL DT_0.4.tar.gz

rm shiny_1.0.5.tar.gz shinyBS_0.61.tar.gz shinyFiles_0.6.2.tar.gz ggplot2_2.2.1.tar.gz gridExtra_2.3.tar.gz corrplot_0.84.tar.gz DT_0.4.tar.gz

