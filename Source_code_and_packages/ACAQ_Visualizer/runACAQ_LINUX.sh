#!/bin/bash
Rscript -e "library(shiny);library(ggplot2);library(gridExtra);library(corrplot);library(DT);library(shinyFiles);library(shinyBS);shiny::runApp('./ACAQ_Shiny', launch.browser = TRUE)"

