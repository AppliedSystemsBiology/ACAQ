cat(".Rprofile: Setting R repository")
r = getOption("repos") # hard code the repo for CRAN
r["CRAN"] = "http://cloud.r-project.org/"
options(repos = r)
rm(r)

install.packages("shiny")
install.packages("ggplot2")
install.packages("gridExtra")
install.packages("corrplot")
install.packages("DT")
install.packages("shinyFiles")
install.packages("shinyBS")


