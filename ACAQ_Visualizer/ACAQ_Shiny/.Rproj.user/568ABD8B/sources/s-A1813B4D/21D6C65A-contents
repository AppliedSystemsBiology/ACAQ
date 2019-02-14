# Copyright (C) 2018 Research Group Applied Systems Biology,
# Leibniz Institute for Natural Product Research and Infection Biology – 
# Hans Knöll Institute (HKI)
# All Rights Reserved
# Author: Naim Al-Zaben
# You may use, distribute and modify this code under the terms of the GPL-3 license.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GPL-3 license with this file. 
# If not, please visit https://www.gnu.org/licenses/gpl-3.0.en.html 

library(shiny)
library(shinyBS)
library(shinyjs)
library(shinyFiles)
# Define UI for application that plots ACAQ-v4 plots 
jscode <- "shinyjs.closeWindow = function() { window.close(); }"

shinyUI(fluidPage( 
  useShinyjs(),
  extendShinyjs(text = jscode, functions = c("closeWindow")),
  # Favicon
  list(tags$head(HTML('<link rel="icon", href="logo.png",
                                   type="image/png" />'))),
  # Application title
  titlePanel(title= tags$div(img(src="logo.png", width="150", height="150"), "ACAQ Visualizer"), windowTitle ="ACAQ Visualizer" ),
  
  # Sidebar with functional controls
  sidebarPanel(
    # Input: Select a file ----
    fileInput("file1", "Choose input file",
              multiple = FALSE,
              accept = c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")),
    
    
    # Horizontal line ----
    tags$hr(),
    radioButtons(inputId = "fileType", label = "Select output files", choices = list("png plots", "pdf plots")),
    checkboxInput("checkboxCSV", label = "csv data per plot", value = FALSE),
    shinyDirButton('directory', 'Select download folder', 'Please select a folder'),
    verbatimTextOutput('OutputDir', placeholder = FALSE),
    actionButton("downData", "Create results"),
    tags$hr(),
    fluidRow(column(12, align="center", actionButton("exitApp", "Exit")))
  
  ),
  
  
  # Show generated plots according to their type
  mainPanel("",
    tabsetPanel(
      tabPanel("Scatter Plots", uiOutput("scatterPlot")), 
      tabPanel("Histograms", uiOutput("histogram")), 
      tabPanel("Notched box plots", uiOutput("notchedBox")),
      tabPanel("Correlation matrix", fluidRow(column( 12,align="center",plotOutput("correlationMatrix")))), 
  # Show the csv file with search feature enabled
      tabPanel("Raw data", DT::dataTableOutput("table")),
  # Explain how to use tool (External HTML)
      tabPanel("How to use", fluidRow(column(12, includeHTML("howToUse.html")))),
  # About Us page
      tabPanel("About", fluidRow(column(12, includeHTML("aboutUs.html"))))
    )
  )
)
)
