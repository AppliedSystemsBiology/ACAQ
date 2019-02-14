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
library(ggplot2)
library(gridExtra)
library(corrplot)
library(DT)
library(shinyFiles)
library(stringi)
library(shinyjs)


cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

#Notification
progressNotification <- function(message, id = NULL) {
  
  if(is.null(id)) {
    id <- stringi::stri_rand_strings(1, 16)
  }
  
  showNotification(
    ui = tags$span(icon("circle-o-notch", class = "fa-spin"), message),
    id = id,
    closeButton = F,
    type = "default",
    duration = NULL
  )
  
  return(id)
}


# Define server logic required to generate and plot a random distribution
shinyServer(function(input, output, session) {
  # Expression that generates a plot of the distribution. The expression
  # is wrapped in a call to renderPlot to indicate that:
  #
  #  1) It is "reactive" and therefore should be automatically 
  #     re-executed when inputs change
  #  2) Its output type is a plot 
  #
  #observeEvent(input$file1, {js$reset()}) 
  output$scatterPlot <- renderUI({
    req(input$file1)
    #pathOut <- input$file1$datapath
    AllTable <- read.delim(header = TRUE,input$file1$datapath)
    #AllTable <- read.delim(header = TRUE,"/home/nzaben/Projects/Mohammad_Zoltan/test_large/example_data.csv")

    even_indexes <- rev(seq(2,ncol(AllTable),2))
    
    for (i in even_indexes) {
      AllTable[i] = NULL
    }
    Names <- colnames(AllTable)
    #for(t in 1:length(Names)){
    #  Names[t] <- gsub('\\.', '_', Names[t])
    #}
    j<-2
    k<-3
    max_plots <- 0
    plot_output_list <- list()
    temp <- ncol(AllTable)-1
    
    for (j in 2:temp) {
      temp2 <- j+1
      for (k in temp2:ncol(AllTable)) {
        max_plots <- max_plots + 1
        plotname <- paste0(Names[j],"_",Names[k])
        plot_output_list<- c(plot_output_list,list(plotOutput(plotname, height = 600, width = 800)))
        #dev.off()
      }
    }
    
    do.call(tagList,plot_output_list)
    
  })
  
  
  
    observe({
      
    #reactive({
    req(input$file1)
    
    notification.id <- progressNotification(paste0("Reading table..."))
    AllTable <- read.delim(header = TRUE,input$file1$datapath)
    
    even_indexes <- rev(seq(2,ncol(AllTable),2))
    
    for (i in even_indexes) {
      AllTable[i] = NULL
    }
    
    AllTable[[1]] <- as.character(AllTable[[1]]) 
    j = 1
    for (i in 1:nrow(AllTable)) {
      AllTable[[1]][i] <- gsub('\\','/',as.character(AllTable[[j]][i]), fixed=TRUE)
    }
    FileDirSep <- "/"
    
    
    pos = gregexpr(FileDirSep, AllTable[[1]][j])
    pos[[1]] <- as.numeric(pos[[1]])
    ExName <- substring(AllTable[[1]][j],1,pos[[1]][length(pos[[1]])-1])
    
    Names <- colnames(AllTable)
    j<-2
    k<-3
    temp <- ncol(AllTable)-1
    
    for (j in 2:temp) {
      temp2 <- j+1
      for (k in temp2:ncol(AllTable)) {
        local({
          
          # Need local so that each item gets its own number. Without it, the value
          # of i in the renderPlot() will be the same across all instances, because
          # of when the expression is evaluated.
          plotname <- paste0(Names[j],"_",Names[k])
          
          exp1 <- paste0("abc <- ggplot(AllTable, aes(x=",Names[j],", y=",Names[k],")) + ggtitle(paste0(ExName,\"\n\",Names[j],\"_\",Names[k]))+ theme(plot.title = element_text(hjust = 0.5),axis.text=element_text(size=16,face=\"bold\"), axis.title=element_text(size=18,face=\"bold\"))  + geom_point(size = 3,shape=19, color=\"blue\")+  labs(x=Names[j], y=Names[k])")
          eval(parse(text=exp1))
          abc <- abc + theme(axis.text=element_text(size=12),
                             axis.title=element_text(size=14,face="bold"))
          output[[plotname]] <- renderPlot({
            #StartNote
            notification.id <- progressNotification(paste0("Loading ",plotname," scatterplot"))
            print(abc)
            #EndNote
            removeNotification(notification.id) 
          })
          
          
        })
      }
    }
    #EndNote
    removeNotification(notification.id) 
  })
    
    
    
    output$histogram <- renderUI({
      
      req(input$file1)
      
      AllTable <- read.delim(header = TRUE,input$file1$datapath)
      even_indexes <- rev(seq(2,ncol(AllTable),2))
      
      for (i in even_indexes) {
        AllTable[i] = NULL
      }
      Names <- colnames(AllTable)
      
      
      
      plot_output_list_hist <- list()
      temp <-ncol(AllTable)
      for (j in 2:temp) {
        plotnamehist <- paste0("Histogram_",Names[j])
        plot_output_list_hist<- c(plot_output_list_hist,list(plotOutput(plotnamehist, height = 600, width = 800)))
      }
      do.call(tagList,plot_output_list_hist)
    })
    
    observe({
      req(input$file1)
      
      AllTable <- read.delim(header = TRUE,input$file1$datapath)
      even_indexes <- rev(seq(2,ncol(AllTable),2))
      
      for (i in even_indexes) {
        AllTable[i] = NULL
      }
      Names <- colnames(AllTable)
      
      AllTable[[1]] <- as.character(AllTable[[1]]) 
      j = 1
      for (i in 1:nrow(AllTable)) {
        AllTable[[1]][i] <- gsub('\\','/',as.character(AllTable[[j]][i]), fixed=TRUE)
      }
      FileDirSep <- "/"
      
      
      pos = gregexpr(FileDirSep, AllTable[[1]][j])
      pos[[1]] <- as.numeric(pos[[1]])
      ExName <- substring(AllTable[[1]][j],1,pos[[1]][length(pos[[1]])-1])
      
      j = 2
      
      temp <-ncol(AllTable)
      for (j in 2:temp) {
        local({
          plotnamehist <- paste0("Histogram_",Names[j])
           # histabc <- ggplot(AllTable, aes(AllTable[,j])) +
           #   geom_histogram( color="darkblue", fill="lightblue")  +
           #   ggtitle(paste0(ExName,"\n","Histogram_",Names[j]))+
           #   labs(x=Names[j], y="Count") +
           #   theme(plot.title = element_text(size = 10,hjust = 0.5))
           # histabc <- histabc + theme(axis.text=element_text(size=12),
           #                    axis.title=element_text(size=14,face="bold"))
          histabc <- hist(AllTable[,j],
                          main=paste0(ExName,"\n","Histogram_",Names[j]), 
                          xlab=Names[j], 
                          border="black", 
                          col="light blue",
                          breaks = 30,
                          freq=TRUE,
                          yaxt ='n',
                          cex.axis=1.5,
                          cex.lab=1.5)                    #histabc <- hist(AllTable[,j])
          histabc$density <- histabc$counts/sum(histabc$counts)*100
          
          
          #histabc <-  ggplot(AllTable, aes(AllTable[,j])) +
          #  geom_histogram(aes(y=..count../sum(..count..)) , color="darkblue", fill="lightblue")  +
          #  ggtitle(paste0(ExName,"\n","Histogram_",Names[j])) +
          #  labs(x=Names[j], y="frequency") +
          #  theme(plot.title = element_text(size = 18,hjust = 0.5)) +
          #histabc <- histabc + theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))
          
          recordedPlotHist = recordPlot(histabc)
           #dev.off()
          output[[plotnamehist]] <- renderPlot({
            #StartNote
            notification.id <- progressNotification(paste0("Loading ",plotnamehist," histogram"))
            print(recordedPlotHist)
            
            #EndNote
            removeNotification(notification.id) 
            
          })
        })
        
      }
    })
  
    
    output$notchedBox <- renderUI({
      
      req(input$file1)
      
      AllTable <- read.delim(header = TRUE,input$file1$datapath)
      
      even_indexes <- rev(seq(2,ncol(AllTable),2))
      
      for (i in even_indexes) {
        AllTable[i] = NULL
      }
      Names <- colnames(AllTable)
      
      AllTable[[1]] <- as.character(AllTable[[1]])
      j = 1
      for (i in 1:nrow(AllTable)) {
        AllTable[[1]][i] <- gsub('\\','/',as.character(AllTable[[j]][i]), fixed=TRUE)
      }
      FileDirSep <- "/"
      
      pos = gregexpr(FileDirSep, AllTable[[1]][j])
      pos[[1]] <- as.numeric(pos[[1]])
      ExName <- substring(AllTable[[1]][j],1,pos[[1]][length(pos[[1]])-1])
      
      pos = gregexpr(FileDirSep, AllTable[[1]][1])
      pos[[1]] <- as.numeric(pos[[1]])
      pos[[1]][length(pos[[1]])-1]
      
      for(j in 1:nrow(AllTable)){
        pos = gregexpr(FileDirSep, AllTable[[1]][j])
        pos[[1]] <- as.numeric(pos[[1]])
        AllTable[[1]][j] <- substring(AllTable[[1]][j],pos[[1]][length(pos[[1]])-1]+1)
      }
      pos = gregexpr(FileDirSep, AllTable[[1]][1])
      pos[[1]] <- as.numeric(pos[[1]])
      pos[[1]][length(pos[[1]])]
      
      for(j in 1:nrow(AllTable)){
        pos = gregexpr(FileDirSep, AllTable[[1]][j])
        pos[[1]] <- as.numeric(pos[[1]])
        AllTable[[1]][j] <- substring(AllTable[[1]][j],1,pos[[1]][length(pos[[1]])]-1)
      }
      
      AllTable[[1]] <- gsub(' ', '_', AllTable[[1]])
      plot_output_list_box <- list()
      
      i <- 2
      temp <- ncol(AllTable)
      for(i in 2:ncol(AllTable)){
        plotnamebox <- paste0("Box_",Names[i])
        plot_output_list_box<- c(plot_output_list_box,list(plotOutput(plotnamebox, height = 600, width = 800)))
        
      }
      do.call(tagList,plot_output_list_box)
    })
    
    observe({
      
      req(input$file1)
      
      AllTable <- read.delim(header = TRUE,input$file1$datapath)
      
      even_indexes <- rev(seq(2,ncol(AllTable),2))
      
      for (i in even_indexes) {
        AllTable[i] = NULL
      }
      Names <- colnames(AllTable)
      
      AllTable[[1]] <- as.character(AllTable[[1]])
      
      j = 1
      for (i in 1:nrow(AllTable)) {
        AllTable[[1]][i] <- gsub('\\','/',as.character(AllTable[[j]][i]), fixed=TRUE)
      }
      FileDirSep <- "/"
        
        
      pos = gregexpr(FileDirSep, AllTable[[1]][j])
      pos[[1]] <- as.numeric(pos[[1]])
      ExName <- substring(AllTable[[1]][j],1,pos[[1]][length(pos[[1]])-1])
      
      pos = gregexpr(FileDirSep, AllTable[[1]][1])
      pos[[1]] <- as.numeric(pos[[1]])
      pos[[1]][length(pos[[1]])-1]
      
      for(j in 1:nrow(AllTable)){
        pos = gregexpr(FileDirSep, AllTable[[1]][j])
        pos[[1]] <- as.numeric(pos[[1]])
        AllTable[[1]][j] <- substring(AllTable[[1]][j],pos[[1]][length(pos[[1]])-1]+1)
      }
      pos = gregexpr(FileDirSep, AllTable[[1]][1])
      pos[[1]] <- as.numeric(pos[[1]])
      pos[[1]][length(pos[[1]])]
      
      for(j in 1:nrow(AllTable)){
        pos = gregexpr(FileDirSep, AllTable[[1]][j])
        pos[[1]] <- as.numeric(pos[[1]])
        AllTable[[1]][j] <- substring(AllTable[[1]][j],1,pos[[1]][length(pos[[1]])]-1)
      }
      
      AllTable[[1]] <- gsub(' ', '_', AllTable[[1]])
      #AllTable[[1]][1] <- substring(AllTable[[1]][1],pos[[1]][length(pos[[1]])-1]+1)
      i <- 2
      temp <- ncol(AllTable)
      for(i in 2:ncol(AllTable)){
          local({
            # Need local so that each item gets its own number. Without it, the value
            # of i in the renderPlot() will be the same across all instances, because
            # of when the expression is evaluated.
            #eval(parse(text=exp1)) 
            plotnamebox <- paste0("Box_",Names[i])
            
            #exp1 <- paste0("boxplot(AllTable$",Names[i],"~Image, data=AllTable, notch=TRUE, col=(c(\"gold\",\"lightblue\")), main=paste0(ExName,\"\n\",Names[i]), ylab=\"count\")")
            exp1 <- paste0("boxplot(AllTable$",Names[i],"~Image, data=AllTable, cex.axis=1.5, cex.lab=1.5, notch=TRUE, col=(c(\"gold\",\"lightblue\")), main=paste0(ExName,\"\n\",Names[i]), ylab=\"",Names[i],"\")")
            eval(parse(text=exp1))
            means <- eval(parse(text=paste0("tapply(AllTable$",Names[i],",AllTable$Image,mean)")))
            points(means,col="red",pch=18)
            recordedPlot = recordPlot(eval(parse(text=exp1)))#recordPlot()
            dev.off()
            
            #plot.new()
            output[[plotnamebox]] <- renderPlot({
              #StartNote
              notification.id <- progressNotification(paste0("Loading ",plotnamebox," boxplot"))
                print(recordedPlot)
              #EndNote
              removeNotification(notification.id) 
              
            })
            
          })
      }
    }) 
    
  
  output$correlationMatrix <- renderPlot({
    req(input$file1)
    
    AllTable <- read.delim(header = TRUE,input$file1$datapath)
    even_indexes <- rev(seq(2,ncol(AllTable),2))
    AllTable[[1]] <- as.character(AllTable[[1]]) 
    
    j = 1
    for (i in 1:nrow(AllTable)) {
      AllTable[[1]][i] <- gsub('\\','/',as.character(AllTable[[j]][i]), fixed=TRUE)
    }
    FileDirSep <- "/"
    
    
    pos = gregexpr(FileDirSep, AllTable[[1]][j])
    pos[[1]] <- as.numeric(pos[[1]])
    ExName <- substring(AllTable[[1]][j],1,pos[[1]][length(pos[[1]])-1])
    
    for (i in even_indexes) {
      AllTable[i] = NULL
    }
    AllTable[1] = NULL
    
    for(i in ncol(AllTable):1){
      if(sum(AllTable[i]) == 0 || is.nan(AllTable[[i]][1])){
        AllTable[i] = NULL
      }
    }
    
    M <- cor(AllTable)
    p.mat <- cor.mtest(AllTable)
    title <-paste0(ExName,"\n","Correlation Matrix")
    
    corrplot(M,title = title ,type="upper", p.mat = p.mat, sig.level = 0.01, mar=c(0,0,3,0))
  })
  
  
  
  output$table <- DT::renderDataTable({
    req(input$file1)
    
    if (is.null(input$file1))
      return(NULL)
    AllTable <- read.delim(header = TRUE,input$file1$datapath)
    even_indexes <- rev(seq(2,ncol(AllTable),2))
    
    for (i in even_indexes) {
      AllTable[i] = NULL
    }
    #return(AllTable)
    DT::datatable(AllTable, options = list(lengthMenu = c(5, 30, 50), pageLength = 10))
  })
  
  # load disk drives
  volumes <- getVolumes()
  shinyDirChoose(input, 'directory', roots=volumes, session=session)
  # 
  path1 <- reactive({
    return(print(parseDirPath(volumes, input$directory)))
    
  })
  # write plots to the choosen path on disk
  dataruw1 <- eventReactive(input$directory, {
    #StartNote
    notification.id <- progressNotification("Downloading...")
    
    
    #showNotification("Downloading", type = "message")
    #The directory for plots output
    pathIn <- parseDirPath(volumes, input$directory)
    print(paste0("Sys.info",' ',Sys.info()['sysname']))
    if(as.character(Sys.info()['sysname']) == "Windows")
      DirSep <- "\\"
    else
      DirSep <- "/"
    
    
    #Table preperation
    MainAllTable <- read.delim(header = TRUE,input$file1$datapath)
    even_indexes <- rev(seq(2,ncol(MainAllTable),2))
    for (i in even_indexes) {
      MainAllTable[i] = NULL
    }
    MainNames <- colnames(MainAllTable)
    j = 1
    MainAllTable[[1]] <- as.character(MainAllTable[[1]]) 
    
    for (i in 1:nrow(MainAllTable)) {
      MainAllTable[[1]][i] <- gsub('\\','/',as.character(MainAllTable[[1]][i]), fixed=TRUE)
    }
    FileDirSep <- "/"
    
    pos = gregexpr(FileDirSep, MainAllTable[[1]][j],fixed = TRUE)
    pos[[1]] <- as.numeric(pos[[1]])
    ExName <- substring(MainAllTable[[1]][j],1,pos[[1]][length(pos[[1]])-1])
    
    ##########################################################################
    ########################Scatter Plots#####################################
    ##########################################################################
    showNotification("Scatter plots started", type = "message")
    AllTable <- MainAllTable
    Names <- MainNames
    folderName <- "Scatter_Plots"
    pathOut <- paste0(pathIn,DirSep,folderName,DirSep)
    #print(pathOut)
    if(!dir.exists(pathOut))
      dir.create(pathOut)
    
    j<-2
    k<-3
    ptlist <- list()
    temp <- ncol(AllTable)-1
    for (j in 2:temp) {
      temp2 <- j+1
      for (k in temp2:ncol(AllTable)) {
        #png(filename=paste(pathOut,Names[j],"_",Names[k],".png",sep = ''),width = 800,height = 600)
        ptitle <- paste0(Names[j],"_",Names[k])
        
        exp1 <- paste0("abc <- ggplot(AllTable, aes(x=",Names[j],", y=",Names[k],")) + ggtitle(paste0(ExName,\"\n\",Names[j],\"_\",Names[k]))+ theme(plot.title = element_text(hjust = 0.5),axis.text=element_text(size=16,face=\"bold\"), axis.title=element_text(size=18,face=\"bold\")) + geom_point(size = 3,shape=19, color=\"blue\")+  labs(x=Names[j], y=Names[k])")
        eval(parse(text=exp1))
        if(input$fileType == "png plots")
          png(paste0(pathOut,ptitle,".png"),width = 800,height = 600)
        else
          pdf(paste0(pathOut,ptitle,".pdf"),width = 12,height = 9)
        print(abc)
        dev.off()
        if(input$checkboxCSV){
          expDT<- paste0("data.frame(",Names[1]," = MainAllTable[1],",Names[j]," =AllTable[j] ,",Names[k] ,"=AllTable[k])")# )) ggplot(AllTable, aes(x=",Names[j],", y=",Names[k],")) + ggtitle(paste0(Names[j],\"_\",Names[k]))+ theme(plot.title = element_text(hjust = 0.5)) + geom_point(shape=18, color=\"blue\")+ geom_smooth(method=lm,  linetype=\"dashed\", color=\"darkred\", fill=\"blue\") + labs(x=Names[j], y=Names[k])")
          df <- eval(parse(text=expDT))
          write.csv(df, file = paste0(pathOut,Names[j],"_",Names[k],'.csv'))
        }
        
      }
    }
    
    showNotification("Scatter plots done", type = "message")  
        
      
      
    #########################################################################
    ########################Histograms#######################################
    #########################################################################
    showNotification("Histograms started", type = "message")
    AllTable <- MainAllTable
    Names <- MainNames
      folderName <- "Histograms"
      pathOut <- paste0(pathIn,DirSep,folderName,DirSep)
      print(pathOut)
      if(!dir.exists(pathOut))
        dir.create(pathOut)
      
      for(t in 1:length(Names)){
        Names[t] <- gsub('\\.', '_', Names[t])
      }
      
      
      #j = 18
      temp <-ncol(AllTable)
      for (j in 2:temp) {
        ptitle <- paste0("Histogram_",Names[j])
        if(input$fileType == "png plots")
          png(paste0(pathOut,ptitle,".png"),width = 800,height = 600)
        else
          pdf(paste0(pathOut,ptitle,".pdf"),width = 12,height = 9)
        abc <- ggplot(AllTable, aes(AllTable[,j])) +
          geom_histogram(aes(y=..count../sum(..count..)) , color="darkblue", fill="lightblue")  +
          ggtitle(paste0(ExName,"\n","Histogram_",Names[j]))+
          labs(x=Names[j], y="Frequency") +
          theme(plot.title = element_text(size = 18,hjust = 0.5),axis.text=element_text(size=16,face="bold"), axis.title=element_text(size=18,face="bold"))
        
        #abc <- abc + theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))
        print(abc)
        dev.off()
        if(input$checkboxCSV){
          expDT<- paste0("data.frame(",Names[1]," = MainAllTable[1],",Names[j]," =AllTable[j] )")# )) ggplot(AllTable, aes(x=",Names[j],", y=",Names[k],")) + ggtitle(paste0(Names[j],\"_\",Names[k]))+ theme(plot.title = element_text(hjust = 0.5)) + geom_point(shape=18, color=\"blue\")+ geom_smooth(method=lm,  linetype=\"dashed\", color=\"darkred\", fill=\"blue\") + labs(x=Names[j], y=Names[k])")
          df <- eval(parse(text=expDT))
          write.csv(df, file = paste0(pathOut,"Histogram_",Names[j],'.csv'))
        }
      }
      
      
      showNotification("Histograms done", type = "message")
      
    #########################################################################
    ########################Box Plots########################################
    #########################################################################
      showNotification("Box plots started", type = "message")
      AllTable <- MainAllTable
      Names <- MainNames
      folderName <- "Box_Plots"
      pathOut <- paste0(pathIn,DirSep,folderName,DirSep)
      if(!dir.exists(pathOut))
        dir.create(pathOut)
      
      AllTable[[1]] <- as.character(AllTable[[1]])
      
      ##########################################################################
      j=1
      for (i in 1:nrow(AllTable)) {
        AllTable[[1]][i] <- gsub('\\','/',as.character(AllTable[[j]][i]), fixed=TRUE)
      }
      FileDirSep <- "/"
      
      pos = gregexpr(FileDirSep, AllTable[[1]][1])
      pos[[1]] <- as.numeric(pos[[1]])
      pos[[1]][length(pos[[1]])-1]
      
      for(j in 1:nrow(AllTable)){
        pos = gregexpr(FileDirSep, AllTable[[1]][j])
        pos[[1]] <- as.numeric(pos[[1]])
        AllTable[[1]][j] <- substring(AllTable[[1]][j],pos[[1]][length(pos[[1]])-1]+1)
      }
      pos = gregexpr(FileDirSep, AllTable[[1]][1])
      pos[[1]] <- as.numeric(pos[[1]])
      pos[[1]][length(pos[[1]])]
      
      for(j in 1:nrow(AllTable)){
        pos = gregexpr(FileDirSep, AllTable[[1]][j])
        pos[[1]] <- as.numeric(pos[[1]])
        AllTable[[1]][j] <- substring(AllTable[[1]][j],1,pos[[1]][length(pos[[1]])]-1)
      }
      ##########################################################################
      
      AllTable[[1]] <- gsub(' ', '_', AllTable[[1]])
      
      #AllTable[[1]][1] <- substring(AllTable[[1]][1],pos[[1]][length(pos[[1]])-1]+1)
      i <- 2
      temp <- ncol(AllTable)
      for(i in 2:ncol(AllTable)){
        ptitle <- paste0("Box_",Names[i])
        if(input$fileType == "png plots")
          png(paste0(pathOut,ptitle,".png"),width = 800,height = 600)
        else
          pdf(paste0(pathOut,ptitle,".pdf"),width = 12,height = 9)
        exp1 <- paste0("boxplot(AllTable$",Names[i],"~Image, data=AllTable, notch=TRUE, cex.lab=1.5, cex.axis=1.5, col=(c(\"gold\",\"lightblue\")), main=paste0(ExName,\"\n\",Names[i]), ylab=\"",Names[i],"\")")
        eval(parse(text=exp1))
        means <- eval(parse(text=paste0("tapply(AllTable$",Names[i],",AllTable$Image,mean)")))
        points(means,col="red",pch=18)
        #par(cex.axis=2)
        dev.off()
        
        if(input$checkboxCSV){
          expDT<- paste0("data.frame(",Names[1]," =MainAllTable[1],",Names[i]," =AllTable[i] )")
          df <- eval(parse(text=expDT))
          write.csv(df, file = paste0(pathOut,"Box_",Names[i],'.csv'))
        }
        
      }
      
      
      showNotification("Box plots done", type = "message")
      
    #########################################################################
    #######################Correlation Matrix################################
    #########################################################################
      showNotification("Correlations matrix plot started", type = "message")
      AllTable <- MainAllTable
      Names <- MainNames
      folderName <- "Correlation_Matrix"
      pathOut <- paste0(pathIn,DirSep,folderName,DirSep)
      print(pathOut)
      if(!dir.exists(pathOut))
        dir.create(pathOut)
      
      AllTable[1] = NULL
      #AllTable[20] = NULL
      #AllTable[19] = NULL
      #t <- ncol(AllTable)
      for(i in ncol(AllTable):1){
        if(sum(AllTable[i]) == 0 || is.nan(AllTable[[i]][1])){
          AllTable[i] = NULL
        }
      }
      
      M <- cor(AllTable)
      p.mat <- cor.mtest(AllTable)
      #png(filename=paste(pathOut,"Correlation_Matrix",".png",sep = ''),width = 800,height = 600)
      if(input$fileType == "png plots")
        png(paste0(pathOut,"Correlation_Matrix",".png"),width = 800,height = 600)
      else
        pdf(paste0(pathOut,"Correlation_Matrix",".pdf"),width = 12,height = 9)
      corrplot(M,title = paste0(ExName,"\n","Correlation Matrix"),type="upper", p.mat = p.mat, sig.level = 0.01, 
               # hide correlation coefficient on the principal diagonal
               mar=c(0,0,3,0))
      dev.off()
      if(input$checkboxCSV){
        write.csv(M, file = paste0(pathOut,"Correlation_Matrix",'.csv'))
      }
      
      showNotification("Correlations matrix plot done", type = "message")
      #EndNote
      removeNotification(notification.id) 
      showNotification("Download Completed", type = "message")
  })
  
  output$summary <- renderText({
    
  })
  
  observeEvent(input$directory, {
    output$OutputDir <- renderText({
      print(parseDirPath(volumes, input$directory))
    })
  })
  
  # triger to download the plots
  observeEvent(input$downData, {
    if(length(input$file1$datapath) > 0){
      if(length(path1()) > 0){
        dataruw1()
        #path1()
      }else
      {
        showNotification("Please select a download directory", type = "error")
      }
    }else{
      showNotification("Please load input file first", type = "error")
    }
  })
  # Exits the application
  observeEvent(input$exitApp, {
    js$closeWindow()
    stopApp()
  })
  # resets download path text output
  observeEvent(input$file1, {
    output$OutputDir <- renderText({
      print("")
    })
  })
  
  
  
})
