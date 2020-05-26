#source('directoryInput.R')
library(shiny)
library(DT)
library(shinyjs)
library(lme4)
library(msqrob2)
library(MSnbase)
library(grDevices)
library(limma)
library(graphics)
library(ggplot2)
source("utilities.R")
source("helper.R")
library("ExploreModelMatrix")
variables <- reactiveValues(pe=NULL)
#Max file size: 500 Mb
options(shiny.maxRequestSize=500*1024^2)
# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {

  ###########################################
  #Input Tab
  ###########################################

  #Change the default for protein groups: if MaxQuant: set to true, if not, set to FALSE

  #counter to avoid app to be killed when running on shiny server.
  output$keepAlive <- renderText({
        req(input$count)
        paste("keep alive ", input$count)
      })



      anotationFileToDownload <- eventReactive(input$goAnnotation, {
      data.frame(run = colnames(variables$pe[[1]]))
      })
        output$DownloadAnnot <- downloadHandler(
            filename = function() {
              paste0(input$project_name, gsub(" |:","-",Sys.time()),"_experimental_annotation.xlsx")
            },
            content = function(file) {
              openxlsx::write.xlsx(anotationFileToDownload(), file, colNames = TRUE)
            }
          )





  ########################################################
  #Clear datapaths of backslashes (Needed on Windows only)
  ########################################################

  featuresDatapath <- reactive({getDataPath(input$peptides$datapath)})
  annotationDatapath <- reactive({getDataPath(input$annotation$datapath)})

  ########################################
  #Set Filter option
  ########################################
  filterOptions <- reactive({
    if(is.null(input$peptides)){
      NULL
    } else{
      make.names(as.vector(as.matrix(read.table(featuresDatapath(), nrows=1, sep="\t", quote="", comment.char = ""))))
    }
  })

  selectedFilter <- reactive({
    if(!any(c("Reverse", "Contaminant", "Potential contaminant", "Potential.contaminant") %in% filterOptions())) {
      NULL
    } else{c("Reverse", "Contaminant", "Potential contaminant", "Potential.contaminant")[c("Reverse", "Contaminant", "Potential contaminant", "Potential.contaminant") %in% filterOptions()]}
  })

  output$selectFilters <- renderUI({
    selectInput("filter", NULL, filterOptions(), multiple=TRUE, selected=selectedFilter(), width = '100%')})

  selNormalisation <- reactive({
      selNormalisation <- NULL
  })

  output$selectNormalisation <- renderUI({
    selectInput("normalisation", NULL, c("quantiles","vsn","quantiles.robust", "center.median", "center.mean", "div.mean","div.median","diff.median" ,"max", "sum", "none"), selected=selNormalisation(), width = '100%') #"loess.affy" and "loess.pairs" left out on purpose because they remove everything with at least 1 NA!
    })

  ########################################
  #Generate options for fixed effect variables
  ########################################

  colClasses <- reactive({
    if(isTRUE(input$asis_numeric)){
      colClasses <- "keep"
    } else{colClasses <- "factor"}
    return(colClasses)
  })


observe({
  shinyjs::onclick("button_project_name",
                   shinyjs::toggle(id = "tooltip_project_name", anim = TRUE))
  shinyjs::onclick("button_input_type",
                   shinyjs::toggle(id = "tooltip_input_type", anim = TRUE))
  shinyjs::onclick("button_features",
                   shinyjs::toggle(id = "tooltip_features", anim = TRUE))
  shinyjs::onclick("button_annotation",
                   shinyjs::toggle(id = "tooltip_annotation", anim = TRUE))
  shinyjs::onclick("button_asis_numeric",
                   shinyjs::toggle(id = "tooltip_asis_numeric", anim = TRUE))
  shinyjs::onclick("button_newExpAnnText",
                   shinyjs::toggle(id = "tooltip_newExpAnnText", anim = TRUE))
  shinyjs::onclick("button_cite",
                   shinyjs::toggle(id = "tooltip_cite", anim = TRUE))
  shinyjs::onclick("button_notinlist",
                   shinyjs::toggle(id = "tooltip_notinlist", anim = TRUE))
  shinyjs::onclick("button_logtransform",
                   shinyjs::toggle(id = "tooltip_logtransform", anim = TRUE))
  shinyjs::onclick("button_log_base",
                   shinyjs::toggle(id = "tooltip_log_base", anim = TRUE))
  shinyjs::onclick("button_normalisation",
                   shinyjs::toggle(id = "tooltip_normalisation", anim = TRUE))
  shinyjs::onclick("button_onlysite",
                   shinyjs::toggle(id = "tooltip_onlysite", anim = TRUE))
  shinyjs::onclick("button_smallestUniqueGroups",
                   shinyjs::toggle(id = "tooltip_smallestUniqueGroups", anim = TRUE))
  shinyjs::onclick("button_minIdentified",
                   shinyjs::toggle(id = "tooltip_minIdentified", anim = TRUE))
  shinyjs::onclick("button_filter",
                   shinyjs::toggle(id = "tooltip_filter", anim = TRUE))
  shinyjs::onclick("button_selColPlotNorm",
                   shinyjs::toggle(id = "tooltip_selColPlotNorm", anim = TRUE))
  shinyjs::onclick("button_preprocessing_extension",
                   shinyjs::toggle(id = "tooltip_preprocessing_extension", anim = TRUE))
  shinyjs::onclick("button_h4_int_transformation",
                   shinyjs::toggle(id = "tooltip_h4_int_transformation", anim = TRUE))
  shinyjs::onclick("button_h4_full_preprocessing",
                   shinyjs::toggle(id = "tooltip_h4_full_preprocessing", anim = TRUE))
  shinyjs::onclick("button_h4_MDS_full_preprocessing",
                   shinyjs::toggle(id = "tooltip_h4_MDS_full_preprocessing", anim = TRUE))

  shinyjs::onclick("button_proteins",
                   shinyjs::toggle(id = "tooltip_proteins", anim = TRUE))
  })



  ##############################################
  #Normalization tab
  #############################################
  #Check waarom er nog altijd "Select accession" staat en niet "prot"!!!!
  selectedProteins <- reactive({
    if(!any(c("Proteins","prot","Accession","accession") %in% filterOptions())) {
      ""
    }else{
      c("Proteins","prot","Accession","accession")[which(c("Proteins","prot","Accession","accession") %in% filterOptions())[1]]
      }
  })

  output$selectProteins <- renderUI({
    selectInput("proteins", label=NULL, filterOptions(), multiple=FALSE, selected=selectedProteins() )
  #!Be careful with "selectizeInput" -> fixed and random all of a sudden might not work anymore...
  })

  ###Drop down menu for plot normalization Plot###
  plotDependentVars <- reactive({
    if (!is.null(variables$pe))
    as.list(c("none",colnames(colData(variables$pe))))
  })


  output$selectColPlotNorm1 <- renderUI({
    div(class="MSqRob_input_container",
        list(
          tags$label("Color variable", `for`="selColPlotNorm", class="MSqRob_label"),
          tags$button(id="button_selColPlotNorm", tags$sup("[?]"), class="MSqRob_tooltip"),
          selectInput("selColPlotNorm1", label=NULL,  plotDependentVars()),
          hidden(helpText(id="tooltip_selColPlotNorm","Select the variable by which the densities should be colored."))
        )
    )
  })

  output$selectColPlotProt <- renderUI({
    div(class="MSqRob_input_container",
        list(
          tags$label("Color variable", `for`="selColPlotProt", class="MSqRob_label"),
          tags$button(id="button_selColPlotProt", tags$sup("[?]"), class="MSqRob_tooltip"),
          selectInput("selColPlotProt", label=NULL,  plotDependentVars()),
          hidden(helpText(id="tooltip_selColPlotPlot","Select the variable by which the densities should be colored."))
        )
    )
  })

  output$nfeaturesNormalized = renderText(NULL)

  ####Raw feature density plot####
  colorsNorm <- reactive({
    colors <- 1
    try(
      {colordata <- colData(variables$pe)[,input$selColPlotNorm1]
      colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8,"Spectral"))(length(unique(colordata)))
      colors <- colors[as.numeric(droplevels(as.factor(colordata)))]
      },silent=TRUE)
    return(colors)
  })


  ####NEW######

  observeEvent({input$peptides
  input$annotation},{
      if (!is.null(input$peptides$name)&!is.null(input$annotation$name)){
      removeNotification(id="noPeptideFile")
      if(input$input_type=="MaxQuant"){
          ecols <- MSnbase::grepEcols(featuresDatapath(), "Intensity ", split = "\t")

      }
      if(input$input_type=="moFF"){
          ecols <- MSnbase::grepEcols(featuresDatapath(), "sumIntensity_", split = "\t")
      }
      if(input$input_type=="mzTab"){
          ecols <- MSnbase::grepEcols(featuresDatapath(), "sumIntensity_", split = "\t")
      }
      peOut <- try(readFeatures(table = featuresDatapath(), fnames = 1, ecol = ecols,
              name = "featureRaw", sep="\t"))
      if (class(peOut)[1]=="try-error"|length(ecols)<1) {
          showNotification("Upload proper peptide file that is compliant with the input type",id="noProperPeptideFile",type="error",duration=NULL,closeButton=FALSE)
      } else {
      removeNotification(id="noProperPeptideFile")
      rowData(peOut[["featureRaw"]])$nNonZero <- rowSums(assay(peOut[["featureRaw"]]) > 0)
      peOut <- zeroIsNA(peOut,"featureRaw")
      if (!is.null(input$annotation$name)) {
          removeNotification(id="noAnnotation")
          if(isTRUE(as.logical(grep(".xlsx[/\\]*$",input$annotation$name)))){
              exp_annotation <- as.data.frame(unclass(openxlsx::read.xlsx(annotationDatapath())))
              }
              else {
              exp_annotation <- read.table(annotationDatapath(), sep="\t", header=TRUE, row.names = NULL, quote="", stringsAsFactors = TRUE)
              }
              runName <- which(vapply(exp_annotation, function(x) return(identical(sort(as.character(colnames(peOut[[1]]))),
                    sort(as.character(x)))), FUN.VALUE = TRUE))
              if (length(runName)==0) {
                  showNotification("Adjust the annotation file and reupload it. Make sure that exactly one column in the annotation file has elements equal to the mass spec run names",type="error",duration=NULL,id="annotationError",closeButton=FALSE)
                  } else {
                  removeNotification(id="annotationError")

                  if (length(runName)>0) runName <- runName[1]
                  rownames(exp_annotation) <- as.character(exp_annotation[,runName])
                  exp_annotation <- exp_annotation[colnames(peOut[[1]]),]
                  for (j in colnames(exp_annotation))
                      colData(peOut)[[j]] <- exp_annotation[,j]
                  }
              } else {showNotification("Upload Annotation File",type="error",duration=NULL,id="noAnnotation",closeButton=FALSE)}
      if (input$logtransform) peOut <- logTransform(peOut, base = 2,i="featureRaw",name="featureLog")
      variables$pe <- peOut
      output$nfeaturesRaw <- renderText(nrow(peOut[["featureRaw"]]))
      }
      } else showNotification("Upload peptide and annotation file",id="noPeptideFile",closeButton=FALSE,type="error",duration=NULL)
  })


  observeEvent(input$goNorm,{
    if (!is.null(variables$pe)){
      if ("featureLog" %in% names(variables$pe) & input$normalisation!="vsn"){
        i <- "featureLog"
        } else {i <- "featureRaw"}

      peOut <- variables$pe[,,which(names(variables$pe)%in%c("featureRaw","featureLog"))]
      if(input$input_type=="MaxQuant"){
      for (j in input$selectedFilters)
         peOut[[i]] <- peOut[[i]][rowData(peOut[[i]])[[j]]!= "+", ]
      }

      if (input$smallestUniqueGroups) peOut[[i]] <- peOut[[i]][rowData(peOut[[i]])[[selectedProteins()]] %in% smallestUniqueGroups(rowData(peOut[[i]])[[selectedProteins()]]),]

      peOut[[i]] <- peOut[[i]][rowData(peOut[[i]])$nNonZero >= input$minIdentified, ]
      if (input$normalisation=="none"){
          peOut <- addAssay(peOut,peOut[[i]],"featureNorm")
      } else{
          peOut <- try(normalize(peOut,i=i,name="featureNorm",method=input$normalisation))
          }
      if (class(peOut)!="try-error") {
         variables$pe <- peOut
         output$nfeaturesNormalized <- renderText(nrow(peOut[["featureNorm"]]))
        }
      }
  })


  ###Possibilities for zooming
  rangesRaw <- reactiveValues(x = NULL, y = NULL)
  rangesNorm1 <- reactiveValues(x = NULL, y = NULL)
  rangesMDS <- reactiveValues(x = NULL, y = NULL)
  rangesMDSProt <- reactiveValues(x = NULL, y = NULL)

  ####Raw density plot with zoom####
  observeEvent(input$plotRaw_dblclick, {
    brush <- input$plotRaw_brush
    if (!is.null(brush)) {
      rangesRaw$x <- c(brush$xmin, brush$xmax)
      rangesRaw$y <- c(brush$ymin, brush$ymax)

    } else {
      rangesRaw$x <- NULL
      rangesRaw$y <- NULL
    }
  })

  ####Normal density plot with zoom####
  observeEvent(input$plotNorm1_dblclick, {
    brush <- input$plotNorm1_brush
    if (!is.null(brush)) {
      rangesNorm1$x <- c(brush$xmin, brush$xmax)
      rangesNorm1$y <- c(brush$ymin, brush$ymax)

    } else {
      rangesNorm1$x <- NULL
      rangesNorm1$y <- NULL
    }
  })

  ####MDS plot with zoom####
  observeEvent(input$plotMDS_dblclick, {
    brush <- input$plotMDS_brush
    if (!is.null(brush)) {
      rangesMDS$x <- c(brush$xmin, brush$xmax)
      rangesMDS$y <- c(brush$ymin, brush$ymax)

    } else {
      rangesMDS$x <- NULL
      rangesMDS$y <- NULL
    }
  })

  ####MDS plot proteins with zoom####
  observeEvent(input$plotMDSProt_dblclick, {
    brush <- input$plotMDSProt_brush
    if (!is.null(brush)) {
      rangesMDSProt$x <- c(brush$xmin, brush$xmax)
      rangesMDSProt$y <- c(brush$ymin, brush$ymax)

    } else {
      rangesMDSProt$x <- NULL
      rangesMDSProt$y <- NULL
    }
  })


  #Render the preprocessing plots

  output$plotRaw<- renderPlot({
    if(!is.null(variables$pe)){
      if (input$logtransform){
          se <- variables$pe[["featureLog"]]
          } else{
          se <- variables$pe[["featureRaw"]]
         }
      plotDens(se,rangesRaw$x,rangesRaw$y,colorsNorm())
      }
  })

  output$plotNorm1<- renderPlot({
  if(!is.null(variables$pe)){
        if ("featureNorm" %in% names(variables$pe)){
            plotDens(variables$pe[["featureNorm"]],rangesNorm1$x,rangesNorm1$y,colorsNorm())
        }
        }
  })

  output$plotMDS <- renderPlot({
  if(!is.null(variables$pe)){
        if ("featureNorm" %in% names(variables$pe)){
            if (length(colorsNorm()) > 1)  {
                colors <- colorsNorm()
                } else {
                colors <- 1
                }
            if (input$plotMDSPoints) {
                labels <- NULL
                pch <- 1
                } else{
                labels <- colnames(variables$pe[["featureNorm"]])
                pch <- NULL
                }
            limma::plotMDS(assay(variables$pe[["featureNorm"]]),col=colors,xlim=rangesMDS$x,ylim=rangesMDS$y,pch=pch,labels=labels)
            }
        }
  })

  ######################
  #### Summarisation tab
  ######################

  observeEvent(input$goSum,{
  withBusyIndicatorServer("goSum", {
    if (!is.null(variables$pe)){
      if ("featureNorm" %in% names(variables$pe)){
          if (input$summarisation=="robust") fun <-MsCoreUtils::robustSummary
          if (input$summarisation=="sum") fun <- base::colSums
          if (input$summarisation=="mean") fun <- base::colMeans
          if (input$summarisation=="median") fun <- matrixStats::colMedians
          if (input$summarisation=="medpolish") fun <- stats::medpolish
          hlpNames<-which(names(variables$pe) %in% c("featureRaw","featureLog","featureNorm"))
          peOut <- variables$pe[,,hlpNames]
          if (input$summarisation!="none"){
              peOut <- aggregateFeatures(peOut, i="featureNorm", fcol = selectedProteins(), na.rm = TRUE, name="summarized",fun=fun)
              } else {
              peOut <- addAssay(peOut,peOut[["featureNorm"]],name="summarized")
              }
          variables$pe <- peOut
          } else (showNotification("Run normalisation first",type="error"))
      }
  })
  })


  ####MDS plot proteins with zoom####
  observeEvent(input$plotMDSProt_dblclick, {
    brush <- input$plotMDSProt_brush
    if (!is.null(brush)) {
      rangesMDSProt$x <- c(brush$xmin, brush$xmax)
      rangesMDSProt$y <- c(brush$ymin, brush$ymax)

    } else {
      rangesMDSProt$x <- NULL
      rangesMDSProt$y <- NULL
    }
  })


    colorsProt <- reactive({
        colors <- 1
        try(
          {colordata <- colData(variables$pe)[,input$selColPlotProt]
          colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8,"Spectral"))(length(unique(colordata)))
          colors <- colors[as.numeric(droplevels(as.factor(colordata)))]
          },silent=TRUE)
        return(colors)
    })

    output$plotMDSProt <- renderPlot({
    if(!is.null(variables$pe)){
          if ("summarized" %in% names(variables$pe)) {
              se <- assay(variables$pe[["summarized"]])
              }  else {se <- NULL
              }
          if (!is.null(se)){
              if (length(colorsProt()) > 1)  {
                  colors <- colorsProt()
                  } else {
                  colors <- 1
                  }
              if (input$plotMDSPointsProt) {
                  labels <- NULL
                  pch <- 1
                  } else{
                  labels <- colnames(se)
                  pch <- NULL
                  }
              plotMDS(se,col=colors,xlim=rangesMDSProt$x,ylim=rangesMDSProt$y,pch=pch,labels=labels)
              } else {
              plot(0,0,xaxt="none",yaxt="none",col=0,xlab="",ylab="")
              text(0,0,"Run summarisation",col="red",cex=3)
              }
          }
    })

  #######################
  #Model Tab
  #######################

  output$selectFixed <- renderUI({
      h4(paste(colnames(colData(variables$pe)),collapse="\n"), class="MSqRob_topheader")
      })

  visDesign <- reactive({
      out <- VisualizeDesign(colData(variables$pe),input$designformula)
      rank <- qr(out$designmatrix)$rank
      if (rank==nrow(out$designmatrix)) {
        showNotification(paste0("The model is overparameterized. ",
                        "The residual degrees of freedom is 0. ",
                        "Variances and standard errors can not be ",
                        "estimated from data with this design."),
                        type="error",
                        closeButton=FALSE,
                        duration=NULL,
                        id="modError"
                        )
        } else {
          removeNotification(id="modError")
        }
      if (rank<ncol(out$designmatrix)) {
          showNotification(paste0("The model is poorly parametrised. ",
                          "Not all model parameters can be estimated ",
                          "and the model will not be fitted."
                          ),
                          type="error",
                          closeButton=FALSE,
                          duration=NULL,
                          id="modError2"
                          )
          } else {
            removeNotification(id="modError2")
          }
      return(out)
      })
      output$fitted_values_plot <- renderUI({
          lapply(1:length(visDesign()[[2]]),
                 function(j){
                 renderPlot(visDesign()[[2]][[j]])
                 })
        })
      observeEvent(input$fitModel,{
      withBusyIndicatorServer("fitModel", {
          peOut <- variables$pe
          peOut <- try(msqrob(object=peOut,i="summarized", formula=stats::as.formula(input$designformula),overwrite=TRUE, ridge=input$doRidge==1))

          if (class(peOut)=="Features") {
          variables$pe <- peOut
          }
      })
      })
      output$annotationDataMatrix <- DT::renderDataTable(as.data.frame(colData(variables$pe)))
  ###########
  #Inference tab
  ###########
  ranges <- reactiveValues(x = NULL, y = NULL)

  output$modelParams <- renderUI({paste(colnames(visDesign()[[3]]),collapse=" \n")})
  dataAll <- reactive({
      data <- topFeatures(rowData(variables$pe[["summarized"]])$msqrobModels,makeContrast(input$contrast,parameterNames=colnames(visDesign()[[3]])))
      data$minusLog10Pval <- - log10(data$pval)
      return(data)
      }
  )

  output$volcanoPlot <- renderPlot(
      makeVolcanoPlot2(dataAll(),clickInfo,input,ranges)
  )

    #When a double-click happens, check if there's a brush on the plot.
    #If so, zoom to the brush bounds; if not, reset the zoom.
    observeEvent(input$plotVolcano_dblclick, {
      brush <- input$plotVolcano_brush
      if (!is.null(brush)) {
        ranges$x <- c(brush$xmin, brush$xmax)
        ranges$y <- c(brush$ymin, brush$ymax)

      } else {
        ranges$x <- NULL
        ranges$y <- NULL
      }

      #Set selection to zero: happens already if ranges change, but should also happen on normal double click
      proxy %>% DT::selectRows(NULL)

    })

    #for zooming
    clickInfo <- reactive({
      # Because it's a ggplot2, we don't need to supply xvar or yvar; if this
      # were a base graphics plot, we'd need those.
      if(!is.null(ranges$x) && !is.null(ranges$y)){clickInfo <- subset(dataAll(), (dataAll()$logFC>ranges$x[1] & dataAll()$logFC<ranges$x[2] & (-log10(dataAll()$pval))>ranges$y[1] & (-log10(dataAll()$pval))<ranges$y[2]))
      } else if(is.null(ranges$x) && is.null(ranges$y)){clickInfo <- dataAll()}
      return(clickInfo)
    })

    data <- reactive({
        data <- clickInfo()
        if (input$sigOnly){
                          return(na.exclude(data[data$adjPval<input$alpha,]))
                          } else {
                          return(data)
                          }
   })

output$table<-DT::renderDT(
  formatSignif(datatable(data()[,1:6]),columns=1:ncol(data()),digits=3)
  )

#Set table Proxy so as to reduce the table according to the zoom in the plot and to highlight points
proxy <- dataTableProxy('table')

observeEvent(input$plotVolcano_click, {

  selected <- nearPoints(clickInfo(), input$plotVolcano_click, addDist = TRUE,maxpoints=1, xvar="logFC", yvar="minusLog10Pval")
  sel_rows <- which(rownames(clickInfo()) %in% rownames(selected))
  #Rows which were selected and selected again are removed, rows which were already selected but not selected again are retained
  #Don't sort this! Otherwise reacalculated.
  new_rows <- c(sel_rows[!sel_rows%in%input$table_rows_selected], input$table_rows_selected[!input$table_rows_selected%in%sel_rows])

  proxy %>% DT::selectRows(new_rows)
})

#Enable or disable add brush to selection and remove brush from selection buttons
observe({
  if (is.null(input$plotVolcano_brush)) {
    shinyjs::disable("add_area_selection")
    shinyjs::disable("remove_area_selection")
  } else {
    shinyjs::enable("add_area_selection")
    shinyjs::enable("remove_area_selection")
  }
})

observeEvent(input$add_area_selection, {

  selected <- brushedPoints(clickInfo(), input$plotVolcano_brush, xvar="logFC", yvar="minusLog10Pval")
  sel_rows <- which(rownames(clickInfo) %in% rownames(selected))
  #Rows which were selected and selected again are removed, rows which were already selected but not selected again are retained
  #Don't sort this! Otherwise reacalculated.
  new_rows <- unique(c(input$table_rows_selected,sel_rows))

  proxy %>% DT::selectRows(new_rows)
})

observeEvent(input$remove_area_selection, {

  selected <- brushedPoints(clickInfo(), input$plot1_brush, xvar="logFC", yvar="minusLog10Pval")
  sel_rows <- which(rownames(clickInfo) %in% rownames(selected))
  #Rows which were selected and selected again are removed, rows which were already selected but not selected again are retained
  #Don't sort this! Otherwise reacalculated.
  new_rows <- input$table_rows_selected[!(input$table_rows_selected %in% sel_rows)]

  proxy %>% DT::selectRows(new_rows)
})

observeEvent(input$remove_all_selection, {
  proxy %>% DT::selectRows(NULL)
})





  detailPlotHlp <- reactive(
      makeDetailPlots(variables$pe,clickInfo,input)
      )

  output$detailPlots <- renderUI(
      lapply(1:length(detailPlotHlp()),function(j) renderPlot(detailPlotHlp()[[j]]))
  )

  output$selectColDetailPlot2 <- renderUI({
    div(class="MSqRob_input_container",
        list(
          tags$label("Color variable", `for`="selColDetailPlot2", class="MSqRob_label"),
          tags$button(id="button_selColDetailPlot2", tags$sup("[?]"), class="MSqRob_tooltip"),
          selectInput("selColDetailPlot2", label=NULL,  plotDependentVars()),
          hidden(helpText(id="tooltip_selColDetailPlot2","Select the variable by which the densities should be colored."))
        )
    )
  })

  output$selectHorizontalDetailPlot2 <- renderUI({
    div(class="MSqRob_input_container",
        list(
          tags$label("Horizontal split", `for`="selHorDetailPlot2", class="MSqRob_label"),
          tags$button(id="button_selHorDetailPlot2", tags$sup("[?]"), class="MSqRob_tooltip"),
          selectInput("selHorDetailPlot2", label=NULL,  plotDependentVars()),
          hidden(helpText(id="tooltip_selHorDetailPlot2","Select the variable by which the plots should be stratified horizontally. Select \'none\' to show all data in one plot."))
        )
    )
  })

  output$selectVerticalDetailPlot2 <- renderUI({
      div(class="MSqRob_input_container",
          list(
            tags$label("Vertical split", `for`="selVertDetailPlot2", class="MSqRob_label"),
            tags$button(id="button_selVertDetailPlot2", tags$sup("[?]"), class="MSqRob_tooltip"),
            selectInput("selVertDetailPlot2", label=NULL,  plotDependentVars()),
            hidden(helpText(id="tooltip_selVertDetailPlot2","Select the variable by which the plots should be stratified vertically. Select \'none\' to show all data in one plot."))
          )
      )
    })
  #Stop the App when closing the browser or ending the session
  session$onSessionEnded(stopApp)
})
