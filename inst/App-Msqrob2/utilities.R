
#A function to prompt a folder chooser under Mac OS X,
#normal choose.dir sometimes gives just NA because system has no access to finder
choose.dir2 <- function(default = NA, caption = NA) {
  command = 'osascript'

  #Dit was het:
  #args = "-e 'tell application \"Finder\"' -e 'activate' -e 'POSIX path of (choose folder{{prompt}}{{default}})' -e 'end tell'"
  #'-e "POSIX path of (choose folder{{prompt}}{{default}})"'

  #osascript -e 'tell application "Terminal" to return POSIX path of (choose file)'

  #Find App that is in front
  args1 = "-e 'tell application \"System Events\"' -e 'set frontApp to name of first application process whose frontmost is true' -e 'end tell'"
  suppressWarnings({
    frontApp = system2(command, args = args1, stderr = TRUE)
  })

  #The application that is in front should open the  file choose window
  args = paste0("-e 'tell application \"",frontApp,"\" to return POSIX path of (choose folder{{prompt}}{{default}})'")

  #system2('osascript', args = "-e 'tell application \"Finder\"' -e 'POSIX path of (choose folder{{prompt}}{{default}})' -e 'end tell'", stderr = TRUE)

  if (!is.null(caption) && !is.na(caption) && nzchar(caption)) {
    prompt = sprintf(' with prompt \\"%s\\"', caption)
  } else {
    prompt = ''
  }
  args = sub('{{prompt}}', prompt, args, fixed = TRUE)

  if (!is.null(default) && !is.na(default) && nzchar(default)) {
    default = sprintf(' default location \\"%s\\"', path.expand(default))
  } else {
    default = ''
  }
  args = sub('{{default}}', default, args, fixed = TRUE)

  suppressWarnings({
    path = system2(command, args = args, stderr = TRUE)
  })
  if (!is.null(attr(path, 'status')) && attr(path, 'status')) {
    # user canceled
    path = NA
  }

  return(path)
}

choose.dir_Linux <- function(default = NA, caption = NA) {

  command = "zenity"

  args1 = "--file-selection --directory"

  if (!is.null(caption) && !is.na(caption) && nzchar(caption)) {
    prompt = sprintf(' with prompt \\"%s\\"', caption)
  } else {
    prompt = ''
  }

  if (!is.null(default) && !is.na(default) && nzchar(default)) {
    default = path.expand(default) #sprintf(' default location \\"%s\\"', path.expand(default))
  } else {
    default = "/dev/null"
  }

  suppressWarnings({
    path = system2(command, args = args1, stderr = default, stdout=TRUE)
  })

  if (!is.null(attr(path, 'status')) && attr(path, 'status')) {
    # user canceled
    path = NA
  }

  return(path)
}


#Function to convert data paths if on windows
getDataPath <- function(datapath){
  if(Sys.info()['sysname']=="Windows"){
    datapath <- gsub("\\","/",datapath, fixed=TRUE)
  }
  return(datapath)
}


#Function to plot densities
getDensXlimYlim <- function(se){
  densAll=apply(assay(se),2,density,na.rm=TRUE)
  ymax=max(vapply(densAll,function(d) max(d$y),1))
  rangematrix <- vapply(densAll,function(d) range(d$x, na.rm=TRUE), c(1,1)) #no longer range(eset), but range of d$x!
  xlim=range(rangematrix,na.rm=TRUE)
  ylim=c(0,ymax)
  return(list(densAll=densAll, xlim=xlim, ylim=ylim))
}

plotDens <- function(se, xlim=NULL, ylim=NULL, colors=1, las=1, frame.plot=FALSE, ...)
{

      hlp <- getDensXlimYlim(se)
      if (is.null(xlim)) xlim<-hlp$xlim
      if (is.null(ylim)) ylim <- hlp$ylim
      if (length(colors)>1) {
        plot(hlp$densAll[[1]],col=colors[1],xlim=xlim,ylim=ylim, las=las, frame.plot=frame.plot, main="",...)
        for (i in 2:ncol(se)) lines(hlp$densAll[[i]],col=colors[i])
      } else {
        plot(hlp$densAll[[1]],xlim=xlim,ylim=ylim, las=las, frame.plot=frame.plot, main="",...)
        for (i in 2:ncol(se)) lines(hlp$densAll[[i]])
      }
}


makeVolcanoPlot <- function(dataset,clickInfo,input,ranges){
  if (!is.null(dataset)){
  volcano <- ggplot(dataset,
                    aes(x = logFC,
                    y = -log10(pval),
                    color = adjPval < input$alpha)) +
            geom_point(cex = 2.5) +
            scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
            theme_minimal() +
            coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = TRUE)


  s <- input$table_rows_selected
  if (length(s)) {
    subdataset <- clickInfo()[s, , drop = FALSE]
    return(volcano +
           geom_point(data=subdataset[subdataset$adjPval<input$alpha,],col="purple",size=5,show.legend=FALSE) +
           geom_point(data=subdataset[subdataset$adjPval>=input$alpha,],col="grey",size=5,show.legend=FALSE)
           )
    } else {
    return(return(volcano))
    }
  }
  }


  makeDetailPlots <- function(pe,clickInfo,input){
      if (!is.null(pe)){
        s <- input$table_rows_selected
        if (length(s)==1)
        {
          featureName <- rownames(clickInfo())[s]
          pePlot <- pe[featureName,,c("featureNorm","summarized")]
          pePlotDf <- data.frame(longFormat(pePlot))
          pePlotDf$assay <- factor(pePlotDf$assay,
                                  levels = c("featureNorm", "summarized"))
          if (input$selColDetailPlot2!="none"){
              if(class(colData(pePlot)[[input$selColDetailPlot2]])=="factor"){
              pePlotDf[,input$selColDetailPlot2] <- as.factor(as.character(colData(pePlot)[pePlotDf$colname,input$selColDetailPlot2]))
              }
          }
          if (input$selVertDetailPlot2!="none"){
              if(class(colData(pePlot)[[input$selVertDetailPlot2]])=="factor"){
              pePlotDf[,input$selVertDetailPlot2] <- as.factor(as.character(colData(pePlot)[pePlotDf$colname,input$selVertDetailPlot2]))
              }
          }
          if (input$selHorDetailPlot2!="none"){
              if(class(colData(pePlot)[[input$selHorDetailPlot2]])=="factor"){
              pePlotDf[,input$selHorDetailPlot2] <- as.factor(as.character(colData(pePlot)[pePlotDf$colname,input$selHorDetailPlot2]))
              }
          }
          if (input$logtransform) {
             ylab <- "feature intensity (log2)"
             } else {
             ylab <- "feature intensity"
             }
          p1 <- ggplot(data = pePlotDf,
                 aes(x = colname,
                     y = value,
                     group = rowname)) +
              geom_line() +
              geom_point() +
              facet_grid(~ assay) +
              labs(title = featureName, x = "sample", y = ylab)
          if (input$selColDetailPlot2!="none") {
              if (class(colData(pePlot)[[input$selColDetailPlot2]])=="factor") {
                p2 <- ggplot(pePlotDf, aes(x = colname, y = value,fill=pePlotDf[,input$selColDetailPlot2]))
                }
              } else {
              p2 <- ggplot(pePlotDf, aes(x = colname, y = value))
              }
          p2 <- p2 +
              geom_boxplot(outlier.shape = NA) +
              geom_point(position = position_jitter(width = .1), aes(shape = rowname)) +
              scale_shape_manual(values = 1:nrow(pePlotDf)) +
              labs(title = featureName, x = "sample", y = ylab)
          if (input$selVertDetailPlot2!="none"|input$selHorDetailPlot2!="none"){
              if (input$selVertDetailPlot2=="none") {
                  if (class(colData(pePlot)[[input$selHorDetailPlot2]])=="factor")
                      p2 <- p2 + facet_grid(~pePlotDf[,input$selHorDetailPlot2])
                  } else {
                  if (input$selHorDetailPlot2=="none"){
                      p2 <- p2 + facet_grid(pePlotDf[,input$selVertDetailPlot2]~.)
                      } else {
                      p2 <- p2 + facet_grid(pePlotDf[,input$selVertDetailPlot2]~pePlotDf[,input$selHorDetailPlot2])
                      }
                  }
          }
          return(list(p2,p1))
          }
        }
      }



#' Folder Upload Control
#'
#' Create a folder upload control that can be used to upload one or more filepaths pointing to folders. Strongly based on Shiny's File Upload Control.
#'
#' Whenever a folder upload completes, the corresponding input variable is set
#' to a character path.
#'
#' @family input elements
#'
#' @param inputId	The \code{input} slot that will be used to access the value.
#' @param label	Display label for the control, or \code{NULL} for no label.
#' @param value	Initial value.
#' @param width	The width of the input, e.g. \code{'400px'}, or \code{'100%'}; see \code{\link{validateCssUnit}}.
#' @param placeholder	A character string giving the user a hint as to what can be entered into the control. Internet Explorer 8 and 9 do not support this option.
#' @param multiple Whether the user should be allowed to select and upload
#'   multiple folders at once. \bold{Does not work on older browsers, including
#'   Internet Explorer 9 and earlier.}
#' @param accept A character vector of MIME types; gives the browser a hint of
#'   what kind of folders the server is expecting.
#' @param style The style attribute for the HTML tag. Used to hide/unhide the progress bar.
#'
#' @examples
#' ## Only run examples in interactive R sessions
#' if (interactive()) {
#'
#' ui <- fluidPage(
#'   sidebarLayout(
#'     sidebarPanel(
#'       fileInput("file1", "Choose CSV File",
#'         accept = c(
#'           "text/csv",
#'           "text/comma-separated-values,text/plain",
#'           ".csv")
#'         ),
#'       tags$hr(),
#'       checkboxInput("header", "Header", TRUE)
#'     ),
#'     mainPanel(
#'       tableOutput("contents")
#'     )
#'   )
#' )
#'
#' server <- function(input, output) {
#'   output$contents <- renderTable({
#'     # input$file1 will be NULL initially. After the user selects
#'     # and uploads a file, it will be a data frame with 'name',
#'     # 'size', 'type', and 'datapath' columns. The 'datapath'
#'     # column will contain the local filenames where the data can
#'     # be found.
#'     inFile <- input$file1
#'
#'     if (is.null(inFile))
#'       return(NULL)
#'
#'     read.csv(inFile$datapath, header = input$header)
#'   })
#' }
#'
#' shinyApp(ui, server)
#' }
#' @export
folderInput <- function(inputId, label, value = NA, multiple = FALSE, accept = NULL,
                        width = NULL, style="") {

  restoredValue <- restoreInput(id = inputId, default = NULL)

  # Catch potential edge case - ensure that it's either NULL or a data frame.
  if (!is.null(restoredValue) && !dir.exists(restoredValue)) {
    warning("Restored value for ", inputId, " has incorrect format.")
    restoredValue <- NULL
  }

  if (!is.null(restoredValue)) {
    restoredValue <- toJSON(restoredValue, strict_atomic = FALSE)
  }

  inputTag <- tags$input(
    id = inputId,
    name = inputId,
    type = "button",
    style = "display: none;",
    `data-restore` = restoredValue,
    class = "btn action-button"
    # webkitdirectory = NA,
    # directory = NA
  )

  if (multiple)
    inputTag$attribs$multiple <- "multiple"
  if (length(accept) > 0)
    inputTag$attribs$accept <- paste(accept, collapse=',')


  div(class = "form-group shiny-input-container",
      style = if (!is.null(width)) paste0("width: ", validateCssUnit(width), ";"),
      shiny:::`%AND%`(label, tags$label(label)),

      div(class = "input-group",
          tags$label(class = "input-group-btn",
                     span(id = paste(inputId, "_label", sep = ""), class = "btn btn-default btn-file",
                          "Browse...",
                          inputTag
                     )
          ),
          tags$input(type = "text", class = "form-control", value=value,
                     placeholder = "No folder selected", readonly = "readonly"
          )
      ),

      tags$div(
        id=paste(inputId, "_progress", sep=""),
        class="progress progress-striped active shiny-file-input-progress", style=style, #"visibility: visible;"
        tags$div(class="progress-bar", style="width: 100%;","Folder selected")
      )
  )
}





#'@export
fileInput <- function (inputId, label, multiple = FALSE, accept = NULL, width = NULL)
{
  restoredValue <- restoreInput(id = inputId, default = NULL)
  if (!is.null(restoredValue) && !is.data.frame(restoredValue)) {
    warning("Restored value for ", inputId, " has incorrect format.")
    restoredValue <- NULL
  }
  if (!is.null(restoredValue)) {
    restoredValue <- toJSON(restoredValue, strict_atomic = FALSE)
  }
  inputTag <- tags$input(id = inputId, name = inputId, type = "file",
                         style = "display: none;", `data-restore` = restoredValue)
  if (multiple)
    inputTag$attribs$multiple <- "multiple"
  if (length(accept) > 0)
    inputTag$attribs$accept <- paste(accept, collapse = ",")
  div(class = "form-group shiny-input-container", style = if (!is.null(width))
    paste0("width: ", validateCssUnit(width), ";"),
    shiny:::`%AND%`(label, tags$label(label)),
    div(class = "input-group", tags$label(class = "input-group-btn",
                                                               span(id = paste(inputId, "_label", sep = ""), class = "btn btn-default btn-file", "Browse...",
                                                                    inputTag)), tags$input(type = "text", class = "form-control",
                                                                                           placeholder = "No file selected", readonly = "readonly")),
    tags$div(id = paste(inputId, "_progress", sep = ""),
             class = "progress progress-striped active shiny-file-input-progress",
             tags$div(class = "progress-bar")))
}
