library(shiny)
library(DT)
library(shinyjs)
library(shinythemes)
source("utilities.R")
source("helper.R")


#runApp("pathToApp/App-MsqRob-devel")

####################################
###User Interface
####################################

shinyUI(fluidPage(theme = "MSqRob.css",

  #Use shinyjs package
  shinyjs::useShinyjs(),

  #CSS styles
  # tags$head(
  #   tags$style(HTML("
  # [disabled] {
  # color: pink;
  # }
  #   "))
  # ),
  tags$head(
         HTML(
           "
           <script>
           var socket_timeout_interval
           var n = 0
           $(document).on('shiny:connected', function(event) {
           socket_timeout_interval = setInterval(function(){
           Shiny.onInputChange('count', n++)
           }, 15000)
           });
           $(document).on('shiny:disconnected', function(event) {
           clearInterval(socket_timeout_interval)
           });
           </script>
           "
         )
         ),
    textOutput("keepAlive"),
############################################################################
#Navigation bar with 4 panel:Input, preprocessing, summarisation, quantification
############################################################################
   navbarPage("MSqRob Shiny App v 0.7.6", inverse=TRUE,


    ####################################
    #input tab
    ####################################

    tabPanel('Input',
    #sidebar
    sidebarLayout(
    #Project name
  	sidebarPanel(

  	  h3("Settings", class="MSqRob_topheader"),

      div(class="MSqRob_input_container",
  	      list(
  	  tags$label("Project name", `for`="project_name", class="MSqRob_label"),
  	  tags$button(id="button_project_name", tags$sup("[?]"), class="MSqRob_tooltip"),
  	  textInput("project_name", NULL, value = "project", width = '100%', placeholder = NULL),
  	hidden(helpText(id="tooltip_project_name",
  	       "Give your project a meaningful name.
  	       This name will be given to your results files.
  	       A time stamp will be automatically appended to name."))
  	      )
  	  ),

      div(class="MSqRob_input_container",
        	    list(
        	      tags$label("Input type", `for`="input_type", class="MSqRob_label"),
        	      tags$button(id="button_input_type", tags$sup("[?]"), class="MSqRob_tooltip"),
        	      selectInput("input_type", NULL, c("MaxQuant", "moFF", "mzTab"), width = '100%'),
        	      hidden(helpText(id="tooltip_input_type",
        	                      "Select the type of input.
        	                      "))
        	           )
        	     ),

  	#Peptides.txt files
  	div(class="MSqRob_input_container",
  	    list(
  	tags$label("Peptides file", `for`="peptides", class="MSqRob_label"),
  	tags$button(id="button_peptides", tags$sup("[?]"), class="MSqRob_tooltip"),
  	fileInput(inputId="peptides", label=NULL, multiple = FALSE, accept = NULL, width = NULL),
  	hidden(helpText(id="tooltip_peptides","Specify the location of the file that contains
  	                the peptide-specific intensities.
  	                When analyzing a MaxQuant shotgun proteomics experiment, this should the peptides.txt file.
  	                When using moFF, this file should start with \"peptide_summary_intensity\" and end with \".tab\".
			When using mzTab, this file should be a tab-delimited file with data summarized at the peptide level (\".tsv\" output file).
  	                When using Progenesis, this should be a \".csv\" file with data summarized at the peptide level.
			"))
  	    )
  	),

  	#Annotation file
  	div(class="MSqRob_input_container",
  	    list(
  	tags$label("Annotation file", `for`="annotation", class="MSqRob_label"),
  	tags$button(id="button_annotation", tags$sup("[?]"), class="MSqRob_tooltip"),
		fileInput(inputId="annotation", label=NULL, multiple = FALSE, accept = NULL, width = NULL),
		hidden(helpText(id="tooltip_annotation","Specify the location of your experimental annotation file."))
		)
		)
  	),

		#Main panel with number of output and plots
    mainPanel(width = 5,
            h3("Frequently asked questions", class="MSqRob_topheader"),
            htmlOutput("folderError"),
            div(class="MSqRob_h4_container",
            list(
            h4("What is an annotation file?"),
            tags$button(id="button_newExpAnnText",tags$sup("[show]"), class="MSqRob_tooltip"),
            actionButton(inputId="goAnnotation", label="Generate Annotation File!", class="MSqRob_button_space"),
            htmlOutput("downloadButtonDownloadAnnot"),
            hidden(helpText(id="tooltip_newExpAnnText",
              "An experimental annotation file contains the description of your experiment.
              Indeed, each mass spec run corresponds to e.g. a certain treatment, biological repeat, etc.
              This should be told to MSqRob via an Excel file or a tab delimited file wherein the first column contains all run names
              and the other columns contain all predictors of interest.
              Examples of experimental annotation files for the Francisella and CPTAC experiments can be found ",
              a("here", href="https://github.com/statOmics/MSqRobData/blob/master/inst/extdata/Francisella/label-free_Francisella_annotation.xlsx"),
              "and",
              a("here.", href="https://github.com/statOmics/MSqRobData/blob/master/inst/extdata/CPTAC/label-free_CPTAC_annotation.xlsx"),
              "Click the button to initialize an Excel file with a \"run\" column (works only if peptides.txt is already uploaded!).
              The annotation file will be saved in the output location.
              You still need to add other relevant columns (treatments, biological repeats, technical repeat, etc.) manually!"))
            )
            ),

            div(class="MSqRob_h4_container",
            list(
            h4("How do I cite MSqRob?"),
            tags$button(id="button_cite",tags$sup("[show]"), class="MSqRob_tooltip")
            )
            ),
            hidden(helpText(id="tooltip_cite",
            "MSqRob is free for you to use and completely open source.
            When making use of MSqRob, we would appreciate it if you could cite our two published articles.",
            br(),
            span("(1) The MSqRob algorithm: ", class="bold"),
            br(),
            "
            Goeminne, L. J. E., Gevaert, K., and Clement, L. (2016) Peptide-level Robust Ridge Regression Improves Estimation, Sensitivity,
            and Specificity in Data-dependent Quantitative Label-free Shotgun Proteomics. Molecular & Cellular Proteomics 15(2), pp 657-668.",
            br(),
            span("(2) The MSqRob GUI tutorial article:", class="bold"),
            br(),
            "
            Goeminne, L. J. E., Gevaert, K. and Clement, L. (2017).
            Experimental design and data-analysis in label-free quantitative LC/MS proteomics:
            A tutorial with MSqRob. Journal of Proteomics (in press).")),
            div(class="MSqRob_h4_container",
            list(
            h4("My question is not in this list!"),
            tags$button(id="button_notinlist",tags$sup("[show]"), class="MSqRob_tooltip")
            )
            ),
            hidden(helpText(id="tooltip_notinlist",
            "We are always ready to help you with any kind of issue that you might encounter!
            If for some reason, using MSqRob is hard or counter-intuitive to you,
            or you encounter some weird and unexpected results,
            please do not hesitate to contact us at", a("ludger.goeminne@vib-ugent.be.", href="mailto:ludger.goeminne@vib-ugent.be"),
            strong("User feedback is very important to us in order to improve MSqRob's user-friendliness.")
            ))
  )
   )
    )



    ############################
    #Preprocessing tab
    ###########################
    ,tabPanel('Preprocessing',
    sidebarLayout(
	#Sidebar with input
        sidebarPanel(

        h3("Settings", class="MSqRob_topheader"),

        div(class="MSqRob_input_container",
            list(
              tags$label("Group by", `for`="proteins", class="MSqRob_label"),
              tags$button(id="button_proteins", tags$sup("[?]"), class="MSqRob_tooltip"),
              htmlOutput("selectProteins"),
              hidden(helpText(id="tooltip_proteins","
                              Select the level on which the data should be grouped.
                              This is mostly the column that contains the protein identifiers (\"Proteins\" for MaxQuant data), as for a traditional shotgun experiment, one is mostly interested in which proteins are differentially abundant.
                              However, sometimes, one would for example like to do inference on the peptides.
                              In these more advanced cases, select the appropriate grouping level.
                              "))
              )
              ),

        h4("Transformation", class=c("MSqRob_sidebar")),

        div(class="MSqRob_input_container",
            list(
        checkboxInput("logtransform", label="Log-transform data", value=TRUE),
        tags$button(id="button_logtransform", tags$sup("[?]"), class="MSqRob_tooltip"),
        hidden(helpText(id="tooltip_logtransform",
                        "Leave this box ticked to log-transform the data.
                        Log-transformation is almost always performed to make the data less skewed.
                        Only when the data has already been log-transformed, this box can be unticked."
        ))
            )
        ),

        div(class="MSqRob_input_container",
            list(
          conditionalPanel(
            condition = "input.logtransform == true",
            tags$label("Base", `for`="log_base", class="MSqRob_label"),
            tags$button(id="button_log_base", tags$sup("[?]"), class="MSqRob_tooltip"),
            numericInput("log_base", NULL, value=2, min = 1, max = NA, step = NA, width = '100%'),
            hidden(helpText(id="tooltip_log_base",
                            "The base of the logarithm when log-transformation is performed.
                            Often, a base of 2 is chosen, because the results will then have interpretations
                            in terms of log2 fold changes.
                            "))
          )
          )
        ),

      	h4("Filtering", class="MSqRob_sidebar"),


        div(class="MSqRob_input_container",
            list(
        checkboxInput("smallestUniqueGroups", "Remove compromising protein groups", value=TRUE),
        tags$button(id="button_smallestUniqueGroups", tags$sup("[?]"), class="MSqRob_tooltip"),
        hidden(helpText(id="tooltip_smallestUniqueGroups",
                        "Remove protein groups for which any of its member proteins is present in a smaller protein group.
                        This might be done to remove any overlap of proteins in different protein groups.
                        "))
            )
        ),

        #Filter on features number of occurances
        div(class="MSqRob_input_container",
            list(
              tags$label("Minimum number of features", `for`="minIdentified", class="MSqRob_label"),
              tags$button(id="button_minIdentified", tags$sup("[?]"), class="MSqRob_tooltip"),
      	numericInput("minIdentified", label=NULL, value=2, min = 1, max = NA, step = 1, width = '100%'),
      	hidden(helpText(id="tooltip_minIdentified","
      	                The minimal number of times a feature sequence should be identified over all samples.
      	                feature sequences that are identified less than this number will be removed from the dataset.
      	                The default of 2 has the rationale that it is impossible to discern between the feature-specific effect and
      	                any other effects for a feature that has only been identified once.
      	                "))
      	)
        ),

        div(class="MSqRob_input_container",
            list(
              tags$label("Filter columns", `for`="filter", class="MSqRob_label"),
              tags$button(id="button_filter", tags$sup("[?]"), class="MSqRob_tooltip"),
      htmlOutput("selectFilters"),
      hidden(helpText(id="tooltip_filter","
                      Indicate the columns on which filtering should be done.
                      Features for which a \"+\" is present in these columns will be removed from the dataset.
                      This kind of filtering is typically done for common contaminants (e.g. operator's keratin)
                      and reversed sequences from the identificiation step that are still present in the data.
                      "))
            )
        ),
        h4("Normalization", class="MSqRob_sidebar"),

        div(class="MSqRob_input_container",
            list(
        tags$label("Normalization", `for`="normalisation", class="MSqRob_label"),
        tags$button(id="button_normalisation", tags$sup("[?]"), class="MSqRob_tooltip"),
        #selectInput("normalisation", NULL, c("quantiles","vsn","quantiles.robust", "center.median", "center.mean", "div.mean","div.median","diff.median" ,"max", "sum", "none"), width = '100%'), #"loess.affy" and "loess.pairs" left out on purpose because they remove everything with at least 1 NA!
        htmlOutput("selectNormalisation"),
        hidden(helpText(id="tooltip_normalisation",
                        "Select the type of normalisation from the dropdown menu.
                        Choose \"none\" if no normalisation should be performed
                        or if the data has already been normalised.
                        Note that with Progenesis data, we try to import the Normalized abundance.
                        Therefore, the default normalisation for Progenesis data is set to \"none\".
                        "))
            )
        ),
        actionButton(inputId="goNorm", label="Start Normalization!", class="MSqRob_button_space")

	),


	#Main panel with number of output and plots
        mainPanel(width = 5,

                  h3("Diagnostic plots", class="MSqRob_topheader"),

        strong('Number of features before preprocessing:'),textOutput('nfeaturesRaw',container = span),div(),
        strong('Number of features after preprocessing:'),textOutput('nfeaturesNormalized',container = span),div(),
        htmlOutput("selectColPlotNorm1"),


        div(class="MSqRob_h4_container",
            list(
              h4("Intensities after transformation"),
              tags$button(id="button_h4_int_transformation",tags$sup("[?]"), class="MSqRob_tooltip")
            )
        ),
        hidden(helpText(id="tooltip_h4_int_transformation","
                        A density plot showing the distribution of the feature intensities when only
                        the transformation is executed.
                        Transformation is included because a density plot of untransformed intensities is often uninformative
                        due to a strong skew to the right.
                        Brush and double-click on a selected area to zoom in.
                        Double click outside a selected area to zoom out.")),

        plotOutput('plotRaw',
                   click = "plotRaw_click",
                   dblclick = "plotRaw_dblclick",
                   brush = brushOpts(
                     id = "plotRaw_brush",
                     resetOnNew = TRUE
                   )
                   ),

        div(class="MSqRob_h4_container",
            list(
              h4("feature intensities after normalisation"),
              tags$button(id="button_h4_normalisation",tags$sup("[?]"), class="MSqRob_tooltip")
            )
        ),
        hidden(helpText(id="tooltip_h4_normalisation","
                        A density plot showing the distribution of the feature intensities
                        after execution of all preprocessing steps.
                        This allows you to evaluate the effect of the preprocessing.
                        Brush and double-click on a selected area to zoom in.
                        Double click outside a selected area to zoom out.")),

        plotOutput('plotNorm1',
                   click = "plotNorm1_click",
                   dblclick = "plotNorm1_dblclick",
                   brush = brushOpts(
                     id = "plotNorm1_brush",
                     resetOnNew = TRUE
                   )
                   ),

        div(class="MSqRob_h4_container",
            list(
              h4("MDS plot based on normalized feature intensities"),
              tags$button(id="button_h4_MDS_normalisation",tags$sup("[?]"), class="MSqRob_tooltip")
            )
        ),
        hidden(helpText(id="tooltip_h4_MDS_normalisation","A multidimensional scaling plot. This plot shows a two-dimensional scatterplot
                      so that distances on the plot approximate the typical log2 fold changes between the samples based on a pairwise comparison
                        of the 500 most different features.
                        Brush and double-click on a selected area to zoom in.
                        Double click outside a selected area to zoom out.")),


        div(checkboxInput("plotMDSPoints", "Plot MDS points", value=FALSE)),
        plotOutput('plotMDS',
				click = "plotMDS_click",
                                 dblclick = "plotMDS_dblclick",
                                 brush = brushOpts(
                                  id = "plotMDS_brush",
                                   resetOnNew = TRUE
                                 )
        )

				)
       )
)

###
#Summarisation tab
########
,tabPanel("Summarization",
sidebarLayout(
 sidebarPanel(
 h3("Settings", class="MSqRob_topheader"),
 h4("Summarisation", class="MSqRob_sidebar"),

 div(class="MSqRob_input_container",
     list(
 tags$label("Summarisation", `for`="Summarisation", class="MSqRob_label"),
 tags$button(id="button_summarisation", tags$sup("[?]"), class="MSqRob_tooltip"),
 selectInput("summarisation", NULL, c("none","robust","medpolish","mean","median","sum"), width = '100%'),
 htmlOutput("selectSummarisation"),
 hidden(helpText(id="tooltip_summarisation",
                 "Select the type of summarization from the dropdown menu.
                 ")),
              #    downloadButton("downloadProtSum", "Download protein intensities")
              withBusyIndicatorUI(actionButton(inputId="goSum", label="Start Summarisation!", class="MSqRob_button_space")) ,
              htmlOutput("downloadButtonProtSum")
     )
 )
 ),
 mainPanel(
   h3("Diagnostic Plots", class="MSqRob_topheader"),
   htmlOutput("selectColPlotProt"),
   div(class="MSqRob_h4_container",
       list(
         h4("MDS plot after full preprocessing"),
         tags$button(id="button_h4_MDS_full_preprocessing",tags$sup("[?]"), class="MSqRob_tooltip")
       )
   ),
   hidden(helpText(id="tooltip_h4_MDS_full_preprocessing","A multidimensional scaling plot. This plot shows a two-dimensional scatterplot
                 so that distances on the plot approximate the typical log2 fold changes between the samples based on a pairwise comparison
                   of the 500 most different features.
                   Brush and double-click on a selected area to zoom in.
                   Double click outside a selected area to zoom out.")),
                   div(checkboxInput("plotMDSPointsProt", "Plot MDS points", value=FALSE)),
   plotOutput('plotMDSProt',
   click = "plotMDSProt_click",
                            dblclick = "plotMDSProt_dblclick",
                            brush = brushOpts(
                             id = "plotMDSProt_brush",
                              resetOnNew = TRUE
                            )
   )
   )
)
)



    ###########################
    #Build Model
    ###########################
    ,tabPanel("Model",
      sidebarLayout(
       sidebarPanel(
           h3("Build Model", class="MSqRob_topheader"),
           h4("Following variables can be selected to build the model: "),
  	       h4(htmlOutput("selectFixed")),
           div(class="MSqRob_input_container",
               list(
                    tags$label("Design formula", class="MSqRob_label"),
                    tags$button(id="button_fixed", tags$sup("[?]"), class="MSqRob_tooltip"),
                    shiny::textInput("designformula", label=NULL,"~1"),
                    hidden(helpText(id="tooltip_fixed","
                            Make Formula"))
                    )
                ),
           div(class="MSqRob_input_container",
           list(
                tags$label("Ridge regression for fixed effects?", `for`="doRidge", class="MSqRob_label"),
                tags$button(id="button_doRidge", tags$sup("[?]"), class="MSqRob_tooltip"),
                radioButtons("doRidge", label=NULL,c("No"=0,"Yes" = 1)),
                hidden(helpText(id="tooltip_doRidge","
                    When \"Yes\" is selected the fixed effects are estimated using ridge regression. This shrinks the estimates with low evidence for differential abundance towards zero and improves the performance.
                    But, the method is computationally much more demanding.
                    We therefore suggest to switch ridge regression off \"No\" until you want to perform the final analysis. "
                    )
                    )
                )
            ),
            actionButton(inputId="fitModel", label="Fit Model!", class="MSqRob_button_space")
           ),
      mainPanel(
      fluidRow(column(width = 12, h3("Design Variables"),DT::dataTableOutput('annotationDataMatrix'))),
        h3("Visualize Design"),
        uiOutput('fitted_values_plot')
        )
      )
    )

#############
#Hypothesis tests
#############
,tabPanel("Inference",
    sidebarLayout(
      sidebarPanel(
      h3("Specify contrast"),
      h4("Following parameters can be used in contrasts for hypothesis tests: "),
      h4(htmlOutput("modelParams")),
      div(class="MSqRob_input_container",
          list(
               tags$label("Contrast", class="MSqRob_label"),
               tags$button(id="button_fixed", tags$sup("[?]"), class="MSqRob_tooltip"),
               shiny::textInput("contrast", label=NULL,""),
               hidden(helpText(id="tooltip_fixed","
                       Make Formula"))
               )
           ),
           div(class="MSqRob_input_container",
     		      list(
     		          tags$label("Significance level", `for`="alpha", class="MSqRob_label"),
     		          tags$button(id="button_alpha", tags$sup("[?]"), class="MSqRob_tooltip"),
     		          numericInput("alpha", label=NULL, value=.05, min = 0, max = 1, step = 0.01, width = NULL),
     		          hidden(helpText(id="tooltip_alpha","Select the significance level (alpha) at which the type I error needs to be performed.
     		                      Tests are traditionally performed at the 5% false discovery rate (FDR) level, but more stringent control (e.g. 1% FDR or even less) is sometimes adopted in experiments where false positives are highly unwanted (e.g. clinical settings).
     		                      The lower this level, the more stringent the cut-off and thus the less proteins that will be declared significant, but the higher the statistical certainty of the differential abundance of these proteins.
     		                      An FDR of 5% means that on average an expected 5% of the proteins that are called significant will be in fact false positives."))
     		           )),
            div(class="MSqRob_input_container",
      		      list(
      		          tags$label("Only significant features in table", `for`="sigOnly", class="MSqRob_label"),
      		          tags$button(id="button_sigOnly", tags$sup("[?]"), class="MSqRob_tooltip"),
      		          checkboxInput("sigOnly", label=NULL, value=TRUE),
      		          hidden(helpText(id="tooltip_sigOnly","If sigOnly is checked only the features with an adjusted p-value below the significance level are returned."))
      		           )
      		       )


      ),
      mainPanel(
          plotOutput("volcanoPlot",click = "plotVolcano_click",
                           dblclick = "plotVolcano_dblclick",
                           brush = brushOpts(
                             id = "plotVolcano_brush",
                             resetOnNew = TRUE)),
          fluidRow(column(width = 12, h3("Results table"),DT::dataTableOutput('table')))
      )
   )
)

###############
#Detail plots
###############

,tabPanel("DetailPlots",
sidebarLayout(
  sidebarPanel(
    htmlOutput("selectColDetailPlot2"),
    htmlOutput("selectHorizontalDetailPlot2"),
    htmlOutput("selectVerticalDetailPlot2")
    ),
  mainPanel(
    h4("Select one feature in the volcano plot or in the table of the inference tab to visualize the expression values"),
    uiOutput("detailPlots")
    )
  )
)


#close navbar, page, etc.
)))
