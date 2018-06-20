library(shinythemes)


appCSS <- 
  "#sigSelect ~ .selectize-control.single .selectize-dropdown [data-value='Miller Immune B/P'] { color: blue }
    #sigSelect ~ .selectize-control.single .selectize-dropdown [data-value='Miller Immune T/NK'] { color: blue }
    #sigSelect ~ .selectize-control.single .selectize-dropdown [data-value='Miller Immune M/D'] { color: blue }
    #sigSelect ~ .selectize-control.single .selectize-dropdown [data-value='Th1 Immune Response'] { color: blue }
    #sigSelect ~ .selectize-control.single .selectize-dropdown [data-value='Stoll IFg'] { color: blue }
    #sigSelect ~ .selectize-control.single .selectize-dropdown [data-value='Trastuzumab Immune'] { color: blue }
    #sigSelect ~ .selectize-control.single .selectize-dropdown [data-value='TCR Inhibition'] { color: blue }
    #sigSelect ~ .selectize-control.single .selectize-dropdown [data-value='Neutrophil Hypoxia'] { color: blue }
    #sigSelect ~ .selectize-control.single .selectize-dropdown [data-value='IRIS B cell'] { color: blue }
    #sigSelect ~ .selectize-control.single .selectize-dropdown [data-value='IRIS Dendritic Cell'] { color: blue }
    #sigSelect ~ .selectize-control.single .selectize-dropdown [data-value='IRIS Lymphoid'] { color: blue }
    #sigSelect ~ .selectize-control.single .selectize-dropdown [data-value='IRIS Monocyte'] { color: blue }
    #sigSelect ~ .selectize-control.single .selectize-dropdown [data-value='IRIS Myeloid'] { color: blue }
    #sigSelect ~ .selectize-control.single .selectize-dropdown [data-value='IRIS Neutrophil'] { color: blue }
    #sigSelect ~ .selectize-control.single .selectize-dropdown [data-value='IRIS NK cell'] { color: blue }
    #sigSelect ~ .selectize-control.single .selectize-dropdown [data-value='IRIS T Cell'] { color: blue }
    #sigSelect ~ .selectize-control.single .selectize-dropdown [data-value='List of 4 or more genes'] { color: green }
    #sigSelect ~ .selectize-control.single .selectize-dropdown [data-value='Gene Correlation'] { color: green }"

# Define UI for application 
shinyUI(fluidPage(
  theme=shinytheme("flatly"),
  # Application title
  titlePanel("Integrative gene expression analysis for Breast Cancer",windowTitle = "Viewer v1.0e"),
  h5("Version 1.0e- Now compatible with > R 3.4"),
  sidebarLayout(
    sidebarPanel(

    tags$head(tags$style(HTML(appCSS))),
    uiOutput("signature"),
    conditionalPanel(
      condition = "input.seltab != 'HM'",
      selectizeInput("Gsubtype","Molecular Subtype",
                   choices=c("Basal","Her2","LumA","LumB","Normal","None"),
                   selected=c("Basal","Her2","LumA","LumB","Normal"),
                   multiple=TRUE)
    ),
    conditionalPanel(
      condition = "input.sigSelect == 'List of 4 or more genes'",
      textInput("genes","Gene list", 
                value = "TP53 BCL2 ESR1 TAF1A IDH1")
      ),

    conditionalPanel(
      condition = "input.sigSelect == 'Gene Correlation'",
      textInput("corGene","Gene", 
                value = "TP53"),
      numericInput("corThres","Correlation threshold", 
                value = 0.3, min=0, max=1, step=0.05)
    ),

# conditional display of heatmap controls only if heatmap tab is selected
conditionalPanel(
condition = "input.seltab == 'HM'",   
            h4("Heatmap Presentation Selections"),
    selectizeInput("colors", "Color scheme:",
                list("Green/Red" = "greenred",
                     "Blue/Red" = "bluered",
                     "Blue/Yellow" = "blueyellow",
                     "Matlab" = "matlab"),
                selected="blueyellow"),
            

    

      uiOutput("HM_order"),
      uiOutput("HM_subsort"),
      
       selectInput("colourbar", "Select Colour Bars to Display", selected="No",
                  list("Yes" = "Yes",
						"No" = "No"
                       )),
					   
      conditionalPanel(
        condition = "input.order == 'Clustering'",
        radioButtons("radioDist","Distance metric",
                     c("Euclidean","Correlation"),
                     selected="Euclidean")
		 
      ),
        
		conditionalPanel(
			condition = "input.colourbar == 'Yes'",
				checkboxInput("grade.show",label="Grade",value=TRUE),
				checkboxInput("subtype.show",label="Subtypes",value=TRUE),
				checkboxInput("er.show",label="ER IHC Status",value=TRUE),
			  checkboxInput("gene.show",label="Show Ki67 & ESR Expression", value=FALSE),
				checkboxInput("hoc.show",label="Hallmarks of Cancer",value=FALSE),
				checkboxInput("pc.show",label="Principal Components",value=TRUE),
				checkboxInput("ic.show",label="Independent Components",value=TRUE) 
		),
    
      downloadButton('downloadHMPlot', 'Save Plot'),
  
    
    
		width=5), # end condition panel for tab Heatmap
            

        
        # conditional display of heatmap controls only if survival plot tab is selected
conditionalPanel(
condition = "input.seltab == 'SP'",   

      

 br(),                              
        selectizeInput("radioTFsurv","Use an individual gene instead of a gene set?",
         list("Gene Set","Individual Gene"),
         selected="Gene Set"),                      
conditionalPanel(
condition = "input.radioTFsurv == 'Individual Gene'",

 textInput("genesTocompareS","Select Individual Gene ", 
                value = "BCL2")          
),
                   
conditionalPanel(
			condition="input.radioTFsurv == 'Gene Set'",
			selectInput("order.surv", "Gene set statistic used", selected="HMC",
                  list("PC1" = "PC1",
                       "PC2" = "PC2",
                       "PC3" = "PC3",
                       "PC4" = "PC4",
                       "Centroid" = "Centroid",
                       "IC1" = "IC1",
                       "IC2" = "IC2",
                       "IC3" = "IC3",
                       "IC4" = "IC4"
                       ))
					   ),                       
   
  


 selectizeInput("survplottype","Divide survival plots by...",
                selected="2 groups- Median",
                list("2 groups- Median","2 groups- Quantile","2 groups- Top/Bottom")
                   ),


      conditionalPanel(
        condition = "input.survplottype == '2 groups- Quantile'",
        sliderInput("survSlider","Quantile",0,1,0.5)
      ),
      conditionalPanel(
        condition = "input.survplottype == '2 groups- Top/Bottom'",
        sliderInput("survSlider2","Quantile",0,1,c(0.25,0.75))
),
  downloadButton('dlSurvPlot', 'Save Plot'),

 h4("Data subsets"),       
  selectizeInput("radioGd","Grade",selected="All",
                   list("All","Grade 1","Grade 2","Grade 3")
                   ),
                   
  selectizeInput("radioER","ER IHC status",selected="All",
                   list("All","ER+","ER-")
                   ),
  selectizeInput("radioNZ","Isolate NZ data?",selected="All",
                 list("All","NZ Only","Non-NZ")),
 
  selectInput("survadvCont","Enable Advanced Controls?",selected="FALSE",
              list("Yes" = "TRUE",
                   "No" = "FALSE"
              )),
  conditionalPanel(
    condition="input.survadvCont == 'TRUE'",
    selectizeInput("radioRS","Reverse gene or gene set",
                   selected="No",
                  list("Yes","No")
               ),
    
    selectInput("input.coxph", "Run CoxPH model on geneset as continuous", selected="FALSE",selectize=FALSE,
                list("Yes" = "TRUE",
                     "No" = "FALSE"
                ))
    
    )
                   

 ),# end condition panel for tab SP  
 
conditionalPanel(
  condition = "input.seltab == 'AP'",
  sliderInput("AuCyear","Calculate AuC Between ", value=c(0,5), min=0,max=10,step=0.5),
  uiOutput("Add_Who"),
  actionButton("twoBest","Find 2 Best Predictors",class="btn-success"),
  actionButton("threeBest", "Find 3 Best Predictors",class="btn-info"),
  actionButton("fourBest", "Find 4 Best Predictors",class="btn-primary"),
  br(),
  downloadButton('downloadAuCPlot', 'Save Plot')
              
),
         # conditional display of heatmap controls only if GTPC plot tab is selected
conditionalPanel(
condition = "input.seltab == 'GTPC'",   

        
              selectInput("order.gtpc", "Statistic from gene set to be plotted", selected="HM1",
                  list("PC1" = "PC1",
                       "PC2" = "PC2",
                       "PC3" = "PC3",
                       "PC4" = "PC4",
                       "Centroid" = "Centroid",
			                  "IC1" = "IC1",
                       "IC2" = "IC2",
                       "IC3" = "IC3",
                       "IC4" = "IC4"
                       )),
                 
  textInput("genesTocompare","Individual gene that will be compared to the gene set ", 
                value = "ESR1"),
               
                
  selectizeInput("radioRX","Treatment",
                   list("All tumours","No adjuvant Rx","Endocrine Rx no chemoRx"),
                   selected="All tumours"),
                   
                  
  selectizeInput("radioTF","Type of transformation made to \n both gene and gene set data",
                   list("Z-transform","Rank"),
                   selected="Z-transform"),
                  
  selectizeInput("radioRA","Reverse y axis",
                   list("Yes","No"),
                   selected="No"),

 h4("Data subsets"),       

  
    selectizeInput("radiogtpcGd","Grade",selected="All",
                   list("All","Grade 1","Grade 2","Grade 3")
                   )
 
),# end condition panel for tab GTPC  

conditionalPanel(
  condition = "input.seltab == 'GL'",  
  
  selectizeInput("cormethod","Select Correlation Method",
                 list("Pearson"="pearson","Kendall"="kendall","Spearman"="spearman"),
                 selected="spearman"),
  
  selectizeInput("corABS","Absolute Value?",
                 list("Yes","No"),
                 selected="Yes")
 
    ),
# This needs to be doubled or conditional panel below is showing up
conditionalPanel(
  condition = "input.seltab == 'HOC'",
  
  selectizeInput("HOCcormethod","Select Correlation Method",
                 list("Pearson"="pearson","Kendall"="kendall","Spearman"="spearman"),
                 selected="spearman"),
  
  selectizeInput("HOCcorABS","Absolute Value?",
                 list("Yes","No"),
                 selected="No"),
  downloadButton('dlHOCcolors', 'Save Plot')
),
# This panel should never be visible unless debugging is being used.
conditionalPanel(
  condition = "input.seltab == 'start' && input.hideTab == 'No'", 
  selectInput("hideTab","Hide this Tab?",c("Yes","No"),selected="No",selectize=FALSE)
),

width=2),



   # Dynamic UI logic in Server.R
  mainPanel(
        uiOutput("myTabs")
  ) 
)))
