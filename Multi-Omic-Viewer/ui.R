# To do:
# Add correlation of RNAseq to every gene meth site
# Add cg location to table and plot
#

library(shiny)
library(DT)
library(shinydashboard)
library(shinycssloaders)
# Split out Shiny Pieces to keep editing and updating easier
header<-dashboardHeader(title="Multi Omics Viewer")

sidebar<-dashboardSidebar(
  sidebarMenu(
    menuItem("Tutorial",tabName="tutorial",icon=icon("support")),
    menuItem("Data Set",tabName="datasetSelect",icon=icon("dashboard"),selected=TRUE),
    menuItem("Select Gene",tabName="GeneTable",icon=icon("share-alt")),
    menuItem("Single Gene Plots",tabName="Plots1",icon=icon("area-chart")),
    menuItem("Variant Plots",tabName="Plots2",icon=icon("modx")),
    menuItem("Say Hello :-)", href="#",newtab=F),
    menuItem("LinkedIn", href="https://nz.linkedin.com/in/nsknowlton",icon=icon("linkedin"))

  )
)

body<-dashboardBody(
   tabItems(
    tabItem(tabName="datasetSelect",
           h2("Select Your Dataset"),
          fluidRow(
          box(selectInput('database1','Select From',choices=list("Metastatic Skin Cutaneous Melanoma"='SKCM',"Breast Invasive Carcinoma"='BRCA'),selected='SKCM'),solidHeader = T, width=6, status="primary"),
          box(selectInput('variant','Include Variant Information?', choices=list("Yes- Present/Absent"='Yes',"Yes- Factors"='Yes2',"No"='No'),selected='No'),solidHeader = T, width = 6, status="primary")
          ),
          fluidRow(
            box(selectInput('cna.input','Select Copy Number Flavour', choices=list("Gistic"='Gistic',"Gene Segment Variants"='CNA'),selected='CNA'),solidHeader = T, width=6,status="primary"),
            box(selectInput('lm.type','Select the association model', choices=list("Relative Importance"='relimp',"Correlation"='cor'),selected='relimp'),solidHeader = T, width=6, status="primary")
          ),
          fluidRow(
            box(title="Database Information",solidHeader=T,collapsed=T,status="info",textOutput("db.out"))
           )
    ),
    tabItem(tabName="tutorial",
            includeMarkdown("instructions.md")

    ),
    tabItem(tabName="GeneTable",
            h4("Gene Table"),
          fluidRow(
            box(
          title = "Help", width = 8, solidHeader = TRUE, status = "warning", collapsed = T, collapsible = T,
          "The table below displays the Amount of Variation Explained in % by a Linear model of RNAseq = Copy Number + Methylation + Genetic Variant. Each of these columns can be sorted and a gene of interest can be entered #in the search box. If there were more than one Methylation sites in a given gene, only the site that explained the most variance is below. The other sites can be queried on the Single Gene Plots Tab"
            ),
            box(title="Data Status",textOutput("Test"),status="info",solidHeader=T,collapsed=T, width=4)
          ),

        box(title="Select Gene of Interest", solidHeader=TRUE, width=8,withSpinner(DT::dataTableOutput("resultsGENE"),type=6))


        ),
    tabItem(tabName = "Plots2",
            h4("Variant vs Copy Number Plots"),
            box(uiOutput("VariantPlots1"),width=6),
            box(uiOutput("VariantPlots2"),width=6),
            box(title="RNAseq vs Copy Number", width=6, height=500, plotOutput("ScatterVar1"),status="danger",solidHeader = TRUE),
            box(title="RNAseq vs Copy Number", width=6, height=500, plotOutput("ScatterVar2"),status="primary",solidHeader = TRUE),
            box(title="Variance Explained by Mutation Type",withSpinner(DT::dataTableOutput("VarTable"),type=6), solidHeader=TRUE, width=8)
            ),

    tabItem(tabName="Plots1",
            h4("Single Gene Plots"),
            box(title="Select Methylation Site",withSpinner(DT::dataTableOutput("methTable"),type=6),collapsible = F,collapsed = T,width=6),
           fluidRow(
            box(title="Methylation vs Copy Number coloured by Expression",plotOutput("Scatter", height = 800, width = "auto"),status="primary",solidHeader = T,collapsible=T,width=6)
             ),
           fluidRow(
              box(title="Methylation vs Expression",width=4,plotOutput("ScatterMethGG", height = 400, width = 400),status="danger",solidHeader = T),
              box(title="Expression vs Copy Number", width=4,plotOutput("ScatterCN", height = 400, width = 400), status="warning",solidHeader=T)
            )
          )

    )
)




dashboardPage(header,sidebar,body,skin="green")
