# Load packages

if (!require("pacman")) install.packages("pacman")
pacman::p_load(shiny,shinydashboard, shinycssloaders, shinyjs, sqldf, matrixStats, gplots, ggplot2, colorRamps, rms, devtools, fields, DT)


#
jsResetCode <- "shinyjs.refresh = function() { location.reload(); }"
# Split out Shiny Pieces to keep editing and updating easier
header<-dashboardHeader(title="Multi Omics Viewer")

sidebar<-dashboardSidebar(
  sidebarMenuOutput("menu")
  )

body<-dashboardBody(
  useShinyjs(),
  extendShinyjs(text = jsResetCode, functions = "refresh"),
  uiOutput("bodydyn")
  )

dashboardPage(header,sidebar,body,skin="green")
