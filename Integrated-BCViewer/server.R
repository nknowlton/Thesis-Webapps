library(shiny)
library(shinycssloaders)
if (packageVersion('DT') < '0.1.56') devtools::install_github('rstudio/DT')
if (!require("pacman")) install.packages("pacman")
pacman::p_load(shinydashboard, sqldf, stringi, matrixStats, gplots, ggplot2, colorRamps,rms,devtools,fields,survivalROC,data.table,reshape2,devtools,DT,fastICA,shinythemes,Rmisc,rmarkdown,shinycssloaders)

options(contrasts=c("contr.treatment","contr.treatment")) # to allow Design to work properly with ordered factors


dataObject<-c("GeneData.sqlite")
db<-dbConnect(SQLite(),dbname="GeneData.sqlite")

shinyServer(function(input, output,session) {

  # When shiny loads a blank tab is also loaded to prevent data from being loaded.
  # This tab is immediatly redirects to the HM tab by the dynamic UI element `myTabs'

  	# Reset the color bar selections when reselecting "no"
  observeEvent(input$colourbar,{
    
    if(input$colourbar == 'No') {
     
      updateCheckboxInput(session,"grade.show",value=TRUE)
      updateCheckboxInput(session,"subtype.show",value=TRUE)
      updateCheckboxInput(session,"er.show",value=TRUE)
      updateCheckboxInput(session,"gene.show",value=FALSE)
      updateCheckboxInput(session,"hoc.show",value=FALSE)
      updateCheckboxInput(session,"pc.show",value=TRUE)
      updateCheckboxInput(session,"ic.show",value=TRUE) 
      updateRadioButtons(session,"figure1",selected='No')
      
      
    }
    })
  #Reset Survival plot tab when advanced options are no longer selected
  observeEvent(input$survadvCont,{
    
      updateSelectizeInput(session,"radioRS",selected="No")
      updateSelectizeInput(session,"input.coxph",selected="FALSE")
  })
  # Reset AuC choice to PC1, PC2 as this is a "safe" choice
  observeEvent(input$sigSelect,{
    
    updateSelectizeInput(session,"auc.choice",selected=c("PC1","PC2"))
  })
  # UI element that determines which terms can be selected on AuC tab
  output$Add_Who<-renderUI({
    
    choose.me<-auc.choice()
   
    selectizeInput("auc.choice","Add Which Together",
                   choices=choose.me,
                   selected=c("PC1","PC2"),
                   multiple=TRUE)
    
  })
  #UI element that allows the choice of signature
  output$signature<-renderUI({
    sigs.table<-sig.data()
    choose.me<-c(unique(sigs.table[sigs.table$Tissue=="Breast","Signature"]),"List of 4 or more genes","Gene Correlation")
    #selectInput("sigSelect","Gene Set", choices=choose.me,multiple=FALSE)
    
    selectizeInput("sigSelect","Signature or Gene Set", choices=choose.me,multiple=FALSE)
   
  })
  
  
  #UI element for Ordering of Heatmap
  HMsorting<-reactive({
     choose.me<-c("Cohort","NZ Cohort")
    
    if(input$grade.show) choose.me<-c(choose.me,"Grade")
    if(input$subtype.show) choose.me<-c(choose.me,"Subtype")
    if(input$er.show) choose.me<-c(choose.me,"ER IHC Status")
    if(input$gene.show) choose.me<-c(choose.me,"Ki67 Expression", "ESR1 Expression")
    if(input$hoc.show) choose.me<-c(choose.me,"Proliferation","Growth Suppression","Immune","Immortality","Inflammation","Metastasis","Angiogenesis","Genome Instability","Cell Death","Cell Energy")
    if(input$pc.show) choose.me<-c(choose.me,"Heatmap 1st PC","Heatmap 2nd PC","Heatmap 3rd PC","Heatmap 4th PC")
    choose.me<-c(choose.me,"Heatmap Centroid")
    if(input$ic.show) choose.me<-c(choose.me,"Heatmap 1st IC","Heatmap 2nd IC","Heatmap 3rd IC","Heatmap 4th IC")
    choose.me<-c(choose.me,"Clustering")
    return(choose.me)
  })
  #UI element for selecting HM order
  output$HM_order<-renderUI({
   
    selectizeInput("order", "Heatmap column order:", choices=HMsorting(),selected="Heatmap Centroid")
    
  })
  output$HM_subsort<-renderUI({
    if(is.null(input$order)) return(NULL)
    #print(input$order)
    if(input$order == "Cohort" | input$order == "NZ Cohort" | input$order ==  "Grade" | input$order ==  "Subtype" | input$order ==  "ER IHC Status" | input$order=="Ki67 Expression" | input$order == "ESR1 Expression") {
    remove<-c("Clustering")
    
    if(input$order == "NZ Cohort") {
      remove<-c(remove,"Cohort")
    } else {
       
       remove<-c(remove,"NZ Cohort")
    }
    
    if(input$order== "Cohort") remove<-c(remove,"NZ Cohort")
    choose.me<-setdiff(HMsorting(),remove)
    #print(choose.me)
    selectizeInput("subsort","Heatmap subsort order:",choices=choose.me,selected=input$order)} # providing a selection choice screws up plotng here
  })
  
  NZcohort<-reactive({
    coh<-as.numeric(as.factor(clinical.in()$Cohort..))
    NZc<-ifelse(coh<=2,1,2)
    return(NZc)
  })
  
  # Dynamic output of Tabs for the Top selection
  output$myTabs<-renderUI({
     myTabs<- tabsetPanel(id="seltab",
                          #tabPanel(h4("Tutorial"),includeMarkdown("instructions.md")),
                          tabPanel(h4("Expression Heatmap"),br(), withSpinner(plotOutput("expMap", height = "2000px", width = "1920px"),type=6), value="HM"),
                          tabPanel(h4("Survival Plot"), h4(textOutput("survTitle")), plotOutput("survPlot", height = "1000px", width = "1200px"), value="SP"),
                          tabPanel(h4("AuC Plots"),withSpinner(plotOutput("aucPlot",height="1000px", width="600px"),type=6),value="AP"),
                          tabPanel(h4("HOC Associations"), withSpinner(DT::dataTableOutput("HOCtabtable"),type=6),br(), plotOutput("HOCtabPlot"),value="HOC"),
                          tabPanel(h4("Gene -vs- Gene Set"), h4(textOutput("gtpcTitle")), plotOutput("gtpcPlot", height = "1100px", width = "1920px"), value="GTPC"),
                          tabPanel(h4("Gene Associations"), h4(textOutput("glxx")), withSpinner(DT::dataTableOutput("glDisplay"),type=6), value="GL"),
                          selected="HM"
                          )
                         # }
    
   return(myTabs)
    
  })
  # Read in data from signature database-
  # Currently reads in everything, but can be selective if needed in future
  sig.data<-reactive({
    db2<-dbConnect(SQLite(),dbname="GeneSignatures.sqlite",synchronous=NULL)
    sigs.table<-dbReadTable(db2,"GeneSigs")
  })
  # Pull in clinical data for multiple function use
  clinical.in<-reactive({
    dbGetQuery(conn=db, "SELECT * FROM clinDat")
})
  # Global rows for plotting based upon subtype
  global.subtype<-reactive({
    clinical.in()
    kp<-rep(F,nrow(clinical.in()))
    len.store<-length(input$Gsubtype)
    if(len.store>0) {
      for(i in 1:len.store){
        kp.tmp<-clinical.in()$Subtype==input$Gsubtype[i]
        kp<-kp+kp.tmp
      }}
    return(kp)
  })
  clinicalSubtype<-reactive({
    clinical.in()[as.logical(global.subtype()),]
  })
  # Pulls in data from two genes
  gene.show.data<-reactive({
    genes<-c("MKI67","ESR1")
    x<-dbGetQuery(conn=db, paste0("SELECT * FROM expdat WHERE row_names IN ('",paste0(genes,collapse="', '"),"')")) # Only pull in the data we need. Much faster than loading everything
    x<-as.data.frame(x)
    rownames(x)<-x$row_names
    x<-x[,-1]
    
    z<-x[apply(is.na(x),1,sum)==0,]
    data.out<-t(apply(z,1,function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T)))
   
  })
  
  # AuC tab button that generates best 2 terms for AuC model- exaustive search
  observeEvent(input$twoBest,{
    plot.dat<-auc.data()
    
    out<-auc.search(plot.dat$plot.matrix,plot.dat$kp,plot.dat$dmfsTime,plot.dat$dmfsEvent,input$AuCyear[2],nvars=2)
    
    updateSelectizeInput(session,"auc.choice",selected=out)
  })
  # AuC tab button that generates best 3 terms for AuC model- exaustive search
  observeEvent(input$threeBest,{
    
    plot.dat<-auc.data()
    
    out<-auc.search(plot.dat$plot.matrix,plot.dat$kp,plot.dat$dmfsTime,plot.dat$dmfsEvent,input$AuCyear[2],nvars=3)
    
    updateSelectizeInput(session,"auc.choice",selected=out)
  })
  # AuC tab button that generates best 4 terms for AuC model- exaustive search
  observeEvent(input$fourBest,{
    
    plot.dat<-auc.data()
    
    out<-auc.search(plot.dat$plot.matrix,plot.dat$kp,plot.dat$dmfsTime,plot.dat$dmfsEvent,input$AuCyear[2],nvars=4)
    
    updateSelectizeInput(session,"auc.choice",selected=out)
  })
  signatures.in<-reactive({
    sigs.in<-as.data.frame(fread("signature_fits.csv"))
  })
  # Gets and plots AuC data.  Repeats sections of Auc.choice but also grabs all avaialbe data
  auc.data<-reactive({

    plot.matrix<-eDat()$PcIcOut
    plot.matrix<-cbind(plot.matrix,signatures.in()[,names(signatures.in()) %in% input$sigSelect])
    
    
    if(dim(plot.matrix)[2] == 9)    colnames(plot.matrix)<-c("PC1","PC2","PC3","PC4","Centroid","IC1","IC2","IC3","IC4")
    if(dim(plot.matrix)[2] == 10)   colnames(plot.matrix)<-c("PC1","PC2","PC3","PC4","Centroid","IC1","IC2","IC3","IC4","Score")
    
    # Note that what we are using is Z-transformed data from the outset
    gdat<-eDat()$gdat
    
    dmfsEvent<-clinical.in()$DMFS.EVENT..defined.as.distant.metastasis.or.from.breast.cancer.
    dmfsTime<-clinical.in()$DMFS.TIME 
    kp.out<-global.subtype() & dmfsTime >= input$AuCyear[1]  # Create a vector how who to keep
    out<-list("kp"=kp.out,"plot.matrix"=plot.matrix,"dmfsEvent"=dmfsEvent,"dmfsTime"=dmfsTime)
    return(out)
  })
  
  # Dynamic choice for plotting on AuC tab
  auc.choice<-reactive({
    sigs.in<-as.data.frame(fread("signature_fits.csv"))
    score<-sum(names(sigs.in) %in% input$sigSelect)
    if(score == 0)    names<-c("PC1","PC2","PC3","PC4","Centroid","IC1","IC2","IC3","IC4")
    if(score == 1)   names<-c("PC1","PC2","PC3","PC4","Centroid","IC1","IC2","IC3","IC4","Score")
    return(names)
  })
  # Calculate the Hallmarks of Cancer from the hoc.cent function which returns centroids
  hocCalculate<-reactive({
    
      
        load("hoc symbols.RData") # load HOC data from R datafile
        hoc.symbols<-hoc.symbols[,-1]
        genes<-unique(unlist(hoc.symbols)) # find out who we need to look for
        x<-dbGetQuery(conn=db, paste0("SELECT * FROM expdat WHERE row_names IN ('",paste0(genes,collapse="', '"),"')")) # Only pull in the data we need. Much faster than loading everything
        x<-as.data.frame(x)
        rownames(x)<-x$row_names
        x<-x[,-1]
        
        
        output<-matrix(NA,nrow=10,ncol=dim(x)[2])
        
        for(i in 1:10){
          loc<-match(as.character(hoc.symbols[,i]),rownames(x))
          output[i,]<-colMeans(x[loc,],na.rm=TRUE)}
        
        rownames(output)<-colnames(hoc.symbols)
        colnames(output)<-colnames(x)
        
        outputZ<-t(apply(output,1,function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T)))
        # Keep Z scores from exceeding 3
        # Makes the graphs look nicer
        outputZ[outputZ < -3]<- -3
        outputZ[outputZ > 3]<- 3
    
    return(outputZ)
  })
  #Function that returns the gene symbols used throughout script.  Default values are used for inital startup
  # and then are overwritten otherwise
  genes <- reactive({
    sigs.table<-sig.data()
    g<-c("TP53","ESR1","IGF1","PTN") # Prevents an error on start-up
    
    if(input$sigSelect=="List of 4 or more genes"){
      g<-strsplit(gsub(":"," ",gsub(","," ",gsub(";"," ",input$genes,fixed=T),
                                   fixed=T),fixed=T)," ")[[1]]
      g<-g[g!=""]
      g<-toupper(g)
    }
   
    if(input$sigSelect=="Gene Correlation"){ 
      back<-corMat()
      p<-rownames(back$x)[order(abs(back$corDat),decreasing=T)[sort(abs(back$corDat),decreasing=T)>=input$corThres]]
      g<-unlist(p)
      if(length(g)>300) g<-g[1:300] # only Keep top 300 matches
      if(length(g) == 0) g<-toupper(input$corGene)
    }
    
    if(input$sigSelect != "List of 4 or more genes" & input$sigSelect != "Gene Correlation"){
      g<-sigs.table[sigs.table$Signature==input$sigSelect,"Gene"]
      
    }

	  
    return(g)
  })
  
  #Create Correlation gene list
  corMat<- reactive({
    
    g<-toupper(input$corGene)
    p.gene<-dbGetQuery(conn=db, paste0("SELECT * FROM expdat WHERE row_names ='",g,"'")) # Pull in selected single gene
    x<-dbGetQuery(conn=db, "SELECT * FROM expdat") # Pull in entire dataset
    x<-as.data.frame(x)
    rownames(x)<-x[,1]
    x<-x[,-1]
    
    corDat<-cor(as.numeric(p.gene[1,-1]),t(x)) # 
    out<-list("corDat" = corDat,"x"=x)
    return(out)
  })
  
  # Pull in Data from a list of genes and compute IC PC and Centroid 
  eDat <- reactive({
    
    g<-genes() 
    x<-dbGetQuery(conn=db, paste0("SELECT * FROM expdat WHERE row_names IN ('",paste0(g,collapse="', '"),"')"))
    #<-as.data.frame(x) #rownames now added automatically from matrix
    rownames(x)<-x[,1]
    x<-x[,-1] 
   
    z<-x[apply(is.na(x),1,sum)==0,]
    gdat<-t(apply(z,1,function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T)))
    PcIcOut<-pcicout(gdat)
    out<-list("gdat"=gdat,"PcIcOut"=PcIcOut)
    return(out) # the function eDat returns a matrix of data for all tumours but only the genes you are focusing on
  })

  eDatSubtype <- reactive({
    #print("eDat Subtype")
    g<-genes() 
    x<-dbGetQuery(conn=db, paste0("SELECT * FROM expdat WHERE row_names IN ('",paste0(g,collapse="', '"),"')"))
    x<-as.data.frame(x)
    rownames(x)<-x[,1]
    x<-x[,-1]
    
    z<-x[apply(is.na(x),1,sum)==0,]
    z<-z[,as.logical(global.subtype())]
    gdat<-t(apply(z,1,function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T)))
    PcIcOut<-pcicout(gdat)
    out<-list("gdat"=gdat,"PcIcOut"=PcIcOut)
    return(out) # the function eDat returns a matrix of data for all tumours but only the genes you are focusing on
  })
  
  # Make heatmap data prep reactive so UI can grab elements before plot is submitted
  heatmapProcessing<-reactive({
    
    if(input$hoc.show == TRUE){ # Only load data when requested
      hoc<-hocCalculate() 
      }
      
   
        # Note that what we are using is Z-transformed data from the outset
    ss<-eDat()$PcIcOut[,1]
    ss2<-eDat()$PcIcOut[,2]
    ss3<-eDat()$PcIcOut[,3]
    ss4<-eDat()$PcIcOut[,4]
    ss_centroid<-eDat()$PcIcOut[,5]
    ss5<-eDat()$PcIcOut[,6]
    ss6<-eDat()$PcIcOut[,7]
    ss7<-eDat()$PcIcOut[,8]
    ss8<-eDat()$PcIcOut[,9]
    
    gdat<-eDat()$gdat # Can't left assign into reactive value, so must make a copy to pass to plotting function
    
    # Keep Z scores from exceeding 3
    # Makes the graphs look nicer
    gdat[gdat< -3]<- -3
    gdat[gdat>3]<-3
    
    # Generate colour ramps based upon user input
    if(input$colors=="greenred"){
      ccCOL<-greenred(length(ss))
      ccx<-greenred(50)
      GCOL<-greenred(length(gene.show.data()[1,]))
    } else if(input$colors=="bluered"){
      ccCOL<-bluered(length(ss))
      ccx<-bluered(50)
      GCOL<-bluered(length(gene.show.data()[1,]))
    } else if(input$colors=="blueyellow"){
      ccCOL<-blue2yellow(length(ss))
      ccx<-blue2yellow(50)
      GCOL<-blue2yellow(length(gene.show.data()[1,]))
    } else if(input$colors=="matlab"){
      ccCOL<-matlab.like(length(ss8))
      ccx<-matlab.like(50)
      GCOL<-matlab.like(length(gene.show.data()[1,]))
    }
    cc<-rbind(ccCOL[rank(ss)],
      ccCOL[rank(ss2)],
      ccCOL[rank(ss3)],
      ccCOL[rank(ss4)],
      ccCOL[rank(ss_centroid)],
      ccCOL[rank(ss5)],
      ccCOL[rank(ss6)],
      ccCOL[rank(ss7)],
      ccCOL[rank(ss8)])
    rownames(cc)<-c("Heatmap 1st PC","Heatmap 2nd PC","Heatmap 3rd PC","Heatmap 4th PC","Heatmap Centroid","Heatmap 1st IC","Heatmap 2nd IC","Heatmap 3rd IC","Heatmap 4th IC")
    oo<-order(ss)
    
    # Add info above plot    
    
    ifelse(input$order=="NZ Cohort",chc<-c(matlab.like(length(unique(NZcohort())))),chc<-c(matlab.like(length(unique(clinical.in()$Cohort..))))) # Make chart on top dependent on NZ cohort selection
    
    Grade_1<-c("white","darkgreen")
    Grade_2<-c("white","darkgreen")
    Grade_3<-c("white","darkgreen")
    
    Basal<-c("white","orange")
    Her2<-c("white","orange")
    Normal_like<-c("white","orange")
    Luminal_A<-c("white","orange")
    Luminal_B<-c("white","orange")
    erc<-c("lightblue","darkblue")
    

    # Build top bar of data based upon check boxes.  
    # Not using c[...,drop=FALSE] so the last check is "is.character" 
    
    cbar<-NULL
    
    ifelse(input$order=="NZ Cohort",Cohort<-chc[NZcohort()],Cohort<-chc[as.numeric(as.factor(clinical.in()$Cohort..))]) #Color bar based upon Cohort
    
    cbar<-rbind(cbar,Cohort)
    
    if(input$grade.show == TRUE) {
      grade.tmp<-rbind(Grade_1[as.numeric(as.factor(clinical.in()$G1))],
                       Grade_2[as.numeric(as.factor(clinical.in()$G2))],
                       Grade_3[as.numeric(as.factor(clinical.in()$G3))])
      rownames(grade.tmp)<-c("Grade 1","Grade 2","Grade 3")
      cbar<-rbind(cbar,grade.tmp)}
    if(input$subtype.show == TRUE) {
      subtype.tmp<-rbind(Basal[as.numeric(as.factor(clinical.in()$Basal))],
                         Her2[as.numeric(as.factor(clinical.in()$Her2))],
                         Normal_like[as.numeric(as.factor(clinical.in()$Normal))],
                         Luminal_A[as.numeric(as.factor(clinical.in()$LumA))],
                         Luminal_B[as.numeric(as.factor(clinical.in()$LumB))])
      rownames(subtype.tmp)<-c("Basal","Her2","Normal-like","Luminal A", "Luminal B")
      cbar<-rbind(cbar,subtype.tmp)}

    
    if(input$er.show == TRUE) {
      tmp.names<-rownames(cbar)
      cbar<-rbind(cbar,erc[as.numeric(as.factor(clinical.in()$ER.status))])
      rownames(cbar)<-c(tmp.names,"ER IHC Status")}
    
    if(input$gene.show == TRUE){
      
      tmp.names<-rownames(cbar)
      cbar<-rbind(cbar,GCOL[rank(gene.show.data()[2,])],GCOL[rank(gene.show.data()[1,])])
      rownames(cbar)<-c(tmp.names,"Ki67 Expression", "ESR1 Expression")
      
    }
    
    
    if(input$hoc.show == TRUE) {
      hoc.tmp<-rbind(ccCOL[rank(hoc[1,])],ccCOL[rank(hoc[2,])],ccCOL[rank(hoc[3,])],ccCOL[rank(hoc[4,])],ccCOL[rank(hoc[5,])],ccCOL[rank(hoc[6,])],ccCOL[rank(hoc[7,])],
                     ccCOL[rank(hoc[8,])],ccCOL[rank(hoc[9,])],ccCOL[rank(hoc[10,])])
      rownames(hoc.tmp)<-c("Proliferation","Growth Suppression","Immune","Immortality","Inflammation","Metastasis","Angiogenesis","Genome Instability","Cell Death","Cell Energy")
      cbar<-rbind(cbar,hoc.tmp)}
    
    
    # Since PC and IC are coloured above need to take out instead of add if false
    c.ic<-1
    if(input$pc.show == FALSE){
      cc<-cc[-match(c("Heatmap 1st PC","Heatmap 2nd PC","Heatmap 3rd PC","Heatmap 4th PC"),rownames(cc)),]
      
    }
    if(input$ic.show == FALSE){
      cc<-cc[-match(c("Heatmap 1st IC","Heatmap 2nd IC","Heatmap 3rd IC","Heatmap 4th IC"),rownames(cc)),]
      c.ic<-cc
    }
    cc<-rbind(cbar,cc)  # Note: always have centroid as part of cc
    
    if(is.character(c.ic) ){rownames(cc)[dim(cc)[1]]<-c("Heatmap Centroid")} # if pc and ic are gone rename cc to heatmap.
    
    # set up legend
    ifelse(input$order=="NZ Cohort",nm<-c(paste0("NZ (", table(NZcohort())[[1]],")"),paste0("Non-NZ (", table(NZcohort())[[2]],")")),nm<-c(paste(names(table(clinical.in()$Cohort..))," (",table(clinical.in()$Cohort..),")",sep='')))
    
    # Generate break.ties to allow subsorting the data 
    if(is.null(input$subsort))    break.ties <-rep(1,length(clinical.in()$Grade)) # Just a series of 1's
    if(!is.null(input$subsort))  {
      if(input$subsort=="Cohort") break.ties<-clinical.in()$Cohort..
      if(input$subsort=="NZ Cohort") break.ties<-NZcohort()
      if(input$subsort=="Grade") break.ties<-clinical.in()$Grade
      if(input$subsort=="Subtype") break.ties<-clinical.in()$Subtype
      if(input$subsort=="ER IHC Status") break.ties<-clinical.in()$ER.status
      if(input$subsort=="Ki67 Expression") break.ties<-gene.show.data()[2,]
      if(input$subsort=="ESR1 Expression") break.ties<-gene.show.data()[1,]
      if(input$subsort=="Heatmap 1st PC") break.ties<-ss
      if(input$subsort=="Heatmap 2nd PC") break.ties<-ss2
      if(input$subsort=="Heatmap 3rd PC") break.ties<-ss3
      if(input$subsort=="Heatmap 4th PC") break.ties<-ss4
      if(input$subsort=="Heatmap Centroid") break.ties<-ss_centroid
      if(input$subsort =="Heatmap 1st IC") break.ties<-ss5
      if(input$subsort =="Heatmap 2nd IC") break.ties<-ss6
      if(input$subsort =="Heatmap 3rd IC") break.ties<-ss7
      if(input$subsort =="Heatmap 4th IC") break.ties<-ss8
      if(input$subsort =="Proliferation") break.ties<-hoc[1,]
      if(input$subsort =="Growth Suppression") break.ties<-hoc[2,]
      if(input$subsort =="Immune") break.ties<-hoc[3,]
      if(input$subsort =="Immortality") break.ties<-hoc[4,]
      if(input$subsort =="Inflammation") break.ties<-hoc[5,]
      if(input$subsort =="Metastasis") break.ties<-hoc[6,]
      if(input$subsort =="Angiogenesis") break.ties<-hoc[7,]
      if(input$subsort =="Genome Instability") break.ties<-hoc[8,]
      if(input$subsort =="Cell Death") break.ties<-hoc[9,]
      if(input$subsort =="Cell Energy") break.ties<-hoc[10,]
      
      
    } 
    
    
    # need to clean this up - FIXME
    # order heatmap by ...
    if(input$order=="Cohort") oo<-order(clinical.in()$Cohort..,break.ties)
    if(input$order=="NZ Cohort") oo<-order(NZcohort(),break.ties)
    if(input$order=="Grade") oo<-order(clinical.in()$Grade,break.ties)
    if(input$order=="Subtype") oo<-order(clinical.in()$Subtype,break.ties)
    if(input$order=="ER IHC Status") oo<-order(clinical.in()$ER.status,break.ties)
    if(input$order=="Ki67 Expression") oo<-order(gene.show.data()[2,],break.ties)
    if(input$order=="ESR1 Expression") oo<-order(gene.show.data()[1,],break.ties)
    if(input$order=="Heatmap 1st PC") oo<-order(ss)
    if(input$order=="Heatmap 2nd PC") oo<-order(ss2)
    if(input$order=="Heatmap 3rd PC") oo<-order(ss3)
    if(input$order=="Heatmap 4th PC") oo<-order(ss4)
    if(input$order=="Heatmap Centroid") oo<-order(ss_centroid)
    if(input$order =="Heatmap 1st IC") oo<-order(ss5)
    if(input$order =="Heatmap 2nd IC") oo<-order(ss6)
    if(input$order =="Heatmap 3rd IC") oo<-order(ss7)
    if(input$order =="Heatmap 4th IC") oo<-order(ss8)
    if(input$order =="Proliferation") oo<-order(hoc[1,])
    if(input$order =="Growth Suppression") oo<-order(hoc[2,])
    if(input$order =="Immune") oo<-order(hoc[3,])
    if(input$order =="Immortality") oo<-order(hoc[4,])
    if(input$order =="Inflammation") oo<-order(hoc[5,])
    if(input$order =="Metastasis") oo<-order(hoc[6,])
    if(input$order =="Angiogenesis") oo<-order(hoc[7,])
    if(input$order =="Genome Instability") oo<-order(hoc[8,])
    if(input$order =="Cell Death") oo<-order(hoc[9,])
    if(input$order =="Cell Energy") oo<-order(hoc[10,])
    
    out<-list("gdat"=gdat,"oo"=oo,"cc"=cc,"nm_legend"=nm,"ccx"=ccx,"chc"=chc)
    
  })
  
output$expMap <- renderPlot({    
    
   # Pull in the relevent HeatMap stuff
    heat.in<-heatmapProcessing()
    
    if(input$order!="Clustering"){
        suppressWarnings(
       heatmap.mik(heat.in$gdat[,heat.in$oo], trace='none',
                    col=heat.in$ccx, keysize=1, key=T,
                    ColSideColors=heat.in$cc[,heat.in$oo],mar=c(10,10),
                    labCol="",Colv=F,cexCol=1.2)
        )
        legend(0.4,1.05,heat.in$nm_legend,ncol=6,
               fill=c(heat.in$chc),bty='n')
               
        legend(0.18,1.05,"COHORT", bty='n', cex=1.5)
      }
      if(input$order=="Clustering"){
        if(input$radioDist=="Euclidean"){
          suppressWarnings(
          heatmap.mik(heat.in$gdat, trace='none',
                        col=heat.in$ccx, keysize=1, key=T,
                        ColSideColors=heat.in$cc,mar=c(10,10),
                        labCol="",cexCol=1.2)
          )
        }
        if(input$radioDist=="Correlation"){
          suppressWarnings(
          heatmap.mik(heat.in$gdat, trace='none',
                        col=heat.in$ccx, keysize=1, key=T,
                        ColSideColors=heat.in$cc,mar=c(10,10),
                        labCol="",cexCol=1.2,
                        distfun=function(x) as.dist(1-cor(t(x))))
          )
        }

      }
    
  },height=1200,width=1000)
  
output$downloadHMPlot <- downloadHandler(
  filename = function() { paste0(gsub(":","",input$sigSelect)," HM ",input$order,'.pdf')},
  content = function(file) {
    pdf(file,width=9,height=12)
    # Pull in the relevent HeatMap stuff
    heat.in<-heatmapProcessing()
    
    if(input$order!="Clustering"){
      suppressWarnings(
        heatmap.mik(heat.in$gdat[,heat.in$oo], trace='none',
                          col=heat.in$ccx, keysize=1, key=T,
                          ColSideColors=heat.in$cc[,heat.in$oo],mar=c(10,10),
                          labCol="",Colv=F,cexCol=1.2)
      )
       legend(0.36,1.03,heat.in$nm_legend,ncol=6,
           fill=c(heat.in$chc),bty='n',title="Cohort")
    
       legend(0.22,1.0,"COHORT", bty='n')
      
    }
    if(input$order=="Clustering"){
      if(input$radioDist=="Euclidean"){
        suppressWarnings(
          heatmap.mik(heat.in$gdat, trace='none',
                            col=heat.in$ccx, keysize=1, key=T,
                            ColSideColors=heat.in$cc,mar=c(10,10),
                            labCol="",cexCol=1.2)
        )
      }
      if(input$radioDist=="Correlation"){
        suppressWarnings(
          heatmap.mik(heat.in$gdat, trace='none',
                            col=heat.in$ccx, keysize=1, key=T,
                            ColSideColors=heat.in$cc,mar=c(10,10),
                            labCol="",cexCol=1.2,
                            distfun=function(x) as.dist(1-cor(t(x))))
        )
      }
    }
   
    dev.off()
  },
  contentType = 'application/pdf')

codedOutput<-reactive({
  out<-rep(0,5)
  if(sum(grepl("Basal",input$Gsubtype))==1) out[1]<-1
  if(sum(grepl("Her2", input$Gsubtype))==1) out[2]<-1
  if(sum(grepl("LumA", input$Gsubtype))==1) out[3]<-1
  if(sum(grepl("LumB", input$Gsubtype))==1) out[4]<-1
  if(sum(grepl("Normal",input$Gsubtype))==1) out[5]<-1
  return(out)
})

output$dlSurvPlot <- downloadHandler(
  filename = function() { 
     paste0(gsub(":","",input$sigSelect),"Survial ",input$input.radioTFsurv," ",codedOutput()[1],codedOutput()[2],codedOutput()[3],codedOutput()[4],codedOutput()[5],'.pdf')
    },
  content = function(file) {
    
    pdf(file,width=20,height=12)
    par(mfrow=c(2,4))
    par(cex.axis=1.5, cex.lab=1.5, cex.main=1.5,oma=c(0,0,2,0),mar=c(6,4.5,7,1))
    
    surv.obj<-survPlotData()
    # Survival Debut Break Here <-
    
    surv.plot.shiny(ssx=surv.obj$ssx$ssx1,ss=surv.obj$ss,kp=surv.obj$kp$kp1,dmfsTime=surv.obj$dmfsTime,dmfsEvent=surv.obj$dmfsEvent,surv.obj$main[1],p.cox=input$input.coxph)
    surv.plot.shiny(ssx=surv.obj$ssx$ssx2,ss=surv.obj$ss,kp=surv.obj$kp$kp2,dmfsTime=surv.obj$dmfsTime,dmfsEvent=surv.obj$dmfsEvent,surv.obj$main[2],p.cox=input$input.coxph)
    surv.plot.shiny(ssx=surv.obj$ssx$ssx3,ss=surv.obj$ss,kp=surv.obj$kp$kp3,dmfsTime=surv.obj$dmfsTime,dmfsEvent=surv.obj$dmfsEvent,surv.obj$main[3],p.cox=input$input.coxph)
    surv.plot.shiny(ssx=surv.obj$ssx$ssx4,ss=surv.obj$ss,kp=surv.obj$kp$kp4,dmfsTime=surv.obj$dmfsTime,dmfsEvent=surv.obj$dmfsEvent,surv.obj$main[4],p.cox=input$input.coxph)
    surv.plot.shiny(ssx=surv.obj$ssx$ssx5,ss=surv.obj$ss,kp=surv.obj$kp$kp5,dmfsTime=surv.obj$dmfsTime,dmfsEvent=surv.obj$dmfsEvent,surv.obj$main[5],p.cox=input$input.coxph)
    surv.plot.shiny(ssx=surv.obj$ssx$ssx6,ss=surv.obj$ss,kp=surv.obj$kp$kp6,dmfsTime=surv.obj$dmfsTime,dmfsEvent=surv.obj$dmfsEvent,surv.obj$main[6],p.cox=input$input.coxph)
    surv.plot.shiny(ssx=surv.obj$ssx$ssx7,ss=surv.obj$ss,kp=surv.obj$kp$kp7,dmfsTime=surv.obj$dmfsTime,dmfsEvent=surv.obj$dmfsEvent,surv.obj$main[7],p.cox=input$input.coxph)
    surv.plot.shiny(ssx=surv.obj$ssx$ssx8,ss=surv.obj$ss,kp=surv.obj$kp$kp8,dmfsTime=surv.obj$dmfsTime,dmfsEvent=surv.obj$dmfsEvent,surv.obj$main[8],p.cox=input$input.coxph)
    
    dev.off()
  })

output$downloadAuCPlot <- downloadHandler(
  filename = function() { paste0(gsub(":","",input$sigSelect),"AuC ",input$AuCyear[1]," ",input$AuCyear[2], " ",codedOutput()[1],codedOutput()[2],codedOutput()[3],codedOutput()[4],codedOutput()[5],'.pdf')},
  content = function(file) {
    
    pdf(file,height=20)
    par(mfrow=c(5,2)) 
    plot.dat<-auc.data()

    if(length(input$auc.choice)>1)    pcsum<-rowSums(plot.dat$plot.matrix[,input$auc.choice])
    if(length(input$auc.choice)==1)   pcsum<-plot.dat$plot.matrix[,input$auc.choice]
    if(length(input$auc.choice)==0)   pcsum<-plot.dat$plot.matrix[,"Centroid"]

    if(dim(plot.dat$plot.matrix)[2] == 9) plot.me<-"Centroid"
    if(dim(plot.dat$plot.matrix)[2] == 10) plot.me<-"Score"

    auc.plot.shiny(plot.dat$plot.matrix[,plot.me],plot.dat$kp,plot.dat$dmfsTime,plot.dat$dmfsEvent,plot.me,time=input$AuCyear[2])
    auc.plot.shiny(pcsum,plot.dat$kp,plot.dat$dmfsTime,plot.dat$dmfsEvent,name=paste0(input$auc.choice,collapse=" "),time=input$AuCyear[2]) # plot the results from multi selection
    auc.plot.shiny(plot.dat$plot.matrix[,"PC1"],plot.dat$kp,plot.dat$dmfsTime,plot.dat$dmfsEvent,"PC1",time=input$AuCyear[2])
    auc.plot.shiny(plot.dat$plot.matrix[,"PC2"],plot.dat$kp,plot.dat$dmfsTime,plot.dat$dmfsEvent,"PC2",time=input$AuCyear[2])
    auc.plot.shiny(plot.dat$plot.matrix[,"PC3"],plot.dat$kp,plot.dat$dmfsTime,plot.dat$dmfsEvent,"PC3",time=input$AuCyear[2])
    auc.plot.shiny(plot.dat$plot.matrix[,"PC4"],plot.dat$kp,plot.dat$dmfsTime,plot.dat$dmfsEvent,"PC4",time=input$AuCyear[2])
    auc.plot.shiny(plot.dat$plot.matrix[,"IC1"],plot.dat$kp,plot.dat$dmfsTime,plot.dat$dmfsEvent,"IC1",time=input$AuCyear[2])
    auc.plot.shiny(plot.dat$plot.matrix[,"IC2"],plot.dat$kp,plot.dat$dmfsTime,plot.dat$dmfsEvent,"IC2",time=input$AuCyear[2])
    auc.plot.shiny(plot.dat$plot.matrix[,"IC3"],plot.dat$kp,plot.dat$dmfsTime,plot.dat$dmfsEvent,"IC3",time=input$AuCyear[2])
    auc.plot.shiny(plot.dat$plot.matrix[,"IC4"],plot.dat$kp,plot.dat$dmfsTime,plot.dat$dmfsEvent,"IC4",time=input$AuCyear[2])
    dev.off()
  })

output$dlHOCcolors <- downloadHandler(
  filename = function() { paste0(gsub(":","",input$sigSelect)," HOC ",input$HOCcormethod," ",codedOutput()[1],codedOutput()[2],codedOutput()[3],codedOutput()[4],codedOutput()[5],'.pdf')},
  content = function(file) {
    ggsave(file, HOCtabin(),width=15.4,height=8.41)
  })

survPlotData<-reactive({
  # Create names for plots
  main=c("All treatment groups","No adjuvant Rx","No endocrine Rx, may or may not have chemotherapy","No Chemotherapy, may or may not have endocrine Rx","Endocrine Rx , NO chemotherapy","Chemotherapy, NO endocrine Rx","Endocrine Rx , +/- Chemotherapy","Chemotherapy, +/- Endocrine Rx")
  # Pull in the data
  ss<-eDatSubtype()$PcIcOut[,(colnames(eDatSubtype()$PcIcOut) %in% input$order.surv)]  
  
  
  # If individual genes were chosen for survival analysis then replace the gene set centroid/component with the Z-trasnformed individual gene data vector  
  if(input$radioTFsurv=="Individual Gene") {
    
    
    # first get gene name
    g<-strsplit(gsub(":"," ",gsub(","," ",gsub(";"," ",input$genesTocompareS,fixed=T),
                                  fixed=T),fixed=T)," ")[[1]]
    g<-g[g!=""]
    g<-toupper(g)
    
    # then get gene expression data
    #geneData<-cbind(expDat[match(unlist(annDat[match(g,annDat[,3]),1]),rownames(expDat)),])
    geneData<-dbGetQuery(conn=db, paste0("SELECT * FROM expdat WHERE row_names ='",g,"'"))
    geneData <-as.numeric(geneData[-1])
    
    # Z-transform gene expression data 
    ss<-(geneData-mean(geneData,na.rm=T))/sd(geneData,na.rm=T)
    ss<-ss[as.logical(global.subtype())]
    
  }
  if(input$radioRS=="Yes")    {
    ss<-(-ss) }
  dmfsEvent<-clinicalSubtype()$DMFS.EVENT..defined.as.distant.metastasis.or.from.breast.cancer.
  dmfsTime<-clinicalSubtype()$DMFS.TIME 
  
  
  # Select subtypes based on checkbox 
  kpST<-rep(T,length(ss))
  
  
  kpER<-clinicalSubtype()$ER.status==input$radioER # select subtype
  if(input$radioER=="All") kpER <-rep(T,nrow(clinicalSubtype()))
  
  if(input$radioNZ=="All")      kpNZ<- rep(T,nrow(clinicalSubtype()))
  if(input$radioNZ=="NZ Only")  kpNZ<- NZcohort()[as.logical(global.subtype())]==1
  if(input$radioNZ=="Non-NZ")   kpNZ<- NZcohort()[as.logical(global.subtype())]==2
  
  # Add grade selection for survial analysis
   kpGd<-rep(T,nrow(clinicalSubtype())) # Always run this to avoid errors if there are no samples of a given grade
  if(input$radioGd=="Grade 1") kpGd<-clinicalSubtype()$G1==1
  if(input$radioGd=="Grade 2") kpGd<-clinicalSubtype()$G2==1
  if(input$radioGd=="Grade 3") kpGd<-clinicalSubtype()$G3==1
  
    #plot 1 - All treatment groups
    #define subset of data
    kpTRT<-rep(T,nrow(clinicalSubtype()))
    kp1<-kpTRT & kpST & kpGd & kpER & kpNZ
    kp1[is.na(kp1)]<-FALSE # ensure that kp for all tumours that are not TRUE is False (instead of NA)
    
    if(input$survplottype=="2 groups- Quantile") ssx1<-ifelse(ss[kp1]>quantile(ss[kp1],input$survSlider),"High","Low")
    if(input$survplottype=="2 groups- Median")   ssx1<-ifelse(ss[kp1]>quantile(ss[kp1],0.5),"High","Low")
    if(input$survplottype=="2 groups- Top/Bottom") {
      ssx1<-ss  
      ssx1[ss[kp1]<quantile(ss[kp1],input$survSlider2[1])]<-"Low"
      ssx1[ss[kp1]>quantile(ss[kp1],input$survSlider2[2])]<-"High"
      ssx1[ss[kp1]<quantile(ss[kp1],input$survSlider2[2]) & ss[kp1]>quantile(ss[kp1],input$survSlider2[1])]<-NA
      kp1[is.na(ssx1)]<-FALSE # remove the NAs
      ssx1<-ssx1[kp1]   # Reduce ssx to proper size 
      }
    
    #plot 2 - treatment group No adjuvant Rx
    kpTRT<-clinicalSubtype()$adjuv..treated.==0
    kp2<-kpTRT & kpST & kpGd & kpER & kpNZ# set up subset to plot and to perform survival statistics on
    kp2[is.na(kp2)]<-FALSE # ensure that kp for all tumours that are not TRUE is False (instead of NA)
    if(input$survplottype=="2 groups- Quantile") ssx2<-ifelse(ss[kp2]>quantile(ss[kp2],input$survSlider),"High","Low")
    if(input$survplottype=="2 groups- Median")   ssx2<-ifelse(ss[kp2]>quantile(ss[kp2],0.5),"High","Low")
    if(input$survplottype=="2 groups- Top/Bottom") {
      ssx2<-ss  
      ssx2[ss[kp2]<quantile(ss[kp2],input$survSlider2[1])]<-"Low"
      ssx2[ss[kp2]>quantile(ss[kp2],input$survSlider2[2])]<-"High"
      ssx2[ss[kp2]<quantile(ss[kp2],input$survSlider2[2]) & ss[kp2]>quantile(ss[kp2],input$survSlider2[1])]<-NA
      kp2[is.na(ssx2)]<-FALSE # remove the NAs
      ssx2<-ssx2[kp2]   # Reduce ssx to proper size 
    }
    
    #plot 3 - treatment group No endocrine Rx but may or may not have chemotherapy
    kpTRT<-clinicalSubtype()$NO.endo==1
    kp3<-kpTRT & kpST & kpGd & kpER & kpNZ# set up subset to plot and to perform survival statistics on
    kp3[is.na(kp3)]<-FALSE # ensure that kp for all tumours that are not TRUE is False (instead of NA)
    if(input$survplottype=="2 groups- Quantile") ssx3<-ifelse(ss[kp3]>quantile(ss[kp3],input$survSlider),"High","Low")
    if(input$survplottype=="2 groups- Median")   ssx3<-ifelse(ss[kp3]>quantile(ss[kp3],0.5),"High","Low")
    if(input$survplottype=="2 groups- Top/Bottom") {
      ssx3<-ss  
      ssx3[ss[kp3]<quantile(ss[kp3],input$survSlider2[1])]<-"Low"
      ssx3[ss[kp3]>quantile(ss[kp3],input$survSlider2[2])]<-"High"
      ssx3[ss[kp3]<quantile(ss[kp3],input$survSlider2[2]) & ss[kp3]>quantile(ss[kp3],input$survSlider2[1])]<-NA
      kp3[is.na(ssx3)]<-FALSE # remove the NAs
      ssx3<-ssx3[kp3]   # Reduce ssx to proper size 
    }
    
    #plot 4
    #define subset of data treatment group - No Chemotherapy but may or may not have endocrine Rx
    kpTRT<-clinicalSubtype()$NO.chemo==1
    kp4<-kpTRT & kpST & kpGd & kpER & kpNZ# set up subset to plot and to perform survival statistics on
    kp4[is.na(kp4)]<-FALSE # ensure that kp for all tumours that are not TRUE is False (instead of NA)
    if(input$survplottype=="2 groups- Quantile") ssx4<-ifelse(ss[kp4]>quantile(ss[kp4],input$survSlider),"High","Low")
    if(input$survplottype=="2 groups- Median")   ssx4<-ifelse(ss[kp4]>quantile(ss[kp4],0.5),"High","Low")
    if(input$survplottype=="2 groups- Top/Bottom") {
      ssx4<-ss  
      ssx4[ss[kp4]<quantile(ss[kp4],input$survSlider2[1])]<-"Low"
      ssx4[ss[kp4]>quantile(ss[kp4],input$survSlider2[2])]<-"High"
      ssx4[ss[kp4]<quantile(ss[kp4],input$survSlider2[2]) & ss[kp4]>quantile(ss[kp4],input$survSlider2[1])]<-NA
      kp4[is.na(ssx4)]<-FALSE # remove the NAs
      ssx4<-ssx4[kp4]   # Reduce ssx to proper size 
    }
    
    #plot 5
    #define subset of data treatment group - Endocrine Rx , NO chemotherapy
    kpTRT<-clinicalSubtype()$endo.only==1
    kp5<-kpTRT & kpST & kpGd & kpER & kpNZ# set up subset to plot and to perform survival statistics on
    kp5[is.na(kp5)]<-FALSE # ensure that kp for all tumours that are not TRUE is False (instead of NA)
    if(input$survplottype=="2 groups- Quantile") ssx5<-ifelse(ss[kp5]>quantile(ss[kp5],input$survSlider),"High","Low")
    if(input$survplottype=="2 groups- Median")   ssx5<-ifelse(ss[kp5]>quantile(ss[kp5],0.5),"High","Low")
    if(input$survplottype=="2 groups- Top/Bottom") {
      ssx5<-ss  
      ssx5[ss[kp5]<quantile(ss[kp5],input$survSlider2[1])]<-"Low"
      ssx5[ss[kp5]>quantile(ss[kp5],input$survSlider2[2])]<-"High"
      ssx5[ss[kp5]<quantile(ss[kp5],input$survSlider2[2]) & ss[kp5]>quantile(ss[kp5],input$survSlider2[1])]<-NA
      kp5[is.na(ssx5)]<-FALSE # remove the NAs
      ssx5<-ssx5[kp5]   # Reduce ssx to proper size 
    }
    
    #plot 6
    #define subset of data treatment group - Chemotherapy, NO Endocrine Rx 
    kpTRT<-clinicalSubtype()$chemo.only==1
    kp6<-kpTRT & kpST & kpGd & kpER & kpNZ # set up subset to plot and to perform survival statistics on
    kp6[is.na(kp6)]<-FALSE # ensure that kp for all tumours that are not TRUE is False (instead of NA)
    if(input$survplottype=="2 groups- Quantile") ssx6<-ifelse(ss[kp6]>quantile(ss[kp6],input$survSlider),"High","Low")
    if(input$survplottype=="2 groups- Median")   ssx6<-ifelse(ss[kp6]>quantile(ss[kp6],0.5),"High","Low")
    if(input$survplottype=="2 groups- Top/Bottom") {
      ssx6<-ss  
      ssx6[ss[kp6]<quantile(ss[kp6],input$survSlider2[1])]<-"Low"
      ssx6[ss[kp6]>quantile(ss[kp6],input$survSlider2[2])]<-"High"
      ssx6[ssx6 != "Low" & ssx6 !="High"]<-NA
      kp6[is.na(ssx6)]<-FALSE # remove the NAs
      ssx6<-ssx6[kp6]   # Reduce ssx to proper size 
    }
    
    #plot 7
    #define subset of data treatment group - Endocrine Rx , +/- Chemotherapy
    kpTRT<-clinicalSubtype()$endocrine.==1
    #define subset of data
    kp7<-kpTRT & kpST & kpGd & kpER & kpNZ # set up subset to plot and to perform survival statistics on
    kp7[is.na(kp7)]<-FALSE # ensure that kp for all tumours that are not TRUE is False (instead of NA)
    if(input$survplottype=="2 groups- Quantile") ssx7<-ifelse(ss[kp7]>quantile(ss[kp7],input$survSlider),"High","Low")
    if(input$survplottype=="2 groups- Median")   ssx7<-ifelse(ss[kp7]>quantile(ss[kp7],0.5),"High","Low")
    if(input$survplottype=="2 groups- Top/Bottom") {
      ssx7<-ss  
      ssx7[ss[kp7]<quantile(ss[kp7],input$survSlider2[1])]<-"Low"
      ssx7[ss[kp7]>quantile(ss[kp7],input$survSlider2[2])]<-"High"
      ssx7[ss[kp7]<quantile(ss[kp7],input$survSlider2[2]) & ss[kp7]>quantile(ss[kp7],input$survSlider2[1])]<-NA
      kp7[is.na(ssx7)]<-FALSE # remove the NAs
      ssx7<-ssx7[kp7]   # Reduce ssx to proper size 
    }
    
    #plot 8
    #define subset of data treatment group - Chemotherapy, +/- Endocrine Rx
    kpTRT<-clinicalSubtype()$chemo.==1
    kp8<-kpTRT & kpST & kpGd & kpER & kpNZ # set up subset to plot and to perform survival statistics on
    kp8[is.na(kp8)]<-FALSE # ensure that kp for all tumours that are not TRUE is False (instead of NA)
    if(input$survplottype=="2 groups- Quantile") ssx8<-ifelse(ss[kp8]>quantile(ss[kp8],input$survSlider),"High","Low")
    if(input$survplottype=="2 groups- Median")   ssx8<-ifelse(ss[kp8]>quantile(ss[kp8],0.5),"High","Low")
    if(input$survplottype=="2 groups- Top/Bottom") {
      ssx8<-ss  
      ssx8[ss[kp8]<quantile(ss[kp8],input$survSlider2[1])]<-"Low"
      ssx8[ss[kp8]>quantile(ss[kp8],input$survSlider2[2])]<-"High"
      ssx8[ss[kp8]<quantile(ss[kp8],input$survSlider2[2]) & ss[kp8]>quantile(ss[kp8],input$survSlider2[1])]<-NA
      kp8[is.na(ssx8)]<-FALSE # remove the NAs
      ssx8<-ssx8[kp8]   # Reduce ssx1 to proper size 
    }
    
    kp<-list("kp1"=kp1,"kp2"=kp2,"kp3"=kp3,"kp4"=kp4,"kp5"=kp5,"kp6"=kp6,"kp7"=kp7,"kp8"=kp8)
    ssx<-list("ssx1"=ssx1,"ssx2"=ssx2,"ssx3"=ssx3,"ssx4"=ssx4,"ssx5"=ssx5,"ssx6"=ssx6,"ssx7"=ssx7,"ssx8"=ssx8)
    out<-list("main"=main,"kp"=kp,"ssx"=ssx,"ss"=ss,"dmfsTime"=dmfsTime,"dmfsEvent"=dmfsEvent)
})
 
 #Survial plots
output$survPlot <- renderPlot({  
   par(mfrow=c(2,4))
   par(cex.axis=1.5, cex.lab=1.5, cex.main=1.5,oma=c(0,0,2,0),mar=c(6,4.5,7,1))
  
   surv.obj<-survPlotData()
   surv.plot.shiny(ssx=surv.obj$ssx$ssx1,ss=surv.obj$ss,kp=surv.obj$kp$kp1,dmfsTime=surv.obj$dmfsTime,dmfsEvent=surv.obj$dmfsEvent,surv.obj$main[1],p.cox=input$input.coxph)
   surv.plot.shiny(ssx=surv.obj$ssx$ssx2,ss=surv.obj$ss,kp=surv.obj$kp$kp2,dmfsTime=surv.obj$dmfsTime,dmfsEvent=surv.obj$dmfsEvent,surv.obj$main[2],p.cox=input$input.coxph)
   surv.plot.shiny(ssx=surv.obj$ssx$ssx3,ss=surv.obj$ss,kp=surv.obj$kp$kp3,dmfsTime=surv.obj$dmfsTime,dmfsEvent=surv.obj$dmfsEvent,surv.obj$main[3],p.cox=input$input.coxph)
   surv.plot.shiny(ssx=surv.obj$ssx$ssx4,ss=surv.obj$ss,kp=surv.obj$kp$kp4,dmfsTime=surv.obj$dmfsTime,dmfsEvent=surv.obj$dmfsEvent,surv.obj$main[4],p.cox=input$input.coxph)
   surv.plot.shiny(ssx=surv.obj$ssx$ssx5,ss=surv.obj$ss,kp=surv.obj$kp$kp5,dmfsTime=surv.obj$dmfsTime,dmfsEvent=surv.obj$dmfsEvent,surv.obj$main[5],p.cox=input$input.coxph)
   surv.plot.shiny(ssx=surv.obj$ssx$ssx6,ss=surv.obj$ss,kp=surv.obj$kp$kp6,dmfsTime=surv.obj$dmfsTime,dmfsEvent=surv.obj$dmfsEvent,surv.obj$main[6],p.cox=input$input.coxph)
   surv.plot.shiny(ssx=surv.obj$ssx$ssx7,ss=surv.obj$ss,kp=surv.obj$kp$kp7,dmfsTime=surv.obj$dmfsTime,dmfsEvent=surv.obj$dmfsEvent,surv.obj$main[7],p.cox=input$input.coxph)
   surv.plot.shiny(ssx=surv.obj$ssx$ssx8,ss=surv.obj$ss,kp=surv.obj$kp$kp8,dmfsTime=surv.obj$dmfsTime,dmfsEvent=surv.obj$dmfsEvent,surv.obj$main[8],p.cox=input$input.coxph)
},height=1000,width=1600)

  # Surv Title with dynamic list of who is plotted
  
  
  output$survTitle <- renderText({    
      
      whos<-paste(input$radioST,collapse=", ")
      paste("Below are survival graphs divided by treatment group.  These are for the ",whos,"tumour subtype(s). Data are Z-transformed")
    })
  
    # Static title could me moved to UI.R in future
  output$gtpcTitle <- renderText({    
     paste("Graphs comparing an individual gene/gene set to another individual gene/gene set for the tumour subtype(s) and characteristics selected on the left.")
    })
  
  # separate out to allow downloading of images
  aucPlotIn<-reactive({
    # Input$AuCchoice was removed due to startup issue now will be called separatly

    plot.dat<-auc.data()
    
    if(length(input$auc.choice)>1)    pcsum<-rowSums(plot.dat$plot.matrix[,input$auc.choice])
    if(length(input$auc.choice)==1)   pcsum<-plot.dat$plot.matrix[,input$auc.choice]
    if(length(input$auc.choice)==0)   pcsum<-plot.dat$plot.matrix[,"Centroid"]
    
    if(dim(plot.dat$plot.matrix)[2] == 9) plot.me<-"Centroid"
    if(dim(plot.dat$plot.matrix)[2] == 10) plot.me<-"Score"
    
    auc.plot.shiny(plot.dat$plot.matrix[,plot.me],plot.dat$kp,plot.dat$dmfsTime,plot.dat$dmfsEvent,plot.me,time=input$AuCyear[2])
    auc.plot.shiny(pcsum,plot.dat$kp,plot.dat$dmfsTime,plot.dat$dmfsEvent,name=paste0(input$auc.choice,collapse=" "),time=input$AuCyear[2]) # plot the results from multi selection
    auc.plot.shiny(plot.dat$plot.matrix[,"PC1"],plot.dat$kp,plot.dat$dmfsTime,plot.dat$dmfsEvent,"PC1",time=input$AuCyear[2])
    auc.plot.shiny(plot.dat$plot.matrix[,"PC2"],plot.dat$kp,plot.dat$dmfsTime,plot.dat$dmfsEvent,"PC2",time=input$AuCyear[2])
    auc.plot.shiny(plot.dat$plot.matrix[,"PC3"],plot.dat$kp,plot.dat$dmfsTime,plot.dat$dmfsEvent,"PC3",time=input$AuCyear[2])
    auc.plot.shiny(plot.dat$plot.matrix[,"PC4"],plot.dat$kp,plot.dat$dmfsTime,plot.dat$dmfsEvent,"PC4",time=input$AuCyear[2])
    auc.plot.shiny(plot.dat$plot.matrix[,"IC1"],plot.dat$kp,plot.dat$dmfsTime,plot.dat$dmfsEvent,"IC1",time=input$AuCyear[2])
    auc.plot.shiny(plot.dat$plot.matrix[,"IC2"],plot.dat$kp,plot.dat$dmfsTime,plot.dat$dmfsEvent,"IC2",time=input$AuCyear[2])
    auc.plot.shiny(plot.dat$plot.matrix[,"IC3"],plot.dat$kp,plot.dat$dmfsTime,plot.dat$dmfsEvent,"IC3",time=input$AuCyear[2])
    auc.plot.shiny(plot.dat$plot.matrix[,"IC4"],plot.dat$kp,plot.dat$dmfsTime,plot.dat$dmfsEvent,"IC4",time=input$AuCyear[2])
   
  })
  
  # Create AuC Plots
  output$aucPlot<-renderPlot({
    par(mfrow=c(5,2))
    par(cex.axis=1.5, cex.lab=1.5, cex.main=1.5,mar=c(6,6.1,6.1,6.1))
    aucPlotIn() # call the plots
    
  },height=1200,width=800)
 
  output$gtpcPlot <- renderPlot({    
  	
  	#get gene nam
g<-strsplit(gsub(":"," ",gsub(","," ",gsub(";"," ",input$genesTocompare,fixed=T),
                                   fixed=T),fixed=T)," ")[[1]]
      g<-g[g!=""]
 g<-toupper(g)
 
# then get gene expression data
 geneData<-dbGetQuery(conn=db, paste0("SELECT * FROM expdat WHERE row_names ='",g,"'"))
 geneData <-as.numeric(geneData[1,-1])
 
  # Z-transform gene expression data 
geneData<-(geneData-mean(geneData,na.rm=T))/sd(geneData,na.rm=T)


# Note that what we are using is Z-transformed data from the outset
#Pull in the correct data by matching to column names 
ss<-eDat()$PcIcOut[,(colnames(eDat()$PcIcOut) %in% input$order.gtpc)]   # Match names with %in% which always returns logical

 if(input$radioTF=="Rank")    { # gene expression data and gene set data are already Z-transformed
# rank gene set summaries/components for plotting
 ss<-rank(ss)/length(ss)
# rank gene expression data
geneData <-as.numeric(geneData)
geneData<-rank(geneData)/length(geneData)
}

 if(input$radioRA=="Yes")    {
ss<-(-ss)
}

# select subset if tumours to plot    
	
	kpST<-global.subtype()
	  
   
       #  Add grade selection for survial analysis
        if(input$radiogtpcGd=="All") kpGd<-rep(T,nrow(clinical.in()))
        if(input$radiogtpcGd=="Grade 1") kpGd<-clinical.in()$G1==1
        if(input$radiogtpcGd=="Grade 2") kpGd<-clinical.in()$G2==1
        if(input$radiogtpcGd=="Grade 3") kpGd<-clinical.in()$G3==1
        
       # Add treatment
        if(input$radioRX=="All tumours") kpRX<-rep(T,nrow(clinical.in()))
        if(input$radioRX=="No adjuvant Rx") kpRX<-clinical.in()$adjuv..treated.==0
        if(input$radioRX=="Endocrine Rx no chemoRx") kpRX<-clinical.in()$endo.only==1


 #define subset of tumours
kERp<-clinical.in()$ER.status=="ER+"
kERn<-clinical.in()$ER.status=="ER-"
kpTRT<-rep(T,nrow(clinical.in()))
kp<-kpTRT & kpST & kpGd

kpp<-kpTRT & kpST & kpGd & kERp & kpRX
kpn<-kpTRT & kpST & kpGd & kERn & kpRX
ylim=c(min(ss),max(ss))
xlim=c(min(geneData),max(geneData))

if(anyNA(xlim)) xlim=c(-4,4)

plot(geneData[kpp],ss[kpp], col="red", xlab="", ylab="",main="",cex.axis=1.5, cex.lab=1.5, cex=1.5, cex.main=1.5, xlim=xlim, ylim=ylim)

par(new=TRUE)

mainl=paste("Gene ",g," not found")

if(!is.null(geneData) & !is.na(sum(geneData))) { # only calculate correlation if gene data exists!
  size<-round(cor(geneData[kp],ss[kp],use="complete.obs"),3)
  mainl=paste("correlation = ",size)
}
plot(geneData[kpn],ss[kpn], col="blue", xlab=paste(g,"expression"), ylab="Gene set centroid or component",main=mainl,cex.axis=1.5, cex.lab=1.5, cex.main=1.5, pch=19, xlim=xlim, ylim=ylim, sub="blue= ER-, red=ER+", cex=1)
  },height=600,width=600)
 
output$glxx <- renderText({
    proper.name<-paste0(toupper(substr(input$cormethod, 1, 1)), substr(input$cormethod, 2, nchar(input$cormethod)))
  if(input$corABS == "Yes"){
   paste("Absolute value of ",proper.name, " correlations between individual genes and the gene set's centroid and decomposed components",sep="")
  }else {
    paste(proper.name, " correlations between individual genes and the gene set's centroid and decomposed components",sep="")
  }
  
  
})

# Tab to show correlations of genes in list to centroid and PCs
output$glDisplay <- DT::renderDataTable({

           #generate object of correlations to centroid and PCs - all on Z-transformed data
          gll.t<-data.frame(cbind(t(cor(eDatSubtype()$PcIcOut[,5],t(eDatSubtype()$gdat),method=input$cormethod)),
                                t(cor(eDatSubtype()$PcIcOut[,1],t(eDatSubtype()$gdat),method=input$cormethod)),
                                t(cor(eDatSubtype()$PcIcOut[,2],t(eDatSubtype()$gdat),method=input$cormethod)),
                                t(cor(eDatSubtype()$PcIcOut[,3],t(eDatSubtype()$gdat),method=input$cormethod)),
                                t(cor(eDatSubtype()$PcIcOut[,4],t(eDatSubtype()$gdat),method=input$cormethod)),
                                t(cor(eDatSubtype()$PcIcOut[,6],t(eDatSubtype()$gdat),method=input$cormethod)),
                                t(cor(eDatSubtype()$PcIcOut[,7],t(eDatSubtype()$gdat),method=input$cormethod)),
                                t(cor(eDatSubtype()$PcIcOut[,8],t(eDatSubtype()$gdat),method=input$cormethod)),
                                t(cor(eDatSubtype()$PcIcOut[,9],t(eDatSubtype()$gdat),method=input$cormethod))))
            gll.t<-round(gll.t,digits=3)
          if(input$corABS == "Yes"){ gll.t<-abs(gll.t)}
          gll<-data.frame(cbind(rownames(eDatSubtype()$gdat),gll.t))
          colnames(gll)<-c("Gene Names","Centroid","PC1","PC2","PC3","PC4","IC1","IC2","IC3","IC4")
          gll
        },extensions='Scroller', options=list(
          paging=FALSE,
          info=FALSE,
          caseInsensitive = TRUE,
          searchHighlight = TRUE,
          scrollY=600,
          scrolling=TRUE),
          rownames=FALSE,class=c("display","compact"))

# Generate HOC correlation to PC IC and Centroid
# This info is passed to rendertable and renderplots below based upon UI input
HOCcor<-reactive({
  
  # Now correlate them together and create DataTable of the result
  main.table<-data.frame(cbind(t(cor(eDatSubtype()$PcIcOut[,5],t(hocCalculate()[,as.logical(global.subtype())]),method=input$HOCcormethod)),
                               t(cor(eDatSubtype()$PcIcOut[,1],t(hocCalculate()[,as.logical(global.subtype())]),method=input$HOCcormethod)),
                               t(cor(eDatSubtype()$PcIcOut[,2],t(hocCalculate()[,as.logical(global.subtype())]),method=input$HOCcormethod)),
                               t(cor(eDatSubtype()$PcIcOut[,3],t(hocCalculate()[,as.logical(global.subtype())]),method=input$HOCcormethod)),
                               t(cor(eDatSubtype()$PcIcOut[,4],t(hocCalculate()[,as.logical(global.subtype())]),method=input$HOCcormethod)),
                               t(cor(eDatSubtype()$PcIcOut[,6],t(hocCalculate()[,as.logical(global.subtype())]),method=input$HOCcormethod)),
                               t(cor(eDatSubtype()$PcIcOut[,7],t(hocCalculate()[,as.logical(global.subtype())]),method=input$HOCcormethod)),
                               t(cor(eDatSubtype()$PcIcOut[,8],t(hocCalculate()[,as.logical(global.subtype())]),method=input$HOCcormethod)),
                               t(cor(eDatSubtype()$PcIcOut[,9],t(hocCalculate()[,as.logical(global.subtype())]),method=input$HOCcormethod))))
  main.table<-round(main.table,digits=3)
  if(input$HOCcorABS == "Yes"){ main.table<-abs(main.table)}
  colnames(main.table)<-c("Centroid","PC1","PC2","PC3","PC4","IC1","IC2","IC3","IC4")
  rownames(main.table)<-c("Proliferation","Growth Suppression","Immune","Immortality","Inflammation","Metastasis","Angiogenesis","Genome Instability","Cell Death","Cell Energy")
  return(t(main.table))
})

output$HOCtabtable<-DT::renderDataTable({
  HOCcor()},
  options=list(
    paging=FALSE,
    info=FALSE,
    searching=FALSE
    ),
  selection=list(mode='single')
  )

HOCtabin<-reactive({
  plot.data<-melt(HOCcor())
  colnames(plot.data)<-c("Comp","HOC","Corr")
  ggplot(plot.data,aes(x=HOC,y=Comp,fill=Corr))+ geom_tile()+ labs(x = "Hallmarks of Cancer",fill="Correlation",y="Signature Components")+ggtitle("Correlation of Components vs Hallmarks of Cancer")+
    scale_fill_gradient2(limits=c(-1,1),low="red",mid="white",high="steelblue",breaks=c(-1,0,1),midpoint=0)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank(), plot.title=element_text(size=22),axis.title = element_text(face="bold", size=16))
})
outputOptions(output, "signature", suspendWhenHidden=FALSE)
output$HOCtabPlot<-renderPlot({
  
  HOCtabin()
  
},height=600,width=1100)
  })# end Server Function
  
