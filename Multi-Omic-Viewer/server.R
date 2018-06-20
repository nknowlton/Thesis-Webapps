

if (!require("pacman")) install.packages("pacman")
pacman::p_load(shiny,shinydashboard, sqldf, matrixStats, gplots, ggplot2, colorRamps,rms,devtools,fields)


#Databases
databases<-c("D:/btsync/skcm","D:/skcm","E:/skcm","/data/databases") # Enter in your directory to search here- Only first matched used



# check if dB file is in working directory if not go to where it should be.
if(file.exists("Mar16_BRCA_TCGA.sqlite")) {
  # Nothing
} else {
  catch<-which(dir.exists(databases))
  setwd(databases[catch[1]])
  print(paste0("Database directory assisgned to ",getwd(),"."))
}




# AS OF 3 21 2016 need Development version of DT

# Define server logic
shinyServer(function(input, output,session) {

db1.select<-reactive({

  if(input$database1=='BRCA')    db1<-dbConnect(SQLite(shared.cache = T),dbname="Mar16_BRCA_TCGA.sqlite")
  if(input$database1=='COAD')    db1<-dbConnect(SQLite(shared.cache = T),dbname="COAD_TCGA.sqlite")
  if(input$database1=='SKCM')    db1<-dbConnect(SQLite(shared.cache = T),dbname="20151101SKCM.sqlite")
  return(db1)
})

table.select<-reactive({
  if(input$variant == 'No' & input$cna.input=='Gistic' & input$lm.type=='relimp')   table=list(gene="Relimp_Meth_3",all="Relimp_All_3")
  if(input$cna.input=='Gistic' & input$lm.type=='cor')                              table=list(gene="Cor_Meth",all="Cor_All")
  if(input$cna.input=='CNA' & input$lm.type=='cor')                                 table=list(gene="Cor_Meth_CNA",all="Cor_All_CNA")
  if(input$variant == 'Yes' & input$cna.input=='Gistic' & input$lm.type=='relimp')  table=list(gene="Relimp_Meth_4",all="Relimp_All_4")
  if(input$variant == 'Yes' & input$cna.input=='CNA' & input$lm.type=='relimp')     table=list(gene="Relimp_Meth_4_CNA",all="Relimp_All_4_CNA")
  if(input$variant == 'Yes2' & input$cna.input=='Gistic' & input$lm.type=='relimp')  table=list(gene="Relimp_Meth_4_VAR",all="Relimp_All_4_VAR")
  if(input$variant == 'Yes2' & input$cna.input=='CNA' & input$lm.type=='relimp')     table=list(gene="Relimp_Meth_4_CNA_VAR",all="Relimp_All_4_Var_CNA")
  if(input$variant == 'No' & input$cna.input=='CNA' & input$lm.type=='relimp')      table=list(gene="Relimp_Meth_3_CNA",all="Relimp_All_3_CNA")
  return(table)

})

GT<-reactive({
  genes.table<-dbReadTable(db1.select(),table.select()$gene)

  if(input$lm.type=='cor'){
     genes.table<-genes.table[,-2]

     colnames(genes.table)<-c("Symbol","Meth Site", "Copy Number", "Methylation", "Variants")

     genes.table$`Copy Number`<-round(as.numeric(genes.table$`Copy Number`)*100,2)
     genes.table$Methylation<-round(as.numeric(genes.table$Methylation)*100,2)
     genes.table$Variants<-round(as.numeric(genes.table$Variants)*100,2)

     if(input$variant=="No") return(genes.table[,-6])
     return(genes.table)

   }



  genes.table<-genes.table[,-2]

  if(input$variant == 'No'){
    colnames(genes.table)<-c("Symbol","Meth Site", "Total", "Copy Number", "Methylation")
  }
  if(input$variant == 'Yes'){
  colnames(genes.table)<-c("Symbol","Meth Site", "Total", "Copy Number", "Methylation", "Variants")
  genes.table$Variants<-round(as.numeric(genes.table$Variants)*100,2)
  }

  if(input$variant == 'Yes2'){
  colnames(genes.table)<-c("Symbol","Meth Site", "Total", "Copy Number", "Methylation", "Variants","Missense_Mutation","Silent","Nonsense_Mutation","Splice_Site","RNA","Frame_Shift_Ins","Frame_Shift_Del","In_Frame_Ins","In_Frame_Del","Nonstop_Mutation","Translation_Start_Site")

  # Scale output to make it look nice
  genes.table$Variants<-round(as.numeric(genes.table$Variants)*100,2)
  genes.table$Missense_Mutation<-round(as.numeric(genes.table$Missense_Mutation)*100,2)
  genes.table$Silent<-round(as.numeric(genes.table$Silent)*100,2)
  genes.table$Nonsense_Mutation<-round(as.numeric(genes.table$Nonsense_Mutation)*100,2)
  genes.table$Splice_Site<-round(as.numeric(genes.table$Splice_Site)*100,2)
  genes.table$RNA<-round(as.numeric(genes.table$RNA)*100,2)
  genes.table$Frame_Shift_Ins<-round(as.numeric(genes.table$Frame_Shift_Ins)*100,2)
  genes.table$Frame_Shift_Del<-round(as.numeric(genes.table$Frame_Shift_Del)*100,2)
  genes.table$In_Frame_Ins<-round(as.numeric(genes.table$In_Frame_Ins)*100,2)
  genes.table$In_Frame_Del<-round(as.numeric(genes.table$In_Frame_Del)*100,2)
  genes.table$Nonstop_Mutation<-round(as.numeric(genes.table$Nonstop_Mutation)*100,2)
  genes.table$Translation_Start_Site<-round(as.numeric(genes.table$Translation_Start_Site)*100,2)


  }

 # Convert Data into percentages
  genes.table$Total<-round(as.numeric(genes.table$Total)*100,2)
  genes.table$`Copy Number`<-round(as.numeric(genes.table$`Copy Number`)*100,2)
  genes.table$Methylation<-round(as.numeric(genes.table$Methylation)*100,2)

  return(genes.table)

})

# Function will take GT() and return the proper form for plotting functions
GT.plotting <- reactive({
  genes.table<-GT()
  if(length(names(genes.table))<=6) return(genes.table) # Return if not too long

  genes.table<-genes.table[,c(1:6)]

})

whos.samples<-reactive({
  if(input$variant=='Yes') var=as.logical("T")
  if(input$variant=='No')  var=as.logical("F")
  if(input$variant=='Yes2') var=as.logical("T")
  if(input$cna.input == 'Gistic') cna=as.logical("F")
  if(input$cna.input == 'CNA') cna=as.logical("T")
  u.samples<-unique.samples(db1.select(),cna,var)
})

output$db.out<-renderText({
  out<-paste0("Database successfully loaded with ", length(whos.samples())," unique samples across ",dim(GT())[1], " unique genes.")
})

plot1.data <- reactive({
  name.tmp<-GT()
  gene.selected<-input$resultsGENE_rows_selected[1]
  if(is.null(input$resultsGENE_rows_selected)) gene.selected<-10
  name.gene<-name.tmp[gene.selected,1]
  name.cg<-name.tmp[gene.selected,2]

  u.samples<-whos.samples()
  tmp.RNAseq<-get.RNAseq(u.samples,name.gene,db1.select())

  if(input$cna.input=='Gistic')  tmp.Gistic<-get.Gistic(u.samples,name.gene,db1.select())
  if(input$cna.input=='CNA')     tmp.Gistic<-get.CNA(u.samples,name.gene,db1.select())
  tmp.Meth450<-get.Meth450(u.samples,name.gene,db1.select(),index=T,loc=T)
  meth.table<-dbGetQuery(db1.select(), paste0("SELECT * FROM ",table.select()$all, " WHERE [Searched Symbol] ==  '",name.gene,"'"))

   meth.rows<-meth.table$` Meth Location` %in% tmp.Meth450$cg_name # Keeps the data synced


  if(input$lm.type=='cor'){
  print("I'm in!")
    meth.table<-cbind(meth.table[meth.rows,c(2,3)],tmp.Meth450$location,meth.table[meth.rows,c(4:6)])
    colnames(meth.table)<-c("Symbol(s)","Meth Site","Location", "Copy Number", "Methylation", "Variants")
    meth.table$`Copy Number`<-round(as.numeric(meth.table$`Copy Number`)*100,2)
    meth.table$Methylation<-round(as.numeric(meth.table$Methylation)*100,2)

    if(input$variant=='Yes' | input$variant == 'Yes2'){
    tmp.Variant<-get.Variant(u.samples,name.gene,db1.select())
    meth.table$Variants<-round(as.numeric(meth.table$Variants)*100,2)
    }

    if(input$variant=='No') {
      tmp.Variant=NULL
      meth.table<-meth.table[,-6]
    }
    out<-list("RNAseq"=tmp.RNAseq,"Gistic"=tmp.Gistic,"Meth450"=tmp.Meth450,"Variant"=tmp.Variant,"GeneName"=name.gene,"Meth.Table"=meth.table,"Cg"=name.cg)
    return(out)
  }

  if(input$variant=='Yes' | input$variant == 'Yes2')  {
      meth.table<-cbind(meth.table[meth.rows,c(2,3)],tmp.Meth450$location,meth.table[meth.rows,c(4:7)])
      tmp.Variant<-get.Variant(u.samples,name.gene,db1.select())
      colnames(meth.table)<-c("Symbol(s)","Meth Site","Location", "Total", "Copy Number", "Methylation", "Variants")
      meth.table$Variants<-round(as.numeric(meth.table$Variants)*100,2)
  }

  if(input$variant=='No'){
    tmp.Variant=NULL
    meth.table<-cbind(meth.table[meth.rows,c(2,3)],tmp.Meth450$location,meth.table[meth.rows,c(4:6)])
    colnames(meth.table)<-c("Symbol(s)","Meth Site","Location", "Total", "Copy Number", "Methylation")

    }


  meth.table$Total<-round(as.numeric(meth.table$Total)*100,2)
  meth.table$`Copy Number`<-round(as.numeric(meth.table$`Copy Number`)*100,2)
  meth.table$Methylation<-round(as.numeric(meth.table$Methylation)*100,2)


  out<-list("RNAseq"=tmp.RNAseq,"Gistic"=tmp.Gistic,"Meth450"=tmp.Meth450,"Variant"=tmp.Variant,"GeneName"=name.gene,"Meth.Table"=meth.table,"Cg"=name.cg,"GeneRow"=gene.selected)

  return(out)})

output$resultsGENE <- DT::renderDataTable({
  DT::datatable(GT.plotting(),
                selection=list(mode='single', selected=which.max(GT.plotting()$Methylation)),
                options = list(pageLength = 20, lengthMenu=c(10,20,50,100), searching = T, sorting = T, order=list(3,'desc'))
          )})

output$methTable<-DT::renderDataTable({
  DT::datatable(as.data.frame(plot1.data()$Meth.Table),
                selection=list(mode='single', selected=which.max(plot1.data()$Meth.Table[,5])),
                options = list(pageLength = 20, searching = F, sorting = T, lengthMenu=c(10,20,30),order=list(5,'desc'))
  )})

output$VarTable<-DT::renderDataTable({
  search.me<-c("Symbol","Variants",VPlotChoices())
  data<-GT()[plot1.data()$GeneRow,]
  dataTable<-data[1,names(data) %in% search.me]
  DT::datatable(dataTable,
                selection=list(mode='single'),
                options = list(pageLength = 1, searching = F, sorting = F, lengthMenu=c(10,20,30),order=list(5,'desc'))
  )})

output$Test<-renderText({
  paste0("Plot data generated for ",plot1.data()$GeneName)
})

output$Scatter <- renderPlot({
    set.seed(23939209) # Make plots reproducible
  if(is.null(input$methTable_rows_selected)) return(NULL)
  gistic.plot<-as.numeric(plot1.data()$Gistic)
  if(input$cna.input == 'Gistic') gistic.plot<-gistic.plot + runif(length(gistic.plot),min=-0.25,max=0.25)

scatterhist(gistic.plot, plot1.data()$Meth450[input$methTable_rows_selected,-c(1,2,3)],col=greenred(length(plot1.data()$RNAseq))[rank(plot1.data()$RNAseq)], pch=20, xlab="Copy Number", ylab="Methylation", name=plot1.data()$GeneName)
zr=range(plot1.data()$RNAseq)
image.plot(legend.only=T, col= greenred(length(plot1.data()$RNAseq)) , zlim= zr, horizontal=T, smallplot= c(0,0.9,0,0.03),legend.lab="expression", legend.cex=0.75)

},height=800)


output$ScatterCN <- renderPlot({
  gistic.plot<-as.numeric(plot1.data()$Gistic)
  if(input$cna.input == 'Gistic') gistic.plot<-gistic.plot + runif(length(gistic.plot),min=-0.25,max=0.25)
plot(gistic.plot, as.numeric(plot1.data()$RNAseq),col="blue", xlab="Copy Number", ylab="RNA Expression")
cortext<-paste("Spearman Cor =",round(cor(as.numeric(plot1.data()$Gistic), as.numeric(plot1.data()$RNAseq),use="na.or.complete", method="spearman"),3))
#legend("topleft", cortext, bty="n", cex=1)
},height=400)

output$ScatterMeth <- renderPlot({
scatter.smooth(as.numeric(plot1.data()$Meth450[input$methTable_rows_selected,-c(1,2,3)]), as.numeric(plot1.data()$RNAseq),col="blue", xlab="Methylation", ylab="RNA Expression", span = 2/3, normalize=FALSE)
cortext<-paste("Spearman Cor =",round(cor(as.numeric(plot1.data()$Meth450[input$methTable_rows_selected,-c(1,2,3)]), as.numeric(plot1.data()$RNAseq),use="na.or.complete",method="spearman"),3))
legend("topleft", cortext, bty="n", cex=1)
},height=400)

output$ScatterMethGG<-renderPlot({
  if(is.null(input$methTable_rows_selected)) return(NULL)
  plot.me<-as.data.frame(cbind(as.numeric(plot1.data()$RNAseq),as.numeric(plot1.data()$Meth450[input$methTable_rows_selected,-c(1,2,3)])))
  colnames(plot.me)<-c("RNAseq","Methylation")
  ggplot(plot.me,aes(RNAseq,Methylation)) + geom_point(color="blue") + theme_classic() + stat_smooth(color="red",fill="grey50",span=0.9)

})

output$ScatterVar1<-renderPlot({
  var.tmp<-plot1.data()$Variant
  var.tmp<-var.tmp[rownames(var.tmp) %in% input$VP1 ,,drop=FALSE]
  if(!is.null(dim(var.tmp))) {
      if(dim(var.tmp)[1]>1){
        var.tmp<-colSums(var.tmp)# Nested to check for existance first!
        var.tmp<-t(var.tmp)
      }
    }

  if(length(var.tmp)==0) var.tmp<-matrix(data=rep(0,length(plot1.data()$RNAseq)),nrow=1)
  if(!is.null(input$VP1)) var.tmp[var.tmp>=1]<-2 # Remove absolute counts for indication of presense or absense

  gistic.data<-plot1.data()$Gistic
  set.seed(293939)
  if(input$cna.input == 'Gistic') gistic.data<-runif(length(gistic.data),-0.25,0.25)+gistic.data

  data<-cbind(plot1.data()$RNAseq,gistic.data,t(var.tmp))
  colnames(data)<-c("RNAseq","CN","Variant")
  data<-as.data.frame(data)
  # Plot Here
  g<-CN.plot(data,plot1.data()$GeneName)
  plot(g)
})

output$ScatterVar2<-renderPlot({
  var.tmp<-plot1.data()$Variant
  var.tmp<-var.tmp[rownames(var.tmp) %in% input$VP2 ,,drop=FALSE]
  if(!is.null(dim(var.tmp))) {
    if(dim(var.tmp)[1]>1){
      var.tmp<-colSums(var.tmp)# Nested to check for existance first!
      var.tmp<-t(var.tmp)
    }
  }

  if(length(var.tmp)==0) var.tmp<-matrix(data=rep(0,length(plot1.data()$RNAseq)),nrow=1)
  if(!is.null(input$VP2)) var.tmp[var.tmp>=1]<-2 # Remove absolute counts for indication of presense or absense

  gistic.data<-plot1.data()$Gistic
  set.seed(293939)
  if(input$cna.input == 'Gistic') gistic.data<-runif(length(gistic.data),-0.25,0.25)+gistic.data

  data<-cbind(plot1.data()$RNAseq,gistic.data,t(var.tmp))
  colnames(data)<-c("RNAseq","CN","Variant")
  data<-as.data.frame(data)
  # Plot Here
  g<-CN.plot(data,plot1.data()$GeneName)
  plot(g)
})

output$VariantPlots1<-renderUI({
  checkboxGroupInput("VP1","Select Variants to Plot", choices=VPlotChoices(),inline=TRUE)
})

output$VariantPlots2<-renderUI({
  checkboxGroupInput("VP2","Select Variants to Plot", choices=VPlotChoices(),inline=TRUE)
})
# Function will determine which choices are available for Variant PLots.  Used twice so, a F(x)
 VPlotChoices<-reactive({

   if(is.null(plot1.data()$Variant)) {
     choices=c("Not Available with Current Selections")
     return(choices)} #Short circuit F(x) if no Variant Info Available
  tmp.choice<-plot1.data()$Variant[-12,] # Remove All Variants as this is handled through checkboxes

   choices<- rownames(tmp.choice)[rowSums(tmp.choice)>0]

  if(length(choices)==0) choices=c("Not Available with Current Selections")

  return(choices)
})

#UI element that allows the choice of signature
output$VPchoices<-renderUI({
 # variant.data<-GT()
  choose.me<-c(unique(sigs.table[sigs.table$Tissue=="Breast","Signature"]),"List of 4 or more genes","Gene Correlation")
  selectInput("sigSelect","Gene Set", choices=choose.me,multiple=FALSE)

})

#session$onSessionEnded(stopApp)
})
