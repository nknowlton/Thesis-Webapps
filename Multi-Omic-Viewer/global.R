# Define a set of work functions to abstract out the Database Calls

unique.samples<-function(db,CNA=F,var=T){
  
  #print("Generating Sample List")
  samp.rnaseq<-dbListFields(db,"RNAseq")
  samp.meth450<-dbListFields(db,"Meth450")
  if(!CNA) samp.cn<-dbListFields(db,"Gistic")
  if(CNA) samp.cn<-dbListFields(db,"CNA_SNP")
  
  if(var) {
    samp.var<-dbGetQuery(db,'SELECT DISTINCT barcode FROM Variant')
    unique.samps<-Reduce(intersect,list(samp.rnaseq,samp.meth450,samp.cn,samp.var[,1]))
  }
  
  if(!var) unique.samps<-Reduce(intersect,list(samp.rnaseq,samp.meth450,samp.cn))
  # Pull out only intersection
  
  unique.samps<-unique.samps[unique.samps!="Gene_Symbol"]
  
  print(paste("There are ",length(unique.samps)," unique samples identified in this dataset.",sep=""))
  
  return(unique.samps)
  
}

unique.genes<-function(db){
  gene.var<-dbGetQuery(db,'SELECT DISTINCT Hugo_Symbol FROM Variant')
  gene.rnaseq<-dbGetQuery(db,'SELECT DISTINCT Gene_Symbol FROM RNAseq')
  gene.cn<-dbGetQuery(db,'SELECT DISTINCT "Gene_Symbol"FROM Gistic')
  unique.genes<-Reduce(intersect,list(gene.var[,1],gene.rnaseq[,1],gene.cn[,1]))
  print(paste("There are ",length(unique.genes)," unique genes identified in this dataset.",sep=""))
  return(unique.genes)
}


# start with RNAseq

#test.rnaseq<-exp.dat["ERBB2",]

#test.sql<-dbGetQuery(db,"Select * FROM RNAseq WHERE Gene_Symbol  == 'ERBB2'") # data is returned as un logged un Z transformed.

get.RNAseq<-function(unique.samps,gene.name,db){
  if(length(gene.name)>1) stop("You are trying to pass more than one gene name!")
  
  data<-dbGetQuery(db,paste0("Select [",paste0(unique.samps,collapse="], ["), "] FROM RNAseq WHERE Gene_Symbol  == '",gene.name,"'"))
  
  if(dim(data)[1]>1) {
    data<-data[which.max(apply(data,1,sd)),]
    print("Multiple matches for RNAseq Gene Symbol")
  }
  
  log2.data<-log2(as.numeric(data)+1)
  if(sd(log2.data,na.rm=T)!=0){
    zdata<-(log2.data-mean(log2.data,na.rm=T))/sd(log2.data,na.rm=T)
    data.out<-as.matrix(zdata)
    rownames(data.out)<-unique.samps
    colnames(data.out)<-paste0("RNASeq.",gene.name)
    return(data.out)} else{
      return(NULL)}
}
# Test
#get.RNAseq(unique.samps,"ERBB2",db))

# Methylation 450k
get.Meth450<-function(unique.samps,gene.name,db,index=T,loc=F) {
  if(length(gene.name)>1){ stop("You are trying to pass more than one gene name!")}
  
  if(!exists("meth450_location")) meth450_location<<-dbReadTable(db,"Meth450_location")
  
  
  
  if(index){
    #which cg
    if(!exists("meth450_lookup")) meth450_lookup<<-dbReadTable(conn=db,name="Meth450_lookup")
    
    
    cg.names<-meth450_lookup[meth450_lookup$Gene_Symbol==gene.name,1]
    if(length(cg.names)[1] == 0) index=as.logical("F")
    data<-dbGetQuery(db,paste0("SELECT cg_name, Gene_Symbol, [",paste0(unique.samps,collapse="], ["), "] FROM Meth450 WHERE cg_name IN ('",paste0(as.character(unlist(cg.names)),collapse="', '"),"') and [",unique.samps[1],"] != 'NA'"))
  }
  
  if(!index){
    data<-dbGetQuery(db,paste0("SELECT cg_name, Gene_Symbol, [",paste0(unique.samps,collapse="], ["), "] FROM Meth450 WHERE Gene_Symbol LIKE '%",gene.name,"%' and [",unique.samps[1],"] != 'NA'"))
    
    # Need to Fuzzy match in SQL and now will only keep exact matches
    chk<-strsplit(data$Gene_Symbol,";")
    remove<-unlist(lapply(chk,function(x) sum(gene.name %in% x)))
    if(sum(remove)==0) return(NULL) # nothing matches so get out of here!
    
    index<-1:length(data$Gene_Symbol)
    tmp<-cbind(index,remove)
    keep<-tmp[tmp[,2]==1,]
    
    # check if keep changes class to integer
    if(is.matrix(keep)) {
      data<-data[keep[,1],]
    }
  }
  
  if(dim(data)[1] == 0) return(NULL)
  # end no index
  values<-subset(data,select=-c(cg_name,Gene_Symbol))
  m.data<-log2(values/(1-values)) # Convert Beta to M values
  scale.m.data<-t(scale(t(m.data)))  # transpose for function to work down columns
  cg_name<-data$cg_name
  Gene_Symbol<-data$Gene_Symbol
  # merge the Gene Location
  if(loc & index){
    location<-meth450_location[match(data$cg_name,meth450_location$cg_name),2]
    out.data<-cbind(cg_name,Gene_Symbol,location,as.data.frame(scale.m.data))
  } else {
    out.data<-cbind(cg_name,Gene_Symbol,as.data.frame(scale.m.data))
  }
  
  return(out.data)
}
# Test
# get.Meth450(unique.samps,"ERBB2",db)

# Copy Number

get.Gistic<-function(unique.samps,gene.name,db) {
  if(length(gene.name)>1) stop("You are trying to pass more than one gene name!")
  
  data<-dbGetQuery(db,paste0("SELECT [",paste0(unique.samps,collapse="], ["), "] FROM Gistic WHERE Gene_Symbol  == '",gene.name,"'"))
  out.data<-t(data)
  colnames(out.data)<-gene.name
  return(out.data)
}
#Test
# get.Gistic(unique.samps,"ERBB2",db)

get.Variant<-function(unique.samps,gene.name,db,verbose=F){
  if(length(gene.name)>1) stop("You are trying to pass more than one gene name!")
  
  # Get back only those samples and genes wanted
  data<-dbGetQuery(db,paste0("SELECT Hugo_Symbol, Variant_Classification, barcode FROM Variant WHERE barcode IN ('",paste0(unique.samps,collapse="', '"),"') AND Hugo_Symbol == '",gene.name,"'"))
  # Create a variant matrix for plotting
  r.names<-c("Missense_Mutation","Silent","Nonsense_Mutation","Splice_Site","RNA","Frame_Shift_Ins","Frame_Shift_Del","In_Frame_Ins","In_Frame_Del","Nonstop_Mutation","Translation_Start_Site","Any_Mutation")
  data.out<-matrix(data=0, nrow=12,ncol=length(unique.samps),dimnames=list(r.names,unique.samps))
  for(i in r.names) {
    loc<-grep(i,data$Variant_Classification)
    data.tmp<-data$barcode[loc]
    if(length(loc)>0) {
      for(j in 1:length(loc)) {
        data.out[i,data.tmp[j]]<-data.out[i,data.tmp[j]]+1
      }}
  }
  
  a<-colSums(data.out)
  
  data.out[12,]<-ifelse(a>1,1,a)
  if(verbose==T)  print(paste0("Mutations found in ",sum(data.out[12,])," Samples"))
  return(data.out)
  
}

get.CNA<-function(unique.samps,gene.name,db) {
  if(length(gene.name)>1) stop("You are trying to pass more than one gene name!")
  
  data<-dbGetQuery(db,paste0("SELECT [",paste0(unique.samps,collapse="], ["), "] FROM CNA_SNP WHERE Gene_Symbol  == '",gene.name,"'"))
  out.data<-t(data)
  out.data<-(2^out.data)*2
  colnames(out.data)<-gene.name
  return(out.data)
}


# define plotting function
scatterhist = function(x, y, xlab="", ylab="", col, pch, name){
  x=as.numeric(x)
  y=as.numeric(y)
  zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
  layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
  xhist = hist(x, plot=FALSE, breaks=100)
  yhist = hist(y, plot=FALSE, breaks=100)
  top = max(c(xhist$counts, yhist$counts))
  par(mar=c(3,3,1,1))
  plot(x,y, col=col,pch=pch, cex=2)
  legend("topright", name, bty="n", cex=2)
  par(mar=c(0,3,1,1))
  barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0)
  par(mar=c(3,0,1,1))
  barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE)
  par(oma=c(3,3,0,0))
  mtext(xlab, side=1, line=1, outer=TRUE, adj=0,
        at=0.35, padj=1.55)
  mtext(ylab, side=2, line=1, outer=TRUE, adj=0,
        at=0.35, padj=-1.55)
}

# Copy Number plotting function will take data frame and plot using GGPLOT2
# Used as F(x) as reactives don't allow variable passing
CN.plot<-function(data,name){
  # data has fill column
  range.out<-c(1,3)
  if(max(data$Mutations)==0) range.out<-c(1,1)
  
  
  p<-ggplot(data,aes(y=RNAseq,x=CN))+geom_point(aes(colour = factor(Mutations),size=Mutations),show.legend=FALSE)+scale_size(range=range.out)+labs(y="RNAseq", x="Copy Number")+scale_colour_manual(name="",values=c("2"="darkorchid4","0"="bisque3"))+annotate("text",y=max(data$RNAseq), x=max(data$CN),label=name,hjust=1)+theme_light()+theme(panel.grid.minor=element_blank())
  return(p)
}


