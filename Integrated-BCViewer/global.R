# create a helper script to do come of the data handling

pcicout<-function(gdat) {

	 a <- fastICA(t(gdat), n.comp=4, alg.typ = "deflation", method = "C", row.norm = FALSE, maxit = 200, verbose = FALSE)
	 
	 pca<-a$X %*% a$K
	 output<-cbind(pca,colMeans(gdat),a$S) 
	 colnames(output)<-c("PC1","PC2","PC3","PC4","Centroid","IC1","IC2","IC3","IC4")
	 return(output)}
	 
hoc.cent<-function(expDat){
	load("hoc symbols.RData")
	hoc.symbols<-hoc.symbols[,-1]
	output<-matrix(NA,nrow=10,ncol=dim(expDat)[2])
	
	for(i in 1:10){
	loc<-match(as.character(hoc.symbols[,i]),rownames(expDat))
	output[i,]<-colMeans(expDat[loc,],na.rm=TRUE)}
	
rownames(output)<-colnames(hoc.symbols)
colnames(output)<-colnames(expDat)
return(output) }


auc.search<-function(plot.matrix,kp,dmfsTime,dmfsEvent,time,nvars=2){
  who<-c("PC1","PC2","PC3","PC4","Centroid","IC1","IC2","IC3","IC4","Score")
  check.list<-combn(who[1:dim(plot.matrix)[2]],nvars)
  
  store<-rep(0,dim(check.list)[2])
  
  for(i in 1:length(store)){
    tmp.marker<-rowSums(plot.matrix[,check.list[,i]])
    store[i]<-survivalROC.NSK(Stime=dmfsTime[kp],status=dmfsEvent[kp],marker=tmp.marker[kp],predict.time=time,span=0.25*length(kp)^(-0.2))$AUC
  }
  store[store<0.5]<-1-store[store<0.5]
  return(check.list[,which.max(store)])
}
	
auc.plot.shiny<-function(ss,kp,dmfsTime,dmfsEvent,name,time){
  
  # The numbers in the CoxPH legend need to account for missing survival information,
  # otherwise the group sizes include samples with NA values for either or both of
  # dmfsTime and dmfsEvent.
  kp <- kp & !is.na(dmfsTime) & !is.na(dmfsEvent)
  auc.out<-survivalROC.NSK(Stime=dmfsTime[kp],status=dmfsEvent[kp],marker=ss[kp],predict.time=time,span=0.25*length(kp)^(-0.2))
 
  # plot survival curves
  plot(auc.out$FP,auc.out$TP, type="l",col="blue",xlim=c(0,1),ylim=c(0,1), 
       xlab=paste("FP", "\n", "AUC= ",round(auc.out$AUC,3)," or ",round(1-auc.out$AUC,3)),
       ylab="TP",
       main=paste("AuC Plot of",name))
  abline(0,1)  
} 


surv.plot.shiny<-function(ssx,ss,kp,dmfsTime,dmfsEvent,main,p.cox,auc,break.here=FALSE){
  # The numbers in the CoxPH legend need to account for missing survival information,
  # otherwise the group sizes include samples with NA values for either or both of
  # dmfsTime and dmfsEvent.
  ssx1 <- ssx[!is.na(dmfsTime[kp]) & !is.na(dmfsEvent[kp])] # must run before kp is.na below
  kp <- kp & !is.na(dmfsTime) & !is.na(dmfsEvent) # kp is reduced from here on out

  if(break.here) browser()
   
      # log rank test for specific cutoff
  cpv<-1-pchisq(survdiff(Surv(dmfsTime[kp],dmfsEvent[kp])~ssx1)$chisq,1)
  # Avoid scientific notation kicking in below 0.001
  cpvx<-paste("Logrank p=",formatC(cpv, format = "f", digits=4),sep='')
  # Avoid p-values of zero after rounding
  if(cpv<0.0001) cpvx<-"Logrank p<0.0001"
  

  
  # Cox PH p calculation
  if(p.cox){
    fit.cox<-coxph(Surv(dmfsTime[kp],dmfsEvent[kp])~ss[kp])
    hr<-round(summary(fit.cox)$coef[1,2],4)
    cint<-round(exp(confint(fit.cox)),4)
    coxpv<-round(summary(fit.cox)$coef[1,5],4)
    if(coxpv<0.0001) coxpv<-"<0.0001"
    main.out<-paste(main,"\n\n","CoxPH HR ",hr," 95% CI (",cint[1],", ",cint[2],")", "\nP-",coxpv, "\n")
  }

   if(p.cox==FALSE){
     main.out<-main
   }
    # plot survival curves
    plot(survfit(Surv(dmfsTime[kp],dmfsEvent[kp])~ssx1),col=2:3,lwd=3,xlim=c(0,12),main=main.out,xlab="Time (years)",ylab="Proportion metastasis-free")
    # generate legend

    lg<-paste(names(table(ssx))," (",table(ssx1),")",sep='')
    legend.title<-cpvx 
    legend(0.5,0.2,lg,fill=2:3,title=legend.title,cex=1.3)
    
} 

heatmap.mik<-function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                       distfun = dist, hclustfun = hclust, dendrogram = c("both",
                                                                          "row", "column", "none"), symm = FALSE, scale = c("none",
                                                                                                                            "row", "column"), na.rm = TRUE, revC = identical(Colv,
                                                                                                                                                                             "Rowv"), add.expr, breaks, col = "heat.colors", colsep,
                       rowsep, sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote,
                       notecex = 1, notecol = "cyan", na.color = par("bg"), trace = c("column",
                                                                                      "row", "both", "none"), tracecol = "cyan", hline = median(breaks),
                       vline = median(breaks), linecol = tracecol, margins = c(5,
                                                                               5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr),
                       cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL,
                       key = TRUE, keysize = 1.5, density.info = c("histogram",
                                                                   "density", "none"), denscol = tracecol, symkey = min(x <
                                                                                                                          0, na.rm = TRUE), densadj = 0.25, main = NULL, xlab = NULL,
                       ylab = NULL, ...)
{
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if ((Colv == "Rowv") && (!isTRUE(Rowv) || is.null(Rowv)))
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    x <- sweep(x, 1, rowMeans(x, na.rm = na.rm))
    sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    x <- sweep(x, 2, colMeans(x, na.rm = na.rm))
    sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) <
      1)
    if (missing(col))
      breaks <- 16
  else breaks <- length(col) + 1
  if (length(breaks) == 1) {
    breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                  length = breaks)
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  else if (is.character(col) && length(col) == 1)
    col <- do.call(col, list(ncol))
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[] <- ifelse(x < min.breaks, min.breaks, x)
  x[] <- ifelse(x > max.breaks, max.breaks, x)
  lmat <- rbind(4:3, 2:1)
  lhei <- lwid <- c(keysize, 4)
  if (!missing(ColSideColors)) {
    #        if (!is.character(ColSideColors) || length(ColSideColors) !=
    #            nc)
    #            stop("'ColSideColors' must be a character vector of length ncol(x)")
    if(is.null(nrow(ColSideColors))){
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      lhei <- c(lhei[1], 0.2, lhei[2])
    }
    else{
      lmat <- rbind(lmat[1, ] + nrow(ColSideColors),
                    t(matrix(as.numeric(unlist(strsplit(paste("NA",1:nrow(ColSideColors))," "))),2,nrow(ColSideColors))),
                    lmat[2, ] + nrow(ColSideColors))
      #          lhei <- c(lhei[1], 0.2*nrow(ColSideColors), lhei[2])
      lhei <- c(lhei[1], rep(0.2,nrow(ColSideColors)), lhei[2])
    }
  }
  if (!missing(RowSideColors)) {
    if (!is.character(RowSideColors) || length(RowSideColors) !=
        nr)
      stop("'RowSideColors' must be a character vector of length nrow(x)")
    lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1),
                                   1), lmat[, 2] + 1)
    lwid <- c(lwid[1], 0.2, lwid[2])
  }
  lmat[is.na(lmat)] <- 0
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  if (!missing(RowSideColors)) {
    par(mar = c(margins[1], 0, 0, 0.5))
    image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
  }
  if (!missing(ColSideColors)) {
    if(is.null(nrow(ColSideColors))) {
      par(mar = c(0.5, 0, 0, margins[2]))
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    }
    else{
      for(kk in 1:nrow(ColSideColors)){
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[kk,colInd], axes = FALSE)
        axis(4, 0, labels = rownames(ColSideColors)[kk], las = 2,
             line = -0.5, tick = 0, cex.axis = cexCol)
      }
    }
  }
  par(mar = c(margins[1], 0, 0, margins[2]))
  if (!symm || scale != "none") {
    x <- t(x)
    cellnote <- t(cellnote)
  }
  if (revC) {
    iy <- nr:1
    ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 +
          c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col,
        breaks = breaks, ...)
  #    if (!invalid(na.color) & any(is.na(x))) {
  mmat <- ifelse(is.na(x), 1, NA)
  image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
        col = na.color, add = TRUE)
  #    }
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
       cex.axis = cexCol)
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
       cex.axis = cexRow)
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0,
                                                                length(csep)), xright = csep + 0.5 + sepwidth[1],
                              ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1,
                              col = sepcolor, border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) +
                                                      1 - rsep) - 0.5, xright = ncol(x) + 1, ytop = (ncol(x) +
                                                                                                       1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1,
                              col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    for (i in colInd) {
      if (!is.null(vline)) {
        vline.vals <- scale01(vline, min.scale, max.scale)
        abline(v = i - 0.5 + vline.vals, col = linecol,
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    for (i in rowInd) {
      if (!is.null(hline)) {
        hline.vals <- scale01(hline, min.scale, max.scale)
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    if (symkey) {
      max.raw <- max(abs(x), na.rm = TRUE)
      min.raw <- -max.raw
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = breaks,
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row")
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column")
      mtext(side = 1, "Column Z-Score", line = 2)
    else mtext(side = 1, "Value", line = 2)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    }
    else title("Color Key")
  }
  else plot.new()
  invisible(list(rowInd = rowInd, colInd = colInd))
}

survivalROC.NSK<-function (Stime, status, marker, predict.time, span = 0.05) 
{
  PredictTime = predict.time
  x <- marker
  drop <- is.na(Stime) | is.na(status) | is.na(x)
  if (sum(drop) > 0) {
    #ndrop <- sum(drop)
    #cat(paste("\n cases dropped due to missing values:", 
    #          ndrop, "\n\n"))
    Stime <- Stime[!drop]
    status <- status[!drop]
    x <- x[!drop]
  }
  unique.Stimes <- unique(Stime[status == 1])
  unique.Stimes <- unique.Stimes[order(unique.Stimes)]
  unique.x <- unique(x)
  unique.x <- unique.x[order(unique.x)]
  ooo <- order(x)
  Stime <- Stime[ooo]
  status <- status[ooo]
  x <- x[ooo]
  n <- length(x)
  p <- length(unique.Stimes)
  q <- length(unique.x)
  TP <- rep(0, q)
  FP <- rep(0, q)
  SurvT <- 0
  z <- .C("survivalROC", as.double(Stime), as.double(status), 
          as.double(unique.Stimes), as.double(x), as.double(unique.x), 
          as.double(PredictTime), SurvT = as.double(SurvT), as.double(span), 
          TP = as.double(TP), FP = as.double(FP), as.integer(n), 
          as.integer(p), as.integer(q), PACKAGE = "survivalROC")
  area <- function(x, y) {
    if (NROW(x) != NROW(y)) {
      print("ERROR: variables of different length")
      return(0)
    }
    else {
      nr <- NROW(x)
      nr1 <- nr - 1
      a <- sum(y[1:nr1] * (x[2:nr] - x[1:nr1]))
      b <- sum(y[2:nr] * (x[2:nr] - x[1:nr1]))
      temp <- (a + b)/2
      return(temp)
    }
  }
  TP <- c(1, z[["TP"]])
  FP <- c(1, z[["FP"]])
  AUC = area(FP[NROW(FP):1], TP[NROW(TP):1])
  out <- list(cut.values = c(-Inf, unique.x), TP = c(1, z[["TP"]]), 
              FP = c(1, z[["FP"]]), predict.time = PredictTime, Survival = z[["SurvT"]], 
              AUC = AUC)
  return(out)
}
