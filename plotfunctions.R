lambda <- function(results) {
	median((results$BETA/results$SE)^2,na.r=T)/qchisq(0.5,1)
}

QQplot <- function(results, title="", ...) {
	plot(x=-log10(1:length(na.omit(results$P))/(length(na.omit(results$P))+1)),y=-log10(sort(na.omit(results$P))),
		xlab="-log10(Expected)",ylab="-log10(Observed)",pch=".",main=title,abline(0,1,col=1,lty=3), ...) 
}

manh.plot<-function (results, base = 10, hits=5, cutoffs = seq(0,100,by=2), color = NULL, title="", ...) 
{
    results$Chr <- as.integer(as.character(results$CHR))
    results <- results [order (results$Chr,results$BP),]
    chr <- results$Chr
	labs <- sort(unique(chr))
    pos <- results$BP
    p <- results$P
    affy <- as.vector(table(chr))
    CM <- cumsum(affy)
    n.markers <- sum(affy)
    n.chr <- length(affy)
    id <- 1:n.chr
    # the following were in the original code, but I have taken them out
    # eps <- .Machine$double.eps
    # dp <- seq(eps, 1, length = n.markers)
    # y <- -log(dp, base)  
    y <- -log(p, base)
    par(xaxt = "n", yaxt = "n")
    if (title != "") {
        title <- sub("FF", "FEV[1]/FVC", title)
        title <- sub("FEV1", "FEV[1]", title)
        title <- gsub(" ", "~", title)
    }
    plot(1:max(CM), y+0.1, type = "n", xlab = "chromosome", ylab = "-log(P-value)", axes = FALSE,main=parse(text=title),...)
    if (is.null(color)) 
        color <- rep(c(4,8), 13)
    
   for (i in 1:n.chr) {
        u <- CM[i]
        l <- CM[i] - affy[i] + 1
        cat("Plotting points ", l, "-", u, "\n")
        chr <- l:u
        y <- -log(p[chr],base=base)
	r <- subset(results,results$Chr==i & results$P< base^(-hits), select=c(SNP))
        points(l:u, y, pch= ifelse(y>hits,23,21), bg=ifelse(y>hits,"red",color[i]), col=ifelse(y>hits,"red",color[i]), ...)
	#text(l:u,y+0.15,labels= ifelse(y>hits,r$Markerid,""))
    }
    par(xaxt = "s", yaxt = "s")
    axis(1, at = CM-affy/2, labels = labs, cex.lab=0.3, tick = FALSE)
    axis(2, at = cutoffs, tick = TRUE)
    box()
}

efficient.plot <- function (results, base = 10, min.interesting=2, hits=5, cutoffs = seq(0,100,by=2), color = NULL, title="", ...) 
{
    require(data.table)
    results <- data.table(results)
    results[, CHR:=as.integer(CHR) ]
    results$idx <- 1:nrow(results)
    results <- results[ order(CHR, BP) ]
    chr <- results$CHR
	labs <- sort(unique(chr))
    pos <- results$BP
    p <- results$P
    affy <- as.vector(table(chr))
    CM <- cumsum(affy)
    n.markers <- sum(affy)
    n.chr <- length(affy)
    id <- 1:n.chr
    y <- -log(p, base)
    par(xaxt = "n", yaxt = "n")
    plot(1:max(CM), y+0.1, type = "n", xlab = "chromosome", ylab = "-log(P-value)", axes = FALSE,main=title,...)
    if (is.null(color)) 
        color <- rep(c(4,8), 13)
    
    for (i in 1:n.chr) {
        u <- CM[i]
        l <- CM[i] - affy[i] + 1
        cat("Plotting points ", l, "-", u, "\n")
        rect(l, 0, u, min.interesting, col=color[i], border=NA)
        r <- subset(results, CHR == i & P < base^(-min.interesting))
        x <- r$idx
        y <- -log10(r$P)
        points(x, y, pch= ifelse(y>hits,23,21), bg=ifelse(y>hits,"red",color[i]), col=ifelse(y>hits,"red",color[i]), ...)
    }
    par(xaxt = "s", yaxt = "s")
    labs <- sub("23", "X", labs)
    labs <- sub("24", "Y", labs)
    labs <- sub("26", "MT", labs)
    axis(1, at = CM-affy/2, labels = labs, cex.axis=1.5, tick = FALSE)
    axis(2, at = cutoffs, cex.axis=1.5, tick = TRUE)
    box()
}

getgene <- function(x) {
    sub("\\(.*", "", x)
}

getdist <- function(x) {
    as.integer(sub(".*\\(dist=(\\d+)\\)", "\\1", x))
}

selectgene <- function(annot.dt) {
#    annot.dt <- annot.dt[, 1:4, with=F ]
#    setnames(annot.dt, c("location", "genes", "chrom", "pos"))
#    setkeyv(annot.dt, c("chrom", "pos"))
    genes <- annot.dt[, strsplit(genes, ",") ]
    tmp.dt <- data.table(gene1=t(genes[1,]), gene2=t(genes[2,]))
    suppressWarnings(tmp.dt[, `:=`(gene1=getgene(gene1.V1), gene2=getgene(gene2.V1), dist1=getdist(gene1.V1), dist2=getdist(gene2.V1)) ])
    tmp.dt[, gene:=gene1 ]
    tmp.dt[ dist2 < dist1, gene:=gene2 ]
    annot.dt[, gene:=tmp.dt$gene ]
}

manh.anno <- function(results, annot, previous, phe, range=1e6, ylim=NA, title="", ...) {
    suppressMessages(require(data.table))
    results.dt <- data.table(results)
    setkeyv(results.dt, c("CHR", "BP"))
    if (any(is.na(ylim)))
        ylim <- results.dt[, c(min(-log10(P)), max(-log10(P)) + 2) ]
#    manh.plot(results, ylim=ylim, title=title, ...)
    efficient.plot(results, ylim=ylim, title=title, ...)
    annot.dt <- data.table(annot)
    annot.dt <- selectgene(annot.dt)
    results.dt[, idx:=1:.N ]
#    annot.dt[, plotx:=results.dt[ J(chrom, pos), idx ] ]
#    annot.dt[, ploty:=-log10(results.dt[ plotx, P ]) ]
    setkeyv(annot.dt, c("chrom", "pos"))
    annot.dt <- results.dt[ annot.dt ][ order(CHR, BP, P), .SD[1], c("CHR", "BP") ][, .(chrom=CHR, pos=BP, plotx=idx, ploty=-log10(P), gene) ]

    if (missing(previous)) {
        colour <- 1
    } else {
        annot.dt <- flagprevious(annot.dt, previous, phe, range=range)
        legend(0,ylim[2],legend=c("Novel", "Known, same trait", "Known, different trait"), text.col=c(1,2,4))
    }
    null <- annot.dt[, text(plotx, ploty, paste(" ", gene), srt=90, col=colour, font=3, pos=4, offset=0) ]
}

flagprevious <- function(annot, prevloci, phe, range=1e6) {
    suppressMessages(require(data.table))
    suppressMessages(require(GenomicRanges))

    prevloci.GRanges <- makeGRangesFromDataFrame(prevloci, seqnames.field="CHROM", start.field="START", end.field="END", ignore.strand=T, keep.extra.columns=T)
    annot.GRanges <- makeGRangesFromDataFrame(annot, seqnames.field="chrom", start.field="pos", end.field="pos", ignore.strand=T)
    prev.GRanges <- suppressWarnings(subsetByOverlaps(prevloci.GRanges, annot.GRanges))
    idx <- annot[, suppressWarnings(nearest(GRanges(chrom, IRanges(pos, pos)), prev.GRanges)), c("chrom", "pos") ]$V1
    prev.dt <- data.table(as.data.frame(prev.GRanges))
    annot <- cbind(annot, prev.dt[ idx, list(PHE, CHROM=as.integer(as.character(seqnames)), POS, gene) ])
    annot[ abs(pos - POS) > range, PHE:=NA ]
    write.table(annot, "", row.names=F, quote=F, sep="\t")
    if (!missing(phe)) {
        annot[, colour:=ifelse(is.na(PHE), 1, ifelse(PHE==phe, 2, 4)) ]
    }
    annot
}
