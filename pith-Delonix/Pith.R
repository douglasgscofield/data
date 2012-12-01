###########################################################
#                                                         #
# functions for analysis of cell divisions per meter data #
#                                                         #
###########################################################

reload <- function(doit=FALSE) if (doit) source("Pith.R")

library(nlme)
library(ape)
options(prompt="Pith> ")

all.branches <- c("A","B","C","D")

###########################################################
#                                                         #
# variance components by hand                             #
#                                                         #
###########################################################

vca1 <- function(data)
{
    data$tree <- factor(data$tree)
    attach(data)
    
    detach(data)
}

###########################################################
#                                                         #
# sample size calculations                                #
#                                                         #
###########################################################

pith.figureAppendix <- function(ssdf, do.postscript=FALSE)
{
    fig.file <- "pith_samplesizes.eps"
    fig.width <- 7.5
    fig.height <- 2.5
    if (do.postscript) postscript(file=fig.file, width=fig.width, 
        height=fig.height, onefile=F, horizontal=F, paper="special")
    nnn <- sort(unique(ssdf$n))
    ss <- sort(unique(ssdf$s))
    bb <- sort(unique(ssdf$b))
    par(mfcol=c(1,length(nnn)), oma=c(3.5,2,2.5,1), mar=c(0,1,0,0), ps=9)
    do.first <- TRUE
    for (nn in nnn) {
        sssdf <- subset(ssdf, n == nn)
        mm <- matrix(sssdf$T, nrow=max(bb), ncol=max(ss), byrow=TRUE)
        contour(bb, ss, mm, axes=FALSE, frame.plot=TRUE, 
            xlim=range(bb), ylim=range(ss),
            levels=seq(2000,3000,by=25),
            xlab="", ylab="", main=""
            )
        text(max(bb)*0.9, max(ss)*0.95, cex=1.5, 
            as.expression(substitute(italic(n)==n.n, list(n.n=as.character(nn)))))
        abline(h=2*1:5, v=2*1:5, col="black", lty=2, lwd=0.1)
        if (do.first) {
            axis(2, las=2)
            axis(2, at=ss, tcl=-0.25, labels=FALSE)
            do.first <- FALSE
        }
        axis(1)
        axis(1, at=bb, tcl=-0.25, labels=FALSE)
    }
    mtext(side=1, expression("Number of Branches per Tree "*(italic(b))), outer=TRUE, line=2.5)
    mtext(side=2, expression("Number of Segments per Branch "*(italic(s))), outer=TRUE, line=0.9)
    mtext(side=3, expression("Sample-Size Effect on Expected Among-Tree Deviance "*italic(T)==sqrt(italic(T)^2)*", by Number of Measurements per Segment "*(italic(n))), outer=TRUE, line=1)

    if (do.postscript) dev.off()
}

pith.ssdf <- function(fit, bb=1:10, ss=1:10, nn=c(5,10,15,20))
{
    vc <- varcomp(fit)  # tree branch segment within-segment
    ssdf <- data.frame(b=0, s=0, n=0, T2=0, T=0, V2=0, V=0)[0,]
    for (n in nn) {
        for (b in bb) {
            for (s in ss) {
                v2 <- V2(vc=vc, b, s, n)
                t2 <- T2(vc=vc, b, s, n)
                ssdf <- rbind(ssdf, list(b=b, s=s, n=n, T2=t2, T=sqrt(t2), 
                    V2=v2, V=sqrt(v2)))
            }
        }
    }
    ssdf
}

T2 <- function(vc, b=4, s=2.1, n=9.5) 
{
    ans <- vc[4]/(b*s*n) + vc[3]/(s*n) + vc[2]/n + vc[1]
    names(ans) <- "T2"
    ans
}

V2 <- function(vc, b, s, n) 
{
    ans <- vc[4] + n*vc[3] + s*n*vc[2] + b*s*n*vc[1]
    names(ans) <- "V2"
    ans
}

###########################################################
#                                                         #
# plots                                                   #
#                                                         #
###########################################################

pith.figure2 <- function(data, include.marginal=FALSE, do.postscript=FALSE) 
{
    if (include.marginal) {
        fig.file <- "Delonixtwigs_plot_histogram_all.eps" 
        fig.width <- 3.5; fig.height <- 1.5
    } else {
        fig.file <- "Delonixtwigs_plot_histogram.eps" 
        fig.width <- 3.5; fig.height <- 1.5
    }
    if (do.postscript) postscript(file=fig.file, width=fig.width, 
        height=fig.height, onefile=F, horizontal=F, paper="special")
    xlabel=expression(paste("Pith cells per meter (",'' 
        %*% 10^3,")"))
    cd <- subset(data, !is.na(data$celldm))
    cd$celldm <- cd$celldm / 10^3
    celldm <- split(cd$celldm, cd$file)
    par(mar=c(2.5,3.2,0,0), mgp=c(2,0.5,0), tcl=-0.25, cex=0.7, las=1)
    # medial histogram
    dat <- celldm$medial
    n.dat <- length(dat)
    h <- hist(dat, breaks=25, plot=FALSE,
        xlim=if (include.marginal) range(cd$celldm) else range(dat))
    my.plot.histogram(h, col="gray", main="", xlab="", ylab="", cex=0.7,
        xlim=if (include.marginal) range(cd$celldm) else range(dat))
    text(30,130, as.expression(
        if (include.marginal) 
            substitute("medial "*italic(N)==nm, list(nm=n.dat))
        else substitute(italic(N)==nm, list(nm=n.dat))
        ),
        adj=c(0,0), pos=4, offset=0)
    # marginal histogram
    if (include.marginal) {
        axis(1, at=seq(15,55,by=5),labels=FALSE)
        axis(1, at=seq(15,55,by=5),labels=FALSE) # extra one to darken axis
        dat <- celldm$marginal
        n.dat <- length(dat)
        h <- hist(dat, breaks=25, xlim=range(dat), plot=FALSE)
        lines(h$mids, h$counts)
        #my.plot.histogram(h, col="white", add=TRUE, boxwex=0.6,
        #    main="", xlab="", ylab="")
        text(38,40, as.expression(
            substitute("marginal "*italic(N)==nm, list(nm=n.dat))), 
            adj=c(0,0), pos=4, offset=0)
    }
    mtext(xlabel, side=1, line=1.5, cex=0.7)
    #mtext("Frequency", las=0, side=2, line=2.5)
    mtext("No. occurrences", las=0, side=2, line=2.2, cex=0.7)
    if (do.postscript) dev.off()
    cat(fig.file, " width =", fig.width, " height =", fig.height, "\n")
}

pith.figure3 <- function(data, do.postscript=FALSE) 
{
    fig.file <- "Delonixtwigs_plot_sorted.eps" 
    fig.width <- 7.5; fig.height <- 2.5
    if (do.postscript) postscript(file=fig.file, width=fig.width, 
        height=fig.height, onefile=F, horizontal=F, paper="special")
    cd <- subset(data, file == "medial")
    n.tree <- length(unique(cd$tree))
    sp <- split(cd$celldm, cd$tree)
    sp.mean <- tapply(cd$celldm, cd$tree, mean)
    sp.order <- order(sp.mean)
    sp <- sp[sp.order]
    #sp.mean <- lapply(sp, mean)
    sp.forplot <- lapply(sp, function(x) x/10^3)
    par(mar=c(3.5,2,0,0), mgp=c(2.0,0.6,0), las=2, tcl=-0.3, cex=0.7)
    ylabel=expression(paste("Pith cells per meter (",'' %*% 10^3,")"))
    xlabel <- "Tree and branches within tree, in order of increasing mean pith cells per meter"
    
    box.at <- (1:n.tree)# - 0.3
    boxplot(sp.forplot, ylim=range(sp.forplot), varwidth=FALSE, at=box.at, 
        boxwex=0.025, frame.plot=FALSE,
        names=rep("",n.tree), yaxt="n", xlab="", ylab="", axes=FALSE,
        las=2, border="white")
    sp.quantiles <- quantile(unlist(sp.forplot))
    lines(c(0.3,n.tree+0.7), rep(sp.quantiles["25%"], 2), col="black", lty="16")
    lines(c(0.3,n.tree+0.7), rep(sp.quantiles["50%"], 2), col="black", lty="13")
    lines(c(0.3,n.tree+0.7), rep(sp.quantiles["75%"], 2), col="black", lty="16")
    for (br in seq(along=all.branches)) {
        b.cd <- subset(cd, branch == all.branches[br])
        b.sp <- split(b.cd$celldm, b.cd$tree)
        b.sp <- b.sp[sp.order]
        b.sp.forplot <- lapply(b.sp, function(x) x/10^3)
        box.at <- (1:n.tree) - 0.27 + (br - 1)*0.18
        #box.at <- (1:n.tree) - 0.3 + (br - 1)*0.2
        #box.at <- (1:n.tree) - 0.225 + (br - 1)*0.15
        boxplot(b.sp.forplot, varwidth=FALSE, at=box.at, col="black", 
            boxwex=0.05, range=0.00001,
            outline=TRUE, xaxt="n", yaxt="n", xlab="", ylab="", add=TRUE,
            #outpch=18, outcex=0.5, frame.plot=FALSE, axes=FALSE) # diamonds
            #outpch=20, outcex=0.4, frame.plot=FALSE, axes=FALSE) # dots
            outpch=175, outcex=1.5, frame.plot=FALSE, axes=FALSE) # dashes
        med.b.sp.forplot <- lapply(b.sp.forplot, median)
        points(box.at, med.b.sp.forplot, pch=21, bg="white", cex=1.2)
        mtext(rep(as.character(br),n.tree), side=1, line=0, at=box.at, 
            las=1, padj=0, cex=0.7)
    }
    mtext(LETTERS[1:n.tree], side=1, line=1, at=1:n.tree, las=1, padj=0, cex=0.7)
    mtext(xlabel, side=1, line=2.5, las=1, cex=0.7)
    #axis(1, at=1:n.tree, cex=0.7, label=FALSE) #label=rep("",n.tree))
    axis(2, line=-1, cex=0.7)
    axis(4, line=-1, cex=0.7, label=FALSE)
    mtext(ylabel, side=2, line=1, las=0, cex=0.7)
    if (do.postscript) dev.off()
    cat(fig.file, "width =", fig.width, "height =", fig.height, "\n")
}

pith.ranef.plots <- function(analysis)
{
    ans <- c("min.analysis")#c("all.analysis","min.analysis")
    fls <- c("medial")#,"marginal")
    lvl <- c("tree","branch")#,"segment")
    for (an in ans) {
        for (fl in fls) {
            fn <- paste(sep="",an,"_",fl,"_ranef_",lvl[1],".png")
            png(filename=fn)
            print(qqnorm(analysis[[an]][[fl]]$pith.lme.ij,~ranef(.,level=1)))
            dev.off()
            fn <- paste(sep="",an,"_",fl,"_ranef_",lvl[2],".png")
            png(filename=fn)
            print(qqnorm(analysis[[an]][[fl]]$pith.lme.ij,~ranef(.,level=2)))
            dev.off()
        }
    }
}


pith.valcorr <- function(data) 
{
    cd <- subset(data, file == "medial")
    n.tree <- length(unique(cd$tree))
    sp <- split(cd$celldm, cd$tree)
    sp.mean <- tapply(cd$celldm, cd$tree, mean)
    sp.order <- order(sp.mean)
    sp <- sp[sp.order]
    #sp.mean <- lapply(sp, mean)
    sp.quantiles <- quantile(unlist(sp))
    df <- c()
    for (br in seq(along=all.branches)) {
        b.cd <- subset(cd, branch == all.branches[br])
        b.sp <- split(b.cd$celldm, b.cd$tree)
        b.sp <- b.sp[sp.order]
        med.b.sp <- unlist(lapply(b.sp, median))
        mean.b.sp <- unlist(lapply(b.sp, mean))
        max.b.sp <- unlist(lapply(b.sp, max))
        max2.b.sp <- unlist(lapply(b.sp, function(x){x <- x[x != max(x)];max(x)}))
        min.b.sp <- unlist(lapply(b.sp, min))
        min2.b.sp <- unlist(lapply(b.sp, function(x){x <- x[x != min(x)];min(x)}))
        ddf <- data.frame(median=med.b.sp, mean=mean.b.sp,
            min=min.b.sp, min2=min2.b.sp, max2=max2.b.sp, max=max.b.sp)
        df <- rbind(df, ddf)
    }
    df
}

###########################################################
#                                                         #
# statistical analyses using lme etc.                     #
#                                                         #
###########################################################

pith.complete <- function(file="pith.complete.txt", append=FALSE, split=FALSE)
{
    cat("pith.complete: reading and analyzing data...\n")
    analysis <- pith.analysis.complete()
    analysis <<- analysis
    cat("pith.complete: summarizing and comparing fits, dumping to",
        file, "...\n")
    sink(file, append=append, split=split)
    pith.compare.complete(analysis)
    sink(NULL)
    cat("pith.complete: done!\n")
}

pith.compare.complete <- function(analysis)
{
    pith.compare(analysis$all.analysis$medial)
    pith.compare(analysis$min.analysis$medial)
    pith.compare(analysis$all.analysis$marginal)
    pith.compare(analysis$min.analysis$marginal)
}

pith.compare <- function(lst)
{
    dataname <- deparse(substitute(lst))
    cat("=============================================\n")
    cat("=========== pith.compare:", dataname, "\n\n")
    attach(lst, warn.conflicts=FALSE)
    cat("\n=======",dataname,"summary\n\n"); print(data.summary)
    cat("\n=======",dataname,"pith.finalprocessed\n\n");cat(data.pith.finalprocessed)
    cat("\n=======",dataname,"pith.summary\n\n");cat(data.pith.summary)
    prnt <- function(fit, intrv=TRUE) {
        fitname <- deparse(substitute(fit))
        cat("\n=======================\n")
        cat("=======",dataname,fitname,"summary\n\n")
        print(summary(fit))
        cat("\n=======",dataname,fitname,"anova\n\n")
        print(anova(fit))
        if (class(fit) == "lm") {
            cat("=======",dataname,fitname,"lm fit stats\n\n")
            cat("AIC =", AIC(fit), "\n")
            cat("BIC =", BIC(fit), "\n")
            cat("logLik =", logLik(fit), "\n")
            cat("95% conf. int =", confint(fit), "\n")
        } else {
            if (intrv) {
                cat("\n=======",dataname,fitname,"intervals\n\n")
                print(intervals(fit))
            }
        }
    }
    try(prnt(pith.lm))
    try(prnt(pith.lme.i))
    try(prnt(pith.lme.ij))
    try(prnt(pith.lme.ij.wi))
    try(prnt(pith.lme.ij.wii))
    try(prnt(pith.lme.ij.wp))
    try(prnt(pith.lme.ijk))
    try(prnt(pith.lme.ijk.wi))
    try(prnt(pith.lme.ijk.wii))
    try(prnt(pith.lme.ijk.wp))
    cat("\n=======",dataname,"anovas\n\n");
    print(anova.lme(verbose=TRUE, pith.lm, pith.lme.i,
        pith.lme.ij, pith.lme.ij.wi,pith.lme.ij.wii,pith.lme.ij.wp,
        pith.lme.ijk, pith.lme.ijk.wi,pith.lme.ijk.wii,pith.lme.ijk.wp))
    print(anova(pith.lme.ij, pith.lme.ij.wi))
    print(anova(pith.lme.ij, pith.lme.ij.wii))
    print(anova(pith.lme.ij, pith.lme.ij.wp))
    print(anova(pith.lme.ijk, pith.lme.ijk.wi))
    print(anova(pith.lme.ijk, pith.lme.ijk.wii))
    print(anova(pith.lme.ijk, pith.lme.ijk.wp))
    print(anova(pith.lme.ijk.wi, pith.lme.ijk.wii))
    print(anova(pith.lme.ijk.wi, pith.lme.ijk.wp))
    print(anova(pith.lme.ijk.wii, pith.lme.ijk.wp))
    if ("dbh.tests" %in% names(lst) && length(dbh.tests) > 0) {
        cat("\n=======",dataname,"dbh covariate tests\n\n");
        attach(dbh.tests, warn.conflicts=FALSE)
        cat("\n=======",dataname,"dbh summary\n\n"); print(data.summary)
        cat("\n=======",dataname,"dbh pith.finalprocessed\n\n");
        cat(data.pith.finalprocessed)
        cat("\n=======",dataname,"dbh pith.summary\n\n");cat(data.pith.summary)
        try(prnt(pith.lm.dbh))
        try(prnt(pith.lme.ml.i))
        try(prnt(pith.lme.ml.i.dbh))
        try(prnt(pith.lme.ml.ij))
        try(prnt(pith.lme.ml.ij.dbh))
        try(prnt(pith.lme.ml.ijk))
        try(prnt(pith.lme.ml.ijk.dbh))
        cat("\n=======",dataname,"dbh anovas\n\n");
        print(anova.lme(verbose=TRUE, pith.lm.dbh, 
            pith.lme.ml.i, pith.lme.ml.i.dbh,
            pith.lme.ml.ij, pith.lme.ml.ij.dbh,
            pith.lme.ml.ijk, pith.lme.ml.ijk.dbh))
        detach(dbh.tests)
    } else
        cat("\n=======",dataname,"dbh covariate tests DO NOT EXIST\n\n");
    cat("\n=========== pith.compare: DONE", dataname, "\n")
    cat("=============================================\n\n")
    detach()
}

anova.lm.lme <- function(fit.lm, fit.lme) 
{
    if (class(fit.lm) != "lm" && class(fit.lme) != "lme") stop("wrong classes")
    fits <- list(fit.lm, fit.lme)
    aux <- lapply(fits, logLik, TRUE)
    if (length(unique(unlist(lapply(aux, function(el) attr(el, "nall"))))) > 1) {
        stop("all fitted objects must use the same number of observations")
    }
    dfModel <- unlist(lapply(aux, function(el) attr(el, "df")))
    logLik <- unlist(lapply(aux, function(el) c(el)))
    AIC <- unlist(lapply(aux, AIC))
    BIC <- unlist(lapply(aux, BIC))
    rt <- 2
    aod <- data.frame(fit=c("lm","lme"),Model=1:rt,df=dfModel,AIC=AIC,
        BIC=BIC,logLik=logLik,check.names=FALSE)
    ddf <- diff(dfModel)
    if (sum(abs(ddf)) > 0) {
        effects <- rep("", rt)
        for (i in 2:rt) {
            if (ddf[i-1] != 0) {
                effects[i] <- paste(i-1, i, sep=" vs ")
            }
        }
        pval <- rep(NA, rt-1)
        ldf <- as.logical(ddf)
        lratio <- 2*abs(diff(logLik))
        lratio[!ldf] <- NA
        pval[ldf] <- 1 - pchisq(lratio[ldf], abs(ddf[ldf]))
        aod <- data.frame(aod, Test=effects, L.ratio=c(NA,lratio), 
            "p-value"=c(NA,pval),check.names=FALSE)
    }
    aod
}

pith.analysis.complete <- function()
{
    data <- pith.readdata()
    medial <- pith.analysis(subset(data, file == "medial"))
    marginal <- pith.analysis(subset(data, file == "marginal"))
    all.analysis <- list(data=data, medial=medial, marginal=marginal)
    data <- pith.readdata.min()
    medial <- pith.analysis(subset(data, file == "medial"))
    marginal <- pith.analysis(subset(data, file == "marginal"))
    min.analysis <- list(data=data, medial=medial, marginal=marginal)
    list(all.analysis=all.analysis, min.analysis=min.analysis)
}

pith.analysis <- function(data)
{
    pith.lm <- lm(celldm ~ 1, data=data)
    pith.lme.i <- lme(celldm ~1, random = ~1 | tree, data=data,
        control=lmeControl(maxIter=1000, msMaxIter=1000, niterEM=500))
    pith.lme.i.wi <- update(pith.lme.i, weights=varIdent(form= ~1 | tree))
    pith.lme.i.wii <- update(pith.lme.i, weights=varIdent(form= ~1|tree*branch))
    pith.lme.i.wp <- update(pith.lme.i, weights=varPower())
    pith.lme.ij <- update(pith.lme.i, random = ~1 | tree/branch)
    pith.lme.ij.wi <- update(pith.lme.ij, weights=varIdent(form= ~1 | tree))
    pith.lme.ij.wii <- update(pith.lme.ij, weights=varIdent(form= ~1|tree*branch))
    pith.lme.ij.wp <- update(pith.lme.ij, weights=varPower())
    pith.lme.ijk <- update(pith.lme.ij, random = ~1 | tree/branch/segment)
    pith.lme.ijk.wi <- update(pith.lme.ijk, weights=varIdent(form= ~1 | tree))
    pith.lme.ijk.wii <- update(pith.lme.ijk,weights=varIdent(form= ~1|tree*branch))
    pith.lme.ijk.wp <- update(pith.lme.ijk, weights=varPower())
    if ("dbh" %in% names(data)) {
        sdata <- subset(data, !is.na(dbh))
        pith.lm.dbh <- lm(celldm ~ dbh, data=sdata)
        pith.lme.ml.i <- update(pith.lme.i, method="ML", data=sdata)
        pith.lme.ml.i.dbh <- update(pith.lme.ml.i, fixed=celldm ~ dbh)
        pith.lme.ml.ij <- update(pith.lme.ij, method="ML", data=sdata)
        pith.lme.ml.ij.dbh <- update(pith.lme.ml.ij, fixed=celldm ~ dbh)
        pith.lme.ml.ijk <- update(pith.lme.ijk, method="ML", data=sdata)
        pith.lme.ml.ijk.dbh <- update(pith.lme.ml.ijk, fixed=celldm ~ dbh)
        dbh.tests <- list(data=sdata,
            data.summary=summary(sdata),
            data.pith.finalprocessed=pith.finalprocessed(sdata),
            data.pith.summary=pith.summary(sdata),
            pith.lm.dbh=pith.lm.dbh,
            pith.lme.ml.i=pith.lme.ml.i,
            pith.lme.ml.i.dbh=pith.lme.ml.i.dbh,
            pith.lme.ml.ij=pith.lme.ml.ij,
            pith.lme.ml.ij.dbh=pith.lme.ml.ij.dbh,
            pith.lme.ml.ijk=pith.lme.ml.ijk,
            pith.lme.ml.ijk.dbh=pith.lme.ml.ijk.dbh)
    } else 
        dbh.tests <- list()
    list(data=data,
        data.summary=summary(data),
        data.pith.finalprocessed=pith.finalprocessed(data),
        data.pith.summary=pith.summary(data),
        pith.lm=pith.lm,
        pith.lme.i=pith.lme.i,
        pith.lme.i.wi=pith.lme.i.wi,
        pith.lme.i.wii=pith.lme.i.wii,
        pith.lme.i.wp=pith.lme.i.wp,
        pith.lme.ij=pith.lme.ij,
        pith.lme.ij.wi=pith.lme.ij.wi,
        pith.lme.ij.wii=pith.lme.ij.wii,
        pith.lme.ij.wp=pith.lme.ij.wp,
        pith.lme.ijk=pith.lme.ijk,
        pith.lme.ijk.wi=pith.lme.ijk.wi,
        pith.lme.ijk.wii=pith.lme.ijk.wii,
        pith.lme.ijk.wp=pith.lme.ijk.wp,
        dbh.tests=dbh.tests)
}

###########################################################
#                                                         #
# read and organize data files                            #
#                                                         #
###########################################################

pith.readdata.min <- function(include.dbh=FALSE, medial.min=10, marginal.min=3) {
    alldata <- pith.readdata(include.dbh=include.dbh)
    # include medial entries if a tree has at least 10/branch
    sdata <- subset(alldata, file == "medial")
    attach(sdata, warn.conflicts=FALSE)  # FALSE b/c "file" is a base function
    tab <- table(tree, branch)
    trees.ok <- rownames(tab[tab[,"A"]>=10 &
        tab[,"B"]>=10 &
        tab[,"C"]>=10 &
        tab[,"D"]>=10,])
    detach()
    data <- subset(sdata, tree %in% trees.ok)

    # include marginal entries if a tree has at least 3/branch
    sdata <- subset(alldata, file == "marginal")
    attach(sdata, warn.conflicts=FALSE)  # FALSE b/c "file" is a base function
    tab <- table(tree, branch)
    trees.ok <- rownames(tab[tab[,"A"]>=3 &
        tab[,"B"]>=3 &
        tab[,"C"]>=3 &
        tab[,"D"]>=3,])
    detach()

    # combine the datasets
    data <- rbind(data, subset(sdata, tree %in% trees.ok))
    data$tree <- factor(data$tree)
    data
}

pith.sortbranches <- function(data) {
    data$tree <- factor(data$tree)
    data$branch <- factor(data$branch)
    ans <- data[0,]
    for (tt in levels(data$tree)) {
        td <- subset(data, tree == tt)
        by.branch <- split(td$celldm, td$branch)
        mean.by.branch <- unlist(lapply(by.branch, mean))
        sorted.branch <- names(mean.by.branch)[order(mean.by.branch)]
        tmp.branches <- all.branches
        names(tmp.branches) <- paste(sep="", all.branches, "x")
        td$branch <- as.character(td$branch)
        for (b in seq(along=sorted.branch))
            td$branch[td$branch == sorted.branch[b]] <- names(tmp.branches)[b]
        for (b in seq(along=tmp.branches))
            td$branch[td$branch == names(tmp.branches)[b]] <- tmp.branches[b]
        td <- td[order(td$branch),]
        ans <- rbind(ans, td)
    }
}


pith.readdata <- function(include.dbh=FALSE) {
    cd <- read.table("data_cellsizes.txt", as.is=T, sep="\t", header=T)
    split.segmentname <- function(this) {
        ans <- rep("", 3)
        l.this <- nchar(this)
        ans[1] <- substr(this, 1, l.this-3)  # tree
        ans[2] <- substr(this, l.this-2, l.this-2)  # branch
        ans[3] <- substr(this, l.this, l.this)  # segment
        ans
    }
    cd.segment <- cd.file <- cd.len <- cd.celldm <- c()
    new.cd <- c()
    measure.factor <- c(rep("medial", 10), rep("marginal", 4))
    for (i in 1:length(cd$slide)) {
        segment <- split.segmentname(cd$slide[i])
        if (substr(segment[1],1,1) == "C") next  # skip trials which begin with C
        if (segment[2] == ".") print(segment[1])
        len <- as.numeric(cd[i,-1])
        for (j in 1:14) {
            cd.segment <- rbind(cd.segment, t(segment))
            cd.file <- rbind(cd.file, measure.factor[j])
            cd.len <- rbind(cd.len, len[j])
            cd.celldm <- rbind(cd.celldm, 0)
        }
    }
    colnames(cd.segment) <- c("tree", "branch", "segment")
    new.cd <- data.frame(cd.segment, file=cd.file, len=cd.len, celldm=cd.celldm)
    # lengths are in the data file in micrometer units, each of which
    # is 6.13 um.  lengths are the length of 10 contiguous cells in a file.
    new.cd <- subset(new.cd, !is.na(new.cd$len))
    # refactor now that we've removed NA
    new.cd$tree <- factor(new.cd$tree)
    new.cd$branch <- factor(new.cd$branch)
    new.cd$segment <- factor(new.cd$segment)
    new.cd$len <- new.cd$len * 0.1 * 6.13 # now is mean length in um of 1 cell
    new.cd$celldm <- 1 / (new.cd$len*(10^-6))  # cell div/m is recip of cell length
    #cd <<- new.cd
    if (include.dbh) {
        # read in dbh data, but center it around 0 so that the intercept estimate
        # is actually phi
        tdbh <- subset(read.table("data_trees.txt", as.is=T, sep="\t", header=T),
            tree %in% as.character(unique(new.cd$tree)))
        is.na(tdbh$dbh) <- is.na(tdbh$circ) | tdbh$dbh == 0
        rownames(tdbh) <- tdbh$tree
        tdbh$dbh.zero <- tdbh$dbh - mean(tdbh$dbh, na.rm=TRUE)
        #new.cd <- cbind(new.cd, dbh=tdbh[new.cd$tree,"dbh"])
        new.cd <- cbind(new.cd, dbh=tdbh[new.cd$tree,"dbh"],
            dbh.zero=tdbh[new.cd$tree,"dbh.zero"])
    }
    return(new.cd)
}

pith.grouped <- function(dataf=cd, sub.file=FALSE)
{
    if (sub.file == TRUE || ! is.logical(sub.file)) {
        if (! is.logical(sub.file)) dataf <- subset(dataf, file == sub.file)
        else dataf <- subset(dataf, file == "medial")
        dataf$file <- factor(dataf$file)
        gcd <- groupedData(celldm ~ 1 | tree/branch/segment, data=dataf, 
            labels=list(x="Count"),units=list(x="(cells/m)"))
    } else {
        gcd <- groupedData(celldm ~ file | tree/branch/segment, 
            data=dataf, labels=list(x="Count"),units=list(x="(cells/m)"))
    }
    return(gcd)
}


###########################################################
#                                                         #
# summarize data frames                                   #
#                                                         #
###########################################################

pith.finalprocessed <- function(data) {
    tt <- levels(factor(data$tree))
    fp.t <- length(tt)
    fp.b <- fp.s <- fp.med <- fp.mar <- c()
    for (t in tt) {
        t.sub.cd <- subset(data, tree == t)
        bb <- levels(factor(t.sub.cd$branch))
        fp.b <- cbind(fp.b, length(bb))
        for (b in bb) {
            b.sub.cd <- subset(t.sub.cd, branch == b)
            ss <- levels(factor(b.sub.cd$segment))
            fp.s <- cbind(fp.s, length(ss))
            for (s in ss) {
                s.sub.cd <- subset(b.sub.cd, segment == s)
                med.sub.cd <- subset(s.sub.cd, file == "medial")
                mar.sub.cd <- subset(s.sub.cd, file == "marginal")
                fp.med <- cbind(fp.med, length(med.sub.cd$celldm))
                fp.mar <- cbind(fp.mar, length(mar.sub.cd$celldm))
            }
        }
    }
    fp.b <- as.vector(fp.b)
    fp.s <- as.vector(fp.s)
    fp.med <- as.vector(fp.med)
    fp.mar <- as.vector(fp.mar)
    paste(sep="",paste("number of trees:", fp.t,"\n"),
        paste("number of branches:", sum(fp.b),"\n"),
        paste("mean branches per tree:", mean(fp.b),"\n"),
        paste("S.D. branches per tree:", sqrt(var(fp.b)),"\n"),
        paste("S.E. branches per tree:", standard.error(fp.b),"\n"),
        paste("number of segments:", sum(fp.s),"\n"),
        paste("mean segments per branch:", mean(fp.s),"\n"),
        paste("S.D. segments per branch:", sqrt(var(fp.s)),"\n"),
        paste("S.E. segments per branch:", standard.error(fp.s),"\n"),
        paste("number of medial files:", sum(fp.med),"\n"),
        paste("mean medial files per segment:", mean(fp.med),"\n"),
        paste("S.D. medial files per segment:", sqrt(var(fp.med)),"\n"),
        paste("S.E. medial files per segment:", standard.error(fp.med),"\n"),
        paste("number of marginal files:", sum(fp.mar),"\n"),
        paste("mean marginal files per segment:", mean(fp.mar),"\n"),
        paste("S.D. marginal files per segment:", sqrt(var(fp.mar)),"\n"),
        paste("S.E. marginal files per segment:", standard.error(fp.mar),"\n"))
}

pith.summary <- function(data) {
    measure.summary <- function(measure, file, 
            measure.name=deparse(substitute(measure))) {
        file.summary <- function(x, file.name=deparse(substitute(x)),measure.name) {
            if (length(x) == 0) 
                paste(sep="",
                    paste("data source:", measure.name,"\n"),
                    paste("file:", file.name,"\n"),
                    paste("data length == 0","\n"))
            else {
                qnt <- quantile(x)
                paste(sep="",
                    paste("data source:", measure.name,"\n"),
                    paste("file:", file.name,"\n"),
                    paste("mean:", mean(x),"\n"),
                    paste("S.D.:", sqrt(var(x)),"\n"),
                    paste("S.E.:", standard.error(x),"\n"),
                    paste("min:", qnt[1],"\n"),
                    paste("25%:", qnt[2],"\n"),
                    paste("median:", qnt[3],"\n"),
                    paste("75%:", qnt[4],"\n"),
                    paste("max:", qnt[5],"\n"))
            }
        }
        sp <- split(measure, file)
        paste(sep="\n", file.summary(sp$medial,measure.name=measure.name), 
            file.summary(sp$marginal,measure.name=measure.name), "")
    }
    paste(sep="\n", measure.summary(data$celldm, data$file), 
        measure.summary(data$len, data$file), "")
}

#pith.plotfour <- function(postscript=FALSE, dataf=cd, 
#    this.tree="T17")
#{
#    #opar <- par()
#    fig.file <- "Delonixtwigs_plot_branches.eps" 
#    fig.width <- 4; fig.height <- 2.5
#    if (postscript) {
#        postscript(file=fig.file, width=fig.width, 
#            height=fig.height, onefile=F, horizontal=F, paper="special")
#    }
#    cd <- subset(dataf, file == "medial")
#    this <- subset(cd, tree == this.tree)
#    sp <- split(this$celldm, this$branch)
#    sp.mean <- tapply(this$celldm, this$branch, mean)
#    sp <- sp[order(sp.mean)]
#    sp.mean <- lapply(sp, mean)
#    sp.forplot <- lapply(sp, function(x) x/10^3)
#    par(mar=c(3.1,5,0,0), las=1, cex=0.8)
#    ylabel.first <- "Medial pith cell divisions"
#    ylabel.second <- expression(paste("per meter (",'' %*% 10^3,")"))
#    ylabel=expression(paste("Medial pith cell divisions per meter (",'' %*% 10^3,")"))
#    ylabel=expression(paste("per meter (",'' %*% 10^3,")"))
#    xlabel <- paste("Branches of tree", this.tree, ",\n in order of increasing mean")
#    boxplot(sp.forplot, boxwex=0.6, varwidth=T, xaxt="n", xlab="", ylab="")
##    boxplot(sp.forplot, boxwex=0.6, varwidth=T, range=0, xaxt="n", xlab="", ylab="")
#    mtext(paste("Branches of tree", this.tree), side=1, line=0.8, las=1, cex=1)
#    mtext("in order of increasing mean", side=1, line=2, las=1, cex=1)
#    mtext(ylabel.first, side=2, line=4, las=0, cex=1)
#    mtext(ylabel.second, side=2, line=2.6, las=0, cex=1)
#    if (postscript) {
#        dev.off()
#    }
#    cat(fig.file, "width =", fig.width, "height =", fig.height, "\n")
#    #par(opar)
#}
#
#pith.plotthree <- function(postscript=FALSE, dataf=cd, 
#    this.tree="T16", this.branch="D") 
#{
#    #opar <- par()
#    fig.file <- "Delonixtwigs_plot_histogram_tree.eps" 
#    fig.width <- 5; fig.height <- 2.5
#    if (postscript) {
#        postscript(file=fig.file, width=fig.width, 
#            height=fig.height, onefile=F, horizontal=F, paper="special")
#    }
#    cd <- subset(subset(dataf, tree == this.tree), branch == this.branch)
#    celldm <- split(cd$celldm, cd$file)
#    celldm$medial <- celldm$medial / 10^3
#    n.medial <- length(celldm$medial)
#    par(mar=c(4,4,0,0), las=1)
#    xlabel.first="T16, mean pith cell divisions"
#    xlabel.second=expression(paste("per meter (",'' %*% 10^3,")"))
#    hist(celldm$medial, xlab="", ylab="Frequency", main="")
##    hist(celldm$medial, breaks=15, xlim=range(celldm$medial), xlab="", ylab="Frequency", main="")
#    mtext(xlabel.first, side=1, line=3, las=1, cex=1)
#    mtext(xlabel.second, side=1, line=1.6, las=1, cex=1)
#    text(38,130, as.expression(substitute(n==nm, list(nm=n.medial))), 
#        adj=c(0,0), pos=4, offset=0)
#    if (postscript) {
#        dev.off()
#    }
#    cat(fig.file, "width =", fig.width, "height =", fig.height, "\n")
#    #par(opar)
#}
#
#pith.plottwo <- function(postscript=FALSE, dataf=cd) 
#{
#    #opar <- par()
#    fig.file <- "Delonixtwigs_plot_sorted.eps" 
#    fig.width <- 5; fig.height <- 2.5
#    if (postscript) {
#        postscript(file=fig.file, width=fig.width, 
#            height=fig.height, onefile=F, horizontal=F, paper="special")
#    }
#    cd <- subset(dataf, file == "medial")
#    sp <- split(cd$celldm, cd$tree)
#    # sp <- split(dataf$celldm, dataf$tree)
#    sp.mean <- tapply(cd$celldm, cd$tree, mean)
#    sp <- sp[order(sp.mean)]
#    sp.mean <- lapply(sp, mean)
#    sp.forplot <- lapply(sp, function(x) x/10^3)
#    sp.order <- order()
#    par(mar=c(4.5,5,0,0), las=2, cex=0.8)
#    ylabel.first <- "Medial pith cell divisions"
#    ylabel.second <- expression(paste("per meter (",'' %*% 10^3,")"))
#    ylabel=expression(paste("Medial pith cell divisions per meter (",'' %*% 10^3,")"))
#    ylabel=expression(paste("per meter (",'' %*% 10^3,")"))
#    xlabel <- "Tree, in order of increasing mean"
#    boxplot(sp.forplot, varwidth=T, yaxt="n", xlab="", ylab="")
##    boxplot(sp.forplot, range=0, varwidth=T, yaxt="n", xlab="", ylab="")
#    mtext(xlabel, side=1, line=3.5, las=1, cex=1)
#    axis(2, cex=1)
#    mtext(ylabel.first, side=2, line=4, las=0, cex=1)
#    mtext(ylabel.second, side=2, line=2.6, las=0, cex=1)
#    #stripchart(sp.forplot, pch="o", xlab="", ylab="", vertical=T)
#    if (postscript) {
#        dev.off()
#    }
#    cat(fig.file, "width =", fig.width, "height =", fig.height, "\n")
#    #par(opar)
#}
#
#pith.plotone <- function(postscript=FALSE, dataf=cd) 
#{
#    #opar <- par()
#    fig.file <- "Delonixtwigs_plot_histogram.eps" 
#    fig.width <- 5; fig.height <- 2.5
#    if (postscript) {
#        postscript(file=fig.file, width=fig.width, 
#            height=fig.height, onefile=F, horizontal=F, paper="special")
#    }
#    xlabel=expression(paste("Mean medial pith cell divisions per meter (",'' %*% 10^3,")"))
#    cd <- subset(dataf, !is.na(dataf$celldm))
#    celldm <- split(cd$celldm, cd$file)
#    celldm$medial <- celldm$medial / 10^3
#    n.medial <- length(celldm$medial)
#    par(mar=c(4,4,0,0), las=1)
#    hist(celldm$medial, breaks=25, xlim=range(celldm$medial), xlab=xlabel, ylab="Frequency", main="", col="gray")
#    text(38,130, as.expression(substitute(italic(N)==nm, list(nm=n.medial))), adj=c(0,0), pos=4, offset=0)
#    if (postscript) {
#        dev.off()
#    }
#    cat(fig.file, "width =", fig.width, "height =", fig.height, "\n")
#    #par(opar)
#}
#

###########################################################
#                                                         #
# utility functions                                       #
#                                                         #
###########################################################

pith.var <- function() {
    # to determine  the relationoship between the mean and variance of cell
    # sizes when cells are examined singly, or in groups of 10
    ntwigs <- 100
    ncells <- 10
    cell.mean <- 6
    cell.sd <- 1
    this <- rnorm(ntwigs*ncells, mean=cell.mean, sd=cell.sd)
    single.mean <- mean(this)
    single.sd <- sqrt(var(this))
    twigfac <- as.factor(rep(1:ntwigs,each=ncells))
    tens <- tapply(this, twigfac, sum)
    tens.mean <- mean(tens) / 10
    tens.sd <- sqrt(var(tens))
    tens.tenth <- tens / 10
    tens.tenth.mean <- mean(tens.tenth)
    tens.tenth.sd <- sqrt(var(tens.tenth))
    cat(paste("single mean=",single.mean," var=",single.sd^2," sd=",single.sd))
    cat(paste("tens mean=",tens.mean," var=",tens.sd^2,"  sd=",tens.sd))
    cat(paste("tens.tenth mean=",tens.tenth.mean," var=",tens.tenth.sd^2,"  sd=",tens.tenth.sd))
    cat("establishes that mean(10 cells)/10 = mean(1 cell)")
    cat("establishes that var(10 cells)/10 = var(1 cell)")
}

my.plot.histogram <- function (x, freq = equidist, density = NULL, angle = 45,
    col = NULL, border = par("fg"), lty = NULL, main = paste("Histogram of",
    x$xname), sub = NULL, xlab = x$xname, ylab, xlim = range(x$breaks), 
    ylim = NULL, axes = TRUE, labels = FALSE, add = FALSE, boxwex = 1, ...)
{
    equidist <- if (is.logical(x$equidist))
        x$equidist
    else {
        h <- diff(x$breaks)
        diff(range(h)) < 1e-07 * mean(h)
    }
    if (freq && !equidist)
        warning("the AREAS in the plot are wrong -- rather use `freq=FALSE'!")
    y <- if (freq)
        x$counts
    else {
        y <- x$density
        if (is.null(y))
            x$intensities
        else y
    }
    nB <- length(x$breaks)
    if (is.null(y) || 0 == nB)
        stop("`x' is wrongly structured")
    if (!add) {
        if (is.null(ylim))
            ylim <- range(y, 0)
        if (missing(ylab))
            ylab <- if (!freq)
                "Density"
            else "Frequency"
        plot.new()
        plot.window(xlim, ylim, "")
        title(main = main, sub = sub, xlab = xlab, ylab = ylab,
            ...)
        if (axes) {
            axis(1, ...)
            axis(2, ...)
        }
    }
    my.rightbreaks <- x$breaks[-nB] + diff(x$breaks) * boxwex
    rect(x$breaks[-nB], 0, my.rightbreaks, y, col = col, border = border,
        angle = angle, density = density, lty = lty)
    # rect(x$breaks[-nB], 0, x$breaks[-1], y, col = col, border = border,         angle = angle, density = density, lty = lty)
    if ((logl <- is.logical(labels) && labels) || is.character(labels))
        text(x$mids, y, labels = if (logl) {
            if (freq)
                x$counts
            else round(x$density, 3)
        }
        else labels, adj = c(0.5, -0.5))
    invisible()
}

standard.error <- function(x, na.rm=TRUE, catch.zero=TRUE) {
   len <- function(x, na.rm) if (na.rm) sum(!is.na(x)) else length(x)
   if (catch.zero && len(x, na.rm) < 1) return(0)
   sqrt(var(x, na.rm=na.rm)/len(x, na.rm))
}

