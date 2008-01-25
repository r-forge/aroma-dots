#plotting Exon-specific (a matrix) scores

setMethodS3("plotUnits","scores",function(scoresList,plm,units,arrays=NULL, 
    plotType=list("image"),scoreNames=list(NULL),breaksScores=list(NULL),ylim=list(NULL),ablineAdd=list(NULL), colSchemeScores=list(colSchemeProbe),
    colArrays=NULL,arrayFactor=NULL,arrayNames=NULL,lwdArrays=1, ltyArrays=1, pchArrays=19,cex=1, highLightProbes=list(NULL),
    additionalPlots=c("residuals","weights","none"),additionalType=c("image","lines"),colSchemeProbe=brewer.pal(11,"RdBu"),breaksProbe=NULL,
    mar=c(2.1,7.1,.5,.2),intersp=2.1,logscale=FALSE,equalScale=FALSE,ylimIntensities=NULL,plotHeights=NULL,...){
    
    require(RColorBrewer)
    additionalPlots<-match.arg(additionalPlots)
    additionalType<-match.arg(additionalType)
    if(length(units)>1) stop("units must be single integer -- only can plot 1 unit") ##could change this...
    cdf<-getCdf(plm)
    ds<-getDataSet(plm)
    if(additionalPlots=="residuals") rsMat<-log2(mergeGroups(getResidualSet(plm),units=units,fields="eps")[[1]][[1]])
    if(additionalPlots=="weights") rsMat<-log2(mergeGroups(getWeightsSet(plm),units=units,fields="eps")[[1]][[1]])
    nbrOfArrays <- nbrOfArrays(ds)
    nProbesPerExon<-readCdfNbrOfCellsPerUnitGroup(getPathname(cdf),unit=units)[[1]]
    nExons<-length(nProbesPerExon)
    nscores<-length(scoresList)
    plotType<-rep(plotType,length=nscores)
    plotType<-lapply(plotType,match.arg,choices=c("image","lines"))
    colSchemeScores<-rep(colSchemeScores,length=nscores)
    breaksScores<-rep(breaksScores,length=nscores)
    scoreNames<-rep(scoreNames,length=nscores)
    ylim<-rep(ylim,length=nscores)
    ablineAdd<-rep(ablineAdd,length=nscores)


    if(!is.null(arrays)){
        if(any(!is.numeric(arrays))) throw("Invalid array: must be numeric")
        else if(any(arrays>nbrOfArrays | arrays<1)) throw("Invalid array")
    }
    else arrays<-seq(ds)
    nbrOfArrays <- length(arrays)
    if(is.null(arrayFactor)){
        arrayFactor<-as.factor(arrays)
        hlines<-FALSE}
    else{
        if(length(arrayFactor)!= nbrOfArrays) stop("'arrayFactor' must be of same length as number of arrays plotted")
        hlines<-TRUE}
    arrayOrder<-order(arrayFactor)
    arrays<-arrays[arrayOrder]
    arrayFactor<-arrayFactor[arrayOrder]
    repArrayValue<-function(x){
            if(length(x)!=nbrOfArrays && length(x)!=1) warning("Plotting parameter not same length as number of arrays")
        x<-rep(x,length=nbrOfArrays)
        x<-x[arrayOrder]
        return(x)
    }
    if(is.null(colArrays)) colArrays<-rep(palette(),length=nlevels(arrayFactor))[arrayFactor]
    else{colArrays<-repArrayValue(colArrays)    }
    lwdArrays<-repArrayValue(lwdArrays)
    pchArrays<-repArrayValue(pchArrays)
    ltyArrays<-repArrayValue(ltyArrays)
    
    if(is.null(breaksProbe)) breaksProbe<-seq(min(rsMat),max(rsMat),length=length(colSchemeProbe)+1)
    if(is.null(arrayNames)) arrayNames<-getNames(ds)[arrays]
    else arrayNames<-arrayNames[arrays]

    plotSwitch<-function(mat,type,colScheme,breaks=NULL,ylim,highlight){
        if(is.null(breaks) && type=="image") breaks<-seq(min(mat,na.rm=TRUE),max(mat,na.rm=TRUE),length=length(colScheme)+1)
        if(nrow(mat)==nExons) level<-"exon"
        else{
            if(nrow(mat)==sum(nProbesPerExon))level<-"probe"
            else stop(paste("Incorrect number of rows in scores matrix: Should have",sum(nProbesPerExon),"or",nExons,"rows"))
        }
        exonBoundaries<-c(0, cumsum(nProbesPerExon)) + 0.5
        midpoints<-(tail(exonBoundaries,-1)+head(exonBoundaries,-1))/2
        
        if(type=="image"){
            x<-switch(level,"probe"=1:sum(nProbesPerExon),"exon"=exonBoundaries) #boundaries -- may not be evenly spaced...
            image(x,y=1:ncol(mat),mat,breaks=breaks,
                col=colScheme,yaxt="n",xlab="",ylab="",xaxt="n")
            mapply(axis,label=arrayNames,col.axis=colArrays,
                at=1:length(arrays),MoreArgs=list(side=2),las=2)
            if(hlines) abline(h=c(0,cumsum(table(arrayFactor)))+.5,col="grey")
        }
        if(type=="lines"){
            x<-switch(level,"probe"=1:sum(nProbesPerExon),"exon"=midpoints) #
            matplot(x,mat,col=colArrays,xaxs="i",xlab="",ylab="",type="b",pch=pchArrays,
                xlim=range(exonBoundaries),xaxt="n",ylim=ylim,lwd=lwdArrays,lty=ltyArrays,cex=cex)
        }
        abline(v=exonBoundaries,col="grey",lwd=1.5)
        axis(1,at=exonBoundaries,labels=FALSE,col="grey")
        box()    
    }

    
    if(equalScale){ylim<-rep(list(range(unlist(scoresList),na.rm=TRUE)),length=nscores)}

    nplots<-nscores+1 #the original data
    if(additionalPlots!="none") nplots<-nplots+1 #also the residuals

#    par(mfrow=c(nplots,1))

    if(is.null(plotHeights)) plotHeights<-rep(1,nplots)
    else if(length(plotHeights)!=nplots) stop("If 'heights' is specified, must be length of total number of plots (not length of scoresList)")
    layout(matrix(1:nplots,ncol=1),heights=plotHeights)

    if(additionalPlots!="none"){
        par(mar=c(intersp,mar[2],mar[3],mar[4]),las=2,mgp=c(3, 0.5, 0))
        mat<-rsMat[,arrays]
        ylimRes<-range(as.vector(mat),na.rm=TRUE)
        plotSwitch(mat,type=additionalType,colScheme=colSchemeProbe,breaks=breaksProbe,ylim=ylimRes)
    }

    #score matrix
    lapply(1:nscores,function(i){
        if(additionalPlots=="none" && i==1) par(mar=c(intersp,mar[2],mar[3],mar[4]),las=2,mgp=c(3, 0.5, 0))
        else par(mar=c(intersp,mar[2],intersp,mar[4]),las=2,mgp=c(3, 0.5, 0))
        mat<-scoresList[[i]]
        if(ncol(mat)==nbrOfArrays) mat<-mat[,arrays,drop=FALSE]
        if(is.null(ylim[[i]])) ylim[[i]]<-range(as.vector(mat),na.rm=TRUE)
        plotSwitch(mat,type=plotType[[i]], colScheme=colSchemeScores[[i]],breaks=breaksScores[[i]],ylim=ylim[[i]])
        if(is.null(scoreNames[[i]])) scoreTitle<-names(scoresList)[i]
        else scoreTitle<-scoreNames[[i]]
        title(main=scoreTitle)
        if(!is.null(ablineAdd[[i]])) abline(h=ablineAdd[[i]],col="grey",lty=2)
    })
    par(mar=c(mar[1],mar[2],intersp,mar[4]),mgp=c(3,1,0))
    if(!is.null(ylimIntensities))plotUnits(plm,unit=units,col=colArrays[order(arrayOrder)],arrays=arrays[order(arrayOrder)],lwd=lwdArrays[order(arrayOrder)],lty=ltyArrays[order(arrayOrder)],pch=pchArrays[order(arrayOrder)],ylim=ylimIntensities,...)
    else plotUnits(plm,unit=units,col=colArrays[order(arrayOrder)],arrays=arrays[order(arrayOrder)],lwd=lwdArrays[order(arrayOrder)],lty=ltyArrays[order(arrayOrder)],pch=pchArrays[order(arrayOrder)],...)
    invisible(list(breaksProbe=breaksProbe,breaksScores=breaksScores))
})
