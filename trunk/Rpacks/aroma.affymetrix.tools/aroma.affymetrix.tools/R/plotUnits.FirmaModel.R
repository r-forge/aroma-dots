setMethodS3("plotUnits","FirmaModel",function(this,units,arrays=NULL,
    colSchemeProbe=brewer.pal(11,"RdBu"),colSchemeFIRMA=c("white",brewer.pal(9,"YlGnBu")),
    breaksProbe=NULL,breaksFIRMA=NULL,colArrays=NULL,arrayFactor=NULL,arrayNames=NULL,
    mar=c(2.1,7.1,.5,.2),intersp=2.1,plotResiduals=T,...){

    require(RColorBrewer)
    if(length(units)>1) stop("units must be single integer -- only can plot 1 unit") ##could change this...
    plm<-getPlm(this) #input PLM object
    cdf<-getCdf(this)
    ds<-getDataSet(plm)
    fs<-getFirmaScores(this)
    fsMat<-log2(mergeGroups(fs,units=units)[[1]][[1]])
    wtMat<-log2(mergeGroups(getWeightsSet(plm),units=units,fields="wts")[[1]][[1]])
    rsMat<-log2(mergeGroups(getResidualSet(plm),units=units,fields="eps")[[1]][[1]])
    nbrOfArrays <- nbrOfArrays(ds)
    nProbesPerExon<-readCdfNbrOfCellsPerUnitGroup(getPathname(cdf),unit=units)[[1]]
    
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
    
    if(is.null(colArrays)) colArrays<-rep(palette(),length=nlevels(arrayFactor))[arrayFactor]
    else{
        if(length(colArrays)!=nbrOfArrays){
            warning("not same number of colors in 'colArrays' as arrays to be plotted")
            colArrays<-rep(colArrays,length=nbrOfArrays)
        }
        colArrays<-colArrays[arrayOrder] #put in new order
    }
    if(is.null(breaksProbe)) breaksProbe<-seq(-1,1,length=length(colSchemeProbe)+1)
    if(is.null(breaksFIRMA)) breaksFIRMA<-seq(0,1,length=length(colSchemeFIRMA)+1)
    if(is.null(arrayNames)) arrayNames<-getNames(fs)[arrays]
    else arrayNames<-arrayNames[arrays]
    if(plotResiduals) par(mfrow=c(3,1))
    else  par(mfrow=c(2,1))
    #probe level weights
    #c(2.1,7.1,.5,.2)
    if(plotResiduals){
        par(mar=c(intersp,mar[2],mar[3],mar[4]),las=2,mgp=c(3, 0.5, 0))
        mat<-((1-wtMat)*sign(rsMat))[,arrays]
        image(1:sum(nProbesPerExon),y=1:ncol(mat),mat,breaks=breaksProbe,
            col=colSchemeProbe,xaxt="n",yaxt="n",xlab="",ylab="")
        abline(v=c(0,cumsum(nProbesPerExon))+.5,col="grey",lwd=1.5)
        axis(1,at=c(0,cumsum(nProbesPerExon))+.5,labels=F,col="grey")
        mapply(axis,label=arrayNames,col.axis=colArrays,
            at=1:length(arrays),MoreArgs=list(side=2),las=2)
        if(hlines) abline(h=c(0,cumsum(table(arrayFactor)))+.5,col="grey")
        box()
    }
    #exon level score
    par(mar=c(intersp,mar[2],intersp,mar[4]),las=2,mgp=c(3, 0.5, 0))
    mat<-fsMat[,arrays]
    image(c(0,cumsum(nProbesPerExon))+.5,y=1:ncol(mat),mat,breaks=breaksFIRMA,
        col=colSchemeFIRMA,xaxt="n",yaxt="n",xlab="",ylab="")
    abline(v=c(0,cumsum(nProbesPerExon))+.5,col="grey",lwd=1.5)
    axis(1,at=c(0,cumsum(nProbesPerExon))+.5,labels=F,col="grey")
    mapply(axis,label=arrayNames,col.axis=colArrays,
        at=1:length(arrayFactor),MoreArgs=list(side=2),las=2)
    if(hlines) abline(h=c(0,cumsum(table(arrayFactor)))+.5,col="grey")
    box()

    par(mar=c(mar[1],mar[2],intersp,mar[4]),mgp=c(3,1,0))
    plotUnits(plm,unit=units,col=colArrays[order(arrayOrder)],arrays=arrays[order(arrayOrder)],...)
})
