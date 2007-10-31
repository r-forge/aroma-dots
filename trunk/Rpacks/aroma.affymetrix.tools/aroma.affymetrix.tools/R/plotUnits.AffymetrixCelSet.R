setMethodS3("plotUnits","AffymetrixCelSet",function(this,units,dataMat=NULL,arrays=NULL,field="intensities",
    markGroups=T,markcol="grey",labelGroups=F,labelcol="black",plotByArray=T,marklwd=1.5,joinGroups=T,plot.it=T, main=NULL,
    type="l",xaxt="n",ylab=NULL,xlab="",logscale=T, col=1:6,...){
    
    #plot the data with the effects
    #if(plotByArray==F) markGroups<-F #could change this so have different colors/lty for different exons?
    if(is.null(ylab)){
     if(logscale) ylab<-"Log2 Intensities"
     else ylab<-"Intensities"}
    cdf<-getCdf(this)
    if(length(units)>1) stop("units must be single integer -- only can plot 1 unit") ##could change this...
    nCellsPerGroup<-.subset2(readCdfNbrOfCellsPerUnitGroup(getPathname(cdf),unit=units),1)
    if(is.null(arrays)) arrays<-1:length(this)

#    cellIndices<-unlist(getCellIndices(cdf,units=units)) ##destroy grouping so get solid matrix
#    dataMat<-log2(getIntensities(this,indices=cellIndices)[,arrays]) #if getIntensities is fixed, then make this more elegant with arrays=arrays

    #do the same with readUnits to speed it up
    if(is.null(dataMat)){
        dataList<-mergeGroups(this,units=units,fields=field)
        dataMat<-.subset2(.subset2(dataList,1),field)
    }
    if(logscale==T) dataMat<-log2(dataMat)
    dataMat<-dataMat[,arrays]
    if(plot.it){
        if(plotByArray==F){
            xlim<-c(.5,length(arrays)+.5) #maybe change so more flexible???
            if(markGroups) colExon<-rep(rep(col,length=length(nCellsPerGroup)),times=nCellsPerGroup)
            else colExon<-col
            matplot(1:length(arrays),t(dataMat),xlim=xlim,xaxs="i",type=type,xlab=xlab,ylab=ylab,xaxt=xaxt, col=colExon,...)    
        }
        else{
            xlim<-c(.5,sum(nCellsPerGroup)+.5) #maybe change so more flexible???
            if(markGroups==T && joinGroups==F){
                matplot(1:sum(nCellsPerGroup),dataMat,xlim=xlim,xaxs="i",type="n",xlab=xlab,ylab=ylab,xaxt=xaxt,col=col,...)
                xstart<-c(1,head(cumsum(nCellsPerGroup)+1,-1))
                xend<-cumsum(nCellsPerGroup) 
                tmp<-lapply(1:length(xstart),function(i){
                    matplot(xstart[i]:xend[i],dataMat[xstart[i]:xend[i],],type=type,add=T,col=col,...)
                    })
            }
            else{
                matplot(1:sum(nCellsPerGroup),dataMat,xlim=xlim,xaxs="i",type=type,xlab=xlab,ylab=ylab,xaxt=xaxt,col=col,...)    
            }
            abline(v=c(0,cumsum(nCellsPerGroup))+.5,col=markcol,lwd=marklwd)
            axis(1,at=c(0,cumsum(nCellsPerGroup))+.5,labels=F,col=markcol)
            if(markGroups==T && labelGroups==T){
                xstart<-c(1,head(cumsum(nCellsPerGroup)+1,-1))
                xend<-cumsum(nCellsPerGroup) 
                mids<-(xstart+xend)/2
                groups<-readCdfGroupNames(getPathname(cdf),unit=units)[[1]]
                if(length(labelcol)==1) axis(1,at=mids,labels=groups,tick=F,las=2,col.axis=labelcol)
                else{
                    labelcol<-rep(labelcol,length=length(groups))
                    for(i in 1:length(groups)){ axis(1,at=mids[i],labels=groups[i],tick=F,las=2,col.axis=labelcol[i])}
                }
            }
        }
        if(is.null(main)) main<- paste("Expression of  ", getUnitNames(cdf,units=units)," (Data=",getFullName(this),")",sep="")
        title(main=main)
    }
    invisible(dataMat)
})
