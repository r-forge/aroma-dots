setMethodS3("plotUnits","ProbeAffinityFile",function(this,units,intercept=0,x.pos=NULL,add=TRUE,type=ifelse(add,"l","h"),
    joinGroups=TRUE,plot.it=TRUE,xlab="Probes",ylab="Probe Affinity",...){
    #add lines for the effects -- note can't add on top of data, because the chip effects/intercept wrong...
    #takes parameter intercept to fix this if know the offset
    if(length(units)>1) stop("units must be single integer -- only can plot 1 unit") ##could change this...
    cdf<-getCdf(this)
    cellIndices<-unlist(getCellIndices(cdf,units=units)) ##destroy grouping so get solid matrix
    probeEffect<-log2(getData(this,cellIndices,field="intensities"))[,1] #could add stdev plot too...
    if(is.null(x.pos)) x.pos<-1:length(probeEffect)
    nCellsPerGroup<-.subset2(readCdfNbrOfCellsPerUnitGroup(getPathname(cdf),unit=units),1)
    if(plot.it){
        if(!joinGroups){
            if(!add){
                 plot(x.pos,probeEffect+intercept,xaxs="i",type="n",xlab=xlab,ylab=ylab,...)
            }
            xstart<-c(1,head(cumsum(nCellsPerGroup)+1,-1))
            xend<-cumsum(nCellsPerGroup) 
            tmp<-lapply(1:length(xstart),function(i){
                    lines(xstart[i]:xend[i],(probeEffect+intercept)[xstart[i]:xend[i]],xaxs="i",...)
                    })
        }       
        else{
            if(!add) plot(x.pos,probeEffect+intercept,xaxs="i",type=type,xlab=xlab,ylab=ylab,...)
            else lines(x.pos,probeEffect+intercept,xaxs="i",...)
        }
    }
    invisible(probeEffect)
    }
)
