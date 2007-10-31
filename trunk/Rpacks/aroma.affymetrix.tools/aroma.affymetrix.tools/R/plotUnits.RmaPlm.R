setMethodS3("plotUnits","RmaPlm",function(this,units,arrays=NULL,joinGroups=T,includeProbeEffect=F,includeChipEffect=T, 
    plotRaw=F,
    plotProbeEffect=F,probePlotArgs=list(lwd=3,lty=1,col="black"),plotByArray=T,
    plotIntercept=F,interceptPlotArgs=list(col="grey",lwd=2,lty=2),plot.it=T,...){
    
    if(length(units)>1) stop("units must be single integer -- only can plot 1 unit") ##could change this...
    if(plotRaw){
        includeProbeEffect<-includeChipEffect<-TRUE    
    }
    ds<-getDataSet(this) #input dataset
    if(plotProbeEffect || !includeProbeEffect || !includeChipEffect){ 
        dataMat<-plotUnits(ds,units=units,arrays=NULL,joinGroups=joinGroups,plot.it=F,...)#just return from plotUnits, not plot
        ces<-getChipEffects(this)
        paf<-getProbeAffinities(this)
        cdf<-getCdf(this)
        cellIndices<-unlist(getCellIndices(cdf,units=units)) ##destroy grouping so get solid matrix
        chipEffect<-log2(readUnits(ces,units=units)[[1]][[1]]$theta)
        intercept<-mean(chipEffect)
        if(!includeProbeEffect){ #subtract out Probe Effect
            probeEffect<-log2(getData(paf,cellIndices,field="intensities"))[,1] #could add stdev plot too...
            dataMat<-sweep(dataMat,1,probeEffect)
        }
        if(!includeChipEffect){ #subtract out Chip Effect
            dataMat<-sweep(dataMat,2,chipEffect-intercept)            
        }
        dsPlot<-plotUnits(ds,dataMat=dataMat,arrays=arrays,logscale=F,ylab="Log Intensities",units=units,joinGroups=joinGroups,plot.it=plot.it,plotByArray=plotByArray,...)
        if(!plotByArray) plot.it<-F
        if(plotIntercept && plot.it) do.call(abline,c(list(h=intercept),interceptPlotArgs))
        pafPlot<-do.call(plotUnits,args=c(list(this=paf,units=units,add=T,joinGroups=joinGroups,plot.it=(plot.it & plotProbeEffect) ,intercept=intercept),probePlotArgs))
    }
    else{ dsPlot<-plotUnits(ds,units=units,arrays=arrays,joinGroups=joinGroups,plot.it=plot.it,plotByArray=plotByArray,...)}
    invisible(dsPlot)#,probeEffect=probeEffect,chipEffect=chipEffect,intercept=intercept))   
})
