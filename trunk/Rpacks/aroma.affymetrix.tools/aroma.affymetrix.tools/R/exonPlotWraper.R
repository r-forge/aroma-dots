exonPlotWrapper<-function(plm,pbsetUgp,units,previousOutput=NULL,dataToPlot=c("intensities","residuals"),includeProbeAffinity=TRUE,includeChipEffect=TRUE,plot.it=TRUE,...){
    require(aroma.affymetrix)
    dataToPlot<-match.arg(dataToPlot)
    if(is.null(previousOutput)){
        if(!exists("findGroupPositions")) stop("need function 'findGroupPositions'")
        if(!inherits(plm,"ExonRmaPlm") || !plm$mergeGroups)stop("Only currently supports output from ExonRmaPlm with option mergeGroups=TRUE")
        cdfUnits<-getCdf(plm)
        exonPositions<-findGroupPositions(cdfUnits,units=units,pbsetUgp=pbsetUgp)
        nProbes<-readCdfNbrOfCellsPerUnitGroup(getPathname(cdfUnits),unit=units)
        if(!exists("readUnits.RmaPlm")){
            warning("need method of 'readUnits' for RmaPlm to use options 'includeProbeEffect', 'includeChipEffect', or to choose to plot residuals; will instead return standard intensities")
            dataList<-readUnits(plm,units=units)
        }
        else{
            if(dataToPlot=="residuals"){
                includeProbeAffinity<-FALSE
                includeChipEffect<-FALSE
            }
            dataList<-readUnits(plm,units=units,includeProbeAffinity=includeProbeAffinity,includeChipEffect=includeChipEffect)
        }
        nunits<-length(units)
    }
    else{
        dataList<-previousOutput$dataList
        nProbes<-previousOutput$nProbes
        exonPositions<-previousOutput$exonPositions
        nunits<-length(exonPositions)
    }
    
    if(plot.it){
        require(biomaRt)
        require(EnsemblGraphs)
        mart<-useMart("ensembl",dataset= "hsapiens_gene_ensembl")
    
        for(kk in 1:nunits){
            exonPlot(exonarray=log2(dataList[[kk]][[1]]$int), probepos=exonPositions[[kk]][,c(1,3,4)],
                nprobes=nProbes[[kk]],gene = names(nProbes)[kk],biomart=mart,...)
        }
    }
    invisible(list(dataList=dataList,exonPositions=exonPositions,nProbes=nProbes))
}
