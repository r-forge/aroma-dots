exonPlotWrapper<-function(plm,pbsetUgp,units,previousOutput=NULL,dataToPlot=c("intensities","residuals"),includeProbeAffinity=TRUE,includeChipEffect=TRUE,plot.it=TRUE,...){
    dataToPlot<-match.arg(dataToPlot)
    if(!exists("readUnits.RmaPlm")) stop("need method of 'readUnits' for RmaPlm")
    if(!exists("findGroupPositions")) stop("need function 'findGroupPositions'")
    if(!inherits(plm,"ExonRmaPlm") || !plm$mergeGroups)stop("Only currently supports output from ExonRmaPlm with option mergeGroups=TRUE")
    if(is.null(previousOutput)){
        cdfUnits<-getCdf(plm)
        exonPositions<-findGroupPositions(cdfUnits,units=units,pbsetUgp=pbsetUgp)
        nProbes<-readCdfNbrOfCellsPerUnitGroup(getPathname(cdfUnits),unit=units)
        if(dataToPlot=="residuals"){
            includeProbeAffinity<-FALSE
            includeChipEffect<-FALSE
        }
        dataList<-readUnits(plm,units=units,includeProbeAffinity=includeProbeAffinity,includeChipEffect=includeChipEffect)
    }
    else{
        dataList<-previousOutput$dataList
        nProbes<-previousOutput$nProbes
        exonPositions<-previousOutput$exonPositions
    }
    mart<-useMart("ensembl",dataset= "hsapiens_gene_ensembl")
    
    if(plot.it){
        for(kk in 1:length(units)){
            exonPlot(exonarray=log2(dataList[[kk]][[1]]$int), probepos=exonPositions[[kk]][,c(1,3,4)],
                nprobes=nProbes[[kk]],gene = names(nProbes)[kk],biomart=mart,...)
        }
    }
    invisible(list(dataList=dataList,exonPositions=exonPositions,nProbes=nProbes))
}
