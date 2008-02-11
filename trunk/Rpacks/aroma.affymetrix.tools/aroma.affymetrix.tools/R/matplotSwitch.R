#function that plots either lines or image

matplotSwitch<-function(mat,type,colScheme,nProbesPerExon,breaks=NULL,
    valuesPlot=c("positive","zeroSym","original"),ylim=NULL,matplotType="l",
    yaxisLabel=FALSE,hlines=yaxisLabel,arrayNames=NULL,arrayFactor=NULL,col=NULL,
    lwd=1, lty=1, pch=19,plotLast=NULL,...
    ){
        nbrOfArrays<-ncol(mat)
        arrays<-1:nbrOfArrays
        valuesPlot<-match.arg(valuesPlot)
        if(valuesPlot=="positive"){ 
            mat<-abs(mat)
         }
        if(is.null(arrayFactor)){
            arrayFactor<-as.factor(arrays)
            hlines<-FALSE}
        else{
            if(length(arrayFactor)!= nbrOfArrays) stop("'arrayFactor' must be of same length as number of arrays plotted")
            hlines<-TRUE
        }
        arrayOrder<-order(arrayFactor)
        arrays<-arrays[arrayOrder]
        arrayFactor<-arrayFactor[arrayOrder]
        if(!is.null(plotLast)) plotLast<-match(plotLast,arrayOrder)
        repArrayValue<-function(x){
            if(length(x)!=nbrOfArrays && length(x)!=1) warning("Plotting parameter not same length as number of arrays")
            x<-rep(x,length=nbrOfArrays)
            x<-x[arrayOrder]
            return(x)
        }
        if(is.null(col) && !is.null(arrayFactor)) col<-rep(palette(),length=nlevels(arrayFactor))[arrayFactor]
        else{col<-repArrayValue(col)    }
        lwd<-repArrayValue(lwd)
        pch<-repArrayValue(pch)
        lty<-repArrayValue(lty)
        arrayNames<-arrayNames[arrayOrder]
        mat<-mat[,arrayOrder]
        
        
        
        if(is.null(breaks) && type=="image"){
            if(valuesPlot=="positive") breaks<-seq(0,max(mat,na.rm=TRUE),length=length(colScheme)+1)
            if(valuesPlot=="zeroSym") breaks<-seq(-max(abs(mat),na.rm=TRUE),max(abs(mat),na.rm=TRUE),length=length(colScheme)+1)
            if(valuesPlot=="original") breaks<-seq(min(mat,na.rm=TRUE),max(mat,na.rm=TRUE),length=length(colScheme)+1)
        }
        nExons<-length(nProbesPerExon)
        if(nrow(mat)==nExons) level<-"exon"
        else{
            if(nrow(mat)==sum(nProbesPerExon))level<-"probe"
            else stop(paste("Incorrect number of rows in scores matrix: Should have",sum(nProbesPerExon),"or",nExons,"rows"))
        }
        
        #check that probesets in right order
        pbsetNames<-names(nProbesPerExon)
        if(!is.null(pbsetNames)){
            reOrderExons<-order(pbsetNames,decreasing=FALSE)
            reOrderProbes<-order(rep(pbsetNames,times=nProbesPerExon))
            if(level=="probe") mat<-mat[reOrderProbes,]
            if(level=="exon") mat<-mat[reOrderExons,]
            nProbesPerExon<-nProbesPerExon[reOrderExons]
            
        }
        exonBoundaries<-c(0, cumsum(nProbesPerExon)) + 0.5
        midpoints<-(tail(exonBoundaries,-1)+head(exonBoundaries,-1))/2
        
        if(type=="image"){
            x<-switch(level,"probe"=1:sum(nProbesPerExon),"exon"=exonBoundaries) #boundaries -- may not be evenly spaced...
            image(x,y=1:ncol(mat),mat,breaks=breaks,
                col=colScheme,yaxt="n",xlab="",ylab="",xaxt="n",...)
            if(yaxisLabel) mapply(axis,label=arrayNames,col.axis=col,at=1:length(arrays),MoreArgs=list(side=2),las=2)
            if(hlines) {
                abline(h=c(0,cumsum(table(arrayFactor)))+.5,col="grey")
                axis(2,at=c(0,cumsum(table(arrayFactor)))+.5,labels=FALSE,tck=-1,col="grey")
            }
                        
        }
        if(type=="lines"){
            x<-switch(level,"probe"=1:sum(nProbesPerExon),"exon"=midpoints) #
            matplot(x,mat,col=col,xaxs="i",xlab="",ylab="",type=matplotType,pch=pch,
                xlim=range(exonBoundaries),xaxt="n",ylim=ylim,lwd=lwd,lty=lty,...)
            if(!is.null(plotLast)){
                matplot(x,mat[,plotLast],col=col[plotLast],type=matplotType,pch=pch[plotLast],
                    lwd=lwd[plotLast],lty=lty[plotLast],add=T,...)
                
            }
        }
        abline(v=exonBoundaries,col="grey",lwd=1.5)
        axis(1,at=exonBoundaries,labels=FALSE,col="grey")
        box()  
        invisible(list(exonBoundaries=exonBoundaries,midpoints=midpoints,arrayOrder=arrayOrder,type=type,reOrderProbes=reOrderProbes,reOrderExons=reOrderExons,breaks=breaks))  
    }
