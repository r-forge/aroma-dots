#create 2 functions: 
#--createProbesetCDF() creates from scratch the original (no transcript cluster) cdf from the probe information
#   use file:HuEx-1_0-st-v2.probe.tab
#--createTranscriptCDF() creates core/extended/full/etc. cdf using the annotation information and the original cdf file
#   use file: HuEx-1_0-st-v2.na23.hg18.probeset.csv

# nbrOfUnits(cdfOrig)
####    [1] 1432154
#    nbrOfCells(cdfOrig)
####    [1] 6553600
#    nbrOfColumns(cdfOrig) #number of columns on array
####    [1] 2560
#    nbrOfQcUnits(cdfOrig) #number of quality control elements on array
####    [1] 0
#    nbrOfRows(cdfOrig) #number of rows on array
####    [1] 2560


createTranscriptCDF<-function(probesetCdf, csvAnnotFile,linesPerRead=1000,type,dirLoc=getPath(probesetCdf),addTag=""){
    type<-match.arg(type,c("core","extended","full","main","control","all"),keepColumns=c(1,7,16,39),ncolumns=39)
    npbset<-nbrOfUnits(probesetCdf)
    if(!is.character(addTag)) stop("'addTag' must be character")
    #if(!file.exists(csvAnnotFile)) stop("invalid probeset annotation file name")
     #   nLinesFile<-length(readLines(csvAnnotFile))
#    nChunks<-npbset%/%linesPerRead+1 #may not be the same size as annotation file, but close approximation
    cat(paste(npbset,"probesets in original cdf\n"))
 
    keepType<-switch(type,"core"="core","extended"=c("core","extended"),"full"=c("core","extended","full"),
        "main"="main","control"=c("control->affx","control->chip","control->bgp->antigenomic","control->bgp->genomic","normgene->exon","normgene->intron"),
        "all"=NULL)
    if(type%in% c("core","extended","full")) columnCompare<-"level"
    if(type%in%c("main","control")) columnCompare<-"probeset_type"
    if(!all(keepColumns %in% 1:ncolumns)) stop(paste("'keepColumns' must be a subset of 1-",ncolumns,sep=""))
    colClasses<-rep("NULL",times=ncolumns) #don't keep most of the columns
    colClasses[keepColumns]<-"character" #default corresponds to 1=probeset_id,7=transcript_cluster_id,16=level,39=probeset_type; keep as characters
    
    #get first value of first 500 lines to find how many to skip
    firstLines<-substr(scan(csvAnnotFile,nmax=500,what="",quote="",sep="\n"),start=1,stop=1) #get first character of each line
    startRead<-match(F,sapply(firstLines,function(x){x=="#"})) #first time not start with #

    #take out nrow in final version
    pbsetFinalTab<-read.table(csvAnnotFile,sep=",",skip=startRead-1,header=T,stringsAsFactors=F,colClasses=colClasses,row.names=NULL)
    indexTemp<-indexOf(probesetCdf,names=pbsetFinalTab[,"probeset_id"])
    pbsetcount<-nrow(pbsetFinalTab)
    if(pbsetcount > npbset) warning(paste("Annotation has",pbsetcount,"lines (probesets); original cdf only has",npbset))
    
    if(any(is.na(indexTemp))) {
      warning(paste(length(which(is.na(indexTemp))),"probesets in annotation are not in original cdf; will be discarded"))
      out<-pbsetFinalTab[is.na(indexTemp),"probeset_id"]

      pbsetFinalTab<-pbsetFinalTab[!is.na(indexTemp),]
    }
    else out<-NULL
    pbsetcount<-nrow(pbsetFinalTab)
    if(pbsetcount < npbset) {
      cat(paste("Fewer (valid) probesets in annotation than orginal cdf (",pbsetcount,"compared to",npbset,")\n"))
      
    }
    if(!is.null(keepType)){
        keepPbSet<-which(pbsetFinalTab[,columnCompare]%in%keepType)
        pbsetFinalTab<-pbsetFinalTab[keepPbSet,]
        cat(paste(length(keepPbSet),"probesets match requirement (out of",pbsetcount,"valid probesets in annotation\n"))
        if(length(keepPbSet)==0) invisible(NULL)
    }
    rm(indexTemp)
    gc()
                   
#    pbsetFinalTab<-matrix(nrow=0,ncol=length(which(colClasses!="NULL")))
#    header<-unlist(read.table(csvAnnotFile,sep=",",header=F,skip=startRead-1,nrow=1,stringsAsFactors=F,colClasses=colClasses),use.names=F)
#    colnames(pbsetFinalTab)<-header
#    pbsetcount<-0
#per chunk, add only those in type want; considerable memory savings if type=core!;
#    for(i in 1:nChunks){ #for(i in 1:3){ # should be 
#      
#        if(i==nChunks){
#          linesPerRead<- -1 #get remaining lines regardless
#          cat("Reading remaining lines in iteration ",i," out of ", nChunks," total iterations\n",sep="")
#        }
#        else   cat("Reading ",linesPerRead," lines in iteration ",i," out of ", nChunks," total iterations\n",sep="")
#
#        pbsetTab<-read.table(csvAnnotFile,sep=",",nrows=linesPerRead,skip=startRead,header=F,stringsAsFactors=F,colClasses=colClasses,row.names=NULL)
#        names(pbsetTab)<-header
#        pbsetcount<-pbsetcount+nrow(pbsetTab)
#        indexTemp<-indexOf(probesetCdf,names=pbsetTab[,"probeset_id"])
#    
#        if(any(is.na(indexTemp))) warning("Probesets in annotation are not all in original cdf")
#        rm(indexTemp)
#        gc()
#    
#        keepPbSet<-which(pbsetTab[,columnCompare]%in%keepType)
#        pbsetFinalTab<-rbind(pbsetFinalTab,pbsetTab[keepPbSet,])
#        if(i != nChunks) startRead<-startRead+linesPerRead
#      }    
#    if(pbsetcount > nbrOfUnits(probesetCdf)) warning("Annotation has more probesets than original cdf") ##put back!!

    getCdfInfo<-function(uname){
        whichUnit<-which(pbsetFinalTab[,"transcript_cluster_id"]==uname)
        pbsetMat<-pbsetFinalTab[whichUnit,]
        groupNames<-pbsetMat[,"probeset_id"]
        groupUnits<-indexOf(probesetCdf,names=groupNames)
        if(any(is.na(groupUnits))){ warning(paste("not all group names belonging to",uname,"given in annotation file were found in the given cdf. These groups have been dropped"))
            groupUnits<-na.omit(groupUnits)}
      
      groupCdf<-readCdf(getPathname(probesetCdf),units=groupUnits)
#        if(type%in% c("core","extended","full","main")) pbsetType<-as.character(unique(pbsetMat[,"level"])) #change from Ken's so that type indicates level if only have main design probes
#        else pbsetType<- as.character(unique(pbsetMat[,"probeset_type"]))
#        if(length(pbsetType)!=1) stop("Programming Error 1")

       unittype <- .subset2(.subset2(groupCdf,1), "unittype")
       unitdirection <- .subset2(.subset2(groupCdf,1), "unitdirection")
         natoms <- sum(unlist(lapply(groupCdf, .subset2, "natoms")))
        ncells <- sum(unlist(lapply(groupCdf, .subset2, "ncells")))
        ncellsperatom <- .subset2(.subset2(groupCdf,1), "ncellsperatom")
        groups <- lapply(lapply(groupCdf, .subset2, "groups"), .subset2, 1)
        return(list(groups=groups, unittype=unittype, unitdirection=unitdirection, natoms=natoms, ncells=ncells, ncellsperatom=ncellsperatom, unitnumber=unitnumber))
    }
    outList<-list()
    nprobesetsAdded<-0
    
    unitNames<-unique(pbsetFinalTab[,"transcript_cluster_id"])
    unitnumber<-0
    cat(paste("Starting iteration over",length(unitNames),"units\nDone so far:"))
    for(i in 1:length(unitNames)){#can't get the unitnumber to increment within lapply...Ken did though
      unitnumber <- unitnumber + 1
      if(unitnumber %% 100 == 0) cat("\t",unitnumber)
      if(unitnumber %% 2000 == 0) cat("\n")
      outList[[i]]<-getCdfInfo(unitNames[[i]])
      nprobesetsAdded<-nprobesetsAdded+length(outList[[i]]$groups)
    }

    names(outList)<-unitNames
    nunits <- length(outList)
    nprobesets <- nrow(pbsetFinalTab)
    if(nprobesetsAdded!=nprobesets) stop("Programming error: not all probesets in matrix added to new cdf")
    cat(paste("New cdf consists of",nunits,"units and",nprobesets,"probesets\n"))
    rm(pbsetFinalTab)    
    gc()

#Copy header and qc info from original cdf
    qc <- readCdfQc(getPathname(probesetCdf))
    hdr <- readCdfHeader(getPathname(probesetCdf))
    filename<-paste(getChipType(probesetCdf), paste(type,addTag,".cdf",sep=""), sep=",")
    cat(paste("New cdf has filename",filename,"\n"))
    outfile<-paste(dirLoc,"/",filename,sep="")
    hdr$chiptype <- paste(getChipType(probesetCdf),type,sep=",")
    hdr$filename <- outfile
    hdr$probesets <- nprobesets
    hdr$nunits <- nunits

    if (length(qc)==0) { #i.e. no qc elements;
      cat("No inherited qc items. Creating fake qc portion for writeCdf()\n")
      qc <- list(x=0, y=0, indices=1, length=25, pm=FALSE, background=FALSE, type="unknown", ncells=1);
      qc <- list(qc);
      hdr$nqcunits <- 1;
    }
    #return(list(outList,qc,hdr))
    writeCdf(outfile, cdfheader=hdr, cdf=outList, cdfqc=qc, overwrite=TRUE, verbose=10)
    return(out)
}
