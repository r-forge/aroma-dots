#need create 2 functions: 
#--createProbesetCDF() creates from scratch the original (no transcript cluster) cdf from the probe information
#   use file:HuEx-1_0-st-v2.probe.tab
#--createTranscriptCDF() creates core/extended/full/etc. cdf using the annotation information and the original cdf file
#   use file: HuEx-1_0-st-v2.na23.hg18.probeset.csv

createTranscriptCDF<-function(probesetCdf, csvAnnotFile,type="all",dirLoc=getPath(probesetCdf),tag=type,
    keepColumns=c(1,7,16,19,39),ncolumns=39,keepType=NULL,columnCompare=NULL,nrowRead=NULL,writeCdf=TRUE)
{
#if user defined, keepType is the (charcater) value to compare entries of columnCompare for probesets to keep
#nrowRead allows set only read certain number of lines -- for testing
#writeCdf whether to call function to write cdf, otherwise returns list -- for testing
#keepColumns are the columns in the annot. file to keep when do comparisons
#ncolumns are the number of (total) columns in annot. file -- could find it, but safer to just have user give
#tag added to name: <chiptype>,tag.cdf (can contain commas itself)
#example (a test run for keeping only probesets with 3 or 4 probes):
#x<-createTranscriptCDF(cdfProbeset, csvAnnotFile, keepType=c("3","4"),columnCompare="probe_count",
#  dirLoc=getwd(),tag="temp",nrow=500, keepColumns=c(1,6,7,16,19,39),writeCdf=FALSE)

    type<-match.arg(type,c("core","extended","full","main","control","all","cds"))
    npbset<-nbrOfUnits(probesetCdf)
    if(!is.character(tag)) stop("'tag' must be character")
    if(!file.exists(csvAnnotFile)) stop("invalid probeset annotation file name")
    if(length(columnCompare)>1 || !is.character(columnCompare)) stop("Invalid format for 'columnCompare' variable")
    cat(paste(npbset,"probesets in original cdf\n"))
 
    #set up the comparisons to do for paring down
    if(is.null(keepType)) keepType<-switch(type,"core"="core","extended"=c("core","extended"),"full"=c("core","extended","full"),
        "main"="main","control"=c("control->affx","control->chip","control->bgp->antigenomic","control->bgp->genomic","normgene->exon","normgene->intron"),
        "all"=NULL,"cds"="1")
    if(is.null(columnCompare)){
        if(type%in% c("core","extended","full")) columnCompare<-"level"
        if(type%in%c("main","control")) columnCompare<-"probeset_type"
        if(type%in% "cds") columnCompare<-"has_cds"
    }

    #designate which columns in file to keep
    if(!all(keepColumns %in% 1:ncolumns)) stop(paste("'keepColumns' must be a subset of 1-",ncolumns,sep=""))
    colClasses<-rep("NULL",times=ncolumns)
    colClasses[keepColumns]<-"character" #default corresponds to 1=probeset_id,7=transcript_cluster_id,16=level,19=has_cds,39=probeset_type; keep as characters
    
    #get first value of first 500 lines to find how many to skip -- don't really need because R comments out those rows anyway...
    firstLines<-substr(scan(csvAnnotFile,nmax=500,what="",quote="",sep="\n"),start=1,stop=1) #get first character of each line
    startRead<-match(FALSE,sapply(firstLines,function(x){x=="#"})) #first time not start with # 
    header<-unlist(read.table(csvAnnotFile,sep=",",header=FALSE,skip=startRead-1,nrow=1,stringsAsFactors=FALSE,colClasses=colClasses),use.names=FALSE)
    print(header)
    if(! columnCompare %in% header) stop(paste(columnCompare,"is not the name of a column in",csvAnnotFile))
    if(!all(c("probeset_id","transcript_cluster_id") %in% header) ) stop(paste("'csvAnnotFile' must contain headers 'probeset_id' and 'transcript_cluster_id'"))


    #take out nrow in final version
    if(!is.null(nrowRead)) pbsetFinalTab<-read.table(csvAnnotFile,sep=",",skip=startRead-1,nrow=nrowRead,header=TRUE,stringsAsFactors=FALSE,colClasses=colClasses,row.names=NULL)
    else pbsetFinalTab<-read.table(csvAnnotFile,sep=",",skip=startRead-1,header=TRUE,stringsAsFactors=FALSE,colClasses=colClasses,row.names=NULL)
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
                   
    getCdfInfo<-function(uname){
        whichUnit<-which(pbsetFinalTab[,"transcript_cluster_id"]==uname)
        pbsetMat<-pbsetFinalTab[whichUnit,]
        groupNames<-pbsetMat[,"probeset_id"]
        groupNames<-sort(groupNames)
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
        out<-list(groups=groups, unittype=unittype, unitdirection=unitdirection, natoms=natoms, ncells=ncells, ncellsperatom=ncellsperatom, unitnumber=unitnumber)
        return(out)
    }
    outList<-list()
    nprobesetsAdded<-0
    
    unitNames<-unique(pbsetFinalTab[,"transcript_cluster_id"])
    unitnumber<-0
    cat(paste("Starting iteration over",length(unitNames),"units\n"))
    if(length(unitNames)> 100) cat("Done so far:\n")

#    outList<-sapply(1:length(unitNames),function(i){
#      unitnumber<-i
#      if(unitnumber %% 100 == 0) cat("\t",unitnumber)
#      if(unitnumber %% 2000 == 0) cat("\n")
#      outList[[i]]<-getCdfInfo(unitNames[[i]])
#      nprobesetsAdded[[i]]<-length(outList[[i]]$groups)
#    })
#    nprobesetsAdded<-sum(unlist(nprobesetsAdded))
    for(i in 1:length(unitNames)){#can't get the unitnumber to increment within lapply...Ken did though?
      unitnumber <- unitnumber + 1
      if(unitnumber %% 100 == 0) cat("\t",unitnumber)
      if(unitnumber %% 2000 == 0) cat("\n")
      outList[[i]]<-getCdfInfo(unitNames[[i]])
      nprobesetsAdded<-nprobesetsAdded+length(outList[[i]]$groups)
    }
    cat("\n")
    names(outList)<-unitNames
    nunits <- length(outList)
    nprobesets <- nrow(pbsetFinalTab)
    if(nprobesetsAdded != nprobesets) stop("Programming error: not all probesets in matrix added to new cdf")
    cat(paste("New cdf consists of",nunits,"units and",nprobesets,"probesets\n"))
    rm(pbsetFinalTab)    
    gc()

    if(writeCdf){
#Copy header and qc info from original cdf
      qc <- readCdfQc(getPathname(probesetCdf))
      hdr <- readCdfHeader(getPathname(probesetCdf))
      filename<-paste(getChipType(probesetCdf), paste(tag,".cdf",sep=""), sep=",")
      cat(paste("New cdf has filename",filename,"\n"))
      outfile<-paste(dirLoc,"/",filename,sep="")
      hdr$chiptype <- paste(getChipType(probesetCdf),tag,sep=",")
      hdr$filename <- outfile
      hdr$probesets <- nprobesets
      hdr$nunits <- nunits

      if (length(qc)==0) { #i.e. no qc elements;
        cat("No inherited qc items. Creating fake qc portion for writeCdf()\n")
        qc <- list(x=0, y=0, indices=1, length=25, pm=FALSE, background=FALSE, type="unknown", ncells=1);
        qc <- list(qc);
        hdr$nqcunits <- 1;
      }
      writeCdf(outfile, cdfheader=hdr, cdf=outList, cdfqc=qc, overwrite=TRUE, verbose=10)
      return(out)
    }
    else return(outList)
}
