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

#csvAnnotFile<-"Annotation//CSVFiles//HuEx-1_0-st-v2.na23.hg18.probeset.csv"
#linesPerRead<-4
#type<-"extended"
#probesetCdf<-AffymetrixCdfFile$fromChipType("HuEx-1_0-st-v2" )
#chipType<-"HuEx-1_0-st-v2"
#createTranscriptCDF(AffymetrixCdfFile$fromChipType("HuEx-1_0-st-v2" ), "Annotation//CSVFiles//HuEx-1_0-st-v2.na23.hg18.probeset.csv",linesPerRead=3,type="extended")

#   main                       part of the main design
#    control->affx              a standard AFFX control
#    control->chip              a chip control
#    control->bgp->antigenomic  contains antigenomic background probes
#    control->bgp->genomic      contains genomic background probes
#    normgene->exon             from an exonic region of a
#                               normalization control gene
#    normgene->intron           from an intronic region of a
#                               normalization control gene
#    rescue->FLmRNA->unmapped   contains probes tiled across
#                               an mRNA transcript which either did not
#                               align to the genome, or aligned poorly

createTranscriptCDF<-function(probesetCdf, csvAnnotFile,linesPerRead=1000,type){
    type<-match.arg(type,c("core","extended","full","main","control"))
    npbset<-nbrOfUnits(probesetCdf)
    if(!file.exists(csvAnnotFile)) stop("invalid probeset annotation file name")
    
    nChunks<-npbset%/%linesPerRead+1

    keepType<-switch(type,"core"="core","extended"=c("core","extended"),"full"=c("core","extended","full"),
        "main"="main","control"=c("control->affx","control->chip","control->bgp->antigenomic","control->bgp->genomic","normgene->exon","normgene->intron"))
    if(type%in% c("core","extended","full")) columnCompare<-"level"
    if(type%in%c("main","control")) columnCompare<-"probeset_type"
    colClasses<-rep("NULL",times=39) #don't keep most of the columns
    colClasses[c(1,7,16,39)]<-"character" #corresponds to 1=probeset_id,7=transcript_cluster_id,16=level,39=probeset_type; keep as characters
    
    #get first value of first 1000 lines to find how many to skip
    pbsetFinalTab<-matrix(nrow=0,ncol=length(which(colClasses!="NULL")))
    firstLines<-substr(scan(csvAnnotFile,nmax=500,what="",quote="",sep="\n"),start=1,stop=1) #get first character of each line
    startRead<-match(F,sapply(firstLines,function(x){x=="#"})) #first time not start with #
    header<-unlist(read.table(csvAnnotFile,sep=",",header=F,skip=startRead-1,nrow=1,stringsAsFactors=F,colClasses=colClasses),use.names=F)
    print(header)
    print(sapply(colClasses,is.null))
    colnames(pbsetFinalTab)<-header
    
    pbsetcount<-0
    #per chunk, add only those in type want; considerable memory savings if type=core!;
    for(i in 1:3){ # should be for(i in 1:nChunks)
        cat("Reading ",i,"th chunk of ", linesPerRead,"\n",sep="")
        pbsetTab<-read.table(csvAnnotFile,sep=",",nrows=linesPerRead,skip=startRead,header=F,stringsAsFactors=F,colClasses=colClasses,row.names=NULL)
        names(pbsetTab)<-header
        pbsetcount<-pbsetcount+nrow(pbsetTab)
        indexTemp<-indexOf(probesetCdf,names=pbsetTab[,"probeset_id"])
        if(any(is.na(indexTemp))) stop("Probesets in annotation are not all in original cdf")
        rm(indexTemp)
        gc()
        keepPbSet<-which(pbsetTab[,columnCompare]%in%keepType)
        pbsetFinalTab<-rbind(pbsetFinalTab,pbsetTab[keepPbSet,])
        if(i != nChunks) startRead<-startRead+linesPerRead
    }    

    if(pbsetcount != nbrOfUnits(probesetCdf)) stop("Original cdf and annotation have different numbers of probesets") ##put back!!
    unitNames<-unique(pbsetFinalTab[,"transcript_cluster_id"])
    getCdfInfo<-function(uname){
        whichUnit<-which(pbsetFinalTab[,"transcript_cluster_id"]==uname)
        pbsetMat<-pbsetFinalTab[whichUnit,]
        groups<-pbsetMat[,"probeset_id"]
        groupUnits<-indexOf(probesetCdf,names=groups)
        groupCdf<-readCdfUnits(getPathname(probesetCdf),units=groupUnits)
        if(type%in% c("core","extended","full","main")) pbsetType<-as.character(unique(pbsetMat[,"level"])) #change from Ken's so that type indicates level if only have main design probes
        else pbsetType<- as.character(unique(pbsetMat[,"probeset_type"]))
        if(length(pbsetType)!=1) stop("Programming Error 1")
        direction<-unique(sapply(groupCdf,function(x){x$direction}))
        if(length(direction)!=1) stop("Programming Error 2")
        groups<-lapply(groupCdf,function(x){x$groups[[1]]})
        return(list(type=pbsetType,direction=direction,groups=groups))
    }
    outList<-lapply(unitNames,getCdfInfo)
    names(outList)<-unitNames
    nunits <- length(outList)
    nprobesets <- nrow(pbsetFinalTab)
    rm(pbsetFinalTab)    
    gc()

#Copy header and qc info from original cdf
    qc <- readCdfQc(getPathname(probesetCdf))
    hdr <- readCdfHeader(getPathname(probesetCdf))
    filename<-paste(getChipType(probesetCdf), paste(type,".cdf",sep=""), sep=",")
    outfile<-paste(getPath(probesetCdf),"/",filename,sep="")
    hdr$chiptype <- paste(getChipType(probesetCdf),type,sep=",")
    hdr$filename <- outfile
    hdr$probesets <- nprobesets
    hdr$nunits <- nunits

    if (length(qc)==0) { #i.e. no qc elements;
      qc <- list(x=0, y=0, indices=1, length=25, pm=FALSE, background=FALSE, type="unknown", ncells=1);
      qc <- list(qc);
      hdr$nqcunits <- 1;
    }
    #return(list(outList,qc,hdr))
    writeCdf(outfile, cdfheader=hdr, cdf=outList, cdfqc=qc, overwrite=TRUE, verbose=10)
}
