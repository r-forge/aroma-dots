##Note, for 'component names' to work, must have group structure...

#could really just call 'applyScores' and unlist it?

flatScores<-function(filename,path=NULL,FUN=function(x){x},componentNames=NULL){ #apply FUN to each score and returns result as long vector; default just returns long vector of original scores ...
    
    #get all chunk file names
    if(!is.null(path)) chunkfiles<-list.files(path)[grep(paste("^",filename,",",sep=""),list.files(path))]
    else chunkfiles<-list.files()[grep(paste("^",filename,",",sep=""),list.files())]
    nbrOfChunks<-length(chunkfiles)
    fullpathnames<-paste(path,chunkfiles,sep="/")
    verbose && enter(verbose, "Extracting unit data")
    
    x<-list()
    for(i in 1:nbrOfChunks) {
        verbose && enter(verbose, paste("Chunk #",
            i,"of",nbrOfChunks))
        if(is.null(componentNames)) x<-c(x,loadObject(file=fullpathnames[i]))
        else x<-c(x,lapply(loadObject(file=fullpathnames[i]),extractLists,names=componentNames)) #take only specific components
        verbose && exit(verbose)
        }
    verbose && exit(verbose)
    gc()
    return(FUN(unlist(x)))
}
