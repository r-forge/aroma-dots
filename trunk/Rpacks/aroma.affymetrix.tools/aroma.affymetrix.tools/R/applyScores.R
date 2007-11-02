#apply FUN to each group element (so FUN needs to take into account sublist structure, if any)
# returns with unit/group structure, but perhaps not intensity sublist depending on FUN
#default just recreates entire list
applyScores<-function(filename,path=NULL,FUN=function(x){x}){ 
#get all chunk file names
    if(!is.null(path)) chunkfiles<-list.files(path)[grep(paste("^",filename,",",sep=""),list.files(path))]
    else chunkfiles<-list.files()[grep(paste("^",filename,",",sep=""),list.files())]
    nbrOfChunks<-length(chunkfiles)
    fullpathnames<-paste(path,chunkfiles,sep="/")
    verbose && enter(verbose, "Extracting unit data and applying function to groups")
    
    x<-list()
    for(i in 1:nbrOfChunks) {
        verbose && enter(verbose, paste("Chunk ",i,"of",nbrOfChunks))
        y<-loadObject(file=fullpathnames[i])
        y<-lapply(y,function(groupList){lapply(groupList,FUN)})
        x<-c(x,y)
        rm(y)
        gc()
        verbose && exit(verbose)    
        }
    verbose && exit(verbose)
    gc()
    return(x)
}
