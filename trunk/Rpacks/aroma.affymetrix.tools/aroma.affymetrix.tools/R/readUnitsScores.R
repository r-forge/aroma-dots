#note: unit names must be unique!
readUnitsScores<-function(names,filename,path=NULL,verbose=0){ #get the scores for specific unit *names*; was called 'findScores'
    #get all chunk file names
    if(!is.null(path)) chunkfiles<-list.files(path)[grep(paste("^",filename,",",sep=""),list.files(path))]
    else chunkfiles<-list.files()[grep(paste("^",filename,",",sep=""),list.files())]
    nbrOfChunks<-length(chunkfiles)
    fullpathnames<-paste(path,chunkfiles,sep="/")
    verbose && enter(verbose, "Scanning unit data")
    x<-list()
    i<-1
    currnames<-names
    while(length(currnames)>0) {
        verbose && enter(verbose, paste("Chunk #",
            i,"of",nbrOfChunks))
        cat(paste("\n\tFilename:",fullpathnames[i],"\n"))
        x <- c(x,loadObject(file=fullpathnames[i]))
        keep<-names(x) %in% currnames
        currnames<-currnames[!(currnames %in% names(x))] #only keep looking for names not found
        x<-x[keep]
        verbose && exit(verbose)
        i<-i+1
    }
    verbose && exit(verbose)
    gc <- gc()
    return(x[names])
}
