#readFUN function that takes as input 'units' and 'verbose'; 
    #creates the list of input data for calculateFUN; static; 
#calculateFUN takes as input a single element of the list created by 'readFun';
    #create list of scores -- should be same structures as from readUnits, etc. 
    #assume that there is a group, and per group information is saved under a list with some name (like $intensities, only score specific)
    #if want single number summary per unit, still need to save as a single group (like monocell)
calculateScores<-function(cdf,units=NULL, calculateFUN, readFUN,nbrOfArrays, filename,path=NULL,ram=1,verbose=0){
    ######################
    ##Divide up units, etc.
    unitsPerChunk <- ram * 1e+05/nbrOfArrays
    unitsPerChunk <- Arguments$getInteger(unitsPerChunk, range = c(1, Inf))
    verbose <- Arguments$getVerbose(verbose)
    if (verbose) {
        pushState(verbose)
        on.exit(popState(verbose))
    }
    #cdf <- getCdf(this)
    if (is.null(units)) {
        nbrOfUnits <- nbrOfUnits(cdf)
        units <- 1:nbrOfUnits
    }
    else {
        nbrOfUnits <- length(units)
    }
    unitsToDo <- units
#    if (force) {
#        unitsToDo <- units
#    }
#    else { #problem...made force=T
#        unitsToDo <- findUnitsTodo(ws, units = units)
#    }
    verbose && printf(verbose, "Number of units: %d\n", nbrOfUnits)
    nbrOfChunks <- ceiling(nbrOfUnits/unitsPerChunk)
    verbose && printf(verbose, "Number of chunks: %d (%d units/chunk)\n", 
        nbrOfChunks, unitsPerChunk)
    ################
    
    
    ##Calculating Loop
    head <- 1:unitsPerChunk
    verbose && enter(verbose, "Extracting unit data")
    count <- 1
    while (length(unitsToDo) > 0) {
        if (length(unitsToDo) < unitsPerChunk) {
            head <- 1:length(unitsToDo)
        }
        units <- unitsToDo[head]
        verbose && printf(verbose, "Chunk #
            %d of %d (%d units)\n", count, nbrOfChunks, length(units))
        dataList <- readFUN(units = units, verbose = less(verbose))#, stratifyBy = "pm")
        verbose && enter(verbose, "Calculating scores")
        scoresList <- lapply(dataList, FUN = calculateFUN)
        names(scoresList)<-getUnitNames(cdf,unit=units)
        verbose && exit(verbose)
        
        if(nbrOfChunks>1) chunkfilename<-paste(filename,",Chunk",count,".RDATA",sep="")
        else chunkfilename<-paste(filename,",Chunk1",".RDATA",sep="")
        verbose && enter(verbose,paste( "Saving scores, file name:",chunkfilename))
        saveObject(scoresList,file=chunkfilename,path=path)
        verbose && exit(verbose)
    
        #verbose && enter(verbose, "Storing scores") #someday for creating CEL files from arbitrary scores???
#        cdf <- getCellIndices(getCdf(this), units = units, ...) #stratifyBy = "pm", ...)
#        for (kk in seq(ds)) {
#            wf <- getFile(ws, kk)
#            verbose && enter(verbose, sprintf("Array #%d ('%s')", 
#                kk, getName(wf)))
#            data <- lapply(weightsList, function(unit) {
#                lapply(unit, function(group) {
#                  nrow <- nrow(group)
#                  list(intensities = 2^group[, kk], stdvs = rep(1, 
#                    nrow), pixels = rep(1, nrow))
#                })
#            })
#            updateCelUnits(getPathname(wf), cdf = cdf, data = data)
#            verbose && exit(verbose)
#        }
#        verbose && exit(verbose)
        rm(dataList)
        rm(scoresList)
        gc()
        unitsToDo <- unitsToDo[-head]
        count <- count + 1
    }
    gc <- gc()
    verbose && print(verbose, gc)
    verbose && exit(verbose)
} 
