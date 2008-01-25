##searches by 'group' whether any of the scores saved for that group meet cutoff. If they do, then that entire group will be returned under that unit name; if no groups in unit match criteria, then unit is dropped. 
##Note: will return entire 'group', regardless if only 1 score associated with group meets criteria
#right now just pulling from how saved, where maybe not logged depending on how the scores were saved. True for searching the .CEL files as well -- should fix
#allows arbitrary logical function... applied per group
#note that default logical function (and also quantity cutoff function) will unlist the group element, so should not be problem regardless of how save it; if multiple items saved (e.g. fitting) then will join these together

searchScores<-function(this=NULL,filename=NULL,nbrOfArrays,units=NULL, cutoff=NULL, quantity=NULL,
    path=NULL,ram=1,verbose=0, type=c("high","low","highlow")){
    
    type<-match.arg(type)
    func<-switch(type, "high"=function(x){x}, "low"=function(x){-x},"highlow"=function(x){abs(x)}) #allows to deal with search for high/low/both extremems
    cutoffFun<-function(x,cutoff){
                x<-func(na.omit(unlist(x))); 
                return(any(x>=cutoff))
            }
    if(!is.null(cutoff)){
        if(is.numeric(cutoff)){ #default logical function
             testfunc<-function(x){cutoffFun(x,cutoff=func(cutoff))}
        }
        else{
            if(!is.function(cutoff)) stop("cutoff must be either numeric cutoff value or function that returns logical based on group vector/list")
            else testfunc<-cutoff
        }
    }
    if(is.null(cutoff) && is.null(quantity))stop("Must specify either 'cutoff' or 'quantity'")
    if(!is.null(filename) && !is.null(this)) stop("must specify either 'this' OR 'filename'")
    if(!is.null(units) && !is.null(filename)) stop("cannot specify units if reading from .RDATA file under 'filename'")
    
    currentList<-NULL
    currentUnits<-NULL
    currentGroups<-NULL
    verbose <- Arguments$getVerbose(verbose)
    if (verbose) {
        pushState(verbose)
        on.exit(popState(verbose))
    }
    if(is.null(filename)){ #read from CEL files
        ######################
        ##Divide up units, etc.
        unitsPerChunk <- ram * 1e+05/nbrOfArrays
        unitsPerChunk <- Arguments$getInteger(unitsPerChunk, range = c(1, Inf))
        
        cdf <- getCdf(this)
        if (is.null(units)) {
            nbrOfUnits <- nbrOfUnits(cdf)
            units <- 1:nbrOfUnits
        }
        else {
            nbrOfUnits <- length(units)
        }
        unitsToDo <- units

        verbose && printf(verbose, "Number of units: %d\n", nbrOfUnits)
        nbrOfChunks <- ceiling(nbrOfUnits/unitsPerChunk)
        verbose && printf(verbose, "Number of chunks: %d (%d units/chunk)\n", 
            nbrOfChunks, unitsPerChunk)

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
            dataList <- readUnits(this,units = units, verbose = less(verbose))#, stratifyBy = "pm")
            
            verbose && enter(verbose, "Finding qualifying scores")
            if(!is.null(cutoff)){
                aboveList <- lapply(dataList, FUN = function(x){ #function per unit
                    out<-lapply(x,FUN=function(group){
                        if(testfunc(group)) return(group)
                        else return(NULL)
                        })
                    keep<-which(unlist(lapply(out,length))!=0)
                    return(out[keep])
                }) 
                keep<-which(unlist(lapply(aboveList,length))!=0)
                aboveList<-aboveList[keep]
             #   aboveUnits<-units[keep]
            }
            else{
                aboveList<-dataList #so all of them are added to current list, and then take top 1000
             #   aboveUnits<-units
            }            
            currentList<-c(currentList,aboveList)
            #currentUnits<-c(currentUnits,aboveUnits)
            rm(dataList)
            rm(aboveList)
            #rm(aboveUnits)
            gc()
            if(!is.null(quantity)){
                quantCutoff<-tail(sort(func(unlist(currentList))),1001)[1] #smallest value of the top 1000
                aboveList <- lapply(dataList, FUN = function(x){
                    out<-lapply(x,FUN=function(group){
                        if(cutoffFun(group,cutoff=quantCutoff)) return(group)
                        else return(NULL)
                        })
                    keep<-which(unlist(lapply(out,length))!=0)
                    return(out[keep])
                }) 
                keep<-which(unlist(lapply(aboveList,length))!=0)
                currentList<-aboveList[keep]
                #currentUnits<-currentUnits[keep]
                rm(aboveList)
                rm(keep)
                gc()
            }
            verbose && exit(verbose)
            
            
            #verbose && enter(verbose, "Storing scores")
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

            unitsToDo <- unitsToDo[-head]
            count <- count + 1
        }
    }
    else{  #read the number of .rdata files
        if(!is.null(path)) chunkfiles<-list.files(path)[grep(paste("^",filename,",",sep=""),list.files(path))]
        else chunkfiles<-list.files()[grep(paste("^",filename,",",sep=""),list.files())]
        nbrOfChunks<-length(chunkfiles)
        fullpathnames<-paste(path,chunkfiles,sep="/")
        verbose && enter(verbose, "Extracting unit data")
        
        for(i in 1:nbrOfChunks) {
            verbose && enter(verbose, paste("Chunk #",
                i,"of",nbrOfChunks))
            cat(paste("\nFilename:",fullpathnames[i],"\n"))
            dataList <- loadObject(file=fullpathnames[i])
            verbose && enter(verbose, "Finding qualifying scores")
            if(!is.null(cutoff)){
                aboveList <- lapply(dataList, FUN = function(x){
                    out<-lapply(x,FUN=function(group){
                        if(testfunc(group)) return(group)
                        else return(NULL)
                        })
                    
                    keep<-which(unlist(lapply(out,length))!=0) #which group/probeset to keep within unit
                    return(out[keep]) #only return those probesets with some element above value
                }) 
                keep<-which(unlist(lapply(aboveList,length))!=0) #only keep the units with something kept; each element of aboveList is list, so =1 if full list, =0 if empty list
                aboveList<-aboveList[keep]
            #    aboveUnits<-names(dataList)[keep]
            }
            else{
                aboveList<-dataList #so all of them are added to current list, and then take top 1000
            #    aboveUnits<-names(dataList)
            }            
            currentList<-c(currentList,aboveList)
            #currentUnits<-c(currentUnits,aboveUnits)
            rm(dataList)
            rm(aboveList)
            #rm(aboveUnits)
            gc()
            if(!is.null(quantity)){
                quantCutoff<-tail(sort(func(unlist(currentList))),quantity)[1] #smallest value of the top 1000
                aboveList <- lapply(currentList, FUN = function(x){
                    out<-lapply(x,FUN=function(group){
                        if(cutoffFun(group,cutoff=quantCutoff)) return(group)
                        else return(NULL)
                        })
                    keep<-which(unlist(lapply(out,length))!=0)
                    return(out[keep])
                }) 
                keep<-which(unlist(lapply(aboveList,length))!=0)
                currentList<-aboveList[keep]
                #currentUnits<-currentUnits[keep]
                rm(aboveList)
                rm(keep)
                gc()
            }
            verbose && exit(verbose)
            
            verbose && exit(verbose)
        }
    }
    gc <- gc()
    verbose && print(verbose, gc)
    verbose && exit(verbose)
    return(currentList)
}
