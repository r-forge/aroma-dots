findUnit<-function(cdf,names,returnNames=TRUE){ #given group (probeset) name(s), returns the unit name or unit number where group is found
    groupNames <- readCdfGroupNames(cdf$.pathname, truncateGroupNames=FALSE);
    if(length(which(unlist(lapply(groupNames, function(x){all(x=="")}))))>0) stop("Programming error related to group names equal to unit names. Contact author of program.")
#    blanks<-which(unlist(lapply(groupNames, function(x){all(x=="")})))
#    groupNames[blanks] <- names(groupNames)[blanks];
    out<-sapply(names,function(gname){
        x<-which(sapply(groupNames,function(x){gname %in% x}))
        if(length(x)>1){ 
            warning(paste("More than 1 unit for", gname,"."))
            #x<-x[1]
        }
        if(length(x)==0) x<-NA
        return(x)
        })
    if(returnNames){
        names(out)<-names
        unitNames<-sapply(out,function(unitNumb){getUnitNames(cdf,units=unitNumb)})
        names(unitNames)<-names
        if(is.list(out)) return(list(unitNumbers=out,unitNames=unitNames))
        else{
            returnObj<-data.frame(unitNumbers=out,unitNames=unitNames)
            rownames(returnObj)<-names
            return(returnObj)
        }
    }
    else{
        names(out)<-names
        return(out)
    }
}
