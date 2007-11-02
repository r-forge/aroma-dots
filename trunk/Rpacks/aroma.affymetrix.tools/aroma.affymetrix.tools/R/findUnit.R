findUnit<-function(cdf,names){ #given group (probeset) name, returns the unit number (cdf dependent) where group is found
    groupNames <- readCdfGroupNames(cdf$.pathname);
    blanks <- which(unlist(lapply(groupNames, function(x){all(x=="")})));
    groupNames[blanks] <- names(groupNames)[blanks];
    out<-sapply(names,function(gname){
        x<-which(sapply(groupNames,function(x){gname %in% x}))
        if(length(x)>1){ 
            warning(paste("More than 1 unit for", gname," Returning only 1st found"))
            x<-x[1]
        }
        if(length(x)==0) x<-NA
        return(x)
        })
    names(out)<-names
    return(out)
}
