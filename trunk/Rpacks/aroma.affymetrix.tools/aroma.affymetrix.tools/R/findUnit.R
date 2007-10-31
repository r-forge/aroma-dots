findUnit<-function(cdf, names,groupList=NULL){
    if(is.null(groupList)) groupList <- readCdfGroupNames(cdf$.pathname)
    listId<-sapply(names,function(name){ 
        units<-which(sapply(groupList,function(x){name %in% x}))
        unitNames<-names(groupList)[units]
        if(length(units)>1) {
            warning(paste("Non-unique group to unit identification for '", name,"'. Returning first match."))
            unitNames<-unitNames[1]
        }
        if(length(units)==0) return("None")
        return(unitNames)
        })
    return(listId)
}
