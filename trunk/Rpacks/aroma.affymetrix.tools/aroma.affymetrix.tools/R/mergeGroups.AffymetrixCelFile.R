setMethodS3("mergeGroups","AffymetrixCelFile",function(this,units,fields=NULL,force=FALSE,logscale=FALSE,...){
    unitFcn<-getUnitMergeGroupsFunction(this,logscale=logscale,...) #a function that takes a groupList, arrays, and fields
    data<-readUnits(this,units=units,force=force)
    if(!is.null(fields)) mergeData<-lapply(data,unitFcn,fields=fields)
    else mergeData<-lapply(data,unitFcn)
    return(mergeData)
})
