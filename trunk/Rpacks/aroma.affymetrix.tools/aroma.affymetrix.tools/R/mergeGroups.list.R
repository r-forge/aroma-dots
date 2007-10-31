setMethodS3("mergeGroups","list",function(this,arrays=NULL,fields=NULL,logscale=F,...){
    unitFcn<-getUnitMergeGroupsFunction(this,logscale=logscale,...) #a function that takes a groupList, arrays, and fields
    #if(is.null(arrays)) arrays<-1:nbrOfArrays(this)
    if(!is.null(fields)) mergeData<-lapply(this,unitFcn,fields=fields,arrays=arrays)
    else mergeData<-lapply(this,unitFcn,arrays=arrays)
    return(mergeData)
})
