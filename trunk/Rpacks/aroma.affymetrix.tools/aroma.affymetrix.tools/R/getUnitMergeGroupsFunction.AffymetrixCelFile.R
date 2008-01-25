setMethodS3("getUnitMergeGroupsFunction","AffymetrixCelFile",function (this,logscale=FALSE,...) 
{
    mergeFcn<-function(groupList,arrays,fields="intensities"){
        perFieldFunction<-function(field){
            fieldList<-lapply(groupList,function(group){return(.subset2(group,field))})
            nRows<-sapply(fieldList,nrow)
            #nCols<-unique(sapply(fieldList,ncol))
            x<-vector()# matrix(nrow=0,ncol=nCols)
            for(i in 1:length(fieldList)){
                x<-c(x,.subset2(fieldList,i))
            }
            if(logscale)x<-log2(x)
            return(x)    
        }
        
        allFieldsList<-lapply(fields,perFieldFunction)
        names(allFieldsList)<-fields
        return(allFieldsList)
        }
    return(mergeFcn)
})
