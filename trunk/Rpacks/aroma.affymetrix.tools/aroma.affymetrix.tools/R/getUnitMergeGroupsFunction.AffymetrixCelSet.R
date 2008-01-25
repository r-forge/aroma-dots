setMethodS3("getUnitMergeGroupsFunction","AffymetrixCelSet",function (this,logscale=FALSE,...) 
{
    mergeFcn<-function(groupList,arrays,fields="intensities"){
        perFieldFunction<-function(field){
            fieldList<-lapply(groupList,function(group){return(.subset2(group,field))})
            nRows<-sapply(fieldList,nrow)
            nCols<-unique(sapply(fieldList,ncol))
            if(length(nCols)>1) stop("Cannot merge groups with data of different sizes")
            x<-matrix(nrow=0,ncol=nCols)
            for(i in 1:length(fieldList)){
                x<-rbind(x,.subset2(fieldList,i))
            }
            #rownames(x)<-names(fieldList)
            if(logscale) x<-log2(x)
            return(x[,arrays,drop=FALSE])    
        }
        
        allFieldsList<-lapply(fields,perFieldFunction)
        names(allFieldsList)<-fields
        return(allFieldsList)
        }
    return(mergeFcn)
})
