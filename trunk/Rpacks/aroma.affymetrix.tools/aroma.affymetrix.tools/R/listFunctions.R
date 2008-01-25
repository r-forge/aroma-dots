##to combine lists together
#Takes listOfLists of structure
#      +- name # 1
#      |     +- toplevel #1 
#      |     |   +- data
#      |     +- toplevel #2
#      |         +- data
#      .
#      .
#      .
#      +- name #2
#      |     +- toplevel #1 
#      |     |   +- data
#      |     +- toplevel #2
#      |         +- data
#      .
#      .
#      .
#      +- etc.
#returns list of structure
#      +- toplevel #1
#      |     +- name #1 
#      |     |   +- data
#      |     +- name #2
#      |         +- data
#      .
#      .
#      .
#      +- toplevel #2
#      |     +- name #1 
#      |     |   +- data
#      |     +- name #2
#      |         +- data
#      .
#      .
#      .
#      +- etc.

combineLists<-function(listOfLists){
    listLength<-lapply(listOfLists,length)
    if(length(unique(listLength))!=1) stop("elements of 'listOfList' must be lists of same length")
    listNames<-names(listOfLists[[1]])
    if(!all(sapply(listOfLists,function(x){all(names(x)==listNames)}))) stop("elements of 'listOfList' must be lists with the same names")
    metanames<-names(listOfLists)
    outList<-lapply(1:listLength[[1]],function(i){
        iList<-lapply(listOfLists,function(x){x[[i]]})
    })
    names(outList)<-listNames
    return(outList)
}



##extract same named element(s) from components of list
#name=name(s) of elements
#listOfLists is of structure
#      +- toplevel #1
#      |     +- name #1 
#      |     |   +- data
#      |     +- name #2
#      |         +- data
#      +- toplevel #2
#      |     +- name #1 
#      |     |   +- data
#      |     +- name #2
#      |         +- data
#      +- etc.
#creates list of structure (if names=name #1)
#      +- toplevel #1
#      |     +- name #1 
#      |     |   +- data
#      +- toplevel #2
#      |     +- name #1 
#      |     |   +- data
#      +- etc.
#for aroma.affymetrix, must lapply this to a unit list...

extractLists<-function(listOfLists,names){
    listLength<-lapply(listOfLists,length)
    if(length(unique(listLength))!=1) stop("each element of 'listOfList' must be a list, each of same length")
    listNames<-names(listOfLists[[1]]) #possible names to pick from
    if(!all(sapply(listOfLists,function(x){all(names(x)==listNames)}))) stop("each element of 'listOfList' must be a named list, each with the same names")
    listComponent<-sapply(listOfLists,function(x){all(names %in% names(x))})    
    if(!any(listComponent)) stop("each element of 'listOfList' must be named list, with names in 'names'")
    metanames<-names(listOfLists)
    outList<-lapply(listOfLists,function(x){x[names]})
    names(outList)<-metanames
    return(outList)
}

#listOfData should have each element of the list as a vector or matrix; returns a matrix of the data 'unlisted'. 
#join is whether you should join by row or column (like rbind or cbind)
mergeDataList<-function(listOfData,join=c("row","column")){
    join<-match.arg(join)
    dim1<-dim(listOfData[[1]])
    if(!is.null(dim1)){ #matrices
        if(join=="row") listdim<-sapply(listOfData,ncol) #columns must agree in number
        else listdim<-sapply(listOfData,nrow)
        if(length(unique(listdim))!=1) stop("Must have same dimension in each list element")
        x<-switch(join, row=matrix(ncol=listdim),columns=matrix(nrow=listdim))
        if(join=="row")colnames(x)<-colnames(listOfData[[1]])
        if(join=="column")rownames(x)<-rownames(listOfData[[1]])
        for(i in 1:length(listOfData)){
            x<-switch(join,row=rbind(x,listOfData[[i]]),column=cbind(x,listOfData[[i]]))

        }
        x<-switch(join,row=x[-1,],column=x[,-1]) #get rid of NA 
        return(x)
    }
    if(is.null(dim1)){
        numVec<-sapply(listOfData,function(x){is.vector(x) && (is.numeric(x) || all(is.na(x)))})
        if(!all(numVec)) stop("Must be list of numeric data")
        listLength<-sapply(listOfData,length)
        if(length(unique(listLength))!=1) stop("elements of 'listOfList' must be lists of same length")
        allData<-unlist(listOfData)
        x<-switch(join,row=matrix(allData,ncol=unique(listLength),byrow=TRUE),column=matrix(allData,nrow=unique(listLength),byrow=FALSE))
        if(join=="row") rownames(x)<-names(listOfData)
        if(join=="column") colnames(x)<-names(listOfData)
        return(x)
    }
}
