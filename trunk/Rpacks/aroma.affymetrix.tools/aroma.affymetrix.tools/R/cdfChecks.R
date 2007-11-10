cdfChecks<-function(cdf){
    #cells point to multiple probesets/groups?
    out<-list()
    check<-1
    #check repeated cells:
    cells <- getCellIndices(cdf, unlist=F, useNames=T) #unit list; each unit element is $groups, then list with indices per group
    tcells <- table(unlist(cells,use.names=F))
    repsCore <- tcells[tcells > 1] 
    if(length(repsCore)>0) {
        cat(paste(check,". Possible Problem:",length(repsCore),"cells assigned to multiple probesets/units\n"))
        out<-c(out,list(repeatCells=as.numeric(names(repsCore))))}
    else{cat(paste(check,". Passed: Cells assigned to unique probesets/units\n"))}
    check<-check+1
    
    #check repeated group names:
    allGroups<-readCdfGroupNames(getPathname(cdf)) #list of groups by unit
    nperUnit<-sapply(allGroups,length) #how many groups per unit
    whichSingle<-which(nperUnit==1)
    whichBlank<-which(sapply(allGroups,function(x){length(x)==1 && x==""}))
    tgroups<-table(unlist(allGroups[-whichBlank],use.names=F))
    tgroups[tgroups>1]
    whichSingleNotBlank<-whichSingle[!(whichSingle %in% whichBlank)]
    whichSingleGpUnitEqual<-which(sapply(names(allGroups[whichSingleNotBlank]),function(unitName){ allGroups[[unitName]]==unitName})) #indices (allGroups) of single with groupname=unit name
    if(length(tgroups[tgroups>1])>0) {
        cat(paste(check,". Possible Problem:",length(tgroups[tgroups>1]),"groups with same name (not including groups with blank names)\n"))
        out<-c(out,list(repeatGroups=names(tgroups[tgroups>1])))}
    else{cat(paste(check,". Passed: Unique group names\n"))}
    check<-check+1
    
    #check unique unit names
    allUnits<-readCdfUnitNames(getPathname(cdf))
    tUnits<-table(allUnits)
    if(length(tUnits[tUnits>1])>0) {
        cat(paste(check,". Possible Problem:",length(tUnits[tUnits>1]),"units with same name\n"))
        out<-c(out,list(repeatUnits=names(tUnits[tUnits>1])))}
    else{cat(paste(check,". Passed: Unique unit names\n"))}
    check<-check+1
    
    #check if group name is also a unit name
    groupUnitOverlap<-unlist(allGroups,use.names=F)[which(unlist(allGroups,use.names=F) %in% allUnits)] #names of groups that also names of units
    if(length(groupUnitOverlap)>0) {
        cat(paste(check,". Possible Problem:",length(groupUnitOverlap),"group names are also unit names\n"))

        subsetOverlap<-sapply(groupUnitOverlap,function(name){sapply(allGroups[[name]],'%in%',x=name)}) #logical for groupUnitOverlap: whether the group in the same unit
        whichSubsetOverlap<-which(names(allGroups) %in% groupUnitOverlap[which(subsetOverlap)]) #indices (allGroups) of which unit names have group names that are the same
        cat("\t",length(whichSubsetOverlap),"are a group within a unit of the same name\n")

        whichUnitsOverlap<-which(names(allGroups) %in% groupUnitOverlap) #indices (allGroups) of which unit names are also names of groups
        cat(paste("\t",length(intersect(whichUnitsOverlap,whichSingle)),"are a unit with a single probeset\n"))
        
        cat(paste("\t",length(intersect(whichSubsetOverlap,whichSingleGpUnitEqual)),"are a unit with a single probeset, and the probeset name is the same as the unit name\n"))
        
        out<-c(out,list(overlapGroupUnit=groupUnitOverlap))}
    else{cat(paste(check,". Passed: Distinct unit and group names\n"))}
    check<-check+1

    #information about the nprobes per group
    nprobesPerUnit<-lapply(cells,function(x){sapply(x$groups,function(y){length(y$indices)})})
    whichBig<-which(sapply(nprobesPerUnit,function(x){any(x>4)}))
    if(length(whichBig)>0){
        cat(paste(check,". Possible Problem:",length(whichBig),"units have probesets with more than 4 probes\n"))
        cat("\tSummary of distribution of number of probes per probesets (for those with >4 probes):\n")
        print(summary(unlist(lapply(nprobesPerUnitFull[whichBigFull],function(x){x[x>4]}))))
        out<-c(out,unitsLongProbesets=names(nprobesPerUnit[whichBig]))
    }
    else cat(paste(check,". Passed: 4 or less probes per probeset\n"))
    check<-check+1

    cat("Summary Information:\n")
    cat(paste("Number of Units",length(cells),"\n"))
    cat(paste("Number of Groups",sum(nperUnit),"\n"))
    cat(paste("Number of Cells (probes)",length(unlist(cells)),"\n"))
    

    cat(paste("Number of Single Group Units:",length(whichSingle),"(",length(whichBlank),"have blank group names,",length(whichSingleGpUnitEqual),"have group name equal to unit name)\n"))
#    gc(cdf) #put this back
    invisible(out)

}
