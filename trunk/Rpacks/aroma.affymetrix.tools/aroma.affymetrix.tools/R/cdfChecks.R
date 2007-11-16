cdfChecks<-function(cdf){
    out<-list()
    check<-1
    
    #check repeated cells:
    cells <- getCellIndices(cdf, unlist=F, useNames=T) #unit list; each unit element is $groups, then list with indices per group
    tcells <- table(unlist(cells,use.names=F))
    repsCore <- tcells[tcells > 1] 
    if(length(repsCore)>0) {
        cat(paste(check,". Possible Problem:",length(repsCore),"cells assigned to multiple probesets/units\n"))
        out<-c(out,list(repeatCellIndices=as.numeric(names(repsCore))))}
    else{cat(paste(check,". Passed: Cells assigned to unique probesets/units\n"))}
    check<-check+1
    
    #check repeated group names:
    allGroups<-readCdfGroupNames(getPathname(cdf)) #list of groups by unit -- note is "" if group name is unit name? Might want to fix and just use info from cells?
    nperUnit<-sapply(allGroups,length) #how many groups per unit
    whichSingle<-which(nperUnit==1)
    whichSingleBlank<-which(sapply(allGroups,function(x){length(x)==1 && x==""}))#singleton and blank...
    if(length(whichSingleBlank)>0) tgroups<-table(unlist(allGroups[-whichSingleBlank],use.names=F)) #problem if length 0
    else  tgroups<-table(unlist(allGroups,use.names=F))
    tgroups[tgroups>1]
    whichSingleNotBlank<-whichSingle[!(whichSingle %in% whichSingleBlank)]
    if(length(tgroups[tgroups>1])>0) {
        cat(paste(check,". Possible Problem:",length(tgroups[tgroups>1]),"group names appear more than once (not including groups with blank names)\n"))
        out<-c(out,list(repeatGroupNames=names(tgroups[tgroups>1])))}
    else{cat(paste(check,". Passed: Unique group names\n"))}
    check<-check+1
    
    #check unique unit names
    allUnits<-readCdfUnitNames(getPathname(cdf))
    tUnits<-table(allUnits)
    if(length(tUnits[tUnits>1])>0) {
        cat(paste(check,". Possible Problem:",length(tUnits[tUnits>1]),"unit names appear more than once\n"))
        out<-c(out,list(repeatUnitNames=names(tUnits[tUnits>1])))}
    else{cat(paste(check,". Passed: Unique unit names\n"))}
    check<-check+1
    
    #check if group name is also a unit name
    groupUnitOverlap<-unlist(allGroups,use.names=F)[which(unlist(allGroups,use.names=F) %in% allUnits)] #names of groups that also names of units
    whichNotSingleBlank<-which(sapply(allGroups,function(x){length(x)>1 && x==""}))
    if(length(groupUnitOverlap)>0 || length(whichNotSingleBlank)>0) {
        cat(paste(check,". Possible Problem:",length(groupUnitOverlap),"unit names are same as a group name in another unit and",length(whichNotSingleBlank),"units with more than 1 group contain a group with the same name.\n"))

        subsetOverlap<-sapply(groupUnitOverlap,function(name){sapply(allGroups[[name]],'%in%',x=name)}) #logical for groupUnitOverlap: whether the group in the same unit
        whichSubsetOverlap<-which(names(allGroups) %in% groupUnitOverlap[which(subsetOverlap)]) #indices (allGroups) of which unit names have group names that are the same
        whichUnitsOverlap<-which(names(allGroups) %in% groupUnitOverlap) #indices (allGroups) of which unit names are also names of groups
        cat(paste("\t",length(intersect(whichUnitsOverlap,whichSingle)),"units that are also a group name somewhere else have only single probeset\n"))
        out<-c(out,list(overlapGroupUnitNames=c(groupUnitOverlap,names(allGroups)[whichNotSingleBlank])))
    }
    else{cat(paste(check,". Passed: Distinct unit and group names (not counting singleton probeset units)\n"))}
    check<-check+1

    #information about the nprobes per group
    nprobesPerUnit<-lapply(cells,function(x){sapply(x$groups,function(y){length(y$indices)})})
    whichBig<-which(sapply(nprobesPerUnit,function(x){any(x>4)}))
    if(length(whichBig)>0){
        cat(paste(check,". Possible Problem:",length(whichBig),"units have probesets with more than 4 probes\n"))
        cat("\tSummary of distribution of number of probes per probesets (for those with >4 probes):\n")
        print(summary(unlist(lapply(nprobesPerUnit[whichBig],function(x){x[x>4]}))))
        out<-c(out,unitsLongProbesets=names(nprobesPerUnit[whichBig]))
    }
    else cat(paste(check,". Passed: 4 or less probes per probeset\n"))
    check<-check+1

    cat("\nSummary Information:\n")
    cat(paste("Number of Units",length(cells),"\n"))
    cat(paste("Number of Groups",sum(nperUnit),"\n"))
    cat(paste("Number of Cells (probes)",length(unlist(cells)),"\n"))
    cat(paste("Number of Single Group Units:",length(whichSingle),"(",length(whichSingleBlank),"have group names same as unit)\n"))
    gc(cdf) #put this back
    invisible(out)

}
