createProject<-function(files,projectTag,chipType,coreName=NULL,add=FALSE,searchDir="rawData",targetDir="probeData"){
#if coreName=NULL, then search all subfolders under rawData
#doesn't search probeData; creates project in probeData
    if(.Platform$OS.type=="windows") stop("Sorry, cannot create symbolic links in windows")
    files<-gsub(".[Cc][Ee][Ll]$","",files)
    if(is.null(coreName)) projectFilename<-projectTag
    else projectFilename<-paste(coreName,projectTag,sep=",")
    projectPath<-file.path(targetDir,projectFilename,chipType)
    if(file.exists(projectPath)) {
        if(add) warning("Directory with name",projectFilename,"already exists; adding to the directory")
        else stop("Directory with name",projectFilename,"already exists")
    }
    else dir.create(projectPath,recursive=TRUE)
    allRawDir<-list.files(searchDir,full=TRUE)#gives full pathname rawData/... etc.
    allRawDir<-allRawDir[file.info(allRawDir)$isdir]
#    print(allRawDir)
    if(!is.null(coreName)){#only those directories that match the coreName
        allRawDir<-allRawDir[grep(paste("^",searchDir,"/",coreName,sep=""),allRawDir)] 
    }
#    print(allRawDir)
    allFiles<-list.files(file.path(allRawDir,chipType),full=TRUE)
    allbasename<-basename(allFiles)	
    CELIndex<-grep(".[Cc][Ee][Ll]$",allbasename)
    projectIndex<-match(files,gsub(".[Cc][Ee][Ll]$","",allbasename[CELIndex]))
    if(any(is.na(projectIndex))) warning("not all of files are found in the designated raw data folders")
    projectIndex<-projectIndex[!is.na(projectIndex)]
	projectFiles<-allFiles[CELIndex][projectIndex]
projectFiles<-file.path(getwd(),projectFiles)
    if(length(projectIndex)>0) file.symlink(from=projectFiles,to=projectPath)

	return()

#    for(kk in 1:length(allRawDir)){
#        dataDir<-file.path(allRawDir[kk],chipType)
#        CELfiles<-list.files(dataDir)
#        CELfiles<-CELfiles[grep(".[Cc][Ee][Ll]$",CELfiles)]     
#	if(length(CELfiles)>0){   
#        	projectFiles<-CELfiles[gsub(".[Cc][Ee][Ll]$","",CELfiles) %in% files]
#        	projectFiles<-file.path(getwd(),dataDir,projectFiles)
 #       	print(projectFiles)
#        	if(length(projectFiles)>0) file.symlink(from=projectFiles,to=projectPath)
#	}    
#    }
    
    
}
