#files -- names of files (with or without .cel) to search for and include in new project. 
#project Tag -- name to append for new project
#coreName -- beginning string to look for;if coreName=NULL, then search all subfolders under searchDir
#searchDir -- which directory to look for raw data
#targetDir -- where to put the folder of shortcuts that defines new project
#add -- if true, will continue even if target project folder already exists (with warning).

createProject<-function(files,projectTag,chipType,coreName=NULL,add=FALSE,searchDir="rawData",targetDir="probeData"){
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
    if(!is.null(coreName)){#only those directories that match the coreName
        allRawDir<-allRawDir[grep(paste("^",searchDir,"/",coreName,sep=""),allRawDir)] 
    }
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
}
