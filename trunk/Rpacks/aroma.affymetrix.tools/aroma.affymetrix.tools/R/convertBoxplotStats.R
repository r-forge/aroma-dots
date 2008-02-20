convertBoxplotStats<-function(boxplotStats){
    bxpStats <- list()
    bxpStats[["stats"]] <- do.call("cbind", lapply(boxplotStats, 
        function(x) {
            x[["stats"]]
        }))
    bxpStats[["conf"]] <- do.call("cbind", lapply(boxplotStats, 
        function(x) {
            x[["conf"]]
        }))
    bxpStats[["n"]] <- do.call("c", lapply(boxplotStats, function(x) {
        x[["n"]]
    }))
    bxpStats[["out"]] <- do.call("c", lapply(boxplotStats, function(x) {
        x[["out"]]
    }))
    igroup <- 0
    bxpStats[["group"]] <- do.call("c", lapply(boxplotStats, 
        function(x) {
            igroup <<- igroup + 1
            rep(igroup, length(x[["out"]]))
        }))
    bxpStats[["names"]]<-names(boxplotStats)
    return(bxpStats)
}
