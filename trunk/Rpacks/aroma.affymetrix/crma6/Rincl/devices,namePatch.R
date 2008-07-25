# These requires R.utils v1.0.4
dev.set <-
    function(which = dev.next())
{
    if (is.character(which))
      which <- .devIndexOf(which);
    which <- .Internal(dev.set(as.integer(which)))
    names(which) <- .Devices[[which]]
    which
}

dev.off <-
    function(which = dev.cur())
{
    if (is.character(which))
      which <- .devIndexOf(which);
    if(which == 1)
      stop("cannot shut down device 1 (the null device)")
    .Internal(dev.off(as.integer(which)))
    dev.cur()
} 


############################################################################
# HISTORY: 
# 2008-07-18
# o Created.
############################################################################
