setConstructorS3("GdasAnnotationSet", function(...) {
  extend(AffymetrixFileSet(...), "GdasAnnotationSet"
  )
})

setMethodS3("fromFiles", "GdasAnnotationSet", function(static, path="annotations", ...) {
  path <- Arguments$getReadablePath(path, mustExist=TRUE);

  pattern <- "[.]tsv$";
  pathnames <- findFiles(pattern=pattern, path=path, firstOnly=FALSE);
  if (length(pathnames) == 0)
    throw("No files found: ", path);

  files <- lapply(pathnames, FUN=GdasAnnotationFile);

  newInstance(static, files);
}, static=TRUE)

setMethodS3("fromChipType", "GdasAnnotationSet", function(static, chipType, what, ...) {
  files <- vector("list", length(what));
  for (kk in seq(along=what)) {
    files[[kk]] <- GdasAnnotationFile$fromChipType(chipType=chipType, what=what[kk], ...);
  }

  newInstance(static, files);
}, static=TRUE)

setMethodS3("join", "GdasAnnotationSet", function(this, indexColumn=1, orderBy=NULL, ..., force=FALSE) {
  df <- this$.join;
  if (force || is.null(df)) {
    df <- NULL;
    files <- as.list(this);
    for (file in files) {
      data <- readData(file);
      if (is.null(df)) {
        df <- data;
      } else {
        df <- merge(df, data);
      }
    }
    this$.join <- df;
  }

  if (!is.null(orderBy)) {
    args <- vector("list", length(orderBy));
    for (kk in seq(along=args))
      args[[kk]] <- df[,orderBy[kk]];
    o <- do.call("order", args=args);
    df <- df[o,];
    rownames(df) <- 1:nrow(df);
  }
  this$.join <- df;

  df;
})


setMethodS3("[", "GdasAnnotationSet", function(this, i=NULL, j=NULL, drop=TRUE) {
  data <- join(this);
  names <- data[,1];
  data <- data[,-1];
  rownames(data) <- names;
  if (!is.null(i))
    data <- data[i,,drop=FALSE];
  if (!is.null(j))
    data <- data[,j,drop=FALSE];
  data;
})


############################################################################
# HISTORY:
# 2006-08-27
# o Created.
############################################################################
