#########################################################################/**
# @RdocClass Filter
#
# @title "The Filter class"
#
# \description{
#  @classhierarchy
#
#   A Filter object takes a vector of indices from its input(s) and passes
#   some or all of them through given some criteria. There are two main
#   classes of filters; (1) the @see "SerialFilter" class and (2) the
#   @see "ParallelFilter" class. Shortly, the \emph{serial} filters gets
#   some indices from one single source, whereas the \emph{parallel}
#   filters gets indices from several sources. Example of serial filters
#   are the @see "FieldFilter"s and the @see "NotFilter". Example of
#   parallel filters are the \emph{logical} filters @see "AndFilter" and
#   @see "OrFilter".
# }
#
# @synopsis
#
# \arguments{
#   \item{cex}{The scale factor of symbols that this filter highlights.}
#   \item{col}{The color of symbols that this filter highlights.}
#   \item{pch}{The plot symbols that this filter highlights.}
#   \item{visible}{If @TRUE, the data points filtered out by this filter
#     will be highlighted, otherwise not.}
# }
#
# \section{Fields and Methods}{
#  \bold{Fields}
#  \tabular{rll}{
#    \tab parameter \tab A list of parameters for the filter. \cr
#    \tab visible \tab Specifies if the data through this filter should be
#         displayed in plots etc or not. \cr
#  }
#
#  @allmethods "public"
# }
#
#
# @author
#
# \examples{
#    SMA$loadData("mouse.data")
#    layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#    raw <- RawData(mouse.data, layout=layout)
#
#    ma <- getSignal(raw)
#    normalizeWithinSlide(ma, method="s")
#    normalizeAcrossSlides(ma)
#
#    tma <- as.TMAData(ma)
#
#    # Filter out too weak spots, i.e. A < 8 (<256).
#    fA  <- AFilter(tma, range=c(8, Inf), visible=FALSE)
#
#    # Look at the top 5\% M-values.
#    fM  <- MFilter(tma, top=0.05, col="red")
#
#    # Look at the top 5\% T-values, but not those with to small SE.
#    fT  <- TFilter(tma, top=0.05, col="blue")
#    fNotSE <- SEFilter(tma, range=c(-Inf,0.02), col="yellow")
#    fSE <- NotFilter(fNotSE, visible=FALSE)
#
#    # Selects those spot that passes through all filters.
#    myFilter <- AndFilter(fA, fM, fT, fSE, col="purple", pch=19)
#
#    plot(tma, "TvsSE", xlog=2)
#    highlight(myFilter, recursive=TRUE)
# }
#
# \keyword{manip}
#*/#########################################################################
setConstructorS3("Filter", function(cex=NULL, col=NULL, pch=20, visible=TRUE) {
  extend(Object(), "Filter",
    .parameter = list(cex=cex, col=col, pch=pch),
    visible = visible
  )
}, abstract=TRUE);


setMethodS3("as.character", "Filter", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- paste(data.class(this), ": ", sep="");
  s <- paste(s, "visible=", this$visible, sep="");
  parameter <- novirtual(this)$.parameter;
  names     <- names(parameter);
  for (k in 1:length(parameter)) {
    s <- paste(s, ", ", names[k], "=", parameter[k], sep="");
  }
  s;
})



#########################################################################/**
# @RdocMethod getVisible
#
# @title "Checks if the filter is visible or not"
#
# @synopsis
#
# \description{
#   Checks if the filter is visible or not. A filter can be set to visible
#   or non-visible if the data for the indices outputted by the filter
#   should be highlighted, labeled etc. By default a filter is visible.
# }
#
# \value{
#   Returns @TRUE if the filter is visible, otherwise @FALSE.
# }
#
# @author
#
# \seealso{
#   @seemethod "setVisible".
#   @seeclass.
# }
#*/#########################################################################
setMethodS3("getVisible", "Filter", function(obj) {
  obj$visible;
})




#########################################################################/**
# @RdocMethod setVisible
#
# @title "Sets if the filter should bee visible or not"
#
# @synopsis
#
# \description{
#   Sets if the filter should bee visible or not.
#   By default a filter is visible.
# }
#
# \arguments{
#   \item{visible}{Logical value specifying if the filter should be 
#     visible or not.}
# }
#
# @author
#
# \seealso{
#   @seemethod "getVisible".
#   @seeclass.
# }
#*/#########################################################################
setMethodS3("setVisible", "Filter", function(obj, visible) {
  obj$visible <- visible;
})




#########################################################################/**
# @RdocMethod getParameter
#
# @title "Gets the values for a specific parameter for indices accepted by the filter"
#
# @synopsis
#
# \description{
#   Gets the values for a specific parameter for those indices that the 
#   filter lets through. If the filter's has this parameter set, e.g.
#   through @seemethod "setParameter", all indices accepted by the filter will
#   be having this parameter set to the filters value.
#   If the parameter is not set in the filter, or its value is @NULL
#   or \@NA, the filter asks all its input connections for assign
#   parameters to the indices. If more than one input filter tries to set
#   the parameter of an index, only the \emph{last} assignment will be used.
#   Hence, the order in which the input filters are specified might affect
#   the final assignment of paramters.
# }
#
# @author
#
# \seealso{
#   To set a parameter of a filter see @seemethod "setParameter".
#   To get the indices accepted by a filter see @seemethod "getIndex".
#   @seeclass.
# }
#*/#########################################################################
setMethodS3("getParameter", "Filter", function(this, paramName) {
  parameter <- this$.parameter[[paramName]];
  if (is.null(parameter) || is.na(parameter)) {
    # If the parameter is not set, ask all input filters...
    index <- getIndex(this);
    parameter <- rep(NA, length(index));
    inputs <- getInput(this);
    for (input in inputs) {
      if (inherits(input, "Filter")) {
        # Get the values of the parameter from the current input.
        input.parameter <- getParameter(input, paramName);
        if (!is.null(input.parameter) && !is.na(input.parameter)) {
          # Get the indices with this parameter set.
          valid.parameter <- !is.na(input.parameter);

          # Get the indices from the current input with the parameter set.
          input.index <- getIndex(input);

          input.index <- input.index[valid.parameter];
          input.parameter <- input.parameter[valid.parameter];

          # Among the indices check the ones accepted by this filter.
          map <- match(index, input.index);

          # Of course, some of the indices accepted by this filter, might
          # not be from the current input and therefore we don't look at
          # them either.
          valid.map <- !is.na(map);

          # Is there anything left to look at?
          if (sum(valid.map) > 0) {
            # Make sure the length of the parameter value vector is long enough.
            input.parameter <- rep(input.parameter, length.out=length(valid.map));

            # Assignment...
            parameter[valid.map] <- input.parameter[valid.map];
          }
        }
      } # if (inherits(input, "Filter"))
    } # for (input in inputs)

    # If the parameter only contains one unique value use that to represent
    # the whole vector. Maybe, this is unnecessary. /HB 2001-07-18
    unique.parameter <- na.omit(unique(parameter));
    unique.parameter.length <- length(unique.parameter);
    if (unique.parameter.length == 0)
      parameter <- NA
    else if (unique.parameter.length == 1)
      parameter <- unique.parameter;
  } # if (is.null(parameter) || is.na(parameter))
  parameter;
})


#########################################################################/**
# @RdocMethod setParameter
#
# @title "Sets the values for a specific parameter for indices accepted by the filter"
#
# @synopsis
#
# \description{
#   Sets the values for a specific parameter for those indices that the 
#   filter lets through. For more information see @seemethod "getParameter".
# }
#
# @author
#
# \seealso{
#   To get a parameter of a filter see @seemethod "getParameter".
#   To get the indices accepted by a filter see @seemethod "getIndex".
#   @seeclass.
# }
#*/#########################################################################
setMethodS3("setParameter", "Filter", function(this, paramName, value) {
  if (is.null(this$.parameter))
    this$.parameter <- list();
  this$.parameter[[paramName]] <- value;
})



#########################################################################/**
# @RdocMethod getIndex
#
# @title "Gets the indices accepted by the filter"
#
# @synopsis
#
# \description{
#   Gets the indices that this filter passes through. The filter asks
#   all its input connections for indices and passes through only those
#   indices that are accepted by this filter.
# }
#
# @author
#
# \seealso{
#   @seeclass.
# }
#*/#########################################################################
setMethodS3("getIndex", "Filter", abstract=TRUE);




#########################################################################/**
# @RdocMethod highlight
#
# @title "Highlights the data points accepted by the filter"
#
# @synopsis
#
# \description{
#   Highlights the data points accepted by the filter and optionally its
#   input objects.
#   If the argument \code{recursive == TRUE} the filter will \emph{first}
#   ask its input filters to highlight the data accepted by them. If these
#   filters have filters connected to them they will before doing anything
#   asking there input filters and so on. The last filter to highlight the
#   data will therefore be this filter.
# }
#
# \arguments{
#   \item{cex}{The scaling factor for \emph{all} points, i.e. for
#     \emph{all} indices. If @NULL, the global default value will be
#     used. If \code{"filter"} the scaling factor set by the filter will
#     be used.}
#   \item{col}{The color for \emph{all} points, i.e. for
#     \emph{all} indices. If @NULL, the global default value will be
#     used. If \code{"filter"} the color set by the filter will
#     be used.}
#   \item{pch}{The point style for \emph{all} points, i.e. for
#     \emph{all} indices. If @NULL, the global default value will be
#     used. If \code{"filter"} the point style set by the filter will
#     be used.}
#   \item{...}{Any parameters accepted by most standard plot functions.}
#   \item{recursive}{If @TRUE this filter and all filters up the
#     stream will highlight the data points that passes through them.
#     If @FALSE, only this filter will be used to highlight the data.}
# }
#
# @author
#
# \seealso{
#   To label data points see @seemethod "text".
#   @seeclass.
# }
#*/#########################################################################
setMethodS3("highlight", "Filter", function(this, cex="filter", col="filter", pch="filter", ..., recursive=FALSE) {
  if (recursive) {
    for (input in getInput(this)) {
      if (inherits(input, "Filter"))
        highlight(input, cex=cex, col=col, pch=pch, ..., recursive=recursive);
    }
  }
  if (this$visible) {
    lastPlot <- Device$getPlotParameters();
    object <- lastPlot$object;
    incl <- getIndex(this);

    param <- list();
    for (paramName in c("cex", "col", "pch")) {
      arg <- get(paramName);
      if (length(arg) > 1) {
        arg <- arg[incl];
      } else if (!is.null(arg) && arg == "filter") {
        arg <- getParameter(this, paramName);
      }
      if (is.function(arg)) {
        arg <- arg(incl, ...);
      }

      if (length(arg) == 1 && is.na(arg)) arg <- NULL;
#      cat(paramName, "=", arg, "\n");
      param[[paramName]] <- arg;
    }

    cex <- param$cex;
    col <- param$col;
    pch <- param$pch;

    highlight(object, include=incl, cex=cex, col=col, pch=pch, ...);
  }
})




#########################################################################/**
# @RdocMethod text
#
# @title "Labels the data points accepted by the filter"
#
# @synopsis
#
# \description{
#   Labels the data points accepted by the filter and optionally its
#   input objects.
#   If the argument \code{recursive == TRUE} the filter will \emph{first}
#   ask its input filters to label the data accepted by them. If these
#   filters have filters connected to them they will before doing anything
#   asking there input filters and so on. The last filter to label the
#   data will therefore be this filter.
# }
#
# \arguments{
#   \item{labels}{The labels for \emph{all} points, i.e. for \emph{all}
#     indices. If @NULL, the indices will be used as labels.}
#   \item{cex}{The scaling factor for \emph{all} points, i.e. for
#     \emph{all} indices. If @NULL, the global default value will be
#     used. If \code{"filter"} the scaling factor set by the filter will
#     be used.}
#   \item{col}{The color for \emph{all} points, i.e. for
#     \emph{all} indices. If @NULL, the global default value will be
#     used. If \code{"filter"} the color set by the filter will
#     be used.}
#   \item{...}{Any parameters accepted by standard \code{\link[base]{text}}.}
#   \item{recursive}{If @TRUE this filter and all filters up the
#     stream will label the data points that passes through them.
#     If @FALSE, only this filter will be used to label the data.}
# }
#
# @author
#
# \seealso{
#   To highlight data points see @seemethod "highlight".
#   @seeclass.
# }
#*/#########################################################################
setMethodS3("text", "Filter", function(x, labels=NULL, cex="filter", col=NULL, ..., recursive=FALSE) {
  # To please R CMD check...
  this <- x;

  if (recursive) {
    for (input in getInput(this)) {
      if (inherits(input, "Filter"))
        text(input, labels=labels, cex=cex, col=col, ..., recursive=recursive);
    }
  }
  if (this$visible) {
    lastPlot <- Device$getPlotParameters();
    object <- lastPlot$object;
    incl <- getIndex(this);

    param <- list();
    for (paramName in c("cex", "col")) {
      arg <- get(paramName);
      if (length(arg) > 1) {
        arg <- arg[incl];
      } else if (!is.null(arg) && arg == "filter") {
        arg <- getParameter(this, paramName);
      }
      if (is.function(arg)) {
        arg <- arg(incl, ...);
      }

      if (length(arg) == 1 && is.na(arg)) arg <- NULL;
#      cat(paramName, "=", arg, "\n");
      param[[paramName]] <- arg;
    }

    cex <- param$cex;
    col <- param$col;

    if (length(labels) > 1) {
      labels <- labels[incl];
    } else if (is.function(labels)) {
      labels <- labels(incl, ...);
    }

    text(object, include=incl, labels=labels, cex=cex, col=col, ...);
  }
})




#########################################################################/**
# @RdocMethod getInput
#
# @title "Gets all the input objects connected to the filter"
#
# @synopsis
#
# \description{
#   Gets all the input objects connected to the filter.
# }
#
# @author
#
# \seealso{
#   To change one or more input objects see @seemethod "changeInput".
#   @seeclass.
# }
#*/#########################################################################
setMethodS3("getInput", "Filter", abstract=TRUE);





#########################################################################/**
# @RdocMethod changeInput
#
# @title "Change input(s) on this filter and optionally all filters down the stream"
#
# @synopsis
#
# \description{
#   Changes some input(s) to new input(s). A filter has one or more 
#   input connections that each can be connected to for instance  another
#   filter's output or to some data objects. After having designed a
#   network of filters connected to some data objects it is sometime 
#   desireble to connect the same network of filters to another set of
#   data objects. This can be done easily by calling this method on the
#   very last filter and tell it to recursively update all other filters
#   accordingly.
# }
#
# \arguments{
#   \item{newInput}{All connections to this input will be
#     disconnected and replaced by the input object specified by 
#     this argument.}
#   \item{oldInput}{All connections matching to this input will be replaced.
#     If @NULL all connections will be replaced.}
#   \item{recursive}{If @TRUE this filter and all filters up the
#     stream will update its connections according to \code{oldInput} and
#    \code{newInput}. If @FALSE, only this filter will update its
#    input connection(s).}
# }
#
# @author
#
# \seealso{
#   To the current set of connected inputs see @seemethod "getInput".
#   @seeclass.
# }
#*/#########################################################################
setMethodS3("changeInput", "Filter", abstract=TRUE);




############################################################################
# HISTORY:
# 2003-04-21
# o Updated the Rdocs.
# 2002-02-26
# * Updated code to make use of setMethodS3's.
# 2001-09-28
# * Changed the order of the newInput and oldInput arg's in changeInput().
# 2001-07-18
# * Totally removed the specific fields 'cex', 'col' and 'pch' and replaced
#   the with a general 'parameter' list, which can hold any kind of 
#   parameter. This approach is way more general than getCex(), getCol(), 
#   and getPch()!
# * Improved getCol() to color each single data point. If another filter
#   previously has colored an data point, the following filter can recolor
#   it or leave it as is (by not specifying and pch). If two or more input
#   filters colors the same data point, the last input filter will decide
#   the color. For instance, the order in which filters are added to a
#   ParallelFilter decides which filter "rules".
# * Added some Rdoc comments.
# * Started to use Device$getPlotParameters() which removed the need of
#   an argument specifying the last plotting MicroarrayData object in
#   methods such as highlight and text.
# 2001-07-13
# * Rearrange the class structure to contain Serial- and ParallelFilters.
# * Added some Rdoc comments.
# * Extended. Added plot paramters, recursive plotting etc.
# 2001-07-12
# * Created! Eventually I would like to have a SpeedGroupFilter, 
#   EisenFilter etc.
############################################################################
