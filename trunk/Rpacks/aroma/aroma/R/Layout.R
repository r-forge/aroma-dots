#########################################################################/**
# @RdocClass Layout
#
# @title "The Layout class"
#
# \description{
#  @classhierarchy
#
#  The Layout class describes the layout of a microarray slide, 
#  such as the number of spots and the layout of the grids etc.
#  A microarray slide is layout in a number of grids, where each grid
#  is refered to by its row and its column. Within each grid, each spot
#  is refered to by its row and its column \emph{within that grid}. For
#  example, a microarray slide with 4*3 grids and where each grid has 
#  12*10 spots, has in total 4*3*12*10 = 12*120 = 1440 spots. The spot
#  in the top left corner is located at grid (1,1) and spot (1,1) and
#  we say it has the \emph{location} (1,1,1,1). In the above example, 
#  the spot in the lower right corner of the slide has the location 
#  (4,3,12,10). Further more, each spot on a microarray slide has a
#  unique \emph{index}. This index starts counting starts with grid 1,
#  moves along row 1 from column 1 until
#  the last column, then advances to the next row. After all rows in
#  grid 1 are indexed, counting proceeds to grid 2 and so on. This way of
#  indexing the spots is used by for instance GenePix, Spot and ScanAlyze.
#  Continuing the example above, the spot at location (1,1,1,1) has index
#  1, the spot at location (4,3,12,10) has index 1440. The spot at location
#  (2,2,11,3) has index ((2-1)*3+(2-1))*12*10+(11-1)*10+3 = 483. 
#  To get the \emph{index} from a \emph{location} one do
#  \code{getIndex(layout, 2,2,11,3)} where 
#  \code{layout <- Layout(4,3,12,10)}. To get the \emph{location} from a
#  \emph{index} one do \code{getLocation(layout, 483)}, which
#  gives the vector (2, 2, 11, 3). To get the \emph{position}, i.e. the
#  overall row and column of a spot on the microarray slide one can do
#  \code{getPosition(483)} which give the position (23,13).\cr
#
#  All spots are refered to by their unique \emph{indices}.
# }
#
# @synopsis
#
# \arguments{
#   \item{nspot.r}{The number of rows of spots per grid.
#         This first argument can also be a list containing
#         the fields nspot.r, nspot.c, ngrid.r, and ngrid.c.}
#   \item{nspot.c}{The number of columns of spots per grid.}
#   \item{ngrid.r}{The number of rows of grids per slide.}
#   \item{ngrid.c}{The number of columns of grids per slide.}
#   \item{name}{A @vector of names for the spot. \bold{optional}.}
#   \item{id}{A @vector of ids for the spot. \bold{optional}.}
#   \item{printorder}{A @matrix or a @character string specifying the
#         order that spots where printed. For more details 
#         @seemethod "setPrintorder".}
#   \item{geneSpotMap}{.}
#   \item{plate}{.}
# }
#
# \section{Fields and Methods}{
#  \bold{Fields}
#  \tabular{rll}{
#   \tab \code{ngrid.r} \tab The number of rows of grids per slide. \cr
#   \tab \code{ngrid.c} \tab The number of columns of grids per slide. \cr
#   \tab \code{nspot.r} \tab The number of rows of spots in each grid. \cr
#   \tab \code{nspot.c} \tab The number of columns of spots in each grid. \cr
#   \tab \code{name}    \tab A @vector of strings which specified the name of each gene. \cr
#   \tab \code{id}      \tab A @vector of strings which specified the id of each gene. \cr
#  }
#
#  @allmethods "public"
# }
#
# \note{
#  There are several functions that returns a Layout object.
#  It is only in very special cases that you have to create one yourself.
#
#  In the sma package some functions are related to this class. This,
#  class might be backward compatible with these functions, but the reverse
#  is not true. The following functions are known be related to this class:
#  @see "sma::init.grid", @see "sma::id2image" (and \code{image2id}).
# }
#
# \details{
#   GenBank Accession numbers, SwissProt/TrEMBL Accession numbers or 
#   Entry Names.
# }
#
# @author
#
# \examples{
#   layout <- Layout(ngrid.r=4, ngrid.c=4, nspot.r=18, nspot.c=18)
#
#   SMA$loadData("mouse.setup")
#   layout <- Layout(mouse.setup)
#   # or, equivalent...
#   layout <- as.Layout(mouse.setup)
# }
#
# \seealso{
#  @see "MicroarrayData.getLayout".
# }
#*/#########################################################################
setConstructorS3("Layout", function(ngrid.r=0, ngrid.c=0, nspot.r=0, nspot.c=0, geneSpotMap=NULL, name=NULL, id=NULL, plate=NULL, printorder=NULL) {
  # This is to support old sma style too.
  if (is.list(ngrid.r)) {
    if ( !any(is.na(match(c("ngrid.r", "ngrid.c", "nspot.r", "nspot.c"), 
                                                       names(ngrid.r)))) ) {
      l <- ngrid.r;
      nspot.r <- l$nspot.r;
      nspot.c <- l$nspot.c; 
      ngrid.r <- l$ngrid.r; 
      ngrid.c <- l$ngrid.c;
    }
  }

  this <- extend(Object(), "Layout", 
    nspot.r      = nspot.r, 
    nspot.c      = nspot.c, 
    ngrid.r      = ngrid.r, 
    ngrid.c      = ngrid.c,
    id           = id,
    .name        = name,
    plate        = plate,
    .printorder  = NULL,
    geneSpotMap  = geneSpotMap,
    replicates   = NULL,
    plateGrps    = NULL,
    printdipGrps = NULL,
    printtipGrps = NULL,
    geneGrps     = NULL
  );

  if (inherits(ngrid.r, "Layout"))
    set(this, ngrid.r);
  if (!is.null(printorder))
    setPrintorder(this, printorder);

  this;
});


setMethodS3("as.character", "Layout", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- paste(sep="", data.class(this), ": ",
    "Grids: ", this$ngrid.r, "x", this$ngrid.c, " (=", 
                                          nbrOfGrids(this), "), ",
    "spots in grids: ", this$nspot.r, "x", this$nspot.c, " (=", 
                                          gridSize(this), "), ",
    "total number of spots: ", nbrOfSpots(this), "."
  );
  if (hasNames(this))
    s <- paste(sep="", s, " Spot names are specified.")
  if (hasIds(this))
    s <- paste(sep="", s, " Spot ids are specified.")
  if (!is.null(this$geneGrps))
    s <- paste(sep="", s, " ", this$geneGrps)
  if (hasPlates(this))
    s <- paste(sep="", s, " Spot plates are specified.")

  s;
})



#########################################################################/**
# @RdocMethod getId
#
# @title "Gets the id of one or more spots"
#
# @synopsis
#
# \arguments{
#   \item{index}{A @vector of indices indicating which ids to set. If 
#    @NULL, all ids are set.}
# }
#
# \value{
#   Returns the @vector of ids.
# }
#
# \description{
#  @get "title" given their indices.
# }
#
# @author
#
# \examples{
#   SMA$loadData(c("mouse.setup", "mouse.gnames"))
#   layout <- as.Layout(mouse.setup, id=mouse.gnames)
#
#   # Get the id of spot # 2453 and 2412:2417.
#   getId(layout, c(2453, 2412:2417))
# }
#
# \seealso{
#   @seemethod "indexOf",
#   @seemethod "setId".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getId", "Layout", function(this, index=NULL, ...) {
  if (is.null(index))
    as.character(this$id)
  else {
    n <- nbrOfSpots(this);   # Modulo n just in case
    index <- ((index-1) %% n)+1;
    as.character(this$id[index]);
  }
})


setMethodS3("hasIds", "Layout", function(this, ...) {
  (length(this$id) > 0)
})


setMethodS3("hasNames", "Layout", function(this, ...) {
  (length(this$.name) > 0)
})

setMethodS3("hasPlates", "Layout", function(this, ...) {
  (length(this$plate) > 0)
})


#########################################################################/**
# @RdocMethod setId
#
# @title "Sets the id of one or more spots"
#
# @synopsis
#
# \arguments{
#   \item{id}{A @vector of ids.}
#   \item{index}{A @vector of indices indicating which ids to set. If 
#    @NULL, all ids are set.}
# }
#
# \description{
#  @get "title" given their indices.
# }
#
# @author
#
# \examples{
#   SMA$loadData(c("mouse.setup", "mouse.gnames"))
#   layout <- as.Layout(mouse.setup, id=mouse.gnames)
#
#   setId(layout, c("2412r", "2414r"), c(2412, 2414))
#
#   # Get the id of spot # 2453
#   getId(layout, c(2453, 2412:2417))
# }
#
# \seealso{
#   @seemethod "getId".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("setId", "Layout", function(this, id, index=NULL, ...) {
  if (is.null(index))
    this$id <- as.character(id)
  else {
    n <- nbrOfSpots(this);   # Modulo n just in case
    index <- ((index-1) %% n)+1;
    this$id[index] <- as.character(id);
  }
})




#########################################################################/**
# @RdocMethod getName
#
# @title "Gets the name of one or more spots"
#
# @synopsis
#
# \arguments{
#   \item{index}{A @vector of indices indicating which names to set. If 
#    @NULL, all names are set.}
# }
#
# \value{
#   Returns the @vector of names.
# }
#
# \description{
#  @get "title" given their indices.
# }
#
# @author
#
# \examples{
#   SMA$loadData(c("mouse.setup", "mouse.gnames"))
#   layout <- as.Layout(mouse.setup, name=mouse.gnames)
#
#   # Get the name of spot # 2453 and 2412:2417.
#   getName(layout, c(2453, 2412:2417))
# }
#
# \seealso{
#   @seemethod "indexOf",
#   @seemethod "setName".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getName", "Layout", function(this, index=NULL, ...) {
  if (is.null(index))
    as.character(this$.name)
  else {
    n <- nbrOfSpots(this);   # Modulo n just in case
    index <- ((index-1) %% n)+1;
    as.character(this$.name[index]);
  }
})


#########################################################################/**
# @RdocMethod setName
#
# @title "Sets the name of one or more spots"
#
# @synopsis
#
# \arguments{
#   \item{name}{A @vector of names.}
#   \item{index}{A @vector of indices indicating which names to set. If 
#    @NULL, all names are set.}
# }
#
# \description{
#  @get "title" given their indices.
# }
#
# @author
#
# \examples{
#   SMA$loadData(c("mouse.setup", "mouse.gnames"))
#   layout <- as.Layout(mouse.setup, name=mouse.gnames)
#
#   setName(layout, c("2412r", "2414r"), c(2412, 2414))
#
#   # Get the name of spot # 2453
#   getName(layout, c(2453, 2412:2417))
# }
#
# \seealso{
#   @seemethod "getName".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("setName", "Layout", function(this, name, index=NULL, ...) {
  if (is.null(index))
    this$.name <- as.character(name)
  else {
    n <- nbrOfSpots(this);   # Modulo n just in case
    index <- ((index-1) %% n)+1;
    this$.name[index] <- as.character(name);
  }
})


setMethodS3("rename", "Layout", function(this, from, to, what="name", ...) {
  from <- as.character(from);
  to   <- as.character(to);
  value <- this[[what]];
  if (is.null(value))
    throw("Unknown value of argument 'what': ", what);
  this[[what]] <- gsub(from, to, value);
}) # rename()



setMethodS3("getPlate", "Layout", function(this, index=NULL, ...) {
  if (is.null(index))
    as.character(this$plate)
  else {
    n <- nbrOfSpots(this);   # Modulo n just in case
    index <- ((index-1) %% n)+1;
    as.character(this$plate[index]);
  }
})


setMethodS3("getPlateNumber", "Layout", function(this, index=NULL, ...) {
  plate <- getPlate(this, index=index);
  plates <- unique(plate);
  match(plate, plates);
})


setMethodS3("nbrOfPlates", "Layout", function(this, ...) {
  plate <- getPlate(this);
  if (is.null(plate)) return(NULL);
  length(unique(plate));
})


setMethodS3("setPlate", "Layout", function(this, plate, index=NULL, ...) {
  if (is.null(index))
    this$plate <- as.character(plate)
  else {
    n <- nbrOfSpots(this);   # Modulo n just in case
    index <- ((index-1) %% n)+1;
    this$plate[index] <- as.character(plate);
  }
})


#########################################################################/**
# @RdocMethod getPrintorder
#
# @title "Gets the order of the spots in which they were printed"
#
# @synopsis
#
# \arguments{
#   \item{value}{The elements to be placed in the resulting @matrix.}
# }
#
# \description{
#   @get "title". The spots in column one were printed first, then the 
#   spots in column two and so on. 
#   By default the spot indices are returned.
# }
#
# \value{
#   Returns a @matrix with spot (values) printed at the same time in the
#   same column. The first spots printed are in column one and the last 
#   ones printed in the last column. 
#   Often this means that there are \code{nbrOfGrids(layout)} rows and
#   \code{gridSize(layout)} columns in the matrix. However, if the 
#   slide was printed in say two halfs (first half of the grids are
#   printed and then the second), then this is not true.
# }
#
# \section{Print order}{
#   The printing of a microarray is time consuming and often several
#   microarray slides are printed at the same time, since it is even more
#   time consuming to switch between the trays. When printing several
#   microarrays at the same time, the arrayer prints the first spot in
#   all grids on \emph{all} slides, before moving on to the second spot.
#   For a example, printing a batch of 100 slides with 6384 spots in 4x4
#   grids takes about 15 hours to print including manual work to switch
#   trays etc. Each grid contains 19*21 spots, i.e. the arrayer has to put
#   down the print tips 399 times on each slide, and in total 39900 times.
#   This is about 44 put-downs a minute. It takes about 45-50 minutes
#   to finish one row of spots.
# }
#
# \section{Different directions}{
#   The most common print-order directions are \code{"row-by-row"} and
#   \code{"column-by-column"}. 
#   In both cases, when printing a slide at each print step
#   \code{nbrOfGrids(layout)} spots are printed at the same time. 
#   The arrayer start of spotting the first spot in \emph{each} of the
#   grids. Then it cleans the print-tip heads, dries them, and go back
#   to the trays to get a \emph{new} set of cDNA and prints the second
#   spot in each of the grids. The second spot is to the right to
#   (\code{"row-by-row"}) or below (\code{"column-by-column"}) the
#   first spot. When the array gets to the end of a row (column) it
#   moves on to print the next row (column) and so on until all in
#   all grids have been printed.
# }
#
# \section{Print-order effects}{
#   An important factor for the quality of the printed spots is the
#   temperature and the humidity. Too high temperature and humidity
#   tends to produce too large spots that can even overlap [1]. If there
#   is no automatic control for temperature and humidity, the quality of
#   the spots could vary a lot between the spots printed during a 15 hours
#   printing process. With a varying printing climate we should expect to
#   see a variating of the quality of the spots along the order of which
#   the spots are printed.
#   The variation of temperature and humidity varies approximately in the
#   time scale of hours. As it takes about 45-50 minutes to print a row
#   of spots, we should therefore expect to see such a variation between
#   the rows in the grids.
# }
#
# \references{
#   [1] Microarrays in Three Easy Steps, Priti Hedge, The Institute for
#       Genomic Research, 200?.
# }
#
# \examples{
#   layout <- Layout(2,2, 3,3)
#
#   # No printorder specified - assumes de facto standard "row-by-row"
#   print(getPrintorder(layout))
#
#   #      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#   # [1,]    1    2    3    4    5    6    7    8    9
#   # [2,]   10   11   12   13   14   15   16   17   18
#   # [3,]   19   20   21   22   23   24   25   26   27
#   # [4,]   28   29   30   31   32   33   34   35   36
#
#   # Spots (1,10,19,28) were printed first, then (2,11,20,29), ...
#
#   setPrintorder(layout, "column-by-column")
#   print(getPrintorder(layout))
#
#   #       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#   # [1,]    1    4    7    2    5    8    3    6    9
#   # [2,]   10   13   16   11   14   17   12   15   18
#   # [3,]   19   22   25   20   23   26   21   24   27
#   # [4,]   28   31   34   29   32   35   30   33   36
#
#   # Spots (1,10,19,28) were printed first, then (4,13,22,31) below, ...
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getPrintorder", "Layout", function(this, value=1:nbrOfSpots(this), ...) {
  printorder <- this$.printorder;
  if (is.null(printorder)) {
    warning("Printorder not specified. Assume row-by-row printing.");
    index <- 1:nbrOfSpots(this);
    printorder <- matrix(index, nrow=nbrOfGrids(this), 
                              ncol=gridSize(this), byrow=TRUE);
  }
  res <- value[printorder];
  dim(res) <- dim(printorder);
  res;
})


setMethodS3("setPrintorder", "Layout", function(this, printorder=c("row-by-row","column-by-column"), ...) {
  if (is.numeric(printorder)) {
    if (length(printorder) != nbrOfSpots(this)) {
      throw("Vector 'printorder' should contain ", nbrOfSpots(this), 
                                       " elements: ", length(printorder));
    }
    printorder <- as.matrix(printorder);
  } else if (is.character(printorder)) {
    printorder <- match.arg(printorder);
    if (printorder == "row-by-row") {
      index <- 1:nbrOfSpots(this);
      printorder <- matrix(index, nrow=nbrOfGrids(this), 
                                ncol=gridSize(this), byrow=TRUE);
    } else if (printorder == "column-by-column") {
      index <- 1:nbrOfSpots(this);
      printorder <- matrix(index, nrow=nbrOfGrids(this), 
                                ncol=gridSize(this), byrow=TRUE);
      perm <- as.vector(matrix(1:gridSize(this), 
                        nrow=this$nspot.r, ncol=this$nspot.c, byrow=TRUE));
      printorder <- printorder[,perm];
    }
  } else if (!is.null(printorder)) {
    throw("Unknown value of 'printorder': ", printorder);
  } 
  this$.printorder <- printorder;
})



#########################################################################/**
# @RdocMethod indexOf
#
# @title "Gets the index of one or more spots from their name or id"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \arguments{
#   \item{name}{A string @vector (or a single regular expression) to match
#      against the names.}
#   \item{id}{A string @vector (or a single regular expression) to match
#      against the ids.}
#   \item{plate}{A string @vector (or a single regular expression) to match
#      against the plates. If it is @numeric it will be matched against the
#      plate number as given by \code{getPlateNumber()}.}
#   \item{ignoreCase}{If @TRUE, the matching is not case sensitive,
#      otherwise it is.}
#   \item{regexpr}{If @TRUE, regular expression matching is used,
#      otherwise plain string comparison is used.}
#   At least one of the arguments \code{name} and \code{id} must be given.
#   If both are given, indices of spots that match \emph{either} the 
#   \code{name} search pattern \emph{or} the \code{id} search pattern.
# }
#
# \value{
#   Returns the @vector of indices of matched names or ids. Returns
#   \code{numeric(0)} if no matches were found.
# }
#
# @author
#
# \examples{
#   SMA$loadData(c("mouse.setup", "mouse.gnames"))
#   layout <- as.Layout(mouse.setup, id=mouse.gnames)
#
#   # Get the index of spots with id "54" and "232".
#   indexOf(layout, id=c("54", "232"))
#   # [1] 54 232
#
#   # Get the index of all spots with id beginning with "120" and
#   # having at least four characters.
#   indexOf(layout, id="^120.+", regexpr=TRUE)
#   # [1] 1200 1201 1202 1203 1204 1205 1206 1207 1208 1209
# }
#
# \seealso{
#   For more help on regular expressions see
#   @see "base::grep" and @see "base::apropos".
#   @seemethod "getName".
#   @seemethod "getId".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("indexOf", "Layout", function(this, name=NULL, id=name, plate=NULL, ignoreCase=TRUE, regexpr=FALSE, ...) {
  if (is.null(name) && is.null(id) && is.null(plate))
    throw("At least one of the arguments 'name', 'id' and 'plate' must be specified.");

  # Do only search if there actually are anything to search for.
  name.src <- getName(this);
  if (is.null(name.src))
    name <- NULL;
  
  id.src <- getId(this);
  if (is.null(id.src))
    id <- NULL;
  
  if (!is.numeric(plate))
    plate.src <- getPlate(this)
  else
    plate.src <- getPlateNumber(this);
  if (is.null(plate.src))
    plate <- NULL;

  if (is.null(name) && is.null(id) && is.null(plate))
    return(numeric(0));

  if (ignoreCase == TRUE) {
    if (!is.null(name))      name <- tolower(name);
    if (!is.null(id))        id <- tolower(id);
    if (!is.null(name.src))  name.src <- tolower(name.src);
    if (!is.null(id.src))    id.src <- tolower(id.src);
    if (!is.numeric(plate)) {
      if (!is.null(plate))     plate <- tolower(plate);
      if (!is.null(plate.src)) plate.src <- tolower(plate.src);
    }
  }
  
  # Search for names?
  if (!is.null(name)) {
    if (regexpr == TRUE)
      names <- which(regexpr(name, name.src) != -1)
    else
      names <- which(!is.na(match(name.src, name)))
  }
  
  # Search for ids?
  if (!is.null(id)) {
    if (regexpr == TRUE)
      ids <- which(regexpr(id, id.src) != -1)
    else
      ids <- which(!is.na(match(id.src, id)))
  }

  # Search for plates?
  if (!is.null(plate)) {
    if (regexpr == TRUE)
      plates <- which(regexpr(plate, plate.src) != -1)
    else
      plates <- which(!is.na(match(plate.src, plate)))
  }

  res <- NULL;
  if (!is.null(name))
    res <- union(res, names);
  if (!is.null(id))
    res <- union(res, ids);
  if (!is.null(plate))
    res <- union(res, plates);
  res;
})

setMethodS3("as.Layout", "Layout", function(this, ...) {
  this;
})

setMethodS3("as.Layout", "ANY", function(object, ...) {
  Layout(object, ...)
})


#########################################################################/**
# @RdocMethod getIndex
#
# @title "Gets the index of a spot given its location"
#
# @synopsis
#
# \description{
#  @get "title". The location can either be
#  a vector containing the grid row and the grid column and the spot row
#  and the spot column in that grid, or it can be the same fields as 
#  seperate arguments.
# }
#
# @author
#
# \examples{
#   # Example 1
#   layout <- Layout(4,4, 18,18)
#   idx <- getIndex(layout, 2, 3, 4, 3)     # 2001
#   idx <- getIndex(layout, c(2, 3, 4, 3))  # 2001 (equivalent)
#   loc <- getLocation(layout, idx)         # 2 3 4 3
#
#   # Example 2
#   SMA$loadData(c("mouse.data", "mouse.setup"))
#   raw <- RawData(mouse.data, layout=as.Layout(mouse.setup))
#   ma <- getSignal(raw)
#   layout <- getLayout(ma)
#
#   plotSpatial(ma)
#
#   # Highlights spot number 2462
#   idx <- 2462
#   highlight(ma, idx, col="purple")
#
#   # Highlights the spot at grid (2,3) and its spot (4,3)
#   idx <- getIndex(layout, 2, 3, 4, 3);  # Spot #2460
#   highlight(ma, idx, col="purple")
#
#   # Highlights all spots in grid (1,2)
#   idx <- getIndices(layout, 1,2, NULL,NULL)
#   highlight(ma, idx, col="purple")
# }
#
# \seealso{
#   This method corresponds to image2id (see @see "sma::id2image") in
#   the sma package.
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getIndex", "Layout", function(this, gridRow, gridColumn=NULL, spotRow=NULL, spotColumn=NULL, ...) {
  if (length(gridRow) >= 4) {
    location   <- gridRow;
    gridRow    <- location[1];
    gridColumn <- location[2];
    spotRow    <- location[3];
    spotColumn <- location[4];
  }
  
  if (gridRow < 1 || gridRow > this$ngrid.r)
    throw("gridRow is out of range: ", gridRow);
  if (gridColumn < 1 || gridColumn > this$ngrid.c)
    throw("gridColumn is out of range: ", gridColumn);
  if (spotRow < 1 || spotRow > this$nspot.r)
    throw("spotRow is out of range: ", spotRow);
  if (spotColumn < 1 || spotColumn > this$nspot.c)
    throw("spotColumn is out of range:", spotColumn);

  nbrOfSpotsPerGrid <- this$nspot.r*this$nspot.c;
  return ( ((gridRow-1)*this$ngrid.c+(gridColumn-1))*nbrOfSpotsPerGrid +
            (spotRow-1)*this$nspot.c+(spotColumn-1) + 1);
})


#########################################################################/**
# @RdocMethod getIndices
#
# @title "Gets the indices of the spots at the given locations"
#
# @synopsis
#
# \arguments{
#  \item{gridRows}{The grid rows to be included. If @NULL all are included.}
#  \item{gridColumns}{The grid column to be included. If @NULL all are included.}
#  \item{spotRows}{The spot rows to be included. If @NULL all are included.}
#  \item{spotColumns}{The spot column to be included. If @NULL all are included.}
# }
#
# \description{
#  @get "title".
#  This method is very useful for instance when one would like to look for
#  spatial effects such as edge effects etc.
# }
#
# @author
#
# \examples{
#   SMA$loadData(c("mouse.data", "mouse.setup"))
#   raw <- RawData(mouse.data, layout=as.Layout(mouse.setup))
#   ma <- getSignal(raw)
#   layout <- getLayout(ma)
#
#   plotSpatial(ma)
#
#   # Highlights all spots in grid (1,2)
#   idx <- getIndices(layout, gridRows=1, gridColumns=2);
#   highlight(ma, idx, col="purple")
#
#   # Highlight all spots in column 1, 6 and 12 in 
#   # grid (3,2), (3,3), (4,2) and (4,3).
#   idx <- getIndices(layout, gridRows=3:4, gridColumns=2:3, spotColumns=c(1,6,12));
#   highlight(ma, idx, col="orange")
#
#   # Highlight all "alley spots" of each printtip group, i.e. those spots that
#   # do *not* have eight neighbors and do a background "alley" next to them.
#   alley <- getIndices(layout, spotRows=c(1,layout$nspot.r))
#   alley <- union(alley, getIndices(layout, spotColumns=c(1,layout$nspot.c)));
#   highlight(ma, alley, col="pink")
# }
#*/#########################################################################
setMethodS3("getIndices", "Layout", function(this, gridRows=NULL, gridColumns=NULL, spotRows=NULL, spotColumns=NULL, ...) {
  if (is.null(gridRows))
    gridRows <- 1:this$ngrid.r;
  if (is.null(gridColumns))
    gridColumns <- 1:this$ngrid.c;
  if (is.null(spotRows))
    spotRows <- 1:this$nspot.r;
  if (is.null(spotColumns))
    spotColumns <- 1:this$nspot.c;

  if (any(gridRows < 1 || gridRows > this$ngrid.r))
    throw("gridRows is out of range.");
  if (any(gridColumns < 1 || gridColumns > this$ngrid.c))
    throw("gridColumns is out of range.");
  if (any(spotRows < 1 || spotRows > this$nspot.r))
    throw("spotRows is out of range.");
  if (any(spotColumns < 1 || spotColumns > this$nspot.c))
    throw("spotColumns is out of range.");

  indices <- c();
  for(gridRow in gridRows) {
    for(gridColumn in gridColumns) {
      for(spotRow in spotRows) {
        for(spotColumn in spotColumns) {
          indices <- c(indices, 
                   this$getIndex(gridRow, gridColumn, spotRow, spotColumn));
        }
      }
    }
  }
  
  sort(indices);
})




#########################################################################/**
# @RdocMethod getLocation
#
# @title "Gets the location of a spot given its index"
#
# @synopsis
#
# \arguments{
#  \item{index}{The spot index of one or many spots to be found. All values
#   much be within a valid range otherwise an exception is thrown.}
# }
#
# \description{
#  @get "title". Returns a @vector containing
#  the grid row and the grid column and the spot row and the spot column in
#  that grid.
# }
#
# @author
#
# \examples{
#    layout <- Layout(4,4, 18,18)
#    loc <- getLocation(layout, 2001)      # 2 3 4 3
#    idx <- getIndex(layout, loc)          # 2001
# }
#
# \seealso{
#   This method corresponds to @see "sma::id2image" in the sma
#   package.
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getLocation", "Layout", function(this, index=NULL, ...) {
  if (is.null(index)) {
    index <- 1:nbrOfSpots(this);
  } else if (any(index < 1))
    throw("index is out of range.");
  
  n <- nbrOfSpots(this);   # Modulo n just in case
  index <- ((index-1) %% n)+1;

  nbrOfSpotsPerGrid <- this$nspot.r*this$nspot.c;
  nbrOfSpotsPerGridRow <- this$ngrid.c * nbrOfSpotsPerGrid;

  locations <- c();
  for (idx in index) {
    gridIndex <- (idx-1) %/% nbrOfSpotsPerGrid + 1;
    gridOffset <- (gridIndex-1) * nbrOfSpotsPerGrid;
    gridRow <- (idx-1) %/% nbrOfSpotsPerGridRow + 1;
    gridColumn <- (gridIndex-1) %% this$ngrid.c + 1;

    spotOffset <- idx - gridOffset;
    spotRow <- (spotOffset-1) %/% this$nspot.c + 1;  
    spotColumn <- (spotOffset-1) %% this$nspot.c + 1;  
    
    location <- c(gridRow, gridColumn, spotRow, spotColumn);
    locations <- c(locations, location);
  }

  matrix(locations, ncol=4, byrow=TRUE);
})


#########################################################################/**
# @RdocMethod getPosition
#
# @title "Gets the position of a set of spots given their indices"
#
# @synopsis
#
# \arguments{
#  \item{index}{The spot index of one or many spots to be found. All values
#   much be within a valid range otherwise an exception is thrown.
#   If @NULL, all spots are considered.}
# }
#
# \description{
#  @get "title". The position of
#  a spot is the pair (row, column) where row is the row in the microarray
#  and column is the column in the microarray where the spot is positioned.
# }
#
# \value{
#  Returns a "named" @matrix where each row contains the row and the 
#  column of the spots.
# }
#
# @author
#
# \examples{
#    layout <- Layout(4,4, 18,18)
#    xy <- getPosition(layout, 2001)           # 22 39
#    xy <- getPosition(layout, c(2001,2002))   # 22 39; 22 40
# }
#
# \seealso{
#   @seemethod "getIndex", @seemethod "getIndices", 
#   @seemethod "getLocation".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getPosition", "Layout", function(this, index=NULL, ...) {
  if (is.null(index)) {
    index <- 1:nbrOfSpots(this);
  } else if (any(index < 1))
    throw("index is out of range.");

  n <- nbrOfSpots(this);
  index <- ((index-1) %% n)+1;                     # Modulo n just in case
 
  # Index map for the whole layout. Independent of the argument 'index'!
  idxmap <- toXYMatrix(this, 1:n);
  rows <- c();  cols <- c();
  for (kk in 1:nbrOfRows(this)) rows[idxmap[kk,]] <- kk;
  for (kk in 1:nbrOfColumns(this)) cols[idxmap[,kk]] <- kk;

  rownames <- index;
  colnames <- c("row", "column");
  matrix(c(rows[index],cols[index]), ncol=2, byrow=FALSE, 
                                   dimnames=list(rownames, colnames));
})



#########################################################################/**
# @RdocMethod equals
#
# @title "Checks if a Layout object is equals to some other object"
#
# @synopsis
#
# \arguments{
#  \item{obj}{The other object for which this object should be compared to.}
# }
#
# \description{
#  Checks if a \code{Layout} object is equal to some other object, which
#  normally is another \code{Layout} object. It could however also be a
#  list with the same fields as a \code{Layout} object.
#  \cr
#  A layout is equal to another layout if it has 1) the same number of
#  rows and columns of grids as the other object, 2) the same number of
#  rows and columns of spots as the other object.
# }
#
# @author
#
# \examples{
#    layout1 <- Layout(4,4, 18,18)
#    layout2 <- Layout(4,4, 18,18)
#    layout3 <- Layout(4,4, 16,16)
#    equals(layout1, layout1)   # TRUE (of course)
#    equals(layout1, layout2)   # TRUE
#    equals(layout2, layout1)   # TRUE (symmetric)
#    equals(layout1, layout3)   # FALSE
#    equals(layout2, layout3)   # FALSE
# }
#*/#########################################################################
setMethodS3("equals", "Layout", function(this, obj, ...) {
  if (is.null(obj))
    return (FALSE);

  if (!inherits(obj, "Object") && !is.list(obj))
    return (FALSE);

  if (this$nspot.r != obj$nspot.r) return(FALSE);
  if (this$nspot.c != obj$nspot.c) return(FALSE);
  if (this$ngrid.r != obj$ngrid.r) return(FALSE);
  if (this$ngrid.c != obj$ngrid.c) return(FALSE);

  TRUE;
})
            
            

#########################################################################/**
# @RdocMethod set
#
# @title "Sets the layout"
#
# @synopsis
#
# \arguments{
#   \item{nspot.r}{The number of rows of spots per grid.}
#   \item{nspot.c}{The number of columns of spots per grid.}
#   \item{ngrid.r}{The number of rows of grids per slide.}
#   \item{ngrid.c}{The number of columns of grids per slide.}
#   \item{name}{A @vector if names for the spot.}
#   \item{id}{A @vector if ids for the spot.}
# }
#
# \description{
#  Sets the layout by either 1) explicitly setting the number of rows and
#  columns of grids and the number of rows and columns of spots within each
#  grid, or 2) by giving another \code{Layout} object.
# }
#
# @author
#
# \examples{
#    layout1 <- Layout(4,4, 18,18)
#    layout2 <- Layout()
#    set(layout2, 4,4, 16,16)   # Alternative 1
#    set(layout2, layout1)      # Alternative 2
# }
#*/#########################################################################
setMethodS3("set", "Layout", function(this, ngrid.r=NULL, ngrid.c=NULL, 
                   nspot.r=NULL, nspot.c=NULL, name=NULL, id=NULL, ...) {
  if (inherits(ngrid.r, "Layout")) {
    other <- ngrid.r;
    this$nspot.r <- other$nspot.r;
    this$nspot.c <- other$nspot.c;
    this$ngrid.r <- other$ngrid.r;
    this$ngrid.c <- other$ngrid.c;
    this$id      <- other$id;
  } else {
    if (!is.null(nspot.r)) this$nspot.r <- nspot.r;
    if (!is.null(nspot.c)) this$nspot.c <- nspot.c;
    if (!is.null(ngrid.r)) this$ngrid.r <- ngrid.r;
    if (!is.null(ngrid.c)) this$ngrid.c <- ngrid.c;
    if (!is.null(name))    this$.name   <- name;
    if (!is.null(id))      this$id      <- id;
  }
  invisible(this);
})


#########################################################################/**
# @RdocMethod nbrOfSpots
#
# @title "Gets the size of a microarray"
#
# \description{
#  Calculates the total number of spots on the microarray slide.
# }
#
# @synopsis
#
# @author
#
# \examples{
#   layout <- Layout(ngrid.r=4, ngrid.c=4, nspot.r=18, nspot.c=18)
#   print(nbrOfSpots(layout))  # 5184
# }
#*/#########################################################################
setMethodS3("nbrOfSpots", "Layout", function(this, ...) {
  this$ngrid.r * this$ngrid.c * this$nspot.r * this$nspot.c;
})



#########################################################################/**
# @RdocMethod nbrOfGrids
#
# @title "Gets the number of grids on a microarray"
#
# \description{
#  Calculates the number of grids on the microarray slide.
# }
#
# @synopsis
#
# @author
#
# \examples{
#   layout <- Layout(ngrid.r=4, ngrid.c=4, nspot.r=18, nspot.c=18)
#   print(nbrOfGrids(layout))  # 16
# }
#*/#########################################################################
setMethodS3("nbrOfGrids", "Layout", function(this, ...) {
  this$ngrid.r * this$ngrid.c;
})




#########################################################################/**
# @RdocMethod nbrOfRows
#
# @title "Gets the number of rows on a microarray"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# @author
#
# \examples{
#   layout <- Layout(ngrid.r=4, ngrid.c=4, nspot.r=18, nspot.c=18)
#   print(nbrOfRows(layout))  # 72 = 4*18
# }
#*/#########################################################################
setMethodS3("nbrOfRows", "Layout", function(this, ...) {
  this$ngrid.r * this$nspot.r;
})





#########################################################################/**
# @RdocMethod nbrOfColumns
#
# @title "Gets the number of columns on a microarray"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# @author
#
# \examples{
#   layout <- Layout(ngrid.r=4, ngrid.c=4, nspot.r=18, nspot.c=18)
#   print(nbrOfColumns(layout))  # 72 = 4*18
# }
#*/#########################################################################
setMethodS3("nbrOfColumns", "Layout", function(this, ...) {
  this$ngrid.c * this$nspot.c;
})





#########################################################################/**
# @RdocMethod size
#
# @title "Gets the size of a microarray"
#
# \description{
#  Calculates the total number of spots on the microarray slide.
# }
#
# @synopsis
#
# @author
#
# \examples{
#   layout <- Layout(ngrid.r=4, ngrid.c=4, nspot.r=18, nspot.c=18)
#   print(size(layout))  # 5184
# }
#
# \seealso{
#   @seemethod "nbrOfSpots".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("size", "Layout", function(this, ...) {
  nbrOfSpots(this);
})





#########################################################################/**
# @RdocMethod gridSize
#
# @title "Gets the number of spots in each grid"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# @author
#
# \examples{
#   layout <- Layout(ngrid.r=4, ngrid.c=4, nspot.r=18, nspot.c=18)
#   print(gridSize(layout))   # 324 = 18*18
# }
#*/#########################################################################
setMethodS3("gridSize", "Layout", function(this, ...) {
  this$nspot.r * this$nspot.c;
})




#########################################################################/**
# @RdocMethod toXYMatrix
#
# @title "Layouts out values in the same way as the spots are layout"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{values}{The values to be put into the matrix.}
#   \item{flip}{If @TRUE the resulting matrix is flipped
#      vertically before returned.}
# }
#
# @author
#
# \examples{
#   layout <- Layout(3,2, 2,3)
#
#   # Print the indices of the each of the spots as they appear
#   # on the slide(s).
#   print(toXYMatrix(layout))
#
#   #      [,1] [,2] [,3] [,4] [,5] [,6]
#   # [1,]    1    2    3    7    8    9
#   # [2,]    4    5    6   10   11   12
#   # [3,]   13   14   15   19   20   21
#   # [4,]   16   17   18   22   23   24
#   # [5,]   25   26   27   31   32   33
#   # [6,]   28   29   30   34   35   36
# }
#*/#########################################################################
setMethodS3("toXYMatrix", "Layout", function(this, values=seq(nbrOfSpots(this)), flip=FALSE, ...) {
  gc <- this$ngrid.c;
  gr <- this$ngrid.r;
  sc <- this$nspot.c;
  sr <- this$nspot.r;
  
  n <- nbrOfSpots(this);
  nbrOfGrids <- nbrOfGrids(this);
  gridSize <- gridSize(this);

  # Split the values gridwise
  gridIdx <- rep(1:(nbrOfGrids), rep(gridSize, nbrOfGrids));
  gridMs  <- split(values, gridIdx);

  # For each grid, create a matrix. Result is put into a list
  imgGrids <- lapply(gridMs, FUN=matrix, nrow=sr, ncol=sc, byrow=TRUE);

  # Number of spots per grid row
  nbrOfSpotsPerGridRow <- gridSize*gc;

  # Grid column each spot is located in.
  gridRowIdx <- rep(1:gr, rep(nbrOfSpotsPerGridRow, gr));

  # The grid rows, row by row.
  gridRows <- split(unlist(imgGrids), gridRowIdx);

  # Turn the the each row of grids into grid matrices
  grids <- lapply(gridRows, FUN=matrix, nrow=sr);

  img <- NULL;
  # For each grid row
  for (i in 1:gr)
   img <- rbind(img, grids[[i]]);

  if (flip) {
    # Flip image vertically. Image size = rows x columns
    img <- t(apply(img, MARGIN=2, FUN=rev)); 
  }

  img;
})




############################################################################
############################################################################
## 
##  PLOTTING & GRAPHICAL METHODS
## 
############################################################################
############################################################################

setMethodS3("put", "Layout", function(this, x,y, width="20%", height=width, border="black", labels=NULL, offset=c(0.5,0.5), col=NULL, lty=NULL, ...) {

  height <- height;
  if (is.null(offset))
    throw("Argument 'offset' can't be NULL.");
  offset <- rep(offset, length.out=2);

  plot.area <- par("usr");

  if (is.character(width)) {
    tmp <- as.integer(unlist(strsplit(width, "%"))[1]);
    width <- (plot.area[2]-plot.area[1])*tmp/100;
  } else if (is.null(width)) {
    throw("Argument 'width' can not be NULL.");
  }
  if (is.character(height)) {
    tmp <- as.integer(unlist(strsplit(height, "%"))[1]);
    height <- (plot.area[4]-plot.area[3])*tmp/100;
  } else if (is.null(height)) {
    throw("Argument 'height' can not be NULL.");
  }

  if (is.character(x)) {
    tmp <- as.integer(unlist(strsplit(x, "%"))[1]);
    x <- plot.area[1]+(plot.area[2]-plot.area[1]-width)*tmp/100;
  } else if (is.null(x)) {
    throw("Argument 'x' can not be NULL.");
  }
  if (is.character(y)) {
    tmp <- as.integer(unlist(strsplit(y, "%"))[1]);
    y <- plot.area[3]+(plot.area[4]-plot.area[3]-height)*tmp/100;
  } else if (is.null(y)) {
    throw("Argument 'y' can not be NULL.");
  }

  x0 <- x; x1 <- x+width;
  y0 <- y; y1 <- y+height;

  nrow <- this$ngrid.r;
  ncol <- this$ngrid.c;
  dx <- width/ncol;
  dy <- height/nrow;

  # Plot the grid layout
  for (row in 0:nrow) {
    yy <- y0+row*dy;
    lines(x=c(x0,x1), y=c(yy,yy), col=border, ...);
  }
  for (column in 0:ncol) {
    xx <- x0+column*dx;
    lines(x=c(xx,xx), y=c(y0,y1), col=border, ...);
  }

  if (is.null(labels)) labels <- 1:(ncol*nrow);

  if (!is.null(col) && col=="auto")
    col <- rainbow(ncol)
  else if (is.null(col)) 
    col <- par("col");
  col <- rep(col, length=ncol*nrow);

  if (!is.null(lty) && lty=="auto")
    lty <- matrix(1:nrow, nrow=nrow, ncol=ncol, byrow=TRUE);
  lty <- rep(lty, length.out=ncol*nrow);

  grids <- matrix(1:nbrOfGrids(this), ncol=ncol, byrow=TRUE);
  for (row in 1:nrow) {
    rect.y0 <- y1-(row-1+0.1)*dy;
    rect.y1 <- y1-(row-1+0.9)*dy;
    yy <- y1-(row-1+offset[2])*dy;
    for (column in 1:ncol) {
      rect.x0 <- x0+(column-1+0.1)*dx;
      rect.x1 <- x0+(column-1+0.9)*dx;
      xx <- x0+(column-1+offset[1])*dx;
      idx <- grids[row,column];
      text(xx,yy, labels=labels[idx], col=col[idx], ...);
      if (!is.null(lty[idx]))
        rect(rect.x0,rect.y0,rect.x1,rect.y1, border=col[idx], lty=lty[idx], ...)
    }
  }
})



#########################################################################/**
# @RdocMethod toPrintorderMatrix
#
# @title "Gets a matrix of spot indices in the order they were printed"
#
# @synopsis
#
# \arguments{
#   \item{value}{The elements to be placed in the resulting matrix.}
# }
#
# \description{
#   @get "title". The
#   spots in column one were printed first, then the spots in column two
#   and so on. By default the spot indices are returned.
# }
#
# \value{
#   Returns a @matrix with \code{nbrOfGrids(layout)} rows and
#   \code{gridSize(layout)} columns.
# }
#
# \details{
#   When printing a slide at each print step \code{nbrOfGrids(layout)}
#   spots are printed at the same time. The arrayer start of spotting the
#   first spot in \emph{each} of the grids. The it cleans the print-tip
#   heads, dries them, and go back to the trays to get a \emph{new} set of
#   cDNA and prints the second spot in each of the grids. The second spot
#   is to the right of the first spot. When the array gets to the end of
#   a row it has printed the first line of each grid. The arrayer now moves
#   on to print the spots in line two, three and so on until all lines in
#   all grids have been printed.\cr
#
#   The printing of a microarray is time consuming and often several
#   microarray slides are printed at the same time, since it is even more
#   time consuming to switch between the trays. When printing several
#   microarrays at the same time, the arrayer prints the first spot in
#   all grids on \emph{all} slides, before moving on to the second spot.
#   For a example, printing a batch of 100 slides with 6384 spots in 4x4
#   grids takes about 15 hours to print including manual work to switch
#   trays etc. Each grid contains 19*21 spots, i.e. the arrayer has to put
#   down the print tips 399 times on each slide, and in total 39900 times.
#   This is about 44 put-downs a minute. It takes about 45-50 minutes
#   to finish one row of spots.
#
#   An important factor for the quality of the printed spots is the
#   temperature and the humidity. Too high temperature and humidity
#   tends to produce too large spots that can even overlap [1]. If there
#   is no automatic control for temperature and humidity, the quality of
#   the spots could vary a lot between the spots printed during a 15 hours
#   printing process. With a varying printing climate we should expect to
#   see a variating of the quality of the spots along the order of which
#   the spots are printed.
#   The variation of temperature and humidity varies approximately in the
#   time scale of hours. As it takes about 45-50 minutes to print a row
#   of spots, we should therefore expect to see such a variation between
#   the rows in the grids.
# }
#
# \references{
#   [1] Microarrays in Three Easy Steps, Priti Hedge, The Institute for
#       Genomic Research, 200?.
# }
#
# \examples{
#   layout <- Layout(2,2, 3,3)
#
#   print(toPrintorderMatrix(layout))
#
#   #      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#   # [1,]    1    2    3    4    5    6    7    8    9
#   # [2,]   10   11   12   13   14   15   16   17   18
#   # [3,]   19   20   21   22   23   24   25   26   27
#   # [4,]   28   29   30   31   32   33   34   35   36
#
#   # Spot 1, 10, 19 and 28 were printed first, then 2, 11, 20, and 29...
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("toPrintorderMatrix", "Layout", function(this, value=seq(nbrOfSpots(this)), ...) {
  getPrintorder(this, value=value);
}, private=TRUE)





setMethodS3("as.data.frame", "Layout", function(x, ...) {
  # To please R CMD check...
  this <- x;

  grow <- matrix(1:this$ngrid.r, nrow=this$ngrid.r, ncol=this$ngrid.c*gridSize(this), byrow=FALSE);
  grow <- as.vector(t(grow));
  gcol <- matrix(1:this$ngrid.c, nrow=gridSize(this), ncol=this$ngrid.c*this$ngrid.r, byrow=TRUE);
  gcol <- as.vector(gcol);
  row <- matrix(1:this$nspot.r, nrow=this$nspot.r, ncol=this$nspot.c, byrow=FALSE);
  row <- as.vector(t(row));
  row <- rep(row, times=this$ngrid.r*this$ngrid.c);
  col <- matrix(1:this$nspot.c, nrow=this$nspot.c, ncol=this$nspot.r, byrow=FALSE);
  col <- as.vector(col);
  col <- rep(col, times=this$ngrid.r*this$ngrid.c);
  df <- data.frame(gridRow=grow, gridColumn=gcol, spotRow=row, spotColumn=col);

  fields <- getFields(this, private=TRUE);
  fields <- sub("^[.]", "", fields);
#  fields <- c("name", "id", "plate", "acc", "clid", "type");
  for (field in fields) {
    value <- this[[field]];
    if (inherits(value, "AsIs"))
      value <- unclass(value);
    if (is.vector(value) && length(value) == nrow(df)) {
      names <- c(names(df), field);
      df <- cbind(df, I(value));
      names(df) <- names;
    }
  }
  df;  
});


setMethodS3("fromDataFrame", "Layout", function(this, df, ...) {
  if (!is.data.frame(df))
    throw("Argument 'df' is not a data frame.");
  header <- c("gridRow", "gridColumn", "spotRow", "spotColumn");
  if (!all(is.element(header, names(df)))) {
    headerStr <- paste("'", header, "'", sep="", collapse=", ");
    throw("The data frame must contain all of the fields ", headerStr, ": ", paste(names(df), collapse=", "));
  }
  ngrid.r <- max(df[[header[1]]]);
  ngrid.c <- max(df[[header[2]]]);
  nspot.r <- max(df[[header[3]]]);
  nspot.c <- max(df[[header[4]]]);
  id      <- as.vector(df[["id"]]);
  name    <- as.vector(df[["name"]]);
  plate   <- as.vector(df[["plate"]]);
  layout <- Layout(ngrid.r, ngrid.c, nspot.r, nspot.c, id=id, name=name, plate=plate)
  for (name in setdiff(colnames(df), c("gridRow", "gridColumn", "spotRow", "spotColumn", "id", "name", "plate")))
    layout[[name]] <- as.vector(df[[name]]);
  layout;
}, static=TRUE);



#########################################################################/**
# @RdocMethod read
#
# @title "Reads layout information from a tab-delimited file"
#
# \description{
#  Static method that reads layout information from a tab-delimited file.
# }
#
# @synopsis
#
# \arguments{
#   \item{filename}{The name of the file.}
#   \item{path}{Optional path where the data should be written.}
#   \item{sep}{Separator @character between columns.}
#   \item{header}{If @TRUE column headers are written, otherwise not.}
#   \item{...}{Other arguments accepted by \code{read.table()}.}
# }
#
# @author
#
# \examples{
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   print(layout)
# }
#
# \seealso{
#   For writing a Layout object to a file @seemethod "write".
#   See also \code{read.table()}.
#   @seeclass
# }
#*/#########################################################################
setMethodS3("read", "Layout", function(this, filename, path=NULL, header=TRUE, sep="\t", quote="", ...) {
  filename <- Arguments$getReadablePathname(filename, path);  

  # Support gzip'ed files too.
  if (regexpr("[.]gz$", filename) != -1) {
    tmpname <- tempfile();
    n <- gunzip(filename, tmpname);
    filename <- tmpname;
    on.exit(file.remove(tmpname));
  } 

  df <- read.table(file=filename, sep=sep, header=header, quote=quote, ...);
  for (k in seq(along=length(df))) {
  }

  Layout$fromDataFrame(df);
}, static=TRUE); 




#########################################################################/**
# @RdocMethod write
#
# @title "Writes the layout information to a tab-delimited file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{filename}{The name of the file.}
#   \item{path}{Optional path where the data should be written.}
#   \item{row.names}{If @TRUE row names are written, otherwise not.}
#   \item{sep}{Separator @character between columns.}
#   \item{...}{Other arguments accepted by \code{write.table()}.}
# }
#
# @author
#
# \examples{
#   layout <- Layout(3,2, 2,3)
#   tmpfile <- tempfile()
#   write(layout, tmpfile)
#   file.show(tmpfile)
#   unlink(tmpfile)
# }
#
# \seealso{
#   For read a Layout object from file see @seemethod "read".
#   See also \code{write.table()}.
#   @seeclass
# }
#*/#########################################################################
setMethodS3("write", "Layout", function(this, filename, path=NULL, overwrite=FALSE, row.names=FALSE, sep="\t", quote=FALSE, ...) {
  filename <- Arguments$getWritablePathname(filename, path, mustNotExist=!overwrite);  
  
  df <- as.data.frame(this);
  write.table(df, file=filename, row.names=row.names, sep=sep, quote=quote, ...);
});



setMethodS3("getLayoutGroupsByName", "Layout", function(this, name, ...) {
  if (!is.character(name) && length(name) != 1)
    throw("Argument 'name' must be a single character string.");

  methodName <- paste("get", capitalize(name), "Groups", sep="");

  # Subclasses of Layout that we should look in.
  classes <- class(this);
  last <- which(classes == "Layout");
  classes <- classes[1:last];

  methods <- paste(methodName, classes, sep=".");

  for (method in methods) {
    if (exists(method, mode="function")) {
      fcn <- get(method, mode="function");
      layoutGroup <- fcn(this);
      if (!inherits(layoutGroup, "LayoutGroups"))
        throw("Method did not return a LayoutGroups object as expected: ", data.class(layoutGroup));
      return(layoutGroup);
    }
  }

  throw("No method corresponding to LayoutGroup '", name, "' exists: ", methodName);
})


setMethodS3("getGeneGroups", "Layout", function(object, ...) {
  # Cache it
  if (is.null(object$geneGrps))
    object$geneGrps <- GeneGroups(object)
  object$geneGrps;
})

setMethodS3("setGeneGroups", "Layout", function(object, genes, ...) {
  if (is.null(genes) || inherits(genes, "GeneGroups"))
    object$geneGrps <- genes;
})

setMethodS3("getPlateGroups", "Layout", function(object, ...) {
  # Cache it
  if (is.null(object$plateGrps))
    object$plateGrps <- PlateGroups(object)
  object$plateGrps;
})

setMethodS3("setPlateGroups", "Layout", function(object, plates, ...) {
  if (is.null(plates) || inherits(plates, "PlateGroups"))
    object$plateGrps <- plates;
})

setMethodS3("getPrintdipGroups", "Layout", function(object, ...) {
  # Cache it
  if (is.null(object$printdipGrps))
    object$printdipGrps <- PrintdipGroups(object)
  object$printdipGrps;
})

setMethodS3("setPrintdipGroups", "Layout", function(object, printdips, ...) {
  if (is.null(printdips) || inherits(printdips, "PrintdipGroups"))
    object$printdipGrps <- printdips;
})


setMethodS3("getPrinttipGroups", "Layout", function(object, ...) {
  # Cache it
  if (is.null(object$printtipGrps))
    object$printtipGrps <- PrinttipGroups(object)
  object$printtipGrps;
})

setMethodS3("setPrinttipGroups", "Layout", function(object, printtips, ...) {
  if (is.null(printtips) || inherits(printtips, "PrinttipGroups"))
    object$printtipGrps <- printtips;
})


setMethodS3("getSlideRowGroups", "Layout", function(object, ...) {
  # Cache it
  if (is.null(object$slideRowGrps))
    object$slideRowGrps <- SlideRowGroups(object)
  object$slideRowGrps;
})

setMethodS3("setSlideRowGroups", "Layout", function(object, slideRows, ...) {
  if (is.null(slideRows) || inherits(slideRows, "SlideRowGroups"))
    object$slideRowGrps <- slideRows;
})


setMethodS3("getSlideColumnGroups", "Layout", function(object, ...) {
  # Cache it
  if (is.null(object$slideColumnGrps))
    object$slideColumnGrps <- SlideColumnGroups(object)
  object$slideColumnGrps;
})

setMethodS3("setSlideColumnGroups", "Layout", function(object, slideColumns, ...) {
  if (is.null(slideColumns) || inherits(slideColumns, "SlideColumnGroups"))
    object$slideColumnGrps <- slideColumns;
})



setMethodS3("getBlank", "Layout", function(object, blanks="blank|empty", ...) {
  result <- c();
  count <- 0;
  s <- getName(object);
  if (length(s) > 0) {
    s <- tolower(s);
    result <- (regexpr(blanks, s) != -1);
    count <- count + 1;
  }
  s <- getId(object);
  if (length(s) > 0) {
    s <- tolower(s);
    result <- (regexpr(blanks, s) != -1);
    count <- count + 1;
  }
  if (count == 0)
    warning("No blanks were found because no names nor ids are defined.");
  result
})


setMethodS3("anonymize", "Layout", function(object, method=c("numerate", "randomize", "reshuffle"), blanks="blank|empty", ...) {
  BLANK <- "BLANK";

  method <- tolower(method[1]);
  if (!is.element(method, c("numerate", "randomize", "reshuffle")))
    throw("Unknown anonymizer method: ", method);
  
  n <- nbrOfSpots(object);
  # Keep the blanks...
  blanks <- getBlank(object, blanks=blanks);
#  blanks <- rep(FALSE, n);
  if (hasNames(object)) {
    s <- getName(object);
    s <- tolower(s);
    uniques <- unique(s);
    nu <- length(uniques);
    match <- match(s, uniques);
    if (method == "numerate") {
      map <- 1:nu;
    } else if (method == "randomize") {
      map <- sample(1:nu, nu);
    } else if (method == "reshuffle") {
      map <- sample(uniques, nu);
    }
    ws <- (nchar(s) == 0);
    s[!blanks & !ws] <- map[match][!blanks & !ws];
    s[blanks] <- BLANK;
    setName(object, s);
  }
  if (hasIds(object)) {
    s <- getId(object);
    s <- tolower(s);
    uniques <- unique(s);
    nu <- length(uniques);
    match <- match(s, uniques);
    if (method == "numerate") {
      map <- 1:nu;
    } else if (method == "randomize") {
      map <- sample(1:nu, nu);
    } else if (method == "reshuffle") {
      map <- sample(uniques, nu);
    }
    ws <- (nchar(s) == 0);
    s[!blanks & !ws] <- map[match][!blanks & !ws];
    s[blanks] <- BLANK;
    setId(object, s);
  }
  if (hasPlates(object)) {
    s <- getPlate(object);
    s <- tolower(s);
    uniques <- unique(s);
    nu <- length(uniques);
    match <- match(s, uniques);
    if (method == "numerate") {
      map <- 1:nu;
    } else if (method == "randomize") {
      map <- sample(1:nu, nu);
    } else if (method == "reshuffle") {
      map <- sample(uniques, nu);
    }
    s <- map[match];
    setPlate(object, s);
  }
})



# - Deprecated stuff - #

setMethodS3("getID", "Layout", function(this, ...) {
  getId(this, ...);
}, deprecated=TRUE, private=TRUE)

setMethodS3("hasIDs", "Layout", function(this, ...) {
  hasIds(this, ...);
}, deprecated=TRUE, private=TRUE)

setMethodS3("setID", "Layout", function(this, ...) {
  setId(this, ...);
}, deprecated=TRUE, private=TRUE)


setMethodS3("setByLayout", "Layout", function(this, other, ...) {
  this$nspot.r <- other$nspot.r;
  this$nspot.c <- other$nspot.c;
  this$ngrid.r <- other$ngrid.r;
  this$ngrid.c <- other$ngrid.c;
  this$id      <- other$id;

  this;
}, private=TRUE);
            

setMethodS3("setExplicit", "Layout", function(this, ngrid.r=NULL, 
            ngrid.c=NULL, nspot.r=NULL, nspot.c=NULL, name=NULL, id=NULL, ...) {
  if (!is.null(nspot.r)) this$nspot.r <- nspot.r;
  if (!is.null(nspot.c)) this$nspot.c <- nspot.c;
  if (!is.null(ngrid.r)) this$ngrid.r <- ngrid.r;
  if (!is.null(ngrid.c)) this$ngrid.c <- ngrid.c;
  if (!is.null(name))    this$.name   <- name;
  if (!is.null(id))      this$id      <- id;
  this;
}, private=TRUE);

############################################################################
# HISTORY:
# 2005-10-21
# o Replace 'overwrite' arguments with 'mustNotExist' in calls to Arguments. 
# 2005-07-19
# o Replaced all path="" arguments with path=NULL.
# 2005-06-11
# o Making use of Arguments in R.utils.
# 2005-03-08
# o Added automatic detection and reading of gunzipped files (*.gz).
# 2005-02-12
# o getPosition() now returned a "named" matrix, i.e. with column names
#   "row" and "column", and rownames for the spot indicies.
# 2004-08-15
# o Updates "name's" <- "names", "id's" <- "ids" etc. 
# 2003-09-27
# o Added getLayoutGroupsByName(). Very useful in all plot and normalization
#   methods etc that have a groupBy argument.
# 2003-09-20
# o Added getPrintorder() and setPrintorder(). getPrintorder() is to
#   replace toPrintorderMatrix().
# o Removed obsolete '.return' argument from set().
# o Made setByLayout() and setExplicit() deprecated.
# 2003-07-28
# o Added an example to getIndices() how to get the "alley spots", i.e.
#   those spots in each printtip groups that do not have eight neighbors.
# 2002-12-11
# o Renamed getID(), setID() hasIDs() to getId(), setId() and hasIds()
#   according to RCC. The old ones are kept for backward compatibilities.
# o Added rename().
# 2002-12-05
# o Add getSlideRowGroups() and getSlideColumnGroups().
# 2002-11-27
# o fromDataFrame() now deals with all objects (also private) that are
#   vectors and with the correct length.
# o as.data.frame() now includes data from all objects (also private with
#   available get<Field>() functions defined) such that the returned object
#   is a vector and with the correct length.
# o The above will make write() save more fields.
# o as.data.frame() now make use of I() (the "AsIs" class).
# o BUG FIX: setGeneGroups(), setPrinttipGroups() etc all contained syntax
#   errors.
# 2002-11-20
# o Updated as.character() to report data.class(this) instead of hardcoded 
#   "Layout".
# 2002-05-04
# o Added getPrintdipGroups().
# 2002-05-03
# o Added anonymize() and getBlank().
# o Extracted a lot of methods identified to be obsolete into file
#   Layout.obsolete.R
# 2002-05-02
# o TYPO FIX: Type fix in string returned by as.character().
# 2002-05-01
# * Added field plate, getPlate(), setPlate(), hasPlate(), getPlateNumber()
#   and nbrOfPlates().
# * Added as.data.frame(), read() and write().
# 2002-04-21
# * Added getGeneReplicateIndex() and getGeneSlideReplicateIndex()
#   to MicroarrayData.
# * Removed obsolete getReferenceFields().
# 2002-04-13
# * Added optional arguments ignoreCase=TRUE and regexpr=FALSE to indexOf().
#   Updated the Rdoc accordingly and added a few more examples.
#   This update was trigged by questions from Lei Jiang.
# 2002-04-05
# * Added default values of 'index' in getLocation() and getPosition().
# 2002-03-29
# * BUG FIX: Added as.character() to getID() and getName() just to make sure
#   it is not a factor that is returned. Will fix this in GenePixData.R too.
#   This is a double security so it won't happend again.
# 2002-02-26
# * Updated the Rdoc's.
# 2001-01-24
# * Modified source to make use of setMethodS3 and setClassS3.
# * Added nbrOfReplicates().
# 2001-11-18
# * Freshend the code of setReplicates().
# 2001-08-08
# * Added toPrintorderMatrix().
# 2001-08-06
# * BUG FIX: equals() didn't work correctly.
# 2001-08-01
# * Moving the functionalites of replicates back to Layout. For now only
#   the API, but later also the internal code... maybe.
# 2001-07-18
# * Renamed the fields .name and .id to name and id.
# * Created a new class Replicates and moved all replicate/duplicate 
#   functionalities there.
# 2001-07-15
# * Added findReplicates().
# * Added support to get, set and check replicates of genes. Will have to
#   define the concept of a gene index, but it should be pretty straight-
#   forward. Also, I would like to add a method which automatically
#   generates the replicate matrix by look at the list of Names or IDs.
# 2001-07-12
# * Added get- and setName(). GenePix is using both a Name and a ID field.
# 2001-07-11
# * Updated the Rdoc comments.
# 2001-07-09
# * Improved the speed of getPosition() by 30 times!
# 2001-07-07
# * Moved all duplicates/replicates methods in MicroarrayData here.
# * Made all function take modulo n of all indices.
# 2001-07-04
# * Added toXYMatrix(). Useful for plotSpatial() etc.
# 2001-07-01
# * Updated some of the Rdoc comments.
# 2001-06-30
# * Added id's to the Layout class.
# 2001-06-24
# * Made getLocation() work on vectors too. 
# * Added some Rd comments.
# 2001-05-14
# * Added getInternalReferences() for improving gco() performance.
# 2001-04-08
# * Bug fix: Argument indices in getPositions(indices) was misspelled.
# 2001-04-02
# * Added getIndices(),getLocations(), getPostion(), getPositions().
#   However, these methods need some speed optimization because now they're
#   quite slow. I don't wanna spend more time on this right now!
# 2001-03-25
# * Added Rdoc comments.
# 2001-03-19
# * Created.
############################################################################
