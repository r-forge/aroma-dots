setConstructorS3("SSMatrix", function(x=matrix(), nbrOfSpots=nrow(as.matrix(x)), nbrOfSlides=ncol(as.matrix(x)), nbrOfReplicates=1, geneSpotMap="neighboring", ...) {
  x <- extend(Matrix(x, nrow=nbrOfSpots, ncol=nbrOfSlides, ...), "SSMatrix");

  K <- ncol(x);  # nbrOfSlides;
  S <- nrow(x);  # nbrOfSpots;
  if (!is.null(K) && !is.null(S) && K != 0 && S != 0) {
    R <- nbrOfReplicates;
    G <- S/R;
    if (G %% 1 != 0)
      throw("Incompatible structure of SSMatrix. The number of spots is not a multiple of the number of within-slide replicates.");
  }

  attr(x, "nbrOfReplicates") <- nbrOfReplicates;
  attr(x, "geneSpotMap") <- geneSpotMap;
  # CLASS CONSTANTS
  attr(x, "PER.SPOT") <- 1;
  attr(x, "PER.SLIDE") <- 2;
  x;
});

setMethodS3("asSSMatrix", "ANY", function(object, ...) {
  SSMatrix(object, ...);
})

setMethodS3("asSSMatrix", "matrix", function(object, ...) {
  SSMatrix(object, ...);
})

setMethodS3("asSSMatrix", "SSMatrix", function(object, ...) {
  object;
})


setMethodS3("asSSMatrix", "GSRArray", function(Mnew, geneSpotMap="neighboring") {
  G <- dim(Mnew)[1];  # nbrOfGenes
  K <- dim(Mnew)[2];  # nbrOfSlides
  R <- dim(Mnew)[3];  # nbrOfReplicates
  S <- G*R;           # nbrOfSlides

  M <- rep(NA, S*K);
  dim(M) <- c(S, K);

  if (geneSpotMap == "neighboring") {
    for (k in 1:K) {
      # Get the column of spot measurements from slide k
      Mslide <- Mnew[,k,];
      # Convert it into a vector of length S=G*R
      Mslide <- as.vector(t(Mslide));
      M[,k] <- Mslide;
    }
  } else {
    throw("Unknown value of argument 'geneSpotMap': ", geneSpotMap);
  }

  Matrix(M);
})


setMethodS3("plotXY", "SSMatrix", function(object, object2, rows=NULL, columns=NULL, ...) {
  # Assert correct values on arguments
  if (!inherits(object, "SSMatrix"))
    throw("Second object must be of class SSMatrix: ", data.class(object));

  # Argument 'columns'
  if (is.null(columns)) columns <- seq(ncol(object));
  if (any(columns < 0))
    throw("Argument 'columns' must contain postive integers.");
  if (any(columns > ncol(object)))
    throw("There a fewer columns in 'object' than specified by argument 'columns'.");
  if (any(columns > ncol(object2)))
    throw("There a fewer columns in 'object2' than specified by argument 'columns'.");

  # Argument 'rows'
  if (is.null(rows)) rows <- seq(nrow(object));
  if (any(rows < 0))
    throw("Argument 'rows' must contain postive integers.");
  if (any(rows > nrow(object)))
    throw("There a fewer rows in 'object' than specified by argument 'rows'.");
  if (any(rows > nrow(object2)))
    throw("There a fewer rows in 'object2' than specified by argument 'rows'.");

  plot(object[rows,columns], object2[rows,columns], ...);
})


############################################################################
# HISTORY:
# 2002-09-12
# o Making use of S3 instead of S4 object oriented style. S4 is buggy.
# o Converted from SSMatrix.R in com.braju.sma v0.48.
# 2002-05-21
# * BUG FIX: Moved the setClass() method(s) to "000.R" to make sure they
#   are called in the correct order.
# 2002-01-10
# * Created!
############################################################################
