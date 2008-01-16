setConstructorS3("GSRArray", function(x=matrix(), nbrOfGenes=nrow(as.matrix(x)), nbrOfSlides=ncol(as.matrix(x)), nbrOfReplicates=1, ...) {
  if (missing(nbrOfGenes)) {
    if (!missing(nbrOfReplicates))
      nbrOfGenes <- nbrOfGenes / nbrOfReplicates;
    if (nbrOfGenes %% 1 != 0)
      throw("Number of replicates is not compatible with the observed number of spots.");
    if (length(dim(x)) == 3) {
      # Assuming a 3D array with (G,S,R) structure.
      dim <- dim(x);
    } else {
      dim <- c(nbrOfGenes, nbrOfSlides, nbrOfReplicates);
    }
  } else {
    dim <- c(nbrOfGenes, nbrOfSlides, nbrOfReplicates);
  }

  x <- extend(MultiwayArray(x, dim=dim), "GSRArray");
  # CLASS CONSTANTS
  attr(x, "PER.GENE") <- 1;
  attr(x, "PER.SLIDE") <- 2;
  attr(x, "PER.REPLICATE") <- 3;
  attr(x, "PER.GENE.AND.SLIDE") <- c(1,2)
  attr(x, "PER.GENE.AND.REPLICATE") <- c(1,3)
  attr(x, "PER.SLIDE.AND.REPLICATE") <- c(2,3)
  x;
})


setMethodS3("asGSRArray", "ANY", function(object, ...) {
  GSRArray(object, ...);
})

setMethodS3("asGSRArray", "matrix", function(object, ...) {
  GSRArray(object, ...);
})

setMethodS3("asGSRArray", "GSRArray", function(object, ...) {
  object;
})

setMethodS3("asGSRArray", "SSMatrix", function(object) {
  K <- ncol(object);  # nbrOfSlides;
  S <- nrow(object);  # nbrOfSpots;
  R <- object@nbrOfReplicates;
  G <- S/R;

  Mnew <- rep(NA, S*K);
  dim(Mnew) <- c(G, K, R);

  if (object@geneSpotMap == "neighboring") {
		for (k in 1:K) {
			# Get the column of spot measurements from slide k
			Mslide <- object[,k];
			# Convert it into a matrix KxR matrix
			Mslide <- matrix(Mslide, nrow=G, ncol=R, byrow=TRUE);
			Mnew[,k,] <- Mslide;
		}
  } else {
    throw("Unknown value of argument 'geneSpotMap': ", object@geneSpotMap);
  }

  GSRArray(Mnew);
})


setMethodS3("as.character", "GSRArray", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- paste(data.class(this), ": ", sep="");
})

setMethodS3("dim<-", "GSRArray", function(x, value) {
  # To please R CMD check
  this <- x;

  # Asserts:
  if (length(value) != 3)
    throw("A GSRMatrix must have three and only three dimensions.");
  if (!is.numeric(value) || any(value < 0))
    throw("The specified dimension must be a vector of non-negative integers.");

  class.this <- class(this);
  class(this) <- NULL;
  dim(this) <- value;
  class(this) <- class.this;
  this;
}, appendVarArgs=FALSE)  # 'FALSE' not needed with future versions of R.oo.


# MARGIN = NULL: return the mean of all spot measures.
# MARGIN = 1: returns a vector containing the mean for each gene.
# MARGIN = 2: returns a vector containing the overall mean for each slide.
# MARGIN = 3: returns a vector containing the overall mean for each within-slide replicate. Not really interesting.
# MARGIN = c(1,2): returns a matrix containing the mean for each gene per slide.
# MARGIN = c(1,3): returns a matrix containing the mean for each gene per within-slide replicate.
# MARGIN = c(2,3): returns a matrix containing the mean for each slide per within-slide replicate.
# MARGIN = c(1,2,3): (with na.rm=TRUE) returns a  matrix with the same values and of the same dimensions as the given matrix.
#
#    mean(gsr, MARGIN=c(a,b)) = t(mean(gsr, MARGIN=c(b,a)))
#
setMethodS3("apply", "GSRArray", function(this, MARGIN=c(), FUN, ...) {
  warn <- unlist(options("warn")); options(warn=-1);
  idx <- is.na(as.numeric(MARGIN));
  options(warn=warn);
  if (any(idx)) {
    newMARGIN <- as.numeric(MARGIN);
    newMARGIN[idx] <- unlist(lapply(MARGIN[idx],
                                  FUN=function(x) attr(this, x)));
    MARGIN <- newMARGIN;
    if (any(is.null(MARGIN)))
      throw("MARGIN contains an invalid value.");
  }

  res <- apply.MultiwayArray(this, MARGIN=MARGIN, FUN=FUN, ...);
  if (is.null(MARGIN))
    res
  else if (all(MARGIN == this@PER.GENE.AND.SLIDE))
    GSRArray(res)
  else
    res;
})


setMethodS3("fromSSMatrix", "GSRArray", function(M, nbrOfReplicates=2, geneSpotMap="neighboring") {
  K <- ncol(M);  # nbrOfSlides;
  S <- nrow(M);  # nbrOfSpots;
  R <- nbrOfReplicates;
  G <- S/R;

  Mnew <- rep(NA, S*K);
  dim(Mnew) <- c(G, K, R);

  if (geneSpotMap == "neighboring") {
		for (k in 1:K) {
			# Get the column of spot measurements from slide k
			Mslide <- M[,k];
			# Convert it into a matrix KxR matrix
			Mslide <- matrix(Mslide, nrow=G, ncol=R, byrow=TRUE);
			Mnew[,k,] <- Mslide;
		}
  } else {
    throw("Unknown value of argument 'geneSpotMap': ", geneSpotMap);
  }

  GSRArray(Mnew);
})

           

############################################################################
# HISTORY:
# 2002-09-12
# o Making use of S3 instead of S4 object oriented style. S4 is buggy.
# o Converted from GSRMatrix.R in com.braju.sma v0.48.
# 2002-05-21
# * BUG FIX: Moved the setClass() method(s) to "000.R" to make sure they
#   are called in the correct order.
# 2001-11-15
# * Added dim() and "dim<-"().
# 2001-11-14
# * Created!
############################################################################
