%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.R
% 
% on Tue Jan 15 18:34:02 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{append.MicroarrayData}
\alias{append.MicroarrayData}
\alias{MicroarrayData.append}
\alias{append.MicroarrayData}
\alias{append,MicroarrayData-method}


\title{Appends another MicroarrayData object}

\usage{\method{append}{MicroarrayData}(this, other, ...)}

\arguments{
 \item{...}{Other MicroarrayData object(s) to be appended to this object.}
}

\description{
  Appends another identically structured MicroarrayData object(s).
}

\examples{\dontrun{
  # Example that reads all GenePixData files in the current directory
  # and extracts the log-ratios and the log-intensities into *one*
  # MAData object. This approach only stores one GPR file in memory
  # at the time, which makes it possible to read in approximately
  # 10-20 times more files.
  files <- list.files(pattern="*.gpr");
  ma <- MAData();
  for (file in files) {
    gpr <- GenePixData$read(file);
    raw <- getRawData(gpr);
    rm(gpr); # Delete as quick as possible to optimize memory.
    append(ma, getSignal(raw, bgSubtract=TRUE));
    rm(raw); # Delete as quick as possible to optimize memory.
  }
}}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma:removeSlides.MicroarrayData]{*removeSlides}()} and \code{\link[aroma:keepSlides.MicroarrayData]{*keepSlides}()}.
  For more information see \code{\link{MicroarrayData}}.
}
\keyword{internal}
\keyword{methods}
