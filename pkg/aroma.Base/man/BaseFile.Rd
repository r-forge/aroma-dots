%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BaseFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{BaseFile}
\docType{class}
\alias{BaseFile}


 \title{The BaseFile class}

 \description{
  Package:  aroma.Base \cr
\bold{Class BaseFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{BaseFile}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{BaseFile}\cr
extends \link[R.oo]{Object}\cr


 }

 \usage{BaseFile(sections=NULL, ...)}

 \arguments{
   \item{sections}{A \code{\link[base]{list}} BaseFileSection objects.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.Base:appendSection.BaseFile]{appendSection}} \tab Adds a section to this BASE file structure.\cr
 \tab \code{\link[aroma.Base:as.character.BaseFile]{as.character}} \tab Gets a string description of object.\cr
 \tab \code{\link[aroma.Base:as.list.BaseFile]{as.list}} \tab Gets a list representation of object.\cr
 \tab \code{\link[aroma.Base:equals.BaseFile]{equals}} \tab Checks if this object equals another.\cr
 \tab \code{\link[aroma.Base:getAllDataFiles.BaseFile]{getAllDataFiles}} \tab Gets the names of a external data files used by a BASE file.\cr
 \tab \code{\link[aroma.Base:getSection.BaseFile]{getSection}} \tab Gets one section of this BASE file structure.\cr
 \tab \code{\link[aroma.Base:getSections.BaseFile]{getSections}} \tab Gets a subset or all BASE sections.\cr
 \tab \code{\link[aroma.Base:hasSerialBioAssaySets.BaseFile]{hasSerialBioAssaySets}} \tab Checks if the bioassay sets have data in serial format.\cr
 \tab \code{\link[aroma.Base:isSerial.BaseFile]{isSerial}} \tab Checks if a BASE file is in serial format or not.\cr
 \tab \code{length} \tab  -\cr
 \tab \code{\link[aroma.Base:nbrOfSections.BaseFile]{nbrOfSections}} \tab Gets the number of sections in this BASE file structure.\cr
 \tab \code{\link[aroma.Base:read.BaseFile]{read}} \tab Static method to read a BASE file structure.\cr
 \tab \code{\link[aroma.Base:removeAllDataFiles.BaseFile]{removeAllDataFiles}} \tab Removes all external data files used by a BASE file structure.\cr
 \tab \code{\link[aroma.Base:removeSection.BaseFile]{removeSection}} \tab Removes a section from a BASE file structure.\cr
 \tab \code{\link[aroma.Base:replaceSection.BaseFile]{replaceSection}} \tab Replaces a section in a BASE file structure.\cr
 \tab \code{\link[aroma.Base:seq.BaseFile]{seq}} \tab Gets an index sequence for the sections in this BASE file structure.\cr
 \tab \code{\link[aroma.Base:write.BaseFile]{write}} \tab Writes this BASE file structure.\cr
}


 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save

 }

 \details{
 }

 \examples{
# Read a simple BASE file structure (adopted from [2])
filename <- system.file("data-ex", "baseExample.base", package="aroma.Base")
base <- BaseFile$read(filename)
print(base)

# Write it back to file
tmpfile <- tempfile()
write(base, tmpfile)

# Read the newly create file
base2 <- BaseFile$read(tmpfile)
file.remove(tmpfile)

# Compare with the original structure
stopifnot(equals(base, base2))
}

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
