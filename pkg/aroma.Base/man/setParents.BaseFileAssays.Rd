%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BaseFileAssays.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{setParents.BaseFileAssays}
\alias{setParents.BaseFileAssays}
\alias{BaseFileAssays.setParents}
\alias{setParents.BaseFileAssays}
\alias{setParents,BaseFileAssays-method}


 \title{Sets the ids of the assay parents of one or all bioassays in the bioassay set}

 \description{
   Sets the ids of the assay parents of one or all bioassays in the bioassay set.
 }

 \usage{\method{setParents}{BaseFileAssays}(this, ids, parents, ...)}

 \arguments{
  \item{...}{Arguments with names of assay ids and values that are
     \code{\link[base]{character}} strings of parents, either as a \code{\link[base]{character}} \code{\link[base]{vector}} or
     a singe \code{\link[base]{character}} string where parents are split by \code{/}
     (or \code{,} or \code{;}).}
 }

 \value{
   Returns nothing.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   For more information see \code{\link{BaseFileAssays}}.
 }
\keyword{internal}
\keyword{methods}
