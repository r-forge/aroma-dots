%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  999.package.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{aroma-package}
\alias{aroma-package}
\alias{aroma}
\docType{package}

\title{Package aroma}


 \description{
   Methods for microarray analysis after the images have been scanned. Formerly known as the com.braju.sma package.

   This package was formerly known as the \pkg{com.braju.sma} package.
   On 2003-12-31 it was renamed to the current name \pkg{aroma}, which
   is an achronym for \emph{An R Object-oriented Microarray Analysis}
   package.
 }

 \section{Requirements}{
   This package depends on several other packages, namely
   \pkg{R.oo}, \pkg{R.utils}, \pkg{R.basic}, \pkg{R.graphics},
   \pkg{R.io}, and \pkg{aroma.light}.
 }

 \section{Installation}{
   To install this package and all required package, see instructions at
   \url{http://www.braju.com/R/}.\cr
 }

 \section{To get started}{
   To get started, see:
   \enumerate{
     \item \code{\link{GenePixData}}, \code{\link{QuantArrayData}}, \code{\link{ScanAlyzeData}},
           \code{\link{SpotData}}, or \code{\link{SpotfinderData}}.
     \item \code{\link{calibrateMultiscan.RGData}}().
     \item \code{\link{normalizeAffine.RGData}}().
   }
 }

 \section{Further readings}{
   For an overview of \pkg{aroma}, see [1].
   See also the online webpage of \code{aroma} -
   \url{http://www.maths.lth.se/help/R/aroma/}.
   Several of the normalization methods proposed in [3]-[5] are
   available in this package.
 }

 \section{How to cite this package}{
   Whenever using this package, please cite [1] as\cr

   \preformatted{@TECHREPORT{BengtssonH_2004,
  author = {Bengtsson, Henrik},
  title = {{aroma} - {A}n {R} {O}bject-oriented {M}icroarray
           {A}nalysis environment},
  institution = {Mathematical Statistics, Centre for Mathematical
                 Sciences, Lund University, Sweden},
  year = {2004},
  type = {{Preprint in Mathematical Sciences}},
  number = {2004:18},
}

@ARTICLE{BengtssonHossjer_2006,
  author = {Henrik Bengtsson and Ola H{\"o}ssjer},
  title = {Methodological study of affine transformations of gene
           expressiondata with proposed robust non-parametric
           multi-dimensional normalization	method},
  journal = {BMC Bioinformatics},
  year = {2006},
  volume = {7},
  number = {100},
  month = {mar},
  doi = {10.1186/1471-2105-7-100}
}

@ARTICLE{BengtssonH_etal_2004,
  author = {Henrik Bengtsson and G\"{o}ran J\"{o}nsson and
            Johan Vallon-Christersson},
  title = {Calibration and assessment of channel-specific biases
           in microarray data with extended dynamical range},
  journal = {BMC Bioinformatics},
  year = {2004},
  volume = {5},
  number = {177},
  doi = {10.1186/1471-2105-5-177},
}}
\emph{}

 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \references{
  Some of the reference below can be found at
  \url{http://www.maths.lth.se/bioinformatics/publications/}.\cr

 [1] H. Bengtsson, \emph{aroma - An R Object-oriented Microarray
     Analysis environment}, Preprints in Mathematical Sciences (manuscript
     in preparation), Mathematical Statistics, Centre for Mathematical
     Sciences, Lund University, 2004.\cr

 [2] H. Bengtsson, \emph{The R.oo package - Object-Oriented Programming with References Using Standard R Code}, In Kurt Hornik, Friedrich Leisch and Achim Zeileis, editors, Proceedings of the 3rd International Workshop on Distributed Statistical Computing (DSC 2003), March 20-22, Vienna, Austria. \url{http://www.ci.tuwien.ac.at/Conferences/DSC-2003/Proceedings/}
 \cr

 [3] Henrik Bengtsson and Ola Hössjer, \emph{Methodological Study of Affine Transformations of Gene Expression Data}, Methodological study of affine transformations of gene expression data with proposed robust non-parametric multi-dimensional normalization method, BMC Bioinformatics, 2006, 7:100.
 \cr

 [4] H. Bengtsson, J. Vallon-Christersson and G. Jönsson, \emph{Calibration and assessment of channel-specific biases in microarray data with extended dynamical range}, BMC Bioinformatics, 5:177, 2004.
 \cr

 [5] H. Bengtsson, \emph{Identification and normalization of plate effects
     in cDNA microarray data}, Preprints in Mathematical Sciences,
     2002:28, Mathematical Statistics, Centre for Mathematical Sciences,
     Lund University, 2002.\cr
 }
\keyword{package}
