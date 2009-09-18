%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  readBaseFileSection.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{readBaseFileSection}
\alias{readBaseFileSection.default}
\alias{readBaseFileSection}

\title{Low-level function to read a BASE file section from a connection or a file}

\description{
 Low-level function to read a BASE file section from a connection or a file. This a supportive function to readBaseFile().
}

\usage{\method{readBaseFileSection}{default}(con, suppressSingleAssayLabels=TRUE, extractSpotsData=FALSE, verbose=FALSE, ...)}

\arguments{
  \item{con}{A \code{\link[base:connections]{connection}} or a \code{\link[base]{character}} string filename.}
  \item{suppressSingleAssayLabels}{If \code{\link[base:logical]{TRUE}} and only a single assay is
    read, field names in data table is non suffixed with assay name.}
  \item{extractSpotsData}{If \code{\link[base:logical]{TRUE}}, data tables in 'spots' sections are
    written to seperate files such that each assay has its own file.
    Files are written as tab-delimited files to current directory,
    using filename format "spots\_<assay>.txt".
    Written data is removed from the final list structure and replaced
    with a header \code{dataFiles} of the filenames for each assay.}
  \item{verbose}{Either a \code{\link[base]{logical}}, a \code{\link[base]{numeric}}, or a \code{\link[R.utils]{Verbose}}
    object specifying how much verbose/debug information is written to
    standard output. If a Verbose object, how detailed the information is
    is specified by the threshold level of the object. If a numeric, the
    value is used to set the threshold of a new Verbose object. If \code{\link[base:logical]{TRUE}},
    the threshold is set to -1 (minimal). If \code{\link[base:logical]{FALSE}}, no output is written.
  }
  \item{...}{Not used.}
}

\value{
  Returns a named \code{\link[base]{list}} structure containing the parsed BASE structure.
}

\section{Out of memory?}{
  Some BASE files are very large due to large amount of assay data in
  'spots' sections.  By setting argument \code{extractSpotsData=TRUE},
  such data is written to separate files for each assay immediately after
  being read and thereafter removed from memory.  This means that it is
  in practice possible to read BASE files with almost any number of
  assays, also if they are "passed" via stdin.
}

\section{Slow performance?}{
  The number of lines in the data table of a BASE file data section is
  either not known on before hand or specified by the \code{count} header.
  If \code{count} is available, the data table is read approximately
  eight times faster than if it is not.
  To the best of our knowledge, starting with BASE v1.2.17 (released
  November 2005), header \code{count} is included by default and in
  previous versions it is not available at all.
  If \code{count} is not specified, the table is first read line by line
  to a temporary file and then re-read by \code{read.table()}. (The
  reason for this is that \code{read.table()} otherwise will read to end
  of file, and not to the first empty line. Unfortunately, this means
  that the reading of data is a bit slow.)
}

\seealso{
  \code{\link{readBaseFile}}().
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}



\keyword{file}
\keyword{IO}
