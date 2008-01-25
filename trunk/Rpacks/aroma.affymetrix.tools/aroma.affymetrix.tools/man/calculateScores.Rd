\name{calculateScores}
\alias{calculateScores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate Scores  per unit }
\description{
Replicates the structure of functions in aroma.affymetrix to calculate user-defined functions on user-defined unit data. 
}
\usage{
calculateScores(cdf, units = NULL, calculateFUN, readFUN, nbrOfArrays, filename, path = NULL, ram = 1, verbose = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cdf}{ cdf to use to find unit structure, etc. }
  \item{calculateFUN}{user defined function to calculate scores on unit input, see details }
  \item{readFUN}{ user defined function to read data, see details}
  \item{filename}{ the core part of filename of saved files}
  \item{path}{ (optional) the path for the folder where the scores are saved }
  \item{nbrOfArrays}{ number of arrays (really should get rid of this required option) }
  \item{units}{ units to search (but dicey to specify if going to use other scores functions)}
}
\details{
calculateScores replicates the computational approach of aroma.affymetrix by calculating the scores on groups of units, storing them to the hard drive, and then moving to the next group of units. The difference is that we have not written code to write to CEL files; rather each 'chunk' of units is saved as a standard .rdata file of the name 'filename,Chunkxx.RDATA' in the folder specified with 'path'. 

readFUN is the function defined by the user to 'read' the data for specified units. It must take as arguments 'units' and 'verbose' (and ONLY those arguments). It must return a list, where each element of the list corresponds to a unit and will be the input to calculateFUN. The function 'combineLists' can be useful in creating a unified unit-level list from two separate unit-level lists that are the output of readUnits, for example.

calculateFUN is the function defined by the user to calculate the scores based on the output of readFUN. It must take as an argument an element of the list returned by readFUN. In otherwords, the following code is used in the function:
\code{
dataList <- readFUN(units = units, verbose = less(verbose))
scoresList <- lapply(dataList, FUN = calculateFUN)
}
The output from calculateFUN will be the 'unit' object in the returned list.  There's no writing to CEL files, so any output format is allowed. However, in order to use other functions to manage the scores, calculateFUN must return a 'group list'-- a list with an element corresponding to each group. If your scores are not per group, then you can simply define your output to have only 1 group per unit(like a cel file with a monocell cdf, if we were working with cel files). In otherwords, to return object 'x' that is not a list, instead return 'list(group=x)'. We would recommend that you go even further, and give an additional list structure to your output so that it will match the output from readUnits of an AffymetrixCelSet. This should give a final output of the following format:

     scoreList
     
      +- unit 1
      
      |  +- group 1
      
      |  |   +- Score element 1 
      
      |  |   |   +- data
      
      |  |   +- Score element 2
      
      |  |   |   +- data
      
      .  .   .
 
      |  +-group 2
      
      |  |   +- Score element 1 
      
      |  |   |   +- data
      
      |  |   +- Score element 2
      
      |  |   |   +- data
      
      .  .   .
      
      |  +- group K
      
      |      +- Score element 1 
      
      |      |   +- data
      
      |      +- Score element 2
      
      |      |   +- data
      
      .      .
      
      +- unit 2
      
      .
      
      +- unit J

}
\value{
A list of the length of the number of units, each element corresponding to the output of calculateFUN for that unit.
}
\author{Elizabeth Purdom}
