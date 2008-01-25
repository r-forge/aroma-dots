\name{searchScores}
\alias{searchScores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Search For High Scores}
\description{
Searches saved scores (either .rdata files, or an aroma.affymetrix object) and returns the data corresponding to probesets that have a high enough score. Allows user definitions of what constitutes a high score
}
\usage{
searchScores(this = NULL, filename = NULL, nbrOfArrays, units = NULL, cutoff = NULL, quantity = NULL, path = NULL, ram = 1, verbose = 0, type = c("high", "low", "highlow"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{this}{ an aroma.affymetrix object}
  \item{filename}{ the core part of filename of saved files (given to calculateScores) }
  \item{path}{ (optional) the path for the folder where the scores are saved }
  \item{nbrOfArrays}{ number of arrays (really should get rid of this required option) }
  \item{units}{ units to search (only valid if giving an aroma.affymetrix object, see details)}
  \item{cutoff}{ either a numerical cutoff value to search over OR a logical function, see details }
  \item{quantity}{ number of top scores to return}
  \item{type}{ whether to return high, low, or high and low qualifying scores}
} \details{ The function is designed to return 'high' scores, but note the specific way in which it does it. It 
requires a list in the format of an AffymetrixCelSet list, either from readUnits of an aroma.affymetrix object or a 
saved .rdata file object. The function searches per unit each 'group'; if *any* of the scores saved for that group 
meet the requirement, then all the scores from that group will be returned (not just the scores that meet the 
requirement). The returned list follows the structure of the original searched object, except that groups with no 
scores matching the criteria are dropped, and units with no groups matching the requirement are dropped.  

Two simple search requirements are possible by default. If quantity is specified as 'xxx', the top 'xxx' of the 
scores are returned. Note that if there are more than one score per group, then each of these scores count. If you 
wish to find the top groups, you must create scores that summarize to 1 value per group. If a numeric value is given 
for 'cutoff', then a simple search is done for any scores above that value. Any scores above that level satisfy the 
requirement and their corresponding probesets will be returned.  You can specify both quantity and cutoff and then 
only scores that are both in the top 'xxx' and above the cutoff will be returned. Note that for both quantity and 
cutoff, you can change the option of 'type' to high, low, or highlow to return high valued, low valued or both high 
and low valued scores. 'low' and 'highlow' option is done by taking the function -x or abs(x) of the scores *and* 
the numeric cutoff value, so be careful of the cutoff criterion in this case. 

The user can also define their requirement by giving a value for 'cutoff' that is a function that returns a single 
T/F when applied to the group. Thus the standard cutoff criteria is a function 
    \code{cutoffFun<-function(x,cutoff){
                x<-func(na.omit(unlist(x))); 
                return(any(x>=func(cutoff)))
            }
    }
where 'func' is the function that returns x, -x or abs(x) depending on the value of 'type'. Note that if the saved 
scores follow the AffymetrixCelSet list style, as recommended, then the input x to the function will be a list with 
named elements (e.g. 'intensities' and 'stdevs' for intensities data) 

Currently the function pulls the data exactly as saved, without any logging of the data, even if 
searching .CEL files
} 

\value{ Returns a list of the same form as the searched list object but with only those groups and units that 
contain a score meeting the requirements(see details) } 

\author{Elizabeth Purdom} 
