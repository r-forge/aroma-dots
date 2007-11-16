\name{listFunctions}
\alias{extractLists}
\alias{combineLists}
\alias{mergeDataLists}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ List Helper Functions}
\description{
Tools to manipulate lists.
}
\usage{
extractLists(listOfLists, names)
}
\usage{
combineLists(listOfLists)}
\usage{mergeDataList(listOfData,join=c("row","column"))}

\arguments{
  \item{listOfLists}{a list where each element is a list of the same size}
  \item{names}{ the names of the list elements to extract, see details}
  \item{listOfData}{a list where each element is a matrix or vector}
  \item{join}{whether to join the list elements into a single matrix by rows or columns}
}
\details{
extractLists takes as an argument a list where each element is itself a *named* list. Each named list must have the same set of names and the argument 'names' must be a subset of that. It returns the same length list, only now the element is a subset of the original named list --  the part of the named list corresponding to the argument names. combineLists takes a list of lists and inverts it. mergeDataList merges together the matrix data that compose the element of the list. See examples below. 

}
\value{
Returns named list, where each element is a list of the length of listOfLists
}
\author{ Elizabeth Purdom}
\examples{
lsOfLs<-list(x=list(a=1,b=2),y=list(a=3,b=10))
extractLists(lsOfLs,"a")
combineLists(lsOfLs)
xmat<-matrix(c(1,2,3,4,5,6),nrow=2)
ymat<-matrix(c(11,12,13,14,15,16,17,18),nrow=2)
mergeDataList(list(xmat,ymat),join="column")
mergeDataList(list(t(xmat),t(ymat)),join="row")
}
