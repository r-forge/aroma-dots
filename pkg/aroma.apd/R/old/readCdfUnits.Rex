
filename <- "Mapping50K_Xba240.CDF"

# Read all units in a CDF file [~20s => 0.34ms/unit]
cdf0 <- readCdfUnits(filename, readXY=FALSE, readBases=FALSE)

# Read a subset of units in a CDF file [~6ms => 0.06ms/unit]
units1 <- c(5, 100:109, 34, 50400:50487)
cdf1 <- readCdfUnits(filename, units=units1, readXY=FALSE, readBases=FALSE)
stopifnot(identical(cdf1, cdf0[units1+1]))

# Create a unit name to index map
names <- readCdfUnitNames(filename)
units2 <- match(names(cdf1), names) - 1
stopifnot(identical(units1, units2))
cdf2 <- readCdfUnits(filename, units=units2, readXY=FALSE, readBases=FALSE)

stopifnot(identical(cdf1, cdf2))
