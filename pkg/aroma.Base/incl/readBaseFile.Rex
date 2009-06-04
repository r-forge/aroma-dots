# Read a simple BASE file structure (adopted from [2])
filename <- system.file("data-ex", "baseExample.base", package="aroma.Base")
base <- readBaseFile(filename)
print(base)

# Write it back to file
tmpfile <- tempfile()
writeBaseFile(tmpfile, base)

# Read the newly create file
base2 <- readBaseFile(tmpfile)
file.remove(tmpfile)

# Compare with the original structure
stopifnot(equals(base, base2))
