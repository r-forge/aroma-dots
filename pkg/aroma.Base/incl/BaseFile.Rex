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
