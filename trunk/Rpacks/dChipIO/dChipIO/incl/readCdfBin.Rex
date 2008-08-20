path <- system.file("exData", package="dChipIO")
chipType <- "Test3"
filename <- sprintf("%s.CDF.bin", chipType)
pathname <- file.path(path, filename)

hdr <- readCdfBinHeader(pathname)
print(hdr)

data <- readCdfBin(pathname)
str(data)

