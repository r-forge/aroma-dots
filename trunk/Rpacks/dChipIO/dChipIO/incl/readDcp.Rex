path <- system.file("exData", package="dChipIO")

filename <- "Test3-1-121502.dcp"
pathname <- file.path(path, filename)

hdr <- readDcpHeader(pathname)
print(hdr)

data <- readDcp(pathname)
str(data)

