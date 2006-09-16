# Search for a specific dataset
path <- findCelSet("SlaterH_etal_2004")
print(path)

# Search for a specific dataset of a specific chip type
path <- findCelSet("SlaterH_etal_2004", chipType="Xba")
print(path)

# Find a dataset with a specific chip type
path <- findCelSet(chipType="Sty", firstOnly=TRUE)
print(path)

# Find a dataset with at least 8 arrays
path <- findCelSet(minCount=8, firstOnly=TRUE)
print(path)

