###########################################################################/**
# @RdocDefault extractBaseFileSpotsData
#
# @title "Low-level function to extract assay data and save it to file instead"
#
# \description{
#  @get "title" to save memory.
# }
#
# @synopsis
#
# \arguments{
#   \item{section}{A @list structure containing a 'spots' section.}
#   \item{serialize}{If @TRUE, non-serialized data is serialized before being
#     extracted so that data is saved in seperate files for each assay.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a @list with element \code{section} containing the modified
#   'spots' section. 
# }
#
# @author
#
# \seealso{
#   @see "readBaseFileSection".
# }
#
# @keyword file
# @keyword IO
#*/###########################################################################
setMethodS3("extractBaseFileSpotsData", "default", function(section, serialize=TRUE, ...) { 
  headers <- section$headers;
  type <- headers$section;
  if (type != "spots")
    return(NULL);

  assays <- headers$assays;

  filenames <- c();
  if (length(assays) == 1 || !serialize) {
    assays <- paste(assays, collapse="-");
    filename <- paste(type, "_", assays, ".txt", sep="");
    
    data <- section$data;
    write.table(data, file=filename, quote=FALSE, sep="\t", 
                             row.names=FALSE, col.names=TRUE);

    filenames <- filename;
  } else {
    # "Serialize"  un-serialized data table, i.e. cache in seperate files.
    columns <- setdiff(headers$columns, "assayData");
    assayFields <- headers$assayFields;
    for (assay in assays) {
      filename <- paste(type, "_", assay, ".txt", sep="");
      filenames <- c(filenames, filename);
      fields <- c(columns, paste(assayFields, assay, sep="_of_"));
      data <- section$data[,fields];
      colnames(data) <- c(columns, assayFields);
      write.table(data, file=filename, quote=FALSE, sep="\t", 
                               row.names=FALSE, col.names=TRUE);
    }
  }

  # 1. Add list of new data filenames
  section$headers[["dataFiles"]] <- filenames;

  # 1. Remove data
  section$data <- NULL;

  list(section=section);
}, private=TRUE)


############################################################################
# HISTORY: 
# 2005-12-21
# o Changed Rdoc title.
# 2005-06-19
# o Added Rdoc comments.
# 2005-06-15
# o Created.
############################################################################

