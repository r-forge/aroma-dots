setMethodS3("getFilename", "Arguments", function(static, filename, nchar=c(1,64), class=c("safe"), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'filename':
  filename <- getCharacter(static, filename, nchar=nchar);

  # Argument 'class':
  class <- match.arg(class);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Filter out valid characters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  chars <- filename;

  # Always valid characters
  chars <- gsub("[a-zA-Z0-9_.,]", "", chars);
  chars <- gsub("[-]", "", chars);
  chars <- gsub("[+]", "", chars);

  # Filter out according to classes.
  if ("safe" %in% class) {
    chars <- gsub("[ ]", "", chars);
    chars <- gsub("[\\[\\]]", "", chars);
    chars <- gsub("[#$%&'()`{|}~]", "", chars);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for remaining (=invalid) characters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (nchar(chars) > 0) {
    chars <- unlist(strsplit(chars, split=""));
    chars <- sort(unique(chars));
    chars <- paste(chars, collapse="");
    throw("Not a valid filename. Argument 'filename' contains non-valid filename characters (", chars, "): ", filename);
  }

  filename;
}, static=TRUE, protected=TRUE)


############################################################################
# HISTORY:
# 2006-11-20
# o Added static getFilename() to Arguments to check if a string contains
#   valid filename characters.
############################################################################
