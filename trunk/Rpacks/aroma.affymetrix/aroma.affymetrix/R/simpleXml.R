setMethodS3("listToXml", "list", function(tree, indentStep=" ", collapse="\t", ...) {
##   tree <- list(
##     chipType = "GenomeWideSNP_6",
##     createdOn = "2008-02-13",
##     source = list(
##       url = "http://www.affymetrix.com/",
##       files = ""
##     ),
##     creator = list(
##       name  = "Henrik Bengtsson",
##       email = "hb@stat.berkeley.edu"
##     )
##   )

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  asXml <- function(tree, indentStep=" ", .indent="", collapse="\t") {
    xml <- c();

    names <- names(tree);
    for (name in names) {
      child <- tree[[name]];
      if (is.list(child)) {
        xmlChild <- asXml(child, indentStep=indentStep, 
                               .indent=paste(.indent, indentStep, sep=""));
        xmlChild <- sprintf("%s<%s>\n%s%s</%s>\n", 
                     .indent, name, as.character(xmlChild), .indent, name);
      } else {
        child <- paste(child, collapse=collapse);
        xmlChild <- sprintf("%s<%s>%s</%s>\n",
                                 .indent, name, child, name);
      }
      xml <- c(xml, xmlChild);
    }

    paste(xml, collapse="");
  } # asXml()  


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  asXml(tree, indentStep=indentStep, collapse=collapse, ...);
}) # listToXml()



setMethodS3("xmlToList", "character", function(xml, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  parseXml <- function(xml, ...) {
    xml <- trim(xml);

    res <- list();

    while (nchar(xml) > 0) {
      # Find first tag
      pattern <- "^<([a-zA-Z][^>]*)>(.*)";
      hasStartTag <- (regexpr(pattern, xml) != -1);
      if (hasStartTag) {
        # Find second tag
        tag <- gsub(pattern, "\\1", xml);
        xml <- gsub(pattern, "\\2", xml);
        xml <- trim(xml);
        pattern <- sprintf("(.*)</%s>(.*)", tag);
        hasEndTag <- (regexpr(pattern, xml) != -1);
        if (!hasEndTag) {
          throw("File footer parsing error: Missing </", tag, ">: ", xml);
        }
        body <- gsub(pattern, "\\1", xml);
        body <- trim(body);
        if (is.null(res))
          res <- list();
        res[[tag]] <- parseXml(body);

        xml <- gsub(pattern, "\\2", xml);
        xml <- trim(xml);
      } else {
        res <- xml;
        xml <- "";
      }
    } # while(...)

    if (length(res) == 0)
      res <- "";

    res;
  } # parseXml()

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  xml <- paste(xml, collapse="");
  parseXml(xml);
}) # xmlToList()



############################################################################
# HISTORY:
# 2008-03-05
# o BUG FIX: Regular expression pattern 'a-Z' is illegal on (at least) some
#   OSX systems (where 'A-z' works). Replaced it with the safer 'a-zA-Z'.
# 2008-02-13
# o Validation test: identical(tree, xmlToList(listToXml(tree)))
# o Created to support read/write of footer in AromaTabularBinaryFile.R.
############################################################################
