# Escaping LaTeX symbols
le <- function(s, ...) {
  s <- gsub("_", "\\_", s, fixed=TRUE);
  s;
} # le()
