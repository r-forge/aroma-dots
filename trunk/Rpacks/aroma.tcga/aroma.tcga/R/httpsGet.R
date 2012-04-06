setMethodS3("httpsGet", "character", function(url, username, password, userAgent=getOption("HTTPUserAgent"), ..., method=c("*", "curl", "wget"), ignoreCertificate=FALSE, filename=NULL, path="*", skip=TRUE, overwrite=!skip, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  parseUrl <- function(url, ...) {
    pattern <- "^(http|https)://([^/]+)(.*)";

    # Sanity check
    protocol <- gsub(pattern, "\\1", url);
    protocol <- tolower(protocol);
    host <- gsub(pattern, "\\2", url);
    path <- gsub(pattern, "\\3", url);

    list(url=url, protocol=protocol, host=host, path=path);
  } # parseUrl()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'url':
  url <- Arguments$getCharacter(url);

  # Argument 'username':
  username <- Arguments$getCharacter(username);

  # Argument 'password':
  password <- Arguments$getCharacter(password);

  # Argument 'userAgent':
  if (!is.null(userAgent)) {
    userAgent <- Arguments$getCharacter(userAgent);
  }

  # Argument 'method':
  method <- match.arg(method);

  # Argument 'ignoreCertificate':
  ignoreCertificate <- Arguments$getLogical(ignoreCertificate);

  # Argument 'filename' & 'path':
  if (is.null(filename)) {
    filename <- basename(url);
  }
  if (path == "*") {
    path <- parseUrl(url)$path;
    path <- dirname(path);
    path <- gsub("^/", "", path);
  }
  pathname <- Arguments$getWritablePathname(filename, path=path,
                                                      mustNotExist=FALSE);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Downloading via HTTPS");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Is file already downloaded?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (isFile(pathname)) {
    fi <- file.info(pathname);
    if (!is.na(fi$size) && fi$size == 0) {
      file.remove(pathname);
    } else if (skip) {
      verbose && cat(verbose, "URL: ", url);
      verbose && printf(verbose, "Pathname: %s [%d bytes]\n", pathname, file.info(pathname)$size);
      verbose && cat(verbose, "Already downloaded. Skipping.");
      return(pathname);
      pathname <- Arguments$getWritablePathname(filename, path=path, 
                                                 mustNotExist=!overwrite);
    }
  }

  verbose && cat(verbose, "URL: ", url);
  verbose && cat(verbose, "User agent: ", userAgent);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Parse URL
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Parsing URL");
  res <- parseUrl(url);
  protocol <- res$protocol;
  host <- res$host;
  path <- res$path;
  verbose && cat(verbose, "Protocol: ", protocol);
  verbose && cat(verbose, "Host: ", host);
  verbose && cat(verbose, "Path: ", path);
  stopifnot(is.element(protocol, c("http", "https")));
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Find working method?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (method == "*") {
    verbose && enter(verbose, "Searching for a working download method");

    methods <- eval(formals(httpsGet.character)$method);
    methods <- setdiff(methods, "*");

    verbose && cat(verbose, "Possible methods:");
    verbose && print(verbose, methods);

    for (methodT in methods) {
      if (is.element(methodT, c("wget", "curl"))) {
        cmd <- methodT;
        cmdT <- sprintf("%s --version", cmd);
        res <- system(cmdT, ignore.stdout=TRUE, ignore.stderr=TRUE);
        if (res == 0) {
          method <- methodT;
          break;
        }
      } else {
        throw("Unknown/unsupported method: ", methodT);
      }
    } # for (methodT ...)
    
    if (method == "*") {      
      throw("Failed to identify a working download method: ", hpaste(methods));
    }

    verbose && cat(verbose, "Download method: ", method);
    verbose && exit(verbose);
  }



  # Download to temporary file
  pathnameT <- pushTemporaryFile(pathname);
  verbose && cat(verbose, "Temporary download pathname: ", pathnameT);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # WGET
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (method == "wget") {
    opts <- sprintf("--output-document=\"%s\"", pathnameT);
    if (protocol == "https") {
      opts <- sprintf("%s --http-user=%s --http-passwd=%s", opts, username, password);
      if (ignoreCertificate) {
        opts <- sprintf("%s --no-check-certificate", opts);
      }
    }

    if (!is.null(userAgent)) {
      opts <- sprintf("%s --user-agent=\"%s\"", opts, userAgent);
    }

    cmd <- sprintf("wget %s %s", opts, url);

    verbose && cat(verbose, "System command:");
    verbose && cat(verbose, cmd);

    # Download
    ignore <- !as.logical(less(verbose, 10));
    res <- system(cmd, ignore.stdout=ignore, ignore.stderr=ignore);

    verbose && cat(verbose, "System command results: ", res);

    if (res != 0) {
      throw(sprintf("Failed to download URL (error code %s): %s", res, url));
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # CURL
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (method == "curl") {
    opts <- sprintf("--output \"%s\"", pathnameT);

    if (protocol == "https") {
      opts <- sprintf("%s --user %s:%s", opts, username, password);
      if (ignoreCertificate) {
        opts <- sprintf("%s --insecure", opts);
      }
    }

    if (!is.null(userAgent)) {
      opts <- sprintf("%s --user-agent \"%s\"", opts, userAgent);
    }

    cmd <- sprintf("curl %s %s", opts, url);

    verbose && cat(verbose, "System command:");
    verbose && cat(verbose, cmd);

    # Download
    ignore <- !as.logical(less(verbose, 10));
    res <- system(cmd, ignore.stdout=ignore, ignore.stderr=ignore);

    verbose && cat(verbose, "System command results: ", res);

    if (res != 0) {
      throw(sprintf("Failed to download URL (error code %s): %s", res, url));
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Cleanup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Remove "empty" downloads
  fi <- file.info(pathnameT);
  if (is.na(fi$size) || fi$size == 0) {
    file.remove(pathnameT);
    pathname <- NULL;
  }

  if (!isFile(pathnameT)) {
    throw("Failed to download URL: ", url);
  }

  # Rename temporary file
  pathname <- popTemporaryFile(pathnameT);

  verbose && printf(verbose, "Pathname: %s [%d bytes]\n", pathname, file.info(pathname)$size);

  verbose && exit(verbose);

  invisible(pathname);
}) # httpsGet()


############################################################################
# HISTORY:
# 2012-04-05
# o ROBUSTNESS: Now httpsGet() downloads are "atomic", i.e. they are
#   downloaded to temporary file which is only renamed if the download
#   was complete.
# o Added support for argument path="*" and method="*".
# o Added support for methods "wget" and "curl".
# o Created httpsGet() from private httpGet().
############################################################################ 
