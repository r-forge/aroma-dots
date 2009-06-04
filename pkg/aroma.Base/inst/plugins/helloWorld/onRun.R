onRun <- function(data, welcomeMessage="Hello world!", ..., pluginVersion="0.3.0") {
  # Update progress bar
  setProgress(progress, 0.10);

  #Lc# Output welcome message
  cat("Plugin says: \"", welcomeMessage, "\"\n", sep="");

  #L+# Analyze
  for (kk in 1:10) {
    #Lc# Step ${kk} (of 10)
    increase(progress, 0.05);
  }
  #L-#

  # Update progress bar again
  setProgress(progress, 0.80);

  #Lc# Printing the code for this plugin
  cat("\nThe onParameters() function:\n");
  print(get("onParameters", mode="function", envir=parent.frame()));

  cat("\nThe onRun() function (this function):\n");
  print(get("onRun", mode="function", envir=parent.frame()));

  cat("\nBye! (any messages below are from main() in BasePlugin or from BASE itself)\n");
}
