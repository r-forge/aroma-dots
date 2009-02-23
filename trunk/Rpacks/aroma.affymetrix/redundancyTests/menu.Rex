if (interactive()) savehistory();
library("R.menu");

library("aroma.core");
setOption(aromaSettings, "transitionRules/useVerboseOption", FALSE);
setVerbose(aromaSettings, FALSE);

x <- getVerbose(aromaSettings);
print(getVerbose(aromaSettings));

launchMenu("testScripts");
