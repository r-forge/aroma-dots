var explorer = new ChromosomeExplorer();


window.onresize = function() {
  explorer.update();
}

includeDom("ChromosomeExplorer.onLoad.js");
includeDom("../ChromosomeExplorer.onLoad.js");
includeDom("../samples.js");
includeDom("extras.js");

function onLoad() {
  logAdd("onLoad()...");
  explorer.onLoad();
  logAdd("onLoad()...done");
  explorer.start();
  webcutsOptions['numberLinks'] = false;
}


var nav = null;
var navArea = null;
var navAreaWidth = 0;
var navImage = null;
var navImageOffsetX = 0;
var navImageWidth = 0;

var panel = null;
var panelImage = null;
var panelImageOnLoad = function() {};
var panelImageWidth = 0;
var panelImageOffsetX = 0;
var panelLocator = null;
var panelLocatorTag = null;
var panelWidth = 0;
var panelMaxWidth = 0;


var playSamples = false;
var playDelay = 2000;

function gotoNextSample(step) {
  var nextSampleIdx = sampleIdx + step;
  if (nextSampleIdx >= samples.length) {
    nextSampleIdx = 0;
  } else if (nextSampleIdx < 0) {
    nextSampleIdx = samples.length-1;
  }
  explorer.setSample(nextSampleIdx);
  if (playSamples) {
    cmd = "gotoNextSample(" + step + ");";
    setTimeout(cmd, playDelay);
  }
}

function playAlongSamples(cmd) {
  if (cmd == "start") {
    playSamples = true;
    gotoNextSample(1);
  } else if (cmd == "stop") {
    playSamples = false;
  }
}



/* NOT USED */
var shortcuts = new Array();
var shortcutLabels = new Array();

