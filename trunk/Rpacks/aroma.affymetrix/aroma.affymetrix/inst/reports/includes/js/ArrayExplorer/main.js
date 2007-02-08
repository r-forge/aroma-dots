var explorer = new ArrayExplorer();

window.onresize = function() {
  explorer.update();
}

includeDom("../ArrayExplorer.onLoad.js");
includeDom("extras.js");

function onLoad() {
  explorer.start();
  webcutsOptions['numberLinks'] = false;
}

function changeChipType(chipType) {
  if (explorer.setChipType(chipType))
    explorer.updateImage();
}

function changeSample(sample) {
	if (explorer.setSample(sample))
    explorer.updateImage();
}

function changeColorMap(map) {
  if (explorer.setColorMap(map))
    explorer.updateImage();
}

function changeZoom(scale) {
  if (explorer.setScale(scale))
    explorer.updateImage();
}

