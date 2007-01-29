var explorer = new ArrayExplorer();

window.onresize = function() {
  var y = findXY(explorer.image2d.image).y;
	var h = document.body.clientHeight;
  h = (h - y - 16) + 'px';
	explorer.image2d.container.style.height = h;
  explorer.update();
}

includeDom("../samples.js");
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

