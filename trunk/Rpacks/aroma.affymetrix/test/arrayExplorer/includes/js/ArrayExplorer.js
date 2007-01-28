/****************************************************************
 * ArrayExplorer()
 *
 * Author: Henrik Bengtson, hb@stat.berkeley.edu
 ****************************************************************/
function ArrayExplorer() {
  this.showIndicator = function(state) {
    var statusImage = document.getElementById('statusImage');
    if (state) {
      statusImage.style.visibility = 'visible';
    } else {
      statusImage.style.visibility = 'hidden';
    }
  }

  this.setStatus = function(state) {
    if (state == "") {
      this.showIndicator(false);
      this.image2d.image.style.filter = "alpha(opacity=100)";
      this.image2d.image.style.opacity = 1.0;
    } else if (state == "wait") {
      this.showIndicator(true);
      this.image2d.image.style.filter = "alpha(opacity=50)";
      this.image2d.image.style.opacity = 0.50;
    }
  }

  this.getImageUrl = function() {
    var imgName = this.sample + "," + this.colorMap + ".png";
    var pathname = this.chipType + '/' + imgName;
    return(pathname);
  }

  this.updateImage = function() {
    var pathname = this.getImageUrl();
    this.loadCount = 2;
    this.setStatus('wait');
    this.nav2d.setImage(pathname);
    this.image2d.setImage(pathname);

    /* Update the image link */
    var imageUrl = document.getElementById('imageUrl');
    imageUrl.href = pathname;
    updateText(imageUrl, pathname);

    /* Update the title of the page */
    var title = location.href;
    title = title.substring(0, title.lastIndexOf('\/'));
    title = title.substring(title.lastIndexOf('\/')+1);
    title = title + '/' + pathname;
    document.title = title;
	}

  this.decreaseLoadCount = function() {
    this.loadCount = this.loadCount - 1;
    if (this.loadCount <= 0) {
      this.loadCount = 0;
      this.setStatus("");
    }
  }

  this.setColorMap = function(map) {
    if (this.colorMap == map)
      return(false);

    clearById('colorMap' + this.colorMap);
    highlightById('colorMap' + map);
    this.colorMap = map;
    this.updateImage();
    return(true);
  }

  this.setChipType = function(chipType) {
    if (this.chipType == chipType)
      return(false);

    clearById('chipType' + this.chipType);
    highlightById('chipType' + chipType);
    this.chipType = chipType;
    return(true);
  }

  this.setSample = function(sample) {
    if (this.sample == sample)
      return(false);

    clearById('sample' + this.sample);
    highlightById('sample' + sample);
    updateLabel('sampleLabel', sample);
    this.sample = sample;
    return(true);
  }

  this.setScale = function(scale) {
    if (this.scale == scale)
      return(false);

    clearById('zoom' + this.scale);
    this.scale = scale;
    this.nav2d.setSize(1/scale);
    this.image2d.setSize(scale);
    highlightById('zoom' + scale);
    return(true);
  }

  this.setSamples = function(samples) {
    this.samples = samples;
    if (samples.length > 1) {
      var s = 'Samples: ';
      for (var kk=0; kk < samples.length; kk++) {
        var sample = samples[kk];
        var name = sample;
        if (this.aliases != null)
          name = this.aliases[kk];
        s = s + '[<span id="sample' + sample + '"><a href="javascript:changeSample(\'' + sample + '\');">' + name + '</a></span>]<span style="font-size:1%"> </span>';
      }
      s = s + ' ';
      updateLabel('samplesLabel', s);
    }
  }

  this.setAliases = function(aliases) {
    this.aliases = aliases;
  }

  this.setChipTypes = function(chipTypes) {
    this.chipTypes = chipTypes;

    if (chipTypes.length > 1) {
      var s = 'Chip types: ';
      for (var kk=0; kk < chipTypes.length; kk++) {
        var chipType = chipTypes[kk];
        s = s + '[<span id="chipType' + chipType + '"><a href="javascript:changeChipType(\'' + chipType + '\');">' + chipType + '</a></span>]'; 
      }
      s = s + '<br>';
      updateLabel('chipTypeLabel', s);
    }
  }


  this.setScales = function(scales) {
    function padWidthZeros(x, width) {
      var str = "" + x;
      while (width - str.length > 0)
        str = "0" + str;
      return(str);
    }
   
    this.scales = scales;
    var zWidth = Math.round(Math.log(Math.max(scales)) / Math.log(10) + 0.5);
    var s = 'Zoom: ';
    for (var kk=0; kk < scales.length; kk++) {
      var scale = scales[kk];
      s = s + '[<span id="zoom' + scale + '"><a href="javascript:changeZoom(' + scale + ');">x' + padWidthZeros(scale, zWidth) + '</a></span>]'; 
    }
    s = s + '<br>';
    updateLabel('zoomLabel', s);
  }

  this.setColorMaps = function(colorMaps) {
    this.colorMaps = colorMaps;

    if (colorMaps.length > 1) {
      var s = 'Color map: ';
      for (var kk=0; kk < colorMaps.length; kk++) {
        var colorMap = colorMaps[kk];
        s = s + '[<span id="colorMap' + colorMap + '"><a href="javascript:changeColorMap(\'' + colorMap + '\');">' + colorMap + '</a></span>]'; 
      }
      s = s + '<br>';
      updateLabel('colorMapLabel', s);
    }
  }

  this.samples = new Array();
  this.aliases = null;
  this.chipTypes = new Array();
  this.colorMaps = new Array();
  this.scales = new Array();

  this.loadCount = 0;
  this.scale = 1;
  this.sample = '';
  this.chipType = '';
  this.colorMap = '';

  this.setupEventHandlers = function() {
    var owner = this;

    this.nav2d.onLoad = function() {
      owner.decreaseLoadCount();
    }

    this.image2d.onLoad = function() {
      owner.decreaseLoadCount();
    }

    this.nav2d.onScroll = function() {
	  	owner.image2d.setRelXY(this.x, this.y);
      var pos = owner.image2d.getXY();
      updateLabel('nav2dInfo', '('+Math.round(pos.x)+','+Math.round(pos.y)+')');
		  owner.image2d.update();
      return(false);    
	  }

    this.image2d.onScroll = function() {
      var pos = this.getXY();
      updateLabel('nav2dInfo', '('+Math.round(pos.x)+','+Math.round(pos.y)+')');
      owner.nav2d.setRelXY(this.x, this.y);
      owner.nav2d.update();
      return(false);    
	  }
  }

  this.update = function() {
    this.updateImage();
  }

  this.onLoad = function() { }

  this.start = function() {
    this.nav2d = new Scrollbar2d("nav2d");
    this.image2d = new ScrollImage2d("panel");
    this.setupEventHandlers();

    /* Default settings */
    this.setScales(new Array('1', '2', '4', '8', '16', '32'));
    this.setColorMaps(new Array('gray'));

    this.onLoad();

    this.setSample(this.samples[0]);
    this.setScale(this.scales[0]);
    this.setChipType(this.chipTypes[0]);
    this.setColorMap(this.colorMaps[0]);

    this.update();
  }
} /* ArrayExplorer() */

/****************************************************************
 HISTORY:
 2007-01-27
 o Created.
 ****************************************************************/
