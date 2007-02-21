/****************************************************************
 * ChromosomeExplorer()
 *
 * Author: Henrik Bengtsson, hb@stat.berkeley.edu
 ****************************************************************/
function ChromosomeExplorer() {
  this.loadCount = 0;
  this.imageUrl = null;
  this.bookmarkUrl = null;
  this.cnrUrl = null;

  this.scale = -1;

  this.chromosome = null;
  this.chromosomeIdx = 22;
  this.sample = null;
  this.sampleIdx = 0;
  this.chipType = null;
  this.chipTypeIdx = 0;


  this.showIndicator = function(state) {
    var statusImage = document.getElementById('statusImage');
    if (state) {
      statusImage.style.visibility = 'visible';
    } else {
      statusImage.style.visibility = 'hidden';
    }
  }

  this.setStatus = function(state) {
    navImage = document.getElementById('navigatorImage');
    panelImage = document.getElementById('panelImage');
    if (state == "") {
      this.showIndicator(false);
      navImage.style.filter = "alpha(opacity=50)";
      navImage.style.opacity = 0.50;
      panelImage.style.filter = "alpha(opacity=100)";
      panelImage.style.opacity = 1.0;
      this.updateInfo();
    } else if (state == "wait") {
      this.showIndicator(true);
      navImage.style.filter = "alpha(opacity=20)";
      navImage.style.opacity = 0.20;
      panelImage.style.filter = "alpha(opacity=50)";
      panelImage.style.opacity = 0.50;
    }
  }

  this.setChipType = function(idx) {
    if (this.chipTypeIdx != idx) {
      clearById('chipType' + this.chipType);
      this.loadCount = 2;
      this.setStatus('wait');
      this.chipTypeIdx = idx;
      this.chipType = chipTypes[this.chipTypeIdx];
      highlightById('chipType' + this.chipType);
      this.updatePanel();
      this.updateNavigator();
    }
  }

  this.setChromosome = function(idx) {
    if (this.chromosomeIdx != idx) {
      clearById('chromosome' + this.chromosome);
      this.loadCount = 2;
      this.setStatus('wait');
      this.chromosomeIdx = idx;
      this.chromosome = chromosomes[this.chromosomeIdx];
      highlightById('chromosome' + this.chromosome);
      this.updatePanel();
      this.updateNavigator();
    }
  }

  this.setZoom = function(idx) {
    zoomIdx = idx;
    s = zooms[idx];
    if (this.scale != s) {
      clearById('zoom' + this.scale);
      this.loadCount = 1;
      this.setStatus('wait');
      this.scale = s;
      highlightById('zoom' + this.scale);
      this.updatePanel();
    }
  }

  this.setSample = function(idx) {
    if (this.sampleIdx != idx) {
      clearById('sample' + this.sample);
      this.loadCount = 2;
      this.setStatus('wait');
      this.sampleIdx = idx;
      this.sample = samples[this.sampleIdx];
      highlightById('sample' + this.sample);
      this.updatePanel();
      this.updateNavigator();
    }
  }
  
  this.updateInfo = function() {
    updateLabel('chromosomeLabel', chromosomes[this.chromosomeIdx]);
    var label = samples[this.sampleIdx];
    if (sampleLabels != null) {
      if (sampleLabels[this.sampleIdx] != label) {
        label = sampleLabels[this.sampleIdx] + ' (' + label + ')';
      }
    }
    updateLabel('sampleLabel', label);
  }

  this.getImagePathname = function(chipType, sample, chromosome, zoom) {
    imgName = sample + ",chr" + padWidthZeros(chromosome, 2) + ",x" + padWidthZeros(zoom, 4) + ".png";
    var pathname = chipType + '/glad/' + imgName;
    return(pathname);
  }

  this.resetPositions = function() {
    navAreaX = 0;
    panel.scrollLeft = 0;
    this.updateGlobals();
  }

  this.updateGlobals = function() {
    panelWidth = panel.clientWidth;
    panelMaxWidth = panel.scrollWidth;
    navImageOffsetX = findXY(navImage).x;
    panelImageWidth = panelImage.clientWidth;
    panelImageOffsetX = findXY(panelImage).x;
    navImageWidth = navImage.clientWidth;
    navAreaWidth = navImageWidth * (panelWidth / panelMaxWidth);
  }

  this.locatorUpdated = function() {
    /* Update locator tag */
    var pixelsPerMb = 3; /* /1.0014; */
    var xPx = findXY(panelLocator).x - findXY(panelImage).x + parseFloat(panel.scrollLeft);
    var xMb = (xPx-50)/(this.scale*pixelsPerMb);
    var tag = Math.round(100*xMb)/100 + 'Mb';
    updateText(panelLocatorTag, tag);
  
    var url;
  
    /* Update shortcut link */
    if (this.bookmarkUrl != null) {
      var args = "'" + this.chipType + "', '" + this.sample + "', '" + this.chromosome + "', " + panel.scrollLeft + ", " + this.scale;
      url = 'javascript:explorer.jumpTo(' + args + ');';
      url = 'x:"' + args + '",';
    	/* url = 'javascript:addToFavorites("' + url + '", "sss")';	*/
      this.bookmarkUrl.href = url;
      updateText(this.bookmarkUrl, url);
    }
  
    /* Update CNR link */
    if (this.cnrUrl != null) {
      url = this.chipType + '/glad/' + 'regions.xls';
      this.cnrUrl.href = url;
      updateText(this.cnrUrl, url);
    }
  }

  this.jumpTo = function(newChipType, newSample, newChromosome, newPanelOffset, newZoom) {
    /* Chip type */
    if (this.chipType != '') {
      clearById('chipType' + this.chipType);
      this.chipTypeIdx = chipTypes.indexOf(newChipType);
      this.chipType = chipTypes[this.chipTypeIdx];
      highlightById('chipType' + this.chipType);
    }
  
    /* Sample */
    clearById('sample' + this.sample);
    this.sampleIdx = samples.indexOf(newSample);
    this.sample = samples[this.sampleIdx];
    highlightById('sample' + this.sample);
  
    /* Chromosome */
    clearById('chromosome' + this.chromosome);
    this.chromosomeIdx = chromosomes.indexOf(newChromosome);
    this.chromosome = chromosomes[this.chromosomeIdx];
    highlightById('chromosome' + this.chromosome);
  
    /* Zoom */
    clearById('zoom' + this.scale);
    zoomIdx = -1;
    var kk = 0;
    while (zoomIdx == -1 && kk < zooms.length) {
      if (zooms[kk] == newZoom)
        zoomIdx = kk;
      kk = kk + 1;
    }
    this.scale = zooms[zoomIdx];
    highlightById('zoom' + this.scale);
  
    /* When image is loaded... */
    panelImageOnLoad = function() {
      panel.scrollLeft = newPanelOffset;
      this.panelUpdated();
    }
  
    this.updatePanel();
    this.updateNavigator();
  }

  this.setGlobalCursor = function(status) {
    panel.style.cursor = status;
    panelImage.style.cursor = status;
    nav.style.cursor = status;
    navImage.style.cursor = status;
    navArea.style.cursor = status;
  }
  

  this.setupEventHandlers = function() {
    var owner = this;

    /*******************************************************
     * chromosomePanel
     *******************************************************/
    panel = document.getElementById('panel');
    panelLocator = document.getElementById('panelLocator');
    panelLocatorTag = document.getElementById('panelLocatorTag');
  
    panelOnScroll = function() {
      relOffset = panel.scrollLeft / panelMaxWidth;
      navAreaX = relOffset * navImageWidth;
      owner.navAreaUpdate();
    }
  
    panel.onscroll = panelOnScroll;
  
    panelImage = document.getElementById('panelImage');
  
    var panelLocatorIsLocked = true;
  
    panelImage.ondblclick = function() {
      panelLocatorIsLocked = false;
      panelImage.onmousedown = null;
      var e = arguments[0] || event;
      mouseX = e.clientX;
      panelLocator.style.left = (mouseX-2) + "px";
      owner.locatorUpdated();
    }
    panelLocator.ondblclick = panelImage.ondblclick;
  
    panelImage.onclick = function() {
      panelLocatorIsLocked = true;
      panelImage.onmousedown = panelImageOnMouseDown;
    }
    panelLocator.onclick = panelImage.onclick;
  
    panelImageOnMouseMove = function() {
      if (!panelLocatorIsLocked) {
        var e = arguments[0] || event;
        mouseX = e.clientX;
        panelLocator.style.left = (mouseX-2) + "px";
        owner.locatorUpdated();
      }
      return false;
    }
  
    panelImage.onmousemove = panelImageOnMouseMove;
  
    panelImageOnMouseDown = function() {
      var e = arguments[0] || event;
      var x = panel.scrollLeft + e.clientX;
      owner.setGlobalCursor("move");
      panel.onscroll = null;
      panelImage.onmousemove = null;
  
      document.onmousemove = function() {
        var e = arguments[0] || event;
        isMoving = true;
        panel.scrollLeft = x - e.clientX;
        owner.panelUpdated();
        return false;
      }
  
      document.onmouseup = function() {
        document.onmousemove = null;
        panel.onscroll = panelOnScroll;
        panelImage.onmousemove = panelImageOnMouseMove;
  //      panelLocatorIsLocked = false;
        owner.setGlobalCursor("default");
        return false;
      }
  
      return false;
    }
  
    panelImage.onmousedown = panelImageOnMouseDown;

    /*******************************************************
     * chromosomeNavigator
     *******************************************************/
    nav = document.getElementById('navigator');
    navImage = document.getElementById('navigatorImage');
  
    navArea = document.getElementById('navigatorArea');
    var mouseX = 0;
    var mouseDown = false;
  
    /* Immitate onmousepress, which does not exists */
    navImage.onmousepress = function() {
      if (mouseDown) {
        if (mouseX < navAreaX) {
          owner.navAreaMove(navAreaX - 0.47*navAreaWidth);
        } else if (mouseX > navAreaX + 1*navAreaWidth) {
          owner.navAreaMove(navAreaX + 1.47*navAreaWidth);
        } else {
          owner.navAreaMove(mouseX);
        }
        setTimeout('navImage.onmousepress();', 100);
      }
      return false;
    }
  
    navImage.onmousedown = function() {
      var e = arguments[0] || event;
      mouseDown = true;
      owner.updateGlobals();
      mouseX = (e.clientX - navImageOffsetX);
      if (mouseX < navAreaX) {
        owner.navAreaMove(navAreaX - 0.47*navAreaWidth);
        } else if (mouseX > navAreaX + 1*navAreaWidth) {
        owner.navAreaMove(navAreaX + 1.47*navAreaWidth);
      }
  
      setTimeout('navImage.onmousepress();', 500);
  
      document.onmouseup = function() {
        navImage.onmousemove = null;
        mouseDown = false;
        return false;
      }
  
      navImage.onmousemove = function() {
        var e = arguments[0] || event;
        mouseX = (e.clientX - navImageOffsetX);
        return false;
      }
  
      return false;
    }
  
    navArea.onmousedown = function() {
      var e = arguments[0] || event;
      owner.setGlobalCursor("move");
      panel.onscroll = null;
      owner.updateGlobals();
      mouseX = (e.clientX - navImageOffsetX);
      var dx = navAreaWidth/2 + (navAreaX - mouseX);
      owner.navAreaMove(mouseX + dx);
  
      document.onmousemove = function() {
        var e = arguments[0] || event;
        mouseX = (e.clientX - navImageOffsetX);
        owner.navAreaMove(mouseX + dx);
        return false;
      }
  
      document.onmouseup = function() {
        document.onmousemove = null;
        panel.onscroll = panelOnScroll;
        owner.setGlobalCursor("default");
        return false;
      }
      return false;
    }
  }

  this.start = function() {
    /*******************************************************
     * Set the current sample
     *******************************************************/
    this.sample = samples[this.sampleIdx];
    highlightById('sample' + this.sample);
  
    this.chipType = chipTypes[this.chipTypeIdx];
    highlightById('chipType' + this.chipType);
  
    if (this.chromosome == null)
      this.chromosome = chromosomes[this.chromosomeIdx];
    highlightById('chromosome' + this.chromosome);
  
    this.scale = zooms[zoomIdx];
    highlightById('zoom' + this.scale);
  
    if (navigatorZoom == -1) {
      navigatorZoom = this.scale;
    }
  
    this.imageUrl = document.getElementById('imageUrl');
    this.bookmarkUrl = document.getElementById('bookmarkUrl');
    this.cnrUrl = document.getElementById('cnrUrl');

		this.setupEventHandlers();
  
    this.updateNavigator();
    this.updatePanel();
    this.setStatus('');
    webcutsOptions['numberLinks'] = false;
    setTimeout('explorer.navAreaMoveRel(0.5);', 1000);
  }
    
  this.navAreaUpdate = function() {
    if (navAreaX < 0) {
      navAreaX = 0;
    } else if (navAreaX + navAreaWidth > navImageWidth) {
      navAreaX = navImageWidth - navAreaWidth;
    }
    navArea.style.width = navAreaWidth + "px";
    navArea.style.left = (navAreaX + navImageOffsetX) + "px";
    this.locatorUpdated();
  }

  this.navAreaMove = function(midX) {
    navAreaX = midX - navAreaWidth/2;
    this.navAreaUpdate();
    this.panelMove(navAreaX/navImageWidth);
  }

  this.navAreaMoveRel = function(relX) {
    this.navAreaMove(relX * navImageWidth);
  }

  this.panelMove = function(relOffset) {
    panelX = relOffset*panelMaxWidth;
    if (panelX < 0)
      panelX = 0;
    panel.scrollLeft = panelX;
  }

  this.panelUpdated = function() {
    this.updateGlobals();
    relOffset = panel.scrollLeft / panelMaxWidth;
    navAreaX = relOffset * navImageWidth;
    this.navAreaUpdate();
    this.locatorUpdated();
  }

  this.updatePanel = function() {
    var owner = this;

    var navAreaRelMidX = (navAreaX + navAreaWidth/2) / navImageWidth;
    var pathname = owner.getImagePathname(this.chipType, this.sample, this.chromosomeIdx+1, this.scale);
    panelImage = document.getElementById('panelImage');
    panelImage.onload = function() {
      owner.updateNavigatorWidth();
      owner.updateGlobals();
      owner.navAreaMoveRel(navAreaRelMidX);
      owner.loadCount = owner.loadCount - 1;
      if (owner.loadCount <= 0) {
        owner.loadCount = 0;
        owner.setStatus("");
      }
      panelImageOnLoad();
      panelImageOnLoad = function() {};
    }
    panelImage.src = pathname;
    this.imageUrl.href = pathname;
    updateText(this.imageUrl, pathname);
  
    /* Update the title of the page */
    var title = location.href;
    title = title.substring(0, title.lastIndexOf('\/'));
    title = title.substring(title.lastIndexOf('\/')+1);
    title = title + '/' + pathname;
    document.title = title;
  }
  
  this.updateNavigator = function() {
    var owner = this;

    var pathname = this.getImagePathname(this.chipType, this.sample, this.chromosomeIdx+1, navigatorZoom);
    navImage = document.getElementById("navigatorImage");
    navImage.onload = function() {
      owner.loadCount = owner.loadCount - 1;
      if (owner.loadCount <= 0) {
        owner.loadCount = 0;
        owner.setStatus("");
      }
    }
    navImage.src = pathname;
  } // updateNavigator()


  this.updateNavigatorWidth = function() {
    /* Update the width of the navigator */
    var chromosomeLength = new Array(3840, 3798, 3119, 2993, 2826, 2673, 2482, 2288, 2165, 2117, 2104, 2071, 1785, 1664, 1568, 1387, 1230, 1191, 998, 976, 733, 774, 2417);
    var relWidth = chromosomeLength[this.chromosomeIdx] / chromosomeLength[0];
    navImageWidth = Math.round(relWidth * nav.clientWidth);
    navAreaWidth = Math.round(relWidth * nav.clientWidth);
    navImage.style.width = "" + navImageWidth  + "px";
  }
  
  this.getMouseMb = function(x, chromosome, zoom) {
    return(-1);
  }
         
  this.onLoad = function() { }
} /* ChromosomeExplorer */


/****************************************************************
 HISTORY:
 2007-02-20
 o Updated to <rootPath>/<dataSet>/<tags>/<chipType>/<set>/.
 o Created from old ChromosomeExplorer.js making it more of the
   style of class ArrayExplorer.
 ****************************************************************/
