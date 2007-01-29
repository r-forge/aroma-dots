/****************************************************************
 * ScrollImage2d()
 * Scrollbar2d()
 *
 * Author: Henrik Bengtson, hb@stat.berkeley.edu
 ****************************************************************/
function ScrollImage2d(id) {
  this.setImage = function(url) {
    var owner = this;

    /* Image loader to get the size of the (non-rescaled) image */
    var myImage = new Image();
    myImage.onload = function() {
      owner.imageWidth = this.width;
      owner.imageHeight = this.height;
    }
    myImage.src = url;

    /* Define onload() function */
    this.image.onload = function() {
      owner.onLoad();
      owner.update();
    }

    /* Define onload() function */
    this.image.onerror = function() {
      alert('Image not loaded: ' + url);
    }

    /* Start loading image */
    this.image.src = url;
  }

  this.getImageDimension = function() {
    var res = Object();
    res.width = this.imageWidth;
    res.height = this.imageHeight;
    return(res);
  }

  this.getImageWidth = function() {
    return(this.image.width);
  }

  this.getImageHeight = function() {
  	return(this.image.height);
	}

  this.getAspectRatio = function() {
    return (this.container.clientWidth / this.container.clientHeight);
  }

  this.update = function() {
    this.container.scrollLeft = this.x * this.container.scrollWidth;
    this.container.scrollTop = this.y * this.container.scrollHeight;
    var pos = findXY(this.image);
    this.xOffset = pos.x;
    this.yOffset = pos.y;
    var w = this.getImageWidth();
    var h = this.getImageHeight();
    this.image.style.left = Math.round(this.xOffset + w*this.x) + "px";
    this.image.style.top = Math.round(this.yOffset + h*this.y) + "px";
    this.image.width = this.width * this.container.clientWidth;
    var dim = this.getImageDimension();
    var aspect = this.width / this.height;
    this.image.height = aspect * this.height * this.container.clientWidth;
		/*
    this.container.style.height = winAspect*this.container.clientWidth;
		*/
  }

  this.getXY = function() {
    var w = this.getImageWidth()/this.width;
    var h = this.getImageHeight()/this.height;
    var res = Object();
    res.x = Math.round(w*this.x);
    res.y = Math.round(w*this.y);
    return(res);
  }

  this.getDimension = function() {
    var res = Object();
    var dim = this.getImageDimension();
    res.width = Math.round(dim.width);
    res.height = Math.round(dim.height);
    return(res);
  }

  this.getRegion = function() {
    var xy = this.getXY();
    var dim = this.getDimension();
    var res = Object();
    res.x0 = xy.x;
    res.y0 = xy.y;
    res.x1 = xy.x + dim.width;
    res.y1 = xy.y + dim.height;
    return(res);
  }

  this.setRelXY = function(x,y) {
    x = Math.max(0, x);
    x = Math.min(x, 1);
    y = Math.max(0, y);
    y = Math.min(y, 1);
    this.x = x;
    this.y = y;
  }

  this.setSize = function(scale) {
    this.setRelDimension(scale, scale);
    this.update();
  }

  this.setRelDimension = function(width, height) {
    width = Math.max(0, width);
    height = Math.max(0, height);
    this.width = width;
    this.height = height;
  }

  this.setCursor = function(status) {
    this.container.style.cursor = status;
    this.image.style.cursor = status;
  }

  this.setupEventHandlers = function() {
    var owner = this;

    this.image.onmousedown = function() {
      var e = arguments[0] || event;
      var x0 = owner.x;
      var y0 = owner.y;
      var w = owner.getImageWidth();
      var h = owner.getImageHeight();
      var xStart = e.clientX;
      var yStart = e.clientY;
      var x = owner.container.scrollLeft + e.clientX;
      var y = owner.container.scrollTop + e.clientY;
      owner.setCursor("move");
      owner.image.onmousemove = null;
      owner.onScrollBegin();
  
      document.onmousemove = function() {
        var e = arguments[0] || event;
        var dx = (e.clientX - xStart)/w;
        var dy = (e.clientY - yStart)/h;
        owner.setRelXY(x0-dx, y0-dy);
        owner.update();
        owner.onScroll();
        return false;
      }
  
      document.onmouseup = function() {
        document.onmousemove = null;
        owner.image.onmousemove = null;
        owner.setCursor("default");
        owner.onScrollEnd();
        return false;
      }
  
      return false;
    }
  }

  this.onLoad = function() {}
  this.onScrollBegin = function() {}
  this.onScroll = function() {}
  this.onScrollEnd = function() {}


  /* Initialize */
  this.scale = 1;
  this.x = 0;
  this.y = 0;
  this.xOffset = 0;
  this.yOffset = 0;
  this.width = 1;
  this.height = 1;
  
  this.container = document.getElementById(id);
  this.image = document.getElementById(id + 'Image');
  this.imageWidth = 0;
  this.imageHeight = 0;

  this.setupEventHandlers();
} /* ScrollImage2d() */




function Scrollbar2d(id) {
  this.setImage = function(url) {
    var owner = this;

    /* Define onload() function */
    this.image.onload = function() {
      owner.onLoad();
      owner.update();
    }

    /* Define onload() function */
    this.image.onerror = function() {
      alert('Image not loaded: ' + url);
    }

    /* Start loading image */
    this.image.src = url;
  }

  this.getImageWidth = function() {
    return(this.image.width-6);
  }

  this.getImageHeight = function() {
    return(this.image.height-6);
	}

  this.update = function() {
    var pos = findXY(this.image);
    this.xOffset = pos.x;
    this.yOffset = pos.y;
    var w = this.getImageWidth();
    var h = this.getImageHeight();
    this.marker.style.width = Math.round(w*this.width-0.5) + "px";
    this.marker.style.height = Math.round(h*this.height-0.5) + "px";
    this.marker.style.left = Math.round(this.xOffset + w*this.x) + "px";
    this.marker.style.top = Math.round(this.yOffset + h*this.y) + "px";
    this.marker.style.border = 'solid; black; 2px';
  }

  this.getRegion = function() {
    var res = Object();
    res.x0 = this.x;
    res.y0 = this.y;
    res.x1 = this.x + this.width;
    res.y1 = this.y + this.height;
    return res;
  }

  this.setRelXY = function(x,y) {
    x = Math.max(0, x);
    x = Math.min(x, 1-this.width);
    y = Math.max(0, y);
    y = Math.min(y, 1-this.height);
    this.x = x;
    this.y = y;
  }

  this.setSize = function(scale) {
    this.setRelDimension(scale, scale);
  }

  this.setRelDimension = function(width, height) {
    var xMid = this.x + this.width/2;
    var yMid = this.y + this.height/2;
    width = Math.max(0, width);
    height = Math.max(0, height);
    width = Math.min(width, 1);
    height = Math.min(height, 1);
    this.width = width;
    this.height = height;
    this.setRelXY(xMid - this.width/2, yMid - this.height/2);
  }

  this.setCursor = function(status) {
    this.marker.style.cursor = status;
    this.container.style.cursor = status;
    this.image.style.cursor = status;
  }

  this.setupEventHandlers = function() {
    var owner = this;

    var containerOnClick = function() {
      var w = owner.getImageWidth();
      var h = owner.getImageHeight();
      var e = arguments[0] || event;
      var pos = findXY(owner.marker);
      var dx = (e.clientX - pos.x)/w - owner.width;
  		var dy = (e.clientY - pos.y)/h - owner.height;
			owner.setRelXY(dx,dy);
      owner.update();
    }

		/*
    this.container.onclick = containerOnClick;
    */

    this.marker.onmousedown = function() {
      var x0 = owner.x;
      var y0 = owner.y;
      var w = owner.getImageWidth();
      var h = owner.getImageHeight();
      var e = arguments[0] || event;
      var xStart = e.clientX;
      var yStart = e.clientY;
      owner.setCursor('move');
      owner.onScrollBegin();
  
      document.onmousemove = function() {
        var e = arguments[0] || event;
        var dx = (e.clientX - xStart)/w;
        var dy = (e.clientY - yStart)/h;
        owner.setRelXY(x0+dx, y0+dy);
        owner.update();
        owner.onScroll();
        return false;
      }
  
      document.onmouseup = function() {
        document.onmousemove = null;
        owner.setCursor('default');
        owner.onScrollEnd();
        return false;
      }
  
      return false;
    }
  }

  this.onLoad = function() {}
  this.onScrollBegin = function() {}
  this.onScroll = function() {}
  this.onScrollEnd = function() {}


  /* Initialize */
  this.x = 0;
  this.y = 0;
  this.xOffset = 0;
  this.yOffset = 0;
  this.width = 1;
  this.height = 1;
  
  this.container = document.getElementById(id);
  this.image = document.getElementById(id + 'Image');
  this.marker = document.getElementById(id + 'Marker');

  this.setupEventHandlers();
} /* Scrollbar2d() */


/****************************************************************
 HISTORY:
 2007-01-27
 o Created.
 ****************************************************************/
