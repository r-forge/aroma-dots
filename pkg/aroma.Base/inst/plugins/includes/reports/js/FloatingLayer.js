/*************************************************************************
Title      : FloatingLayer.js

Description: 

Version    : 0.1
Language   : Javascript1.2
Author     : Henrik Bengtsson, hb@maths.lth.se
Date       : September 2005
URL        : 

References:

*************************************************************************/

function FloatingLayer(id) {
  this.id = id;
  this.isStitchy = true;
  this.y = -1;
  this.x = -1;
  this.dx = 0;
  this.dy = 0;
}

function FloatingLayer.prototype.setVisibility(status) {
  if (document.all) {
    element = document.all[this.id];
  } else if (document.layers) {
    element = document[this.id];
  } else if (document.getElementById) {
    element = document[this.id];
  }
  element.style.visibility = (status) ? "visible" : "hidden";
}

function FloatingLayer.prototype.getTop() {
  if (document.all) {
    element = document.all[this.id];
    return(element.style.pixelTop);
  } else if (document.layers) {
    element = document[this.id];
    return(element.top); 
  } else if (document.getElementById) {
    element = document[this.id];
    return(element.style.top);
  }
}

function FloatingLayer.prototype.getLeft() {
  if (document.all) {
    element = document.all[this.id];
    return(element.style.pixelLeft);
  } else if (document.layers) {
    element = document[this.id];
    return(element.left); 
  } else if (document.getElementById) {
    element = document[this.id];
    return(element.style.left);
  }
}


function FloatingLayer.prototype.setTop(top) {
  if (document.all) {
    element = document.all[this.id];
    element.style.pixelTop = top;
  } else if (document.layers) {
    element = document[this.id];
    element.top = top; 
  } else if (document.getElementById) {
    element = document[this.id];
    element.style.top = top;
  }
}

function FloatingLayer.prototype.setLeft(left) {
  if (document.all) {
    element = document.all[this.id];
    element.style.pixelLeft = left;
  } else if (document.layers) {
    element = document[this.id];
    element.left = left; 
  } else if (document.getElementById) {
    element = document[this.id];
    element.style.left = left;
  }
}

function FloatingLayer.prototype.setRight(right) {
  if (document.all) {
    element = document.all[this.id];
    element.style.pixelRight = right;
  } else if (document.layers) {
    element = document[this.id];
    element.right = right; 
  } else if (document.getElementById) {
    element = document[this.id];
    element.style.right = right;
  }
}

function FloatingLayer.prototype.setY(y) {
  if (y >= 0) {
    this.setTop(y);
  } else {
    this.setBottom(y);
  }
}

function FloatingLayer.prototype.setX(x) {
  if (x >= 0) {
    this.setLeft(x);
  } else {
    this.setRight(x);
  }
}

function FloatingLayer.prototype.setTopOffset(offset) {
  this.dy = offset;
}

function FloatingLayer.prototype.setBottomOffset(offset) {
  this.dy = offset;
}

function FloatingLayer.prototype.setRightOffset(offset) {
  this.dx = offset;
}

function FloatingLayer.prototype.setLeftOffset(offset) {
  this.dx = offset;
}



function FloatingLayer.prototype.update() {
  var y = this.y;
  var x = this.x;

  if (this.isStitchy) {
    if (document.all) {
      var element = document.all[this.id];
      if (document.documentElement.scrollTop) {
        y = document.documentElement.scrollTop;
      } else {
        y = 0;
      }
      y = y + document.body.scrollTop;
    } else if (document.layers) {
      var element = document[this.id];
      y = window.pageYOffset;
      y = y + yOffset;
    } else if (document.getElementById) {
      var element = document[this.id];
      y = window.pageYOffset + 'px';
      y = y + yOffset;
    }
    y = y+this.dy;
    x = x+this.dx;
  }

  this.setY(y);
  //  this.setX(x);
}

function FloatingLayer.prototype.toggleStitchy() {
  this.isStitchy = !this.isStitchy;
}





/*************************************************************************
HISTORY:
2005-09-04 [v0.1]
o Created.
*************************************************************************/
