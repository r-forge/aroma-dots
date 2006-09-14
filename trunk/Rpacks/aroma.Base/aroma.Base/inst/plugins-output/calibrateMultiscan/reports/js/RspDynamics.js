/*************************************************************************
Title      : RspDynamics.js

Description: 

Version    : 0.1
Language   : Javascript1.2
Author     : Henrik Bengtsson, hb@maths.lth.se
Date       : September 2005
URL        : 

References:

*************************************************************************/

var toggler = new FloatingLayer('RspTogglerWidget');
toggler.y = 20;
toggler.dy = 20;

updatePos = function() {
  toggler.update();
}

if (document.all) {
  window.onscroll = updatePos;
} else {
  setInterval('updatePos()', 500);
}



function toggleStitchy() {
  var id = "stitchy.toggler";
  var element = document.getElementById(id);
  if (element != null) {
    var isOn = toggler.isStitchy;
    if (isOn) {
      element.style.textDecoration = "underline";
      element.style.color = "black";
    } else {
      element.style.textDecoration = "none";
      element.style.color = "#666666";
    }
  }
}

function expandCollapse(idprefix) {
  var id = idprefix + ".toggler";
  var element = document.getElementById(id);
  var isOn = true;

  if (element != null) {
    isOn = (element.style.textDecoration == "none");
    if (isOn) {
      element.style.textDecoration = "underline";
      element.style.color = "black";
    } else {
      element.style.textDecoration = "none";
      element.style.color = "#666666";
    }
  }

  var id = idprefix;
  var element = document.getElementById(id);
  if (element != null) {
    element.style.display = isOn ? "block" : "none";
  }

  for (var kk=0; kk <= 100; kk++) {
    id = idprefix + "." + kk;
    element = document.getElementById(id);
    if (element != null) {
      element.style.display = isOn ? "block" : "none";
    }
  }
}


/*************************************************************************
HISTORY:
2005-09-04 [v0.1]
o Added collapseExpand().  
  Adopted from http://www.blakems.com/archives/000087.html
*************************************************************************/
