/****************************************************************
 *
 ****************************************************************/
function findXY(obj) {
  var x = 0;
  var y = 0;

  if (obj.offsetParent)  {
    while (obj.offsetParent) {
      x += parseFloat(obj.offsetLeft);
      y += parseFloat(obj.offsetTop);
      obj = obj.offsetParent;
    }
  }  else if (obj.x) {
    x += obj.x;
    y += obj.y;
  }

  pos = new Object();
  pos.x = x;
  pos.y = y;
  return pos;
} /* findXY() */


function updateText(obj, str) {
  obj.innerText = str;
  obj.innerHTML = str;
} /* updateText() */


function updateLabel(id, text) {
  var obj = document.getElementById(id);
  DOM_setInnerText(obj, text); /* From Webcuts.js */
}

// Array.indexOf( value, begin, strict ) - Return index of the first 
// element that matches value.
// Source: http://4umi.com/web/javascript/array.htm#indexof
Array.prototype.indexOf = function(v, b, s) {
  for(var kk=+b || 0, ll=this.length; kk < ll; kk++) {
    if(this[kk] === v || s && this[kk] == v)
      return(kk);
  }
  return(-1);
};


function clearById(id) {
  var obj = document.getElementById(id);
  if (obj != null) {
/*  
    obj.style.visibility = 'hidden'; 
    obj.style.borderBottom = 'none';
*/
    obj.style.background = 'none';
  }
} /* clearById() */


function highlightById(id) {
  var obj = document.getElementById(id);
  if (obj != null) {
/*  
    obj.style.visibility = 'visible'; 
    obj.style.borderBottom = '2px solid black';
*/
    obj.style.background = '#ccccff;';
  }
} /* highlightById() */


/****************************************************************
 HISTORY:
 2007-01-27
 o Created.
 ****************************************************************/
