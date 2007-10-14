var ExplorerSettings = Class.create({
  initialize: function(args) {
    /* Default values */
    this.args = new Hash({
      sample: "NA06985,B5",
      chromosome: 22, 
      zoom: 1
    });

    /* Update by constructor arguments */
    this.importArray(args);
  },

  getArray: function() {
    return this.args;
  },

  get: function(key) {
    if (typeof(key) == "undefined")
      throw "Invalid key: " + key;
    return this.args[key];
  },

  set: function(key, value) {
    if (typeof(key) == "undefined")
      throw "Invalid key: " + key;
    if (typeof(value) == "undefined")
      return false;
    this.args[key] = value;
  },

  /* Update arguments by another Array */
  importArray: function(args) {
    for (key in args) {
      this.set(key, args[key]);
    };
  },

  /* Update arguments by cookies */
  importCookies: function() {
    var that = this;
    this.args.each(function(pair) {
      var value = $jq.cookie(pair.key);
      if (value != null) {
        that.set(pair.key, value);
      }
    });
  },

  /* Update arguments by URL parameters */
  importUrlParameters: function() {
    /* Update arguments by URL parameters */
    var that = this;

    var urlParams = $jq.getUrlParameters();
    this.args.each(function(pair) {
      var value = urlParams[pair.key];
      that.set(pair.key, value);
    });
  },

  load: function() {
    this.importCookies();
    this.importUrlParameters();
  },

  /* Store arguments as cookies */
  save: function() {
    this.args.each(function(pair) {
      $jq.cookie(pair.key, pair.value);
    });
  },

  toQueryString: function() {
    return this.args.toQueryString();
  }
});


var ChromosomeExplorerSettings = Class.create(ExplorerSettings, {
  getImageFilename: function() {
    var fmt = "%s,chr%02d,x%04d.png";
    var a = this.args;
    var filename = sprintf(fmt, a['sample'], a['chromosome'], a['zoom']);
    return filename;
  },

  getImagePathname: function() {
		var path = "imgs";
    return path + "/" + this.getImageFilename();
  },

  getLinkTo: function() {
		var l = window.location;
  	var url = l.protocol + '//' + l.host + l.pathname;
    return url + '?' + this.toQueryString();
  }
});

