Clazz.declarePackage ("com.sparshui.common");
Clazz.load ([], "com.sparshui.common.Location", ["JU.V3"], function () {
var c$ = Clazz.decorateAsClass (function () {
this._x = 0;
this._y = 0;
Clazz.instantialize (this, arguments);
}, com.sparshui.common, "Location", null, java.io.Serializable);
Clazz.makeConstructor (c$, 
function () {
this._x = 0;
this._y = 0;
});
Clazz.makeConstructor (c$, 
function (x, y) {
this._x = x;
this._y = y;
}, "~N,~N");
Clazz.defineMethod (c$, "getX", 
function () {
return this._x;
});
Clazz.defineMethod (c$, "getY", 
function () {
return this._y;
});
Clazz.overrideMethod (c$, "toString", 
function () {
return "x = " + this._x + ", y = " + this._y + (this._x < 1 && this._x > 0 ? "(" + com.sparshui.common.Location.pixelLocation (this).getX () + " " + com.sparshui.common.Location.pixelLocation (this).getY () + ")" : "");
});
Clazz.defineMethod (c$, "getDistance", 
function (location) {
var dx;
var dy;
return Math.sqrt ((dx = this._x - location._x) * dx + (dy = this._y - location._y) * dy);
}, "com.sparshui.common.Location");
Clazz.defineMethod (c$, "getVector", 
function (location) {
return JU.V3.new3 (location._x - this._x, location._y - this._y, 0);
}, "com.sparshui.common.Location");
c$.getCenter = Clazz.defineMethod (c$, "getCenter", 
function (a, b) {
return com.sparshui.common.Location.getCentroid (a, b, 0.5);
}, "com.sparshui.common.Location,com.sparshui.common.Location");
c$.getCentroid = Clazz.defineMethod (c$, "getCentroid", 
function (a, b, w) {
var w1 = 1 - w;
return  new com.sparshui.common.Location (a._x * w1 + b._x * w, a._y * w1 + b._y * w);
}, "com.sparshui.common.Location,com.sparshui.common.Location,~N");
c$.pixelLocation = Clazz.defineMethod (c$, "pixelLocation", 
function (location) {
return location == null ? null :  new com.sparshui.common.Location (location.getX () * com.sparshui.common.Location.screenDim.width, location.getY () * com.sparshui.common.Location.screenDim.height);
}, "com.sparshui.common.Location");
c$.screenLocation = Clazz.defineMethod (c$, "screenLocation", 
function (location) {
return (location == null ? null :  new com.sparshui.common.Location (location.getX () / com.sparshui.common.Location.screenDim.width, location.getY () / com.sparshui.common.Location.screenDim.height));
}, "com.sparshui.common.Location");
c$.screenDim = c$.prototype.screenDim = java.awt.Toolkit.getDefaultToolkit ().getScreenSize ();
});
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
