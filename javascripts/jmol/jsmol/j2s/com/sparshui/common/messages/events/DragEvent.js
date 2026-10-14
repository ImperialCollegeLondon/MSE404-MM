Clazz.declarePackage ("com.sparshui.common.messages.events");
Clazz.load (["com.sparshui.common.Event"], "com.sparshui.common.messages.events.DragEvent", ["com.sparshui.common.utils.Converter"], function () {
var c$ = Clazz.decorateAsClass (function () {
this._dx = 0;
this._dy = 0;
this._nPoints = 1;
this._time = 0;
Clazz.instantialize (this, arguments);
}, com.sparshui.common.messages.events, "DragEvent", null, com.sparshui.common.Event);
Clazz.makeConstructor (c$, 
function () {
this._dx = 0;
this._dy = 0;
});
Clazz.makeConstructor (c$, 
function (dx, dy, nPoints, time) {
this._dx = dx;
this._dy = dy;
this._nPoints = nPoints;
this._time = time;
}, "~N,~N,~N,~N");
Clazz.makeConstructor (c$, 
function (data) {
if (data.length < 17) {
System.err.println ("Error constructing Drag Event.");
this._dx = 0;
this._dy = 0;
} else {
this._dx = com.sparshui.common.utils.Converter.byteArrayToFloat (data, 0);
this._dy = com.sparshui.common.utils.Converter.byteArrayToFloat (data, 4);
this._nPoints = data[8];
this._time = com.sparshui.common.utils.Converter.byteArrayToLong (data, 9);
}}, "~A");
Clazz.defineMethod (c$, "getTime", 
function () {
return this._time;
});
Clazz.defineMethod (c$, "getNPoints", 
function () {
return this._nPoints;
});
Clazz.defineMethod (c$, "getDx", 
function () {
return this._dx;
});
Clazz.defineMethod (c$, "getDy", 
function () {
return this._dy;
});
Clazz.defineMethod (c$, "setDx", 
function (dx) {
this._dx = dx;
}, "~N");
Clazz.defineMethod (c$, "setDy", 
function (dy) {
this._dy = dy;
}, "~N");
Clazz.overrideMethod (c$, "getEventType", 
function () {
return 0;
});
Clazz.overrideMethod (c$, "toString", 
function () {
var ret = "Drag Event: dx = " + this._dx + ", dy = " + this._dy;
return ret;
});
Clazz.overrideMethod (c$, "serialize", 
function () {
var data =  Clazz.newByteArray (21, 0);
com.sparshui.common.utils.Converter.intToByteArray (data, 0, this.getEventType ());
com.sparshui.common.utils.Converter.floatToByteArray (data, 4, this._dx);
com.sparshui.common.utils.Converter.floatToByteArray (data, 8, this._dy);
data[12] = this._nPoints;
com.sparshui.common.utils.Converter.longToByteArray (data, 13, this._time);
return data;
});
});
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
