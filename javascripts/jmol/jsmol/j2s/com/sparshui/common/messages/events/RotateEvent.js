Clazz.declarePackage ("com.sparshui.common.messages.events");
Clazz.load (["com.sparshui.common.Event"], "com.sparshui.common.messages.events.RotateEvent", ["com.sparshui.common.Location", "com.sparshui.common.utils.Converter"], function () {
var c$ = Clazz.decorateAsClass (function () {
this._rotation = 0;
this._center = null;
this._time = 0;
Clazz.instantialize (this, arguments);
}, com.sparshui.common.messages.events, "RotateEvent", null, com.sparshui.common.Event);
Clazz.makeConstructor (c$, 
function () {
this._rotation = 0;
this._center =  new com.sparshui.common.Location ();
});
Clazz.makeConstructor (c$, 
function (rotation, center, time) {
this._rotation = rotation;
this._center = center;
this._time = time;
}, "~N,com.sparshui.common.Location,~N");
Clazz.makeConstructor (c$, 
function (data) {
if (data.length < 12) {
System.err.println ("Error constructing Rotate Event.");
this._rotation = 0;
this._center =  new com.sparshui.common.Location (0, 0);
} else {
this._rotation = com.sparshui.common.utils.Converter.byteArrayToFloat (data, 0);
this._center =  new com.sparshui.common.Location (com.sparshui.common.utils.Converter.byteArrayToFloat (data, 4), com.sparshui.common.utils.Converter.byteArrayToFloat (data, 8));
}}, "~A");
Clazz.overrideMethod (c$, "getEventType", 
function () {
return 1;
});
Clazz.overrideMethod (c$, "serialize", 
function () {
var data =  Clazz.newByteArray (16, 0);
com.sparshui.common.utils.Converter.intToByteArray (data, 0, this.getEventType ());
com.sparshui.common.utils.Converter.floatToByteArray (data, 4, this._rotation);
com.sparshui.common.utils.Converter.floatToByteArray (data, 8, this._center.getX ());
com.sparshui.common.utils.Converter.floatToByteArray (data, 12, this._center.getY ());
return data;
});
Clazz.overrideMethod (c$, "toString", 
function () {
return ("Rotate Event - Rotation: " + this._rotation + ", Center: " + this._center.toString ());
});
Clazz.defineMethod (c$, "getRotation", 
function () {
return this._rotation;
});
Clazz.defineMethod (c$, "getTime", 
function () {
return this._time;
});
Clazz.defineMethod (c$, "getCenter", 
function () {
return this._center;
});
Clazz.defineMethod (c$, "setCenter", 
function (center) {
this._center = center;
}, "com.sparshui.common.Location");
Clazz.defineMethod (c$, "getX", 
function () {
return this._center.getX ();
});
Clazz.defineMethod (c$, "getY", 
function () {
return this._center.getY ();
});
});
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
