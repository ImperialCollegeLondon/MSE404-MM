Clazz.declarePackage ("com.sparshui.common.messages.events");
Clazz.load (["com.sparshui.common.Event"], "com.sparshui.common.messages.events.TouchEvent", ["com.sparshui.common.utils.Converter"], function () {
var c$ = Clazz.decorateAsClass (function () {
this._id = 0;
this._x = 0;
this._y = 0;
this._state = 0;
this._time = 0;
Clazz.instantialize (this, arguments);
}, com.sparshui.common.messages.events, "TouchEvent", null, com.sparshui.common.Event);
Clazz.makeConstructor (c$, 
function () {
this._id = 2147483647;
this._x = 0;
this._y = 0;
this._state = 0;
});
Clazz.makeConstructor (c$, 
function (id, x, y, state) {
this._id = id;
this._x = x;
this._y = y;
this._state = state;
this._time = System.currentTimeMillis ();
}, "~N,~N,~N,~N");
Clazz.makeConstructor (c$, 
function (tp) {
this._id = tp.getID ();
this._x = tp.getLocation ().getX ();
this._y = tp.getLocation ().getY ();
this._state = tp.getState ();
this._time = tp.getTime ();
}, "com.sparshui.server.TouchPoint");
Clazz.defineMethod (c$, "getTouchID", 
function () {
return this._id;
});
Clazz.defineMethod (c$, "getTime", 
function () {
return this._time;
});
Clazz.defineMethod (c$, "getX", 
function () {
return this._x;
});
Clazz.defineMethod (c$, "getY", 
function () {
return this._y;
});
Clazz.defineMethod (c$, "setX", 
function (x) {
this._x = x;
}, "~N");
Clazz.defineMethod (c$, "setY", 
function (y) {
this._y = y;
}, "~N");
Clazz.defineMethod (c$, "getState", 
function () {
return this._state;
});
Clazz.overrideMethod (c$, "getEventType", 
function () {
return 3;
});
Clazz.makeConstructor (c$, 
function (data) {
if (data.length < 24) {
System.err.println ("An error occurred while deserializing a TouchEvent.");
} else {
this._id = com.sparshui.common.utils.Converter.byteArrayToInt (data, 0);
this._x = com.sparshui.common.utils.Converter.byteArrayToFloat (data, 4);
this._y = com.sparshui.common.utils.Converter.byteArrayToFloat (data, 8);
this._state = com.sparshui.common.utils.Converter.byteArrayToInt (data, 12);
this._time = com.sparshui.common.utils.Converter.byteArrayToLong (data, 16);
}}, "~A");
Clazz.overrideMethod (c$, "serialize", 
function () {
var data =  Clazz.newByteArray (28, 0);
com.sparshui.common.utils.Converter.intToByteArray (data, 0, this.getEventType ());
com.sparshui.common.utils.Converter.intToByteArray (data, 4, this._id);
com.sparshui.common.utils.Converter.floatToByteArray (data, 8, this._x);
com.sparshui.common.utils.Converter.floatToByteArray (data, 12, this._y);
com.sparshui.common.utils.Converter.intToByteArray (data, 16, this._state);
com.sparshui.common.utils.Converter.longToByteArray (data, 20, this._time);
return data;
});
Clazz.overrideMethod (c$, "toString", 
function () {
return ("Touch Event: ID: " + this._id + ", X: " + this._x + ", Y: " + this._y + "State: " + this._state);
});
});
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
