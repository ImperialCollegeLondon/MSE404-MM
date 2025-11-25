Clazz.declarePackage ("com.sparshui.common.messages.events");
Clazz.load (["com.sparshui.common.Event"], "com.sparshui.common.messages.events.SpinEvent", ["com.sparshui.common.utils.Converter"], function () {
var c$ = Clazz.decorateAsClass (function () {
this._rotationX = 0;
this._rotationY = 0;
this._rotationZ = 0;
Clazz.instantialize (this, arguments);
}, com.sparshui.common.messages.events, "SpinEvent", null, com.sparshui.common.Event);
Clazz.overrideMethod (c$, "getEventType", 
function () {
return 2;
});
Clazz.makeConstructor (c$, 
function () {
this._rotationX = 0.0;
this._rotationY = 0.0;
this._rotationZ = 0.0;
});
Clazz.makeConstructor (c$, 
function (rotationX, rotationY, rotationZ) {
this._rotationX = rotationX;
this._rotationY = rotationY;
this._rotationZ = rotationZ;
}, "~N,~N,~N");
Clazz.defineMethod (c$, "getRotationX", 
function () {
return this._rotationX;
});
Clazz.defineMethod (c$, "getRotationY", 
function () {
return this._rotationY;
});
Clazz.defineMethod (c$, "getRotationZ", 
function () {
return this._rotationZ;
});
Clazz.defineMethod (c$, "setRotationX", 
function (rotation) {
this._rotationX = rotation;
}, "~N");
Clazz.defineMethod (c$, "setRotationY", 
function (rotation) {
this._rotationY = rotation;
}, "~N");
Clazz.defineMethod (c$, "setRotationZ", 
function (rotation) {
this._rotationZ = rotation;
}, "~N");
Clazz.makeConstructor (c$, 
function (data) {
if (data.length < 12) {
System.err.println ("An error occurred while deserializing a TouchEvent.");
} else {
this._rotationX = com.sparshui.common.utils.Converter.byteArrayToFloat (data, 0);
this._rotationY = com.sparshui.common.utils.Converter.byteArrayToFloat (data, 4);
this._rotationY = com.sparshui.common.utils.Converter.byteArrayToFloat (data, 8);
}}, "~A");
Clazz.overrideMethod (c$, "serialize", 
function () {
var data =  Clazz.newByteArray (16, 0);
com.sparshui.common.utils.Converter.intToByteArray (data, 0, this.getEventType ());
com.sparshui.common.utils.Converter.floatToByteArray (data, 4, this._rotationX);
com.sparshui.common.utils.Converter.floatToByteArray (data, 8, this._rotationY);
com.sparshui.common.utils.Converter.floatToByteArray (data, 12, this._rotationZ);
return data;
});
Clazz.overrideMethod (c$, "toString", 
function () {
return ("Spin Event - rotationX: " + this._rotationX + ", rotationY: " + this._rotationY + ", rotationZ: " + this._rotationZ);
});
});
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
