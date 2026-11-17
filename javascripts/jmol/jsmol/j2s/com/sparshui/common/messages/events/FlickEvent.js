Clazz.declarePackage ("com.sparshui.common.messages.events");
Clazz.load (["com.sparshui.common.Event"], "com.sparshui.common.messages.events.FlickEvent", ["com.sparshui.common.utils.Converter"], function () {
var c$ = Clazz.decorateAsClass (function () {
this.xDirection = 0;
this.yDirection = 0;
this.speedLevel = 0;
Clazz.instantialize (this, arguments);
}, com.sparshui.common.messages.events, "FlickEvent", null, com.sparshui.common.Event);
Clazz.makeConstructor (c$, 
function () {
this.xDirection = 0;
this.yDirection = 0;
this.speedLevel = 0;
});
Clazz.makeConstructor (c$, 
function (absx, absy) {
}, "~N,~N");
Clazz.makeConstructor (c$, 
function (_speedLevel, _xDirection, _yDirection) {
this.speedLevel = _speedLevel;
this.xDirection = _xDirection;
this.yDirection = _yDirection;
}, "~N,~N,~N");
Clazz.makeConstructor (c$, 
function (data) {
if (data.length < 12) {
System.err.println ("Error constructing Flick Event.");
} else {
this.speedLevel = com.sparshui.common.utils.Converter.byteArrayToFloat (data, 0);
this.xDirection = com.sparshui.common.utils.Converter.byteArrayToFloat (data, 4);
this.yDirection = com.sparshui.common.utils.Converter.byteArrayToFloat (data, 8);
}}, "~A");
Clazz.defineMethod (c$, "getSpeedLevel", 
function () {
return this.speedLevel;
});
Clazz.defineMethod (c$, "getXdirection", 
function () {
return this.xDirection;
});
Clazz.defineMethod (c$, "getYdirection", 
function () {
return this.yDirection;
});
Clazz.overrideMethod (c$, "getEventType", 
function () {
return 6;
});
Clazz.overrideMethod (c$, "toString", 
function () {
var ret = "Flick Event";
return ret;
});
Clazz.overrideMethod (c$, "serialize", 
function () {
var data =  Clazz.newByteArray (16, 0);
com.sparshui.common.utils.Converter.intToByteArray (data, 0, this.getEventType ());
com.sparshui.common.utils.Converter.floatToByteArray (data, 4, this.speedLevel);
com.sparshui.common.utils.Converter.floatToByteArray (data, 8, this.xDirection);
com.sparshui.common.utils.Converter.floatToByteArray (data, 12, this.yDirection);
return data;
});
});
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
