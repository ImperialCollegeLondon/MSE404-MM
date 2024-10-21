Clazz.declarePackage ("com.sparshui.common.messages.events");
Clazz.load (["com.sparshui.common.Event"], "com.sparshui.common.messages.events.RelativeDragEvent", ["com.sparshui.common.utils.Converter"], function () {
var c$ = Clazz.decorateAsClass (function () {
this._changeInX = 0;
this._changeInY = 0;
Clazz.instantialize (this, arguments);
}, com.sparshui.common.messages.events, "RelativeDragEvent", null, com.sparshui.common.Event);
Clazz.makeConstructor (c$, 
function () {
this._changeInX = 0;
this._changeInY = 0;
});
Clazz.makeConstructor (c$, 
function (changeInX, changeInY) {
this._changeInX = changeInX;
this._changeInY = changeInY;
}, "~N,~N");
Clazz.defineMethod (c$, "getChangeInX", 
function () {
return this._changeInX;
});
Clazz.defineMethod (c$, "getChangeInY", 
function () {
return this._changeInY;
});
Clazz.makeConstructor (c$, 
function (data) {
if (data.length < 8) {
System.err.println ("Error constructing Drag Event.");
this._changeInX = 0;
this._changeInY = 0;
} else {
this._changeInX = com.sparshui.common.utils.Converter.byteArrayToFloat (data, 0);
this._changeInY = com.sparshui.common.utils.Converter.byteArrayToFloat (data, 4);
}}, "~A");
Clazz.overrideMethod (c$, "getEventType", 
function () {
return 7;
});
Clazz.overrideMethod (c$, "serialize", 
function () {
var data =  Clazz.newByteArray (12, 0);
com.sparshui.common.utils.Converter.intToByteArray (data, 0, this.getEventType ());
com.sparshui.common.utils.Converter.floatToByteArray (data, 4, this._changeInX);
com.sparshui.common.utils.Converter.floatToByteArray (data, 8, this._changeInY);
return data;
});
});
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
