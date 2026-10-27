Clazz.declarePackage ("com.sparshui.common.utils");
Clazz.load (null, "com.sparshui.common.utils.Converter", ["java.lang.Float"], function () {
var c$ = Clazz.declareType (com.sparshui.common.utils, "Converter");
c$.intToByteArray = Clazz.defineMethod (c$, "intToByteArray", 
function (intBits) {
var ret =  Clazz.newByteArray (4, 0);
ret[0] = ((intBits & 0xff000000) >> 24);
ret[1] = ((intBits & 0x00ff0000) >> 16);
ret[2] = ((intBits & 0x0000ff00) >> 8);
ret[3] = ((intBits & 0x000000ff) >> 0);
return ret;
}, "~N");
c$.intToByteArray = Clazz.defineMethod (c$, "intToByteArray", 
function (data, i, idata) {
data[i++] = ((idata & 0xff000000) >> 24);
data[i++] = ((idata & 0x00ff0000) >> 16);
data[i++] = ((idata & 0x0000ff00) >> 8);
data[i] = ((idata & 0x000000ff) >> 0);
}, "~A,~N,~N");
c$.byteArrayToInt = Clazz.defineMethod (c$, "byteArrayToInt", 
function (b) {
return ((b[0] << 24) & 0xFF000000) | ((b[1] << 16) & 0x00FF0000) | ((b[2] << 8) & 0x0000FF00) | (b[3] & 0x000000FF);
}, "~A");
c$.byteArrayToInt = Clazz.defineMethod (c$, "byteArrayToInt", 
function (b, i) {
return ((b[i++] << 24) & 0xFF000000) | ((b[i++] << 16) & 0x00FF0000) | ((b[i++] << 8) & 0x0000FF00) | (b[i] & 0x000000FF);
}, "~A,~N");
c$.floatToByteArray = Clazz.defineMethod (c$, "floatToByteArray", 
function (data, i, fdata) {
com.sparshui.common.utils.Converter.intToByteArray (data, i, Float.floatToIntBits (fdata));
}, "~A,~N,~N");
c$.byteArrayToFloat = Clazz.defineMethod (c$, "byteArrayToFloat", 
function (data, i) {
return Float.intBitsToFloat (com.sparshui.common.utils.Converter.byteArrayToInt (data, i));
}, "~A,~N");
c$.longToByteArray = Clazz.defineMethod (c$, "longToByteArray", 
function (data, i, ldata) {
com.sparshui.common.utils.Converter.intToByteArray (data, i, ((ldata >> 32)));
com.sparshui.common.utils.Converter.intToByteArray (data, i + 4, ((ldata & 0xFFFFFFFF)));
}, "~A,~N,~N");
c$.byteArrayToLong = Clazz.defineMethod (c$, "byteArrayToLong", 
function (data, i) {
return ((com.sparshui.common.utils.Converter.byteArrayToInt (data, i)) << 32) | com.sparshui.common.utils.Converter.byteArrayToInt (data, i + 4) & 0xFFFFFFFF;
}, "~A,~N");
c$.byteArrayToString = Clazz.defineMethod (c$, "byteArrayToString", 
function (bytes) {
var chars =  Clazz.newCharArray (bytes.length, '\0');
for (var i = 0; i < chars.length; i++) chars[i] = String.fromCharCode (bytes[i]);

return  String.instantialize (chars);
}, "~A");
c$.stringToByteArray = Clazz.defineMethod (c$, "stringToByteArray", 
function (s) {
var chars = s.toCharArray ();
var bytes =  Clazz.newByteArray (s.length, 0);
for (var i = 0; i < chars.length; i++) bytes[i] = (chars[i]).charCodeAt (0);

return bytes;
}, "~S");
});
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
