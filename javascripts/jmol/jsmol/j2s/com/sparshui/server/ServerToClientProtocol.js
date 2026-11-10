Clazz.declarePackage ("com.sparshui.server");
Clazz.load (["com.sparshui.common.ClientProtocol"], "com.sparshui.server.ServerToClientProtocol", ["com.sparshui.GestureType", "com.sparshui.common.utils.Converter", "java.io.ByteArrayOutputStream", "$.DataOutputStream", "JU.Lst"], function () {
var c$ = Clazz.decorateAsClass (function () {
this._bufferOut = null;
this._buffer = null;
Clazz.instantialize (this, arguments);
}, com.sparshui.server, "ServerToClientProtocol", com.sparshui.common.ClientProtocol);
Clazz.makeConstructor (c$, 
function (socket) {
Clazz.superConstructor (this, com.sparshui.server.ServerToClientProtocol, [socket]);
this._buffer =  new java.io.ByteArrayOutputStream ();
this._bufferOut =  new java.io.DataOutputStream (this._buffer);
}, "java.net.Socket");
Clazz.defineMethod (c$, "getGestures", 
function (groupID) {
var gestures =  new JU.Lst ();
this._bufferOut.writeInt (groupID);
this.sendBuffer (2);
for (var length = this._in.readInt (); length > 0; length -= 4) {
var gestureID = this._in.readInt ();
if (gestureID < 0) {
var bytes =  Clazz.newByteArray (-gestureID, 0);
this._in.read (bytes);
gestures.addLast ( new com.sparshui.GestureType (com.sparshui.common.utils.Converter.byteArrayToString (bytes)));
length -= bytes.length;
} else {
gestures.addLast ( new com.sparshui.GestureType (gestureID));
}}
return gestures;
}, "~N");
Clazz.defineMethod (c$, "getGroupID", 
function (touchPoint) {
var tempFloat =  Clazz.newByteArray (4, 0);
com.sparshui.common.utils.Converter.floatToByteArray (tempFloat, 0, touchPoint.getLocation ().getX ());
this._bufferOut.write (tempFloat);
com.sparshui.common.utils.Converter.floatToByteArray (tempFloat, 0, touchPoint.getLocation ().getY ());
this._bufferOut.write (tempFloat);
this.sendBuffer (1);
var ret = this._in.readInt ();
return ret;
}, "com.sparshui.server.TouchPoint");
Clazz.defineMethod (c$, "processEvents", 
function (groupID, events) {
for (var i = 0; i < events.size (); i++) {
this._bufferOut.writeInt (groupID);
this._bufferOut.write (events.get (i).serialize ());
this.sendBuffer (0);
}
}, "~N,JU.Lst");
Clazz.defineMethod (c$, "processError", 
function (errCode) {
this._bufferOut.writeInt (-1);
this._bufferOut.writeInt (errCode);
this.sendBuffer (0);
}, "~N");
Clazz.defineMethod (c$, "sendBuffer", 
function (type) {
this._out.writeByte (type);
this._out.writeInt (this._buffer.size ());
this._out.write (this._buffer.toByteArray ());
this._buffer.reset ();
}, "~N");
});
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
