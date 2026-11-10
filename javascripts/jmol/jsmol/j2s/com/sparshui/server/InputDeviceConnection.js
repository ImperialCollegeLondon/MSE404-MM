Clazz.declarePackage ("com.sparshui.server");
Clazz.load (null, "com.sparshui.server.InputDeviceConnection", ["com.sparshui.common.Location", "java.io.DataInputStream", "java.lang.Thread", "java.util.ArrayList", "$.Hashtable"], function () {
var c$ = Clazz.decorateAsClass (function () {
this._gestureServer = null;
this._socket = null;
this._in = null;
this._touchPoints = null;
this._flaggedids = null;
Clazz.instantialize (this, arguments);
}, com.sparshui.server, "InputDeviceConnection", null, Runnable);
Clazz.makeConstructor (c$, 
function (gestureServer, socket) {
this._gestureServer = gestureServer;
this._socket = socket;
this._in =  new java.io.DataInputStream (socket.getInputStream ());
this._touchPoints =  new java.util.Hashtable ();
this._flaggedids =  new java.util.ArrayList ();
this.startListening ();
}, "com.sparshui.server.GestureServer,java.net.Socket");
Clazz.defineMethod (c$, "removeDeadTouchPoints", 
function () {
for (var i = 0; i < this._flaggedids.size (); i++) {
var id = this._flaggedids.get (i);
this._touchPoints.remove (id);
}
this._flaggedids.clear ();
});
Clazz.defineMethod (c$, "flagTouchPointForRemoval", 
function (id) {
this._flaggedids.add (Integer.$valueOf (id));
}, "~N");
Clazz.defineMethod (c$, "receiveData", 
function () {
try {
while (!this._socket.isInputShutdown ()) {
this.readTouchPoints ();
}
} catch (e) {
if (Clazz.exceptionOf(e,"java.io.IOException")){
System.out.println ("[InputDeviceConnection] InputDevice Disconnected");
this._gestureServer.notifyInputLost ();
} else {
throw e;
}
}
});
Clazz.defineMethod (c$, "readTouchPoints", 
function () {
var count = this._in.readInt ();
if (count == 0) {
this._in.close ();
return false;
}var touchPointDataLength;
if (count < 0) {
count = -count;
touchPointDataLength = this._in.readInt ();
} else {
touchPointDataLength = 13;
}var doConsume = false;
for (var i = 0; i < count; i++) doConsume = new Boolean (doConsume | this.readTouchPoint (touchPointDataLength)).valueOf ();

this.removeDeadTouchPoints ();
return doConsume;
});
Clazz.defineMethod (c$, "readTouchPoint", 
function (len) {
var id = this._in.readInt ();
var x = this._in.readFloat ();
var y = this._in.readFloat ();
var state = this._in.readByte ();
var time = (len >= 21 ? this._in.readLong () : System.currentTimeMillis ());
if (len > 21) this._in.read ( Clazz.newByteArray (len - 21, 0));
var location =  new com.sparshui.common.Location (x, y);
var doConsume = this._gestureServer.processTouchPoint (this._touchPoints, id, location, time, state);
if (state == 1) this.flagTouchPointForRemoval (id);
return doConsume;
}, "~N");
Clazz.defineMethod (c$, "startListening", 
function () {
var thread =  new Thread (this);
thread.setName ("SparshUI Server->InputDeviceConnection on port 5947");
thread.start ();
});
Clazz.overrideMethod (c$, "run", 
function () {
this.receiveData ();
});
});
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
