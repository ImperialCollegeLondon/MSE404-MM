Clazz.declarePackage ("com.sparshui.server");
Clazz.load (["J.api.JmolGestureServerInterface", "JU.Lst"], "com.sparshui.server.GestureServer", ["com.sparshui.server.ClientConnection", "$.InputDeviceConnection", "$.TouchPoint", "java.lang.Thread", "java.net.ServerSocket", "JU.Logger"], function () {
var c$ = Clazz.decorateAsClass (function () {
this.clientServer = null;
this.deviceServer = null;
this.main = null;
this.clientThread = null;
this.deviceThread = null;
this._clientSocket = null;
this._deviceSocket = null;
this._mySocket = null;
this._clients = null;
this.port = 0;
this.ic = null;
this.myState = 0;
Clazz.instantialize (this, arguments);
}, com.sparshui.server, "GestureServer", null, [Runnable, J.api.JmolGestureServerInterface]);
Clazz.prepareFields (c$, function () {
this._clients =  new JU.Lst ();
});
Clazz.makeConstructor (c$, 
function () {
JU.Logger.info (this + " constructed");
});
Clazz.overrideMethod (c$, "finalize", 
function () {
if (JU.Logger.debugging) JU.Logger.debug (this + " finalized");
});
Clazz.makeConstructor (c$, 
function (port, main) {
this.port = port;
this.main = main;
}, "~N,com.sparshui.server.GestureServer");
Clazz.overrideMethod (c$, "startGestureServer", 
function () {
this.clientServer =  new com.sparshui.server.GestureServer (5946, this);
this.clientThread =  new Thread (this.clientServer);
this.clientThread.setName ("Jmol SparshUI Client GestureServer on port 5946");
this.clientThread.start ();
this.deviceServer =  new com.sparshui.server.GestureServer (5947, this);
this.deviceThread =  new Thread (this.deviceServer);
this.deviceThread.setName ("Jmol SparshUI Device GestureServer on port 5947");
this.deviceThread.start ();
});
Clazz.overrideMethod (c$, "dispose", 
function () {
try {
this._clientSocket.close ();
} catch (e) {
if (Clazz.exceptionOf(e, Exception)){
} else {
throw e;
}
}
try {
this._deviceSocket.close ();
} catch (e) {
if (Clazz.exceptionOf(e, Exception)){
} else {
throw e;
}
}
try {
this.clientThread.interrupt ();
} catch (e) {
if (Clazz.exceptionOf(e, Exception)){
} else {
throw e;
}
}
try {
this.deviceThread.interrupt ();
} catch (e) {
if (Clazz.exceptionOf(e, Exception)){
} else {
throw e;
}
}
this._clientSocket = null;
this.clientThread = null;
this._deviceSocket = null;
this.deviceThread = null;
this.clientServer = null;
this.deviceServer = null;
});
Clazz.overrideMethod (c$, "run", 
function () {
try {
this.openSocket ();
this.acceptConnections ();
} catch (e) {
if (Clazz.exceptionOf(e, Exception)){
JU.Logger.info ("[GestureServer] connection unavailable");
} else {
throw e;
}
}
});
Clazz.defineMethod (c$, "openSocket", 
function () {
try {
if (this.port == 5946) this._mySocket = this.main._clientSocket =  new java.net.ServerSocket (this.port);
 else this._mySocket = this.main._deviceSocket =  new java.net.ServerSocket (this.port);
JU.Logger.info ("[GestureServer] Socket Open: " + this.port);
this.main.myState = 1;
} catch (e) {
if (Clazz.exceptionOf(e,"java.io.IOException")){
JU.Logger.error ("[GestureServer] Failed to open a server socket.");
e.printStackTrace ();
this.main.myState = 0;
} else {
throw e;
}
}
});
Clazz.defineMethod (c$, "acceptConnections", 
function () {
while (!this._mySocket.isClosed ()) {
try {
if (this.port == 5947) {
JU.Logger.info ("[GestureServer] Accepting device connections");
this.acceptConnection (this._mySocket.accept ());
return;
}JU.Logger.info ("[GestureServer] Accepting client connections");
this.acceptConnection (this._mySocket.accept ());
} catch (e) {
if (Clazz.exceptionOf(e,"java.io.IOException")){
JU.Logger.error ("[GestureServer] Failed to establish connection on port " + this.port);
e.printStackTrace ();
} else {
throw e;
}
}
}
JU.Logger.info ("[GestureServer] Socket Closed on port " + this.port);
});
Clazz.defineMethod (c$, "acceptConnection", 
function (socket) {
var add = socket.getInetAddress ().getAddress ();
if (add[0] != 127 || add[1] != 0 || add[2] != 0 || add[3] != 1) return;
var type = socket.getInputStream ().read ();
if (type == 0) {
JU.Logger.info ("[GestureServer] client connection established on port " + this.port);
this.acceptClientConnection (socket);
} else if (type == 1) {
JU.Logger.info ("[GestureServer] device connection established on port " + this.port);
this.acceptInputDeviceConnection (socket);
}}, "java.net.Socket");
Clazz.defineMethod (c$, "acceptClientConnection", 
function (socket) {
JU.Logger.info ("[GestureServer] Client connection accepted");
var cc =  new com.sparshui.server.ClientConnection (socket);
this.main._clients.addLast (cc);
if (this.main.ic == null) {
cc.processError (-2);
} else {
this.main.myState |= 2;
}}, "java.net.Socket");
Clazz.defineMethod (c$, "acceptInputDeviceConnection", 
function (socket) {
JU.Logger.info ("[GestureServer] Input device connection accepted");
this.main.ic =  new com.sparshui.server.InputDeviceConnection (this, socket);
this.main.myState |= 4;
}, "java.net.Socket");
Clazz.defineMethod (c$, "notifyInputLost", 
function () {
JU.Logger.error ("[GestureServer] sending clients message that input device was lost.");
this.main.ic = null;
this.main.myState &= -5;
this.processBirth (null);
});
Clazz.defineMethod (c$, "processTouchPoint", 
function (inputDeviceTouchPoints, id, location, time, state) {
if (JU.Logger.debugging) {
JU.Logger.debug ("[GestureServer] processTouchPoint id=" + id + " state=" + state + " " + location + " " + time);
}var iid = Integer.$valueOf (id);
if (inputDeviceTouchPoints.containsKey (iid)) {
var touchPoint = inputDeviceTouchPoints.get (iid);
if (!touchPoint.isClaimed ()) return false;
if (JU.Logger.debugging) JU.Logger.debug ("[GestureServer] OK");
{
touchPoint.update (location, time, state);
}return true;
}var touchPoint =  new com.sparshui.server.TouchPoint (id, location, time);
inputDeviceTouchPoints.put (iid, touchPoint);
return this.processBirth (touchPoint);
}, "java.util.Map,~N,com.sparshui.common.Location,~N,~N");
Clazz.defineMethod (c$, "processBirth", 
function (touchPoint) {
var clients_to_remove = null;
var isClaimed = false;
for (var i = 0; i < this.main._clients.size (); i++) {
var client = this.main._clients.get (i);
try {
if (touchPoint == null) client.processError (-2);
 else isClaimed = client.processBirth (touchPoint);
if (isClaimed) break;
} catch (e) {
if (Clazz.exceptionOf(e,"java.io.IOException")){
if (clients_to_remove == null) clients_to_remove =  new JU.Lst ();
clients_to_remove.addLast (client);
} else {
throw e;
}
}
}
if (clients_to_remove != null) for (var i = 0; i < clients_to_remove.size (); i++) {
this.main._clients.removeObj (clients_to_remove.get (i));
JU.Logger.info ("[GestureServer] Client Disconnected");
}
return isClaimed;
}, "com.sparshui.server.TouchPoint");
Clazz.overrideMethod (c$, "getState", 
function () {
return this.myState;
});
});
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
