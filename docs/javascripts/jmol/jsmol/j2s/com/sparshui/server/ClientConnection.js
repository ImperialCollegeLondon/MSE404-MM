Clazz.declarePackage ("com.sparshui.server");
Clazz.load (null, "com.sparshui.server.ClientConnection", ["com.sparshui.server.Group", "$.ServerToClientProtocol", "java.util.Hashtable"], function () {
var c$ = Clazz.decorateAsClass (function () {
this._protocol = null;
this._groups = null;
Clazz.instantialize (this, arguments);
}, com.sparshui.server, "ClientConnection");
Clazz.makeConstructor (c$, 
function (socket) {
this._protocol =  new com.sparshui.server.ServerToClientProtocol (socket);
this._groups =  new java.util.Hashtable ();
}, "java.net.Socket");
Clazz.defineMethod (c$, "processBirth", 
function (touchPoint) {
var groupID = (touchPoint == null ? 0x10000000 : this.getGroupID (touchPoint));
var jmolFlags = (groupID & 0xF0000000);
if (jmolFlags != 0) {
switch (jmolFlags) {
case 0x10000000:
this._groups =  new java.util.Hashtable ();
break;
}
groupID &= ~jmolFlags;
}var group = this.getGroup (groupID);
if (group != null) {
touchPoint.setGroup (group);
return true;
}return false;
}, "com.sparshui.server.TouchPoint");
Clazz.defineMethod (c$, "getGestures", 
function (groupID) {
return this._protocol.getGestures (groupID);
}, "~N");
Clazz.defineMethod (c$, "getGroupID", 
function (touchPoint) {
return this._protocol.getGroupID (touchPoint);
}, "com.sparshui.server.TouchPoint");
Clazz.defineMethod (c$, "getGroup", 
function (groupID) {
if (groupID == 0) return null;
var group = null;
var gid = Integer.$valueOf (groupID);
if (this._groups.containsKey (gid)) {
group = this._groups.get (gid);
} else {
var gestureTypes = this.getGestures (groupID);
group =  new com.sparshui.server.Group (groupID, gestureTypes, this._protocol);
this._groups.put (gid, group);
}return group;
}, "~N");
Clazz.defineMethod (c$, "processError", 
function (errCode) {
try {
this._protocol.processError (errCode);
} catch (e) {
if (Clazz.exceptionOf(e,"java.io.IOException")){
} else {
throw e;
}
}
}, "~N");
});
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
