Clazz.declarePackage ("com.sparshui.client");
Clazz.load (["com.sparshui.common.ClientProtocol"], "com.sparshui.client.ClientToServerProtocol", ["com.sparshui.common.Location", "com.sparshui.common.messages.events.DragEvent", "$.FlickEvent", "$.RelativeDragEvent", "$.RotateEvent", "$.SpinEvent", "$.TouchEvent", "$.ZoomEvent", "com.sparshui.common.utils.Converter"], function () {
var c$ = Clazz.declareType (com.sparshui.client, "ClientToServerProtocol", com.sparshui.common.ClientProtocol);
Clazz.defineMethod (c$, "processRequest", 
function (client) {
try {
var type = this._in.readByte ();
var length = this._in.readInt ();
var data =  Clazz.newByteArray (length, 0);
if (length > 0) this._in.readFully (data);
switch (type) {
case 0:
this.handleEvents (client, data);
break;
case 1:
this.handleGetGroupID (client, data);
break;
case 2:
this.handleGetAllowedGestures (client, data);
break;
}
} catch (e) {
if (Clazz.exceptionOf(e,"java.io.IOException")){
System.err.println ("[Client Protocol] GestureServer Connection Lost");
this.handleEvents (client, null);
return false;
} else {
throw e;
}
}
return true;
}, "com.sparshui.client.SparshClient");
Clazz.defineMethod (c$, "handleEvents", 
function (client, data) {
if (data == null) {
client.processEvent (-1, null);
return;
}if (data.length < 1) {
return;
}var groupID = com.sparshui.common.utils.Converter.byteArrayToInt (data);
var eventType = com.sparshui.common.utils.Converter.byteArrayToInt (data, 4);
var newData =  Clazz.newByteArray (data.length - 8, 0);
System.arraycopy (data, 8, newData, 0, data.length - 8);
var event = null;
switch (eventType) {
case -2:
client.processEvent (-2, null);
return;
case 0:
event =  new com.sparshui.common.messages.events.DragEvent (newData);
break;
case 1:
event =  new com.sparshui.common.messages.events.RotateEvent (newData);
break;
case 2:
event =  new com.sparshui.common.messages.events.SpinEvent ();
break;
case 3:
event =  new com.sparshui.common.messages.events.TouchEvent (newData);
break;
case 4:
event =  new com.sparshui.common.messages.events.ZoomEvent (newData);
break;
case 6:
event =  new com.sparshui.common.messages.events.FlickEvent (newData);
break;
case 7:
event =  new com.sparshui.common.messages.events.RelativeDragEvent (newData);
break;
}
if (event != null) client.processEvent (groupID, event);
}, "com.sparshui.client.SparshClient,~A");
Clazz.defineMethod (c$, "handleGetGroupID", 
function (client, data) {
this._out.writeInt (client.getGroupID ( new com.sparshui.common.Location (com.sparshui.common.utils.Converter.byteArrayToFloat (data, 0), com.sparshui.common.utils.Converter.byteArrayToFloat (data, 4))));
}, "com.sparshui.client.SparshClient,~A");
Clazz.defineMethod (c$, "handleGetAllowedGestures", 
function (client, data) {
var gType;
var gestureTypes = client.getAllowedGestures (com.sparshui.common.utils.Converter.byteArrayToInt (data));
var length = (gestureTypes == null ? 0 : gestureTypes.size ());
var blen = length * 4;
for (var i = 0; i < length; i++) {
gType = gestureTypes.get (i);
if (gType.sType != null) blen += gType.sType.length;
}
this._out.writeInt (blen);
for (var i = 0; i < length; i++) {
gType = gestureTypes.get (i);
if (gType.sType == null) {
this._out.writeInt (gType.iType);
} else {
var len = gType.sType.length;
if (len > 0) {
this._out.writeInt (-len);
this._out.write (com.sparshui.common.utils.Converter.stringToByteArray (gType.sType));
}}}
}, "com.sparshui.client.SparshClient,~A");
});
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
