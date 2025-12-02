Clazz.declarePackage ("com.sparshui.server");
Clazz.load (null, "com.sparshui.server.Group", ["com.sparshui.server.GestureFactory", "JU.Lst"], function () {
var c$ = Clazz.decorateAsClass (function () {
this._id = 0;
this._gestureTypes = null;
this._gestures = null;
this._touchPoints = null;
this._clientProtocol = null;
Clazz.instantialize (this, arguments);
}, com.sparshui.server, "Group");
Clazz.makeConstructor (c$, 
function (id, gestureTypes, clientProtocol) {
this._id = id;
this._gestureTypes = gestureTypes;
this._gestures =  new JU.Lst ();
this._touchPoints =  new JU.Lst ();
this._clientProtocol = clientProtocol;
for (var i = 0; i < this._gestureTypes.size (); i++) {
var gesture = com.sparshui.server.GestureFactory.createGesture (this._gestureTypes.get (i));
if (gesture != null) this._gestures.addLast (gesture);
}
}, "~N,JU.Lst,com.sparshui.server.ServerToClientProtocol");
Clazz.defineMethod (c$, "getID", 
function () {
return this._id;
});
Clazz.defineMethod (c$, "update", 
function (changedPoint) {
var events =  new JU.Lst ();
var state = changedPoint.getState ();
if (state == 0) this._touchPoints.addLast (changedPoint);
for (var i = 0; i < this._gestures.size (); i++) {
var gesture = this._gestures.get (i);
events.addAll (gesture.processChange (this._touchPoints, changedPoint));
}
if (state == 1) this._touchPoints.removeObj (changedPoint);
try {
this._clientProtocol.processEvents (this._id, events);
} catch (e) {
if (Clazz.exceptionOf(e,"java.io.IOException")){
} else {
throw e;
}
}
}, "com.sparshui.server.TouchPoint");
});
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
