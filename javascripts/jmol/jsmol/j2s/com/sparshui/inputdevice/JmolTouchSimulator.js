Clazz.declarePackage ("com.sparshui.inputdevice");
Clazz.load (["java.util.TimerTask", "J.api.JmolTouchSimulatorInterface", "java.util.Hashtable", "$.TreeSet"], "com.sparshui.inputdevice.JmolTouchSimulator", ["java.io.DataOutputStream", "java.lang.Thread", "java.net.Socket", "java.util.Timer", "javax.swing.SwingUtilities", "JU.Logger"], function () {
var c$ = Clazz.decorateAsClass (function () {
this._events = null;
this._active = null;
this._recording = false;
this._touchID = 0;
this._when = 0;
this._timer = null;
this._display = null;
this._out = null;
if (!Clazz.isClassDefined ("com.sparshui.inputdevice.JmolTouchSimulator.TouchData")) {
com.sparshui.inputdevice.JmolTouchSimulator.$JmolTouchSimulator$TouchData$ ();
}
if (!Clazz.isClassDefined ("com.sparshui.inputdevice.JmolTouchSimulator.TouchDataComparator")) {
com.sparshui.inputdevice.JmolTouchSimulator.$JmolTouchSimulator$TouchDataComparator$ ();
}
if (!Clazz.isClassDefined ("com.sparshui.inputdevice.JmolTouchSimulator.TouchTimerTask")) {
com.sparshui.inputdevice.JmolTouchSimulator.$JmolTouchSimulator$TouchTimerTask$ ();
}
Clazz.instantialize (this, arguments);
}, com.sparshui.inputdevice, "JmolTouchSimulator", null, J.api.JmolTouchSimulatorInterface);
Clazz.prepareFields (c$, function () {
this._events =  new java.util.TreeSet (Clazz.innerTypeInstance (com.sparshui.inputdevice.JmolTouchSimulator.TouchDataComparator, this, null));
this._active =  new java.util.Hashtable ();
});
Clazz.makeConstructor (c$, 
function () {
});
Clazz.overrideMethod (c$, "dispose", 
function () {
try {
} catch (e) {
if (Clazz.exceptionOf(e, Exception)){
} else {
throw e;
}
}
try {
this._out.close ();
} catch (e) {
if (Clazz.exceptionOf(e, Exception)){
} else {
throw e;
}
}
try {
this._timer.cancel ();
} catch (e) {
if (Clazz.exceptionOf(e, Exception)){
} else {
throw e;
}
}
});
Clazz.overrideMethod (c$, "startSimulator", 
function (display) {
this._display = display;
var address = "localhost";
this._timer =  new java.util.Timer ();
try {
var socket =  new java.net.Socket (address, 5947);
this._out =  new java.io.DataOutputStream (socket.getOutputStream ());
this._out.writeByte (1);
socket.close ();
return true;
} catch (e$$) {
if (Clazz.exceptionOf(e$$,"java.net.UnknownHostException")){
var e = e$$;
{
JU.Logger.error ("Could not locate a server at " + address);
}
} else if (Clazz.exceptionOf(e$$,"java.io.IOException")){
var e = e$$;
{
JU.Logger.error ("Failed to connect to server at " + address);
}
} else {
throw e$$;
}
}
return false;
}, "~O");
Clazz.overrideMethod (c$, "toggleMode", 
function () {
if (this._recording) {
this.endRecording ();
} else {
this.startRecording ();
}});
Clazz.overrideMethod (c$, "startRecording", 
function () {
this._recording = true;
this._active.clear ();
});
Clazz.overrideMethod (c$, "endRecording", 
function () {
this._recording = false;
this.dispatchTouchEvents ();
});
Clazz.overrideMethod (c$, "mousePressed", 
function (time, x, y) {
this.handleMouseEvent (time, x, y, 0);
}, "~N,~N,~N");
Clazz.overrideMethod (c$, "mouseReleased", 
function (time, x, y) {
this.handleMouseEvent (time, x, y, 1);
}, "~N,~N,~N");
Clazz.overrideMethod (c$, "mouseDragged", 
function (time, x, y) {
this.handleMouseEvent (time, x, y, 2);
}, "~N,~N,~N");
Clazz.defineMethod (c$, "handleMouseEvent", 
function (time, x, y, type) {
var te = Clazz.innerTypeInstance (com.sparshui.inputdevice.JmolTouchSimulator.TouchData, this, null);
te.id = (type == 0) ? ++this._touchID : this._touchID;
var p =  new java.awt.Point (x, y);
try {
javax.swing.SwingUtilities.convertPointToScreen (p, this._display);
} catch (e) {
return;
}
te.x = p.x;
te.y = p.y;
te.type = type;
te.when = time;
if (this._recording) {
if (type == 0) {
te.delay = 0;
this._when = te.when;
} else {
te.delay = te.when - this._when;
}this._events.add (te);
} else {
this.dispatchTouchEvent (te);
if (JU.Logger.debugging) JU.Logger.debug ("[JmolTouchSimulator] dispatchTouchEvent(" + te.id + ", " + te.x + ", " + te.y + ", " + te.type + ")");
}}, "~N,~N,~N,~N");
Clazz.defineMethod (c$, "dispatchTouchEvents", 
function () {
for (var data, $data = this._events.iterator (); $data.hasNext () && ((data = $data.next ()) || true);) {
var task = Clazz.innerTypeInstance (com.sparshui.inputdevice.JmolTouchSimulator.TouchTimerTask, this, null, data);
this._timer.schedule (task, data.delay + 250);
}
this._events.clear ();
this._touchID = 0;
});
Clazz.defineMethod (c$, "dispatchTouchEvent", 
function (data) {
var tk = java.awt.Toolkit.getDefaultToolkit ();
var dim = tk.getScreenSize ();
if (JU.Logger.debugging) JU.Logger.debug ("[JmolTouchSimulator] dispatchTouchEvent(" + data.id + ", " + data.x + ", " + data.y + ", " + data.type + ")");
try {
this._out.writeInt (-1);
this._out.writeInt (21);
this._out.writeInt (data.id);
this._out.writeFloat ((data.x / dim.width));
this._out.writeFloat ((data.y / dim.height));
this._out.writeByte (data.type);
this._out.writeLong (data.when);
} catch (e1) {
if (Clazz.exceptionOf(e1,"java.io.IOException")){
System.err.println ("Failed to send event to server.");
} else {
throw e1;
}
}
}, "com.sparshui.inputdevice.JmolTouchSimulator.TouchData");
c$.$JmolTouchSimulator$TouchData$ = function () {
/*if4*/;(function(){
var c$ = Clazz.decorateAsClass (function () {
Clazz.prepareCallback (this, arguments);
this.type = 0;
this.id = 0;
this.x = 0;
this.y = 0;
this.when = 0;
this.delay = 0;
Clazz.instantialize (this, arguments);
}, com.sparshui.inputdevice.JmolTouchSimulator, "TouchData");
/*eoif4*/})();
};
c$.$JmolTouchSimulator$TouchDataComparator$ = function () {
/*if4*/;(function(){
var c$ = Clazz.decorateAsClass (function () {
Clazz.prepareCallback (this, arguments);
Clazz.instantialize (this, arguments);
}, com.sparshui.inputdevice.JmolTouchSimulator, "TouchDataComparator", null, java.util.Comparator);
Clazz.overrideMethod (c$, "compare", 
function (o1, o2) {
return (o1.delay == o2.delay ? (o1.when < o2.when ? -1 : 1) : o1.delay < o2.delay ? -1 : 1);
}, "com.sparshui.inputdevice.JmolTouchSimulator.TouchData,com.sparshui.inputdevice.JmolTouchSimulator.TouchData");
/*eoif4*/})();
};
c$.$JmolTouchSimulator$TouchTimerTask$ = function () {
/*if4*/;(function(){
var c$ = Clazz.decorateAsClass (function () {
Clazz.prepareCallback (this, arguments);
this.data = null;
Clazz.instantialize (this, arguments);
}, com.sparshui.inputdevice.JmolTouchSimulator, "TouchTimerTask", java.util.TimerTask);
Clazz.makeConstructor (c$, 
function (data) {
Clazz.superConstructor (this, com.sparshui.inputdevice.JmolTouchSimulator.TouchTimerTask, []);
this.data = data;
}, "com.sparshui.inputdevice.JmolTouchSimulator.TouchData");
Clazz.overrideMethod (c$, "run", 
function () {
Thread.currentThread ().setName ("JmolTouchSimulator for type " + this.data.id);
this.b$["com.sparshui.inputdevice.JmolTouchSimulator"].dispatchTouchEvent (this.data);
var iid = Integer.$valueOf (this.data.id);
if (this.data.type == 1) {
this.b$["com.sparshui.inputdevice.JmolTouchSimulator"]._active.remove (iid);
} else {
this.b$["com.sparshui.inputdevice.JmolTouchSimulator"]._active.put (iid, this.data);
}Thread.currentThread ().setName ("JmolTouchSimulator idle");
});
/*eoif4*/})();
};
});
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
