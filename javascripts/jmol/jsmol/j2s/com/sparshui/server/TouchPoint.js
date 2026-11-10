Clazz.declarePackage ("com.sparshui.server");
var c$ = Clazz.decorateAsClass (function () {
this._id = 0;
this._location = null;
this._state = 0;
this._changed = false;
this._time = 0;
this._group = null;
Clazz.instantialize (this, arguments);
}, com.sparshui.server, "TouchPoint");
Clazz.defineMethod (c$, "isClaimed", 
function () {
return (this._group != null);
});
Clazz.makeConstructor (c$, 
function (id, location, time) {
this._id = id;
this._location = location;
this._time = time;
this._state = 0;
}, "~N,com.sparshui.common.Location,~N");
Clazz.makeConstructor (c$, 
function (tp) {
this._id = tp._id;
this._location = tp._location;
this._state = tp._state;
this._time = tp._time;
}, "com.sparshui.server.TouchPoint");
Clazz.defineMethod (c$, "getTime", 
function () {
return this._time;
});
Clazz.defineMethod (c$, "getID", 
function () {
return this._id;
});
Clazz.defineMethod (c$, "getLocation", 
function () {
return this._location;
});
Clazz.defineMethod (c$, "getState", 
function () {
return this._state;
});
Clazz.defineMethod (c$, "setState", 
function (state) {
this._state = state;
}, "~N");
Clazz.defineMethod (c$, "setGroup", 
function (group) {
this._group = group;
this._group.update (this);
}, "com.sparshui.server.Group");
Clazz.defineMethod (c$, "update", 
function (location, time, state) {
this._location = location;
this._state = state;
this._changed = true;
this._time = time;
if (this._group != null) this._group.update (this);
}, "com.sparshui.common.Location,~N,~N");
Clazz.defineMethod (c$, "resetChanged", 
function () {
this._changed = false;
});
Clazz.defineMethod (c$, "isChanged", 
function () {
return this._changed;
});
Clazz.overrideMethod (c$, "clone", 
function () {
return  new com.sparshui.server.TouchPoint (this);
});
Clazz.defineMethod (c$, "isNear", 
function (tp) {
return (Math.abs (this._location.getX () - tp._location.getX ()) < 0.005 && Math.abs (this._location.getY () - tp._location.getY ()) < 0.005);
}, "com.sparshui.server.TouchPoint");
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
