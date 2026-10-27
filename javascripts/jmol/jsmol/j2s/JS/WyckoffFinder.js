Clazz.declarePackage ("JS");
Clazz.load (["java.util.Hashtable", "JU.V3"], "JS.WyckoffFinder", ["JU.JSJSONParser", "$.Measure", "$.P3", "$.P4", "$.PT", "$.Rdr", "JU.SimpleUnitCell", "JV.FileManager"], function () {
var c$ = Clazz.decorateAsClass (function () {
this.positions = null;
this.centerings = null;
Clazz.instantialize (this, arguments);
}, JS, "WyckoffFinder");
Clazz.makeConstructor (c$, 
function () {
});
Clazz.defineMethod (c$, "getWyckoffFinder", 
function (vwr, sgname) {
var helper = JS.WyckoffFinder.helpers.get (sgname);
if (helper == null) {
var itno = JU.PT.parseInt (JU.PT.split (sgname, ":")[0]);
if (itno >= 1 && itno <= 230) {
var resource = this.getResource (vwr, "ita_" + itno + ".json");
if (resource != null) {
var its = resource.get ("its");
if (its != null) {
for (var i = its.size (); --i >= 0; ) {
var map = its.get (i);
if (sgname.equals (map.get ("itaFull"))) {
JS.WyckoffFinder.helpers.put (sgname, helper =  new JS.WyckoffFinder (map));
return helper;
}}
}}}}if (helper == null) {
if (JS.WyckoffFinder.nullHelper == null) JS.WyckoffFinder.nullHelper =  new JS.WyckoffFinder (null);
JS.WyckoffFinder.helpers.put (sgname, JS.WyckoffFinder.nullHelper);
}return helper;
}, "JV.Viewer,~S");
Clazz.makeConstructor (c$, 
function (map) {
if (map != null) {
var wpos = map.get ("wpos");
this.positions = wpos.get ("pos");
var cent = wpos.get ("cent");
if (cent != null) {
this.centerings =  new Array (cent.size ());
for (var i = cent.size (); --i >= 0; ) {
this.centerings[i] = JS.WyckoffFinder.toPoint (cent.get (i));
}
}}}, "java.util.Map");
Clazz.defineMethod (c$, "getWyckoffPosition", 
function (p) {
if (this.positions == null) return "?";
for (var i = this.positions.size (); --i >= 0; ) {
var map = this.positions.get (i);
if (i == 0) {
return map.get ("label");
}var coords = map.get ("coord");
for (var c = 0, n = coords.size (); c < n; c++) {
if (JS.WyckoffFinder.getWyckoffCoord (coords, c).contains (p, this.centerings)) {
return map.get ("label");
}}
}
return "?";
}, "JU.P3");
Clazz.defineMethod (c$, "findPositionFor", 
function (p, letter) {
if (this.positions == null) return null;
for (var i = this.positions.size (); --i >= 0; ) {
var map = this.positions.get (i);
if (map.get ("label").equals (letter)) {
var coords = map.get ("coord");
if (coords != null) JS.WyckoffFinder.getWyckoffCoord (coords, 0).set (p);
return p;
}}
return null;
}, "JU.P3,~S");
c$.getWyckoffCoord = Clazz.defineMethod (c$, "getWyckoffCoord", 
function (coords, c) {
var coord = coords.get (c);
if ((typeof(coord)=='string')) {
coords.set (c, coord =  new JS.WyckoffFinder.WyckoffPos (coord));
}return coord;
}, "JU.Lst,~N");
Clazz.defineMethod (c$, "getResource", 
function (vwr, resource) {
try {
var r = JV.FileManager.getBufferedReaderForResource (vwr, this, "JS/", "sg/json/" + resource);
var data =  new Array (1);
if (JU.Rdr.readAllAsString (r, 2147483647, false, data, 0)) {
return  new JU.JSJSONParser ().parse (data[0], true);
}} catch (e) {
System.err.println (e.getMessage ());
}
return null;
}, "JV.Viewer,~S");
c$.toPoint = Clazz.defineMethod (c$, "toPoint", 
function (xyz) {
var s = JU.PT.split (xyz, ",");
return JU.P3.new3 (JU.PT.parseFloatFraction (s[0]), JU.PT.parseFloatFraction (s[1]), JU.PT.parseFloatFraction (s[2]));
}, "~S");
/*if3*/;(function(){
var c$ = Clazz.decorateAsClass (function () {
this.point = null;
this.line = null;
this.plane = null;
this.type = 0;
this.xyz = null;
Clazz.instantialize (this, arguments);
}, JS.WyckoffFinder, "WyckoffPos");
c$.unitize = Clazz.defineMethod (c$, "unitize", 
function (p) {
JU.SimpleUnitCell.unitizeDim (3, p);
return p;
}, "JU.P3");
Clazz.makeConstructor (c$, 
function (xyz) {
this.create (xyz);
}, "~S");
Clazz.defineMethod (c$, "create", 
function (p) {
var xyz = JU.PT.split (p, ",");
var nxyz = 0;
for (var i = 0; i < 3; i++) {
if (xyz[i].indexOf ('x') >= 0) {
nxyz |= 1;
}if (xyz[i].indexOf ('y') >= 0) {
nxyz |= 2;
}if (xyz[i].indexOf ('z') >= 0) {
nxyz |= 4;
}}
var v1;
var v2;
var v3;
switch (nxyz) {
case 0:
this.type = 1;
this.point = JS.WyckoffFinder.toPoint (p);
break;
case 1:
case 2:
case 4:
this.type = 2;
v1 = JS.WyckoffFinder.WyckoffPos.ptFor (p, 0, 0, 0);
v2 = JS.WyckoffFinder.WyckoffPos.ptFor (p, 1, 1.27, 1.64);
v2.sub2 (v2, v1);
v2.normalize ();
this.point = JS.WyckoffFinder.WyckoffPos.unitize (JU.P3.newP (v1));
this.line = JU.V3.newV (v2);
break;
case 3:
case 5:
case 6:
this.type = 3;
v1 = JS.WyckoffFinder.WyckoffPos.ptFor (p, 0, 0, 0);
v2 = JS.WyckoffFinder.WyckoffPos.ptFor (p, 1.23, 1.47, 1.86);
v3 = JS.WyckoffFinder.WyckoffPos.ptFor (p, 0.1, 0.2, 0.3);
this.plane = JU.Measure.getPlaneThroughPoints (v1, v2, v3, null, null,  new JU.P4 ());
break;
case 7:
break;
}
}, "~S");
c$.ptFor = Clazz.defineMethod (c$, "ptFor", 
function (p, x, y, z) {
var v = JU.PT.split (p, ",");
var a = JS.WyckoffFinder.WyckoffPos.decodeXYZ (v[0], x, y, z);
var b = JS.WyckoffFinder.WyckoffPos.decodeXYZ (v[1], x, y, z);
var c = JS.WyckoffFinder.WyckoffPos.decodeXYZ (v[2], x, y, z);
return JU.P3.new3 (a, b, c);
}, "~S,~N,~N,~N");
c$.decodeXYZ = Clazz.defineMethod (c$, "decodeXYZ", 
function (s, x, y, z) {
s = JU.PT.rep (s, "-", "+-");
s = JU.PT.rep (s, "x", "*x");
s = JU.PT.rep (s, "y", "*y");
s = JU.PT.rep (s, "z", "*z");
s = JU.PT.rep (s, "-*", "-");
s = JU.PT.rep (s, "+*", "+");
var r = 0;
var parts = JU.PT.split (s, "+");
for (var p = parts.length; --p >= 0; ) {
s = parts[p];
if (s.length == 0) continue;
if (s.indexOf ('.') >= 0) {
r += JU.PT.parseFloat (s);
continue;
}var v = 0;
var f2 = 0;
var i0 = 0;
var f = 1;
switch ((s.charAt (0)).charCodeAt(0)) {
case 45:
f = -1;
case 42:
i0++;
break;
}
for (var i = s.length; --i >= i0; ) {
var c = s.charAt (i);
switch ((c).charCodeAt(0)) {
case 120:
v = x;
break;
case 121:
v = y;
break;
case 122:
v = z;
break;
case 47:
v = 1 / v;
case 42:
f *= v;
v = 0;
break;
default:
var u = "0123456789".indexOf (c);
if (u < 0) System.err.println ("WH ????");
if (v == 0) {
v = u;
} else {
f2 = (f2 == 0 ? 10 : f2 * 10);
v += f2 * u;
}break;
}
}
r += f * v;
}
return r;
}, "~S,~N,~N,~N");
Clazz.defineMethod (c$, "contains", 
function (p, centerings) {
if (this.containsPt (p)) return true;
var pc =  new JU.P3 ();
if (centerings != null) for (var i = centerings.length; --i >= 0; ) {
pc.add2 (p, centerings[i]);
JS.WyckoffFinder.WyckoffPos.unitize (pc);
if (this.containsPt (pc)) return true;
}
return false;
}, "JU.P3,~A");
Clazz.defineMethod (c$, "containsPt", 
function (p) {
var d = 1;
switch (this.type) {
case 1:
d = p.distance (this.point);
break;
case 2:
var p1 = JU.P3.newP (p);
JU.Measure.projectOntoAxis (p1, this.point, this.line, JS.WyckoffFinder.WyckoffPos.vtemp1);
d = p.distance (p1);
break;
case 3:
d = JU.Measure.distanceToPlane (this.plane, p);
break;
}
return JS.WyckoffFinder.WyckoffPos.approx (d) == 0;
}, "JU.P3");
Clazz.defineMethod (c$, "set", 
function (p) {
switch (this.type) {
case 1:
p.setT (this.point);
break;
case 2:
JU.Measure.projectOntoAxis (p, this.point, this.line, JS.WyckoffFinder.WyckoffPos.vtemp1);
break;
case 3:
JU.Measure.getPlaneProjection (p, this.plane, JS.WyckoffFinder.WyckoffPos.vtemp1, JS.WyckoffFinder.WyckoffPos.vtemp1);
p.setT (JS.WyckoffFinder.WyckoffPos.vtemp1);
break;
}
}, "JU.P3");
c$.approx = Clazz.defineMethod (c$, "approx", 
function (d) {
return JU.PT.approx (d, 1000);
}, "~N");
c$.vtemp1 = c$.prototype.vtemp1 =  new JU.V3 ();
/*eoif3*/})();
c$.helpers = c$.prototype.helpers =  new java.util.Hashtable ();
Clazz.defineStatics (c$,
"nullHelper", null);
});
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
