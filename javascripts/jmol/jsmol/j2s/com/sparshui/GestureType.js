Clazz.declarePackage ("com.sparshui");
var c$ = Clazz.decorateAsClass (function () {
this.sType = null;
this.iType = 2147483647;
Clazz.instantialize (this, arguments);
}, com.sparshui, "GestureType");
Clazz.makeConstructor (c$, 
function (type) {
this.sType = type;
}, "~S");
Clazz.makeConstructor (c$, 
function (type) {
this.iType = type;
}, "~N");
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
