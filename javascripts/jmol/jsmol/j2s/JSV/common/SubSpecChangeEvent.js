Clazz.declarePackage ("JSV.common");
var c$ = Clazz.decorateAsClass (function () {
this.isub = 0;
this.title = null;
Clazz.instantialize (this, arguments);
}, JSV.common, "SubSpecChangeEvent");
Clazz.makeConstructor (c$, 
function (isub, title) {
this.isub = isub;
this.title = title;
}, "~N,~S");
Clazz.defineMethod (c$, "isValid", 
function () {
return (this.title != null);
});
Clazz.defineMethod (c$, "getSubIndex", 
function () {
return this.isub;
});
Clazz.overrideMethod (c$, "toString", 
function () {
return this.title;
});
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
