Clazz.declarePackage ("JS");
var c$ = Clazz.decorateAsClass (function () {
this.width = 0;
this.height = 0;
Clazz.instantialize (this, arguments);
}, JS, "Dimension");
Clazz.makeConstructor (c$, 
function (w, h) {
this.set (w, h);
}, "~N,~N");
Clazz.defineMethod (c$, "set", 
function (w, h) {
this.width = w;
this.height = h;
return this;
}, "~N,~N");
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
