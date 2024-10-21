Clazz.declarePackage ("JS");
var c$ = Clazz.decorateAsClass (function () {
this.top = 0;
this.left = 0;
this.bottom = 0;
this.right = 0;
Clazz.instantialize (this, arguments);
}, JS, "Insets");
Clazz.makeConstructor (c$, 
function (top, left, bottom, right) {
this.top = top;
this.left = left;
this.bottom = bottom;
this.right = right;
}, "~N,~N,~N,~N");
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
