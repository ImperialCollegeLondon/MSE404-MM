Clazz.declarePackage ("JS");
var c$ = Clazz.decorateAsClass (function () {
this.processName = null;
this.context = null;
Clazz.instantialize (this, arguments);
}, JS, "ScriptProcess");
Clazz.makeConstructor (c$, 
function (name, context) {
this.processName = name;
this.context = context;
}, "~S,JS.ScriptContext");
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
