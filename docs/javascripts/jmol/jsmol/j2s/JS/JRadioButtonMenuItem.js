Clazz.declarePackage ("JS");
Clazz.load (["JS.JMenuItem"], "JS.JRadioButtonMenuItem", null, function () {
var c$ = Clazz.decorateAsClass (function () {
this.isRadio = true;
Clazz.instantialize (this, arguments);
}, JS, "JRadioButtonMenuItem", JS.JMenuItem);
Clazz.makeConstructor (c$, 
function () {
Clazz.superConstructor (this, JS.JRadioButtonMenuItem, ["rad", 3]);
});
});
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
