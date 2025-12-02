Clazz.declarePackage ("JM");
Clazz.load (["JM.MinObject"], "JM.MinTorsion", null, function () {
var c$ = Clazz.declareType (JM, "MinTorsion", JM.MinObject);
Clazz.makeConstructor (c$, 
function (data) {
Clazz.superConstructor (this, JM.MinTorsion, []);
this.data = data;
}, "~A");
});
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
