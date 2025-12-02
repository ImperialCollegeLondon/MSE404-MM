Clazz.declarePackage ("JM");
Clazz.load (["JM.BioPolymer"], "JM.PhosphorusPolymer", null, function () {
var c$ = Clazz.declareType (JM, "PhosphorusPolymer", JM.BioPolymer);
Clazz.makeConstructor (c$, 
function (monomers) {
Clazz.superConstructor (this, JM.PhosphorusPolymer, []);
this.set (monomers);
this.hasStructure = true;
}, "~A");
});
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
