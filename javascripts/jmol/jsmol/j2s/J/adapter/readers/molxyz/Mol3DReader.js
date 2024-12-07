Clazz.declarePackage ("J.adapter.readers.molxyz");
Clazz.load (["J.adapter.readers.molxyz.MolReader"], "J.adapter.readers.molxyz.Mol3DReader", null, function () {
var c$ = Clazz.declareType (J.adapter.readers.molxyz, "Mol3DReader", J.adapter.readers.molxyz.MolReader);
Clazz.overrideMethod (c$, "initializeReader", 
function () {
this.allow2D = false;
});
});
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
