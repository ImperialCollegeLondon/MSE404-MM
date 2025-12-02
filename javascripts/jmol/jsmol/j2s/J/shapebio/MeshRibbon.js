Clazz.declarePackage ("J.shapebio");
Clazz.load (["J.shapebio.Strands"], "J.shapebio.MeshRibbon", null, function () {
var c$ = Clazz.declareType (J.shapebio, "MeshRibbon", J.shapebio.Strands);
Clazz.overrideMethod (c$, "initShape", 
function () {
this.isMesh = true;
});
});
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
