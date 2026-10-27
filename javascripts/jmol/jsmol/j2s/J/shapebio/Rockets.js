Clazz.declarePackage ("J.shapebio");
Clazz.load (["J.shapebio.BioShapeCollection"], "J.shapebio.Rockets", null, function () {
var c$ = Clazz.declareType (J.shapebio, "Rockets", J.shapebio.BioShapeCollection);
Clazz.overrideMethod (c$, "initShape", 
function () {
this.madTurnRandom = 500;
});
});
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
