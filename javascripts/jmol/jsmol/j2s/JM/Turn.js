Clazz.declarePackage ("JM");
Clazz.load (["JM.ProteinStructure"], "JM.Turn", ["J.c.STR"], function () {
var c$ = Clazz.declareType (JM, "Turn", JM.ProteinStructure);
Clazz.overrideConstructor (c$, 
function (apolymer, monomerIndex, monomerCount) {
this.setupPS (apolymer, J.c.STR.TURN, monomerIndex, monomerCount);
this.subtype = J.c.STR.TURN;
}, "JM.AlphaPolymer,~N,~N");
});
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
