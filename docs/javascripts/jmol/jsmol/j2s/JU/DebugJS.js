Clazz.declarePackage ("JU");
var c$ = Clazz.declareType (JU, "DebugJS");
c$._ = Clazz.defineMethod (c$, "_", 
function (msg) {
{
if (Clazz._debugging) {
debugger;
}
}}, "~S");
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
