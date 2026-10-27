Clazz.declarePackage ("JS");
Clazz.load (["JS.ScriptException"], "JS.ScriptInterruption", null, function () {
var c$ = Clazz.declareType (JS, "ScriptInterruption", JS.ScriptException);
Clazz.makeConstructor (c$, 
function (eval, why, millis) {
Clazz.superConstructor (this, JS.ScriptInterruption, [eval, why, "!", millis == -2147483648 || eval.vwr.autoExit]);
if (why.equals ("delay")) eval.delayScript (millis);
}, "JS.ScriptEval,~S,~N");
});
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
