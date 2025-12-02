Clazz.declarePackage ("JS");
Clazz.load (["java.lang.Exception"], "JS.InvalidSmilesException", null, function () {
var c$ = Clazz.declareType (JS, "InvalidSmilesException", Exception);
c$.getLastError = Clazz.defineMethod (c$, "getLastError", 
function () {
return JS.InvalidSmilesException.lastError;
});
c$.clear = Clazz.defineMethod (c$, "clear", 
function () {
JS.InvalidSmilesException.lastError = null;
});
Clazz.overrideMethod (c$, "getMessage", 
function () {
return JS.InvalidSmilesException.lastError;
});
Clazz.makeConstructor (c$, 
function (message) {
Clazz.superConstructor (this, JS.InvalidSmilesException, [message]);
JS.InvalidSmilesException.lastError = (message.startsWith ("Jmol SMILES") ? message : "Jmol SMILES Exception: " + message);
}, "~S");
Clazz.defineStatics (c$,
"lastError", null);
});
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
