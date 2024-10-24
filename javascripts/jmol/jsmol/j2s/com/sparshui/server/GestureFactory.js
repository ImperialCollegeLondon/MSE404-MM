Clazz.declarePackage ("com.sparshui.server");
Clazz.load (null, "com.sparshui.server.GestureFactory", ["JU.Logger"], function () {
var c$ = Clazz.declareType (com.sparshui.server, "GestureFactory");
c$.createGesture = Clazz.defineMethod (c$, "createGesture", 
function (gType) {
if (gType.sType != null) {
try {
return Clazz._4Name (gType.sType).newInstance ();
} catch (e) {
if (Clazz.exceptionOf(e, Exception)){
JU.Logger.error ("[GestureFactory] Error creating instance for " + gType.sType + ": \n" + e.getMessage ());
} else {
throw e;
}
}
return null;
}JU.Logger.error ("[GestureFactory] Gesture not recognized: " + gType.iType);
return null;
}, "com.sparshui.GestureType");
});
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
