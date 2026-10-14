Clazz.declarePackage ("javajs.api");
var c$ = Clazz.declareType (javajs.api, "Interface");
c$.getInterface = Clazz.defineMethod (c$, "getInterface", 
function (name) {
try {
var x = Clazz._4Name (name);
return (x == null ? null : x.newInstance ());
} catch (e) {
if (Clazz.exceptionOf(e, Exception)){
System.out.println ("Interface.java Error creating instance for " + name + ": \n" + e);
return null;
} else {
throw e;
}
}
}, "~S");
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
