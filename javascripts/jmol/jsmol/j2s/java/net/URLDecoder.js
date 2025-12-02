Clazz.declarePackage ("java.net");
var c$ = Clazz.declareType (java.net, "URLDecoder");
c$.decode = Clazz.defineMethod (c$, "decode", 
function (s) {
return decodeURIComponent(s);
}, "~S");
;//5.0.1-v2 Sat Nov 25 17:52:34 CST 2023
