Clazz.declarePackage ("java.net");
var c$ = Clazz.declareType (java.net, "URLEncoder");
c$.encode = Clazz.defineMethod (c$, "encode", 
function (s) {
return encodeURIComponent(s);
}, "~S");
;//5.0.1-v2 Sat Nov 25 17:52:34 CST 2023
