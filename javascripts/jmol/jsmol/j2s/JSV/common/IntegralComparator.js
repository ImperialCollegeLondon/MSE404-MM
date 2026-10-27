Clazz.declarePackage ("JSV.common");
var c$ = Clazz.declareType (JSV.common, "IntegralComparator", null, java.util.Comparator);
Clazz.overrideMethod (c$, "compare", 
function (m1, m2) {
return (m1.getXVal () < m2.getXVal () ? -1 : m1.getXVal () > m2.getXVal () ? 1 : 0);
}, "JSV.common.Measurement,JSV.common.Measurement");
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
