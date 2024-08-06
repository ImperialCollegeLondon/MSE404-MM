Clazz.declarePackage ("JSV.common");
var c$ = Clazz.declareType (JSV.common, "CoordComparator", null, java.util.Comparator);
Clazz.overrideMethod (c$, "compare", 
function (c1, c2) {
return (c1.getXVal () > c2.getXVal () ? 1 : c1.getXVal () < c2.getXVal () ? -1 : 0);
}, "JSV.common.Coordinate,JSV.common.Coordinate");
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
