Clazz.declarePackage ("J.rendersurface");
Clazz.load (["J.rendersurface.IsosurfaceRenderer"], "J.rendersurface.PmeshRenderer", null, function () {
var c$ = Clazz.declareType (J.rendersurface, "PmeshRenderer", J.rendersurface.IsosurfaceRenderer);
Clazz.overrideMethod (c$, "render", 
function () {
return this.renderIso ();
});
});
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
