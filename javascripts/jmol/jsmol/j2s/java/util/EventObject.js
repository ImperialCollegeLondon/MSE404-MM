Clazz.load (null, "java.util.EventObject", ["java.lang.IllegalArgumentException"], function () {
var c$ = Clazz.decorateAsClass (function () {
this.source = null;
Clazz.instantialize (this, arguments);
}, java.util, "EventObject", null, java.io.Serializable);
Clazz.makeConstructor (c$, 
function (source) {
if (source != null) this.source = source;
 else throw  new IllegalArgumentException ();
}, "~O");
Clazz.defineMethod (c$, "getSource", 
function () {
return this.source;
});
Clazz.overrideMethod (c$, "toString", 
function () {
return this.getClass ().getName () + "[source=" + String.valueOf (this.source) + ']';
});
});
;//5.0.1-v2 Sat Nov 25 17:52:34 CST 2023
