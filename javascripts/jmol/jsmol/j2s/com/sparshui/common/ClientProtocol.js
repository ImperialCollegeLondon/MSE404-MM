Clazz.declarePackage ("com.sparshui.common");
Clazz.load (null, "com.sparshui.common.ClientProtocol", ["java.io.DataInputStream", "$.DataOutputStream"], function () {
var c$ = Clazz.decorateAsClass (function () {
this._socket = null;
this._in = null;
this._out = null;
if (!Clazz.isClassDefined ("com.sparshui.common.ClientProtocol.MessageType")) {
com.sparshui.common.ClientProtocol.$ClientProtocol$MessageType$ ();
}
Clazz.instantialize (this, arguments);
}, com.sparshui.common, "ClientProtocol");
Clazz.makeConstructor (c$, 
function (socket) {
this._socket = socket;
this._in =  new java.io.DataInputStream (this._socket.getInputStream ());
this._out =  new java.io.DataOutputStream (this._socket.getOutputStream ());
}, "java.net.Socket");
c$.$ClientProtocol$MessageType$ = function () {
/*if4*/;(function(){
var c$ = Clazz.decorateAsClass (function () {
Clazz.prepareCallback (this, arguments);
Clazz.instantialize (this, arguments);
}, com.sparshui.common.ClientProtocol, "MessageType");
/*eoif4*/})();
};
});
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
