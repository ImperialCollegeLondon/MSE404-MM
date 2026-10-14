Clazz.declarePackage ("com.sparshui.client");
Clazz.load (["java.lang.Thread"], "com.sparshui.client.ClientServerConnection", ["com.sparshui.client.ClientToServerProtocol", "java.io.DataOutputStream", "java.net.Socket"], function () {
var c$ = Clazz.decorateAsClass (function () {
this._client = null;
this._socket = null;
this._protocol = null;
Clazz.instantialize (this, arguments);
}, com.sparshui.client, "ClientServerConnection", Thread);
Clazz.makeConstructor (c$, 
function (address, client) {
Clazz.superConstructor (this, com.sparshui.client.ClientServerConnection, []);
this._client = client;
this._socket =  new java.net.Socket (address, 5946);
var out =  new java.io.DataOutputStream (this._socket.getOutputStream ());
out.writeByte (0);
this._protocol =  new com.sparshui.client.ClientToServerProtocol (this._socket);
this.start ();
}, "~S,com.sparshui.client.SparshClient");
Clazz.overrideMethod (c$, "run", 
function () {
Thread.currentThread ().setName ("SparshUI Client->ServerConnection");
while (this._socket.isConnected ()) {
if (!this._protocol.processRequest (this._client)) break;
}
});
Clazz.defineMethod (c$, "close", 
function () {
try {
this._socket.close ();
} catch (e) {
if (Clazz.exceptionOf(e,"java.io.IOException")){
} else {
throw e;
}
}
});
});
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
