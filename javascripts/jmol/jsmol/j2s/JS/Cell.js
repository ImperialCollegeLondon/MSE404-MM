Clazz.declarePackage ("JS");
var c$ = Clazz.decorateAsClass (function () {
this.component = null;
this.colspan = 0;
this.rowspan = 0;
this.textAlign = 0;
this.c = null;
Clazz.instantialize (this, arguments);
}, JS, "Cell");
Clazz.makeConstructor (c$, 
function (btn, c) {
this.component = btn;
this.colspan = c.gridwidth;
this.rowspan = c.gridheight;
this.c = c;
}, "JS.JComponent,JS.GridBagConstraints");
Clazz.defineMethod (c$, "toHTML", 
function (id) {
var style = this.c.getStyle (false);
return "<td id='" + id + "' " + (this.colspan < 2 ? "" : "colspan='" + this.colspan + "' ") + style + "><span " + this.c.getStyle (true) + ">" + this.component.toHTML () + "</span></td>";
}, "~S");
;//5.0.1-v2 Sat Nov 25 17:51:22 CST 2023
