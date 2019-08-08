"use strict";
(function(){

window.renderers = window.renderers || {};

window.renderers.lines = {};
console.log("Set lines.");

const colorProgram = initShaderProgram(
`
	uniform mat4 mvp_mat4;

	attribute vec4 Position;
	attribute vec4 Color;

	varying vec4 color;

	void main() {
		gl_Position = mvp_mat4 * Position;
		color = Color;
	}
`,`
	varying lowp vec4 color;
	void main() {
		gl_FragColor = color;
	}
`);

const yarnsBuffer = gl.createBuffer();
const yarnStarts = [];

window.renderers.lines.uploadYarns = function lines_uploadYarns() {
	//clear yarnStarts array:
	yarnStarts.splice(0,yarnStarts.length);

	//accmulate Position attribs for all yarns:
	let Positions = [];
	yarns.yarns.forEach(function(yarn){
		yarnStarts.push(Positions.length / 3);
		Positions = Positions.concat(yarn.points);
	});
	yarnStarts.push(Positions.length / 3);

	gl.bindBuffer(gl.ARRAY_BUFFER, yarnsBuffer);
	gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(Positions), gl.STATIC_DRAW);
};

window.renderers.lines.redraw = function lines_redraw() {

	gl.clearColor(0.1, 0.1, 0.1, 1.0);
	gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

	gl.useProgram(colorProgram.program);
	gl.bindBuffer(gl.ARRAY_BUFFER, yarnsBuffer);

	gl.vertexAttrib4f(colorProgram.attribLocations.Color, 1.0, 0.0, 1.0, 1.0);
	gl.disableVertexAttribArray(colorProgram.attribLocations.Color);

	gl.vertexAttribPointer(colorProgram.attribLocations.Position,
		3, //size
		gl.FLOAT, //type
		false, //normalize
		12, //stride
		0 //offset
	);
	gl.enableVertexAttribArray(colorProgram.attribLocations.Position);

	gl.uniformMatrix4fv(
		colorProgram.uniformLocations.mvp_mat4,
		false,
		computeMVP()
	);

	for (let i = 0; i + 1 < yarnStarts.length; ++i) {
		gl.drawArrays(gl.LINE_STRIP, yarnStarts[i], yarnStarts[i+1]-yarnStarts[i]);
	}
};


})();
