"use strict";
(function(){

/*
 * renderer using on-gpu extruded line strips.
 * attribute setup:
 *  - instanced rendering
 *  - attribs that advance one per instance:
 *    - each instance has quadratic control points in attribs (use divisor)
 *    - each instance has length as well (??) (use divisor)
 *  - attribs that advance once per vertex (re-used for all instances):
 *    - parameter position
 *    - left/right indicator
 *  - uniforms
 *    - something about travelling strand min/max/frequency?
 *
 * Questions:
 *  - how to deal with multiple strands?
 *    - instance per strand?
 *    - instance per segment? (many strands in one instance)
 */

window.renderers = window.renderers || {};

let renderer = window.renderers.strips = {};

const stripProgram = initShaderProgram(
`
	uniform mat4 MVP_mat4;
	uniform vec3 EYE_vec3;
	uniform vec3 OUT_vec3;

	attribute float Radius;
	attribute vec4 Start;
	attribute vec4 Middle;
	attribute vec4 End;
	attribute vec2 ParamSide;

	varying vec4 color;
	varying vec3 normal;

	void main() {
		vec3 sm = mix(Start.xyz, Middle.xyz, ParamSide.x);
		vec3 me = mix(Middle.xyz, End.xyz, ParamSide.x);
		vec3 sme = mix(sm, me, ParamSide.x);

		vec3 to_eye = EYE_vec3 - sme;

		vec3 along = mix(Middle.xyz - Start.xyz, End.xyz - Middle.xyz, ParamSide.x);

		sme += Radius * normalize(cross(along, to_eye)) * ParamSide.y;

		//along -= dot(along, OUT_vec3) * OUT_vec3;
		//sme += Radius * normalize(cross(along, OUT_vec3)) * ParamSide.y;

		gl_Position = MVP_mat4 * vec4(sme, 1.0);

		color = vec4(0.5+0.5*normalize(along),1); //Color;
		normal = vec3(1,1,1);
	}
`,`
	varying lowp vec4 color;
	varying highp vec3 normal;
	void main() {
		highp vec3 n = normalize(normal);
		//upper dome light:
		highp vec3 light = mix(vec3(0.1), vec3(1.0), 0.5*n.y+0.5);
		gl_FragColor = vec4(color.rgb * light, color.a);
	}
`);

const perInstanceBuffer = gl.createBuffer();
let instanceCount = 0;
const perVertexBuffer = gl.createBuffer();
let vertexCount = 0;

function quadraticBSpline(points, minDistance, maxAngle, splits) {
	for (let i = 0; i < splits.length; ++i) {
		if (splits[i] < 0 || 3 * splits[i] + 2 >= points.length) {
			throw "split out of range";
		}
		if (i > 0 && !(splits[i-1] < splits[i])) {
			throw "non-monotonic splits";
		}
	}
	if (points.length <= 3) {
		let lengths = [];
		if (splits.length) {
			lengths.push(0.0);
		}
		return {
			points:points.slice(),
			splits:splits.slice(),
			lengths:lengths
		};
	}

	//quadratic b-spline through the points (same as bezier through midpoints)

	let poly = [];

	let threshold = Math.cos(maxAngle * Math.PI / 180.0);
	function curveTo(m, e, mSourceIndex, eSourceIndex) {
		let s = {
			x:poly[poly.length-3],
			y:poly[poly.length-2],
			z:poly[poly.length-1]
		};
		//recursive subdivision version:
		const a = {
			x:(m.x - s.x),
			y:(m.y - s.y),
			z:(m.z - s.z)
		};
		const b = {
			x:(e.x - m.x),
			y:(e.y - m.y),
			z:(e.z - m.z)
		};
		let aLen2 = (a.x * a.x + a.y * a.y + a.z * a.z);
		let bLen2 = (b.x * b.x + b.y * b.y + b.z * b.z);
		let d = (a.x * b.x + a.y * b.y + a.z * b.z);
		if ( //if segment is short:
			(aLen2 <= minDistance * minDistance && bLen2 <= minDistance * minDistance)
			//or segment doesn't bend too much:
			|| (d >= 0 && d * d >= threshold * threshold * aLen2 * bLen2)) {
			//just draw curve:
			if (mSourceIndex >= 0) { //actually include midpoint if source index is associated
				nextPt(m.x, m.y, m.z, mSourceIndex);
			}
			nextPt(e.x, e.y, e.z, eSourceIndex);
		} else {
			//otherwise, subdivide:
			let sm = {
				x:0.5 * (s.x + m.x),
				y:0.5 * (s.y + m.y),
				z:0.5 * (s.z + m.z)
			};
			let me = {
				x:0.5 * (m.x + e.x),
				y:0.5 * (m.y + e.y),
				z:0.5 * (m.z + e.z)
			};
			let sme = {
				x:0.5 * (sm.x + me.x),
				y:0.5 * (sm.y + me.y),
				z:0.5 * (sm.z + me.z)
			};
			curveTo(sm, sme, -1, mSourceIndex);
			curveTo(me, e, -1, eSourceIndex);
		}
	}

	let polySplits = [];
	let polyLengths = [];
	let nextSplit = 0;

	function nextPt(x,y,z, sourceIndex) {
		//accumulate length to this point:
		if (poly.length && polyLengths.length) {
			polyLengths[polyLengths.length-1] += Math.sqrt(
				(x-poly[poly.length-3])*(x-poly[poly.length-3])
				+ (y-poly[poly.length-2])*(y-poly[poly.length-2])
				+ (z-poly[poly.length-1])*(z-poly[poly.length-1])
			);
		}
		//add a new split index if index matches next split index:
		if (nextSplit < splits.length) {
			if (sourceIndex === splits[nextSplit]) {
				polySplits.push(poly.length/3);
				polyLengths.push(0.0);
				++nextSplit;
			} else {
				console.assert(sourceIndex < splits[nextSplit], "must not pass");
			}
		}
		//add this point to the polyline:
		poly.push(x,y,z);
	}

	//straight segment to first midpoint:
	nextPt(points[0], points[1], points[2], 0);

	nextPt(0.5 * (points[0] + points[3]), 0.5 * (points[1] + points[4]), 0.5 * (points[2] + points[5]), -1);
	//curved segments for the rest of the yarn:
	for (let i = 3; i + 5 < points.length; i += 3) {
		//bezier control point: vertex on curve
		let m = {
			x:points[i],
			y:points[i+1],
			z:points[i+2]
		};
		//bezier endpoint: midpoint of line
		let e = {
			x:0.5*(points[i+0]+points[i+3]),
			y:0.5*(points[i+1]+points[i+4]),
			z:0.5*(points[i+2]+points[i+5])
		};
		curveTo(m, e, (i/3), -1);
	}
	//straight segment from last midpoint to end:
	nextPt(
		points[points.length-3],
		points[points.length-2],
		points[points.length-1],
		(points.length-3)/3
	);

	return {
		points:poly,
		splits:polySplits,
		lengths:polyLengths
	};
}

//when using a data view to store Float32 values, do they need to be marked as little endian?
const FLOAT32_LITTLE_ENDIAN = function() {
	let array = new ArrayBuffer(4);
	let f32 = new Float32Array(array);
	let data = new DataView(array);
	f32[0] = 1.0;
	if (data.getFloat32(0, false) === f32[0]) {
		assert(data.getFloat32(0, true) !== f32[0]);
		return false;
	} else if (data.getFloat32(0, true) === f32[0]) {
		assert(data.getFloat32(0, false) !== f32[0]);
		return true;
	}
};

const BYTES_PER_ATTRIB = 3*4 + 3*4 + 4*1;

function addSpline(start, middle, end /*olorOrSplits_, colors_*/) {
	/*
	let colors;
	let splits;

	if (typeof(colors_) !== 'undefined') {
		colors = colors_;
		splits = colorOrSplits_;
		if (colors.length !== splits.length) throw "expecting array of colors + array of splits";
		if (splits[0] !== 0) throw "splits should start at the start";
	} else {
		colors = [colorOrSplits_];
		splits = [0];
	}
	*/

	let tristrip = new ArrayBuffer(attribBytes);
	let data = new DataView(tristrip);
	let dataOffset = 0;


	function attrib(x,y,z, nx,ny,nz, rgba) {
		data.setFloat32(dataOffset, x, FLOAT32_LITTLE_ENDIAN); dataOffset += 4;
		data.setFloat32(dataOffset, y, FLOAT32_LITTLE_ENDIAN); dataOffset += 4;
		data.setFloat32(dataOffset, z, FLOAT32_LITTLE_ENDIAN); dataOffset += 4;
		data.setFloat32(dataOffset, nx, FLOAT32_LITTLE_ENDIAN); dataOffset += 4;
		data.setFloat32(dataOffset, ny, FLOAT32_LITTLE_ENDIAN); dataOffset += 4;
		data.setFloat32(dataOffset, nz, FLOAT32_LITTLE_ENDIAN); dataOffset += 4;
		data.setUint32(dataOffset, rgba, false); dataOffset += 4;
	}
	function prevAttrib(offset) {
		console.assert(offset < 0, "must look before current point");
		offset *= BYTES_PER_ATTRIB;
		console.assert(offset + dataOffset >= 0, "can't look before beginning");
		attrib(
			data.getFloat32(offset + dataOffset + 0, FLOAT32_LITTLE_ENDIAN),
			data.getFloat32(offset + dataOffset + 4, FLOAT32_LITTLE_ENDIAN),
			data.getFloat32(offset + dataOffset + 8, FLOAT32_LITTLE_ENDIAN),
			data.getFloat32(offset + dataOffset + 12, FLOAT32_LITTLE_ENDIAN),
			data.getFloat32(offset + dataOffset + 16, FLOAT32_LITTLE_ENDIAN),
			data.getFloat32(offset + dataOffset + 20, FLOAT32_LITTLE_ENDIAN),
			data.getUint32(offset + dataOffset + 24, false)
		);
	}

	//sweep the cross-section shape along the yarn:
	//..will do this by maintaining a local coordinate system p1,p2 at corners:
	let p1 = { x:NaN, y:NaN, z:NaN};
	let p2 = { x:NaN, y:NaN, z:NaN};
	let pt = { x:spine[0], y:spine[1], z:spine[2] };

	/* DEBUG
	const ColorA = 0xff0000ff;
	const ColorB = 0x00ff00ff;
	const ColorC = 0x0000ffff;
	const ColorD = 0xffffffff;
	*/

	const colorTempArray = new ArrayBuffer(4);
	const colorTempView = new DataView(colorTempArray);
	function rgbaToUint32(color) {
		colorTempView.setUint8(0, Math.max(0x00, Math.min(0xff, Math.round(color.r) )) );
		colorTempView.setUint8(1, Math.max(0x00, Math.min(0xff, Math.round(color.g) )) );
		colorTempView.setUint8(2, Math.max(0x00, Math.min(0xff, Math.round(color.b) )) );
		colorTempView.setUint8(3, Math.max(0x00, Math.min(0xff, Math.round(color.a) )) );
		return colorTempView.getUint32(0, false);
	}
	//console.log(color, rgbaToUint32(color).toString(16)); //DEBUG

	let Color = rgbaToUint32(colors[0]);
	let nextSplit = 1;

	/*
	*/
}

renderer.uploadYarns = function strips_uploadYarns() {
	//accmulate Position attribs for all yarns:
	let checkpointAttribsList = [];
	let yarnAttribsList = [];

	//start-middle-end attributes for every segment:
	// (as groups of 12 floats, I guess?)
	let RSMEAttribs = [];

	console.log("Building yarns...");
	yarns.yarns.forEach(function(yarn){
		//build checkpoints attribs + record split locations:
		const splits = [];
		const lengths = [];
		for (let i = 0; i < yarn.checkpoints.length; ++i) {
			//grab all similar checkpoints:
			let length = 0.0;
			while (true) {
				length += yarn.checkpoints[i].length * yarns.units[yarn.checkpoints[i].unit].length;
				if (!(i + 1 < yarn.checkpoints.length && yarn.checkpoints[i+1].point === yarn.checkpoints[i].point)) break;
				++i;
			}
			const cp = yarn.checkpoints[i];
			splits.push(cp.point);
			lengths.push(length);
		}


		//make sure there is a split point at the start and at the end:
		if (splits.length === 0 || splits[0] !== 0) {
			splits.unshift(0);
			lengths.unshift(NaN); //used as a "don't care" marker
		}
		if (splits[splits.length-1] !== yarn.points.length/3-1) {
			splits.push(yarn.points.length/3-1);
			lengths.push(NaN); //used as a "don't care" marker
		}

		//TODO: anything with splits!


		if (yarn.points.length >= 6) {
			let prev = { x:yarn.points[0], y:yarn.points[1], z:yarn.points[2] };
			let next = { x:yarn.points[3], y:yarn.points[4], z:yarn.points[5] };

			RSMEAttribs.push( yarn.radius,
				prev.x, prev.y, prev.z,
				0.25 * (next.x - prev.x) + prev.x, 0.25 * (next.y - prev.y) + prev.y, 0.25 * (next.z - prev.z) + prev.z,
				0.5 * (prev.x + next.x), 0.5 * (prev.y + next.y), 0.5 * (prev.z + next.z)
			);
		}
		for (let i = 8; i < yarn.points.length; i += 3) {
			let prev = { x:yarn.points[i-8], y:yarn.points[i-7], z:yarn.points[i-6] };
			let cur  = { x:yarn.points[i-5], y:yarn.points[i-4], z:yarn.points[i-3] };
			let next = { x:yarn.points[i-2], y:yarn.points[i-1], z:yarn.points[i] };

			RSMEAttribs.push( yarn.radius,
				0.5 * (prev.x + cur.x), 0.5 * (prev.y + cur.y), 0.5 * (prev.z + cur.z),
				cur.x, cur.y, cur.z,
				0.5 * (next.x + cur.x), 0.5 * (next.y + cur.y), 0.5 * (next.z + cur.z)
			);
		}
		if (yarn.points.length > 6) {
			//last segment
			let prev = { x:yarn.points[yarn.points.length-6], y:yarn.points[yarn.points.length-5], z:yarn.points[yarn.points.length-4] };
			let next = { x:yarn.points[yarn.points.length-3], y:yarn.points[yarn.points.length-2], z:yarn.points[yarn.points.length-1] };

			RSMEAttribs.push( yarn.radius,
				0.5 * (prev.x + next.x), 0.5 * (prev.y + next.y), 0.5 * (prev.z + next.z),
				0.75 * (next.x - prev.x) + prev.x, 0.75 * (next.y - prev.y) + prev.y, 0.75 * (next.z - prev.z) + prev.z,
				next.x, next.y, next.z
			);
		}

		/*
		//compute array of segment colors based on stretch:
		const spineColors = [];
		for (let i = 0; i < spineSplits.length; ++i) {
			let length = lengths[i];
			let spineLength = spineLengths[i];
			if (!(length === length)) {
				spineColors.push(NaN); //marked as "Don't Care"
			} else if (spineLength <= 0.0 || length <= 0.0) {
				spineColors.push(NaN); //invalid -- treat as "Don't Care"
			} else {
				spineColors.push(Math.log2(length / spineLength));
			}
		}
		let minVal = -1.0;
		let maxVal = 1.0;

		//make sure that 0 always ends up centered:
		maxVal = Math.max(maxVal, -minVal);
		minVal = Math.min(minVal, -maxVal);
		console.assert(Math.abs(maxVal) === Math.abs(minVal), "zero is centered");

		//This is BrBG from matplotlib:
		// https://github.com/matplotlib/matplotlib/blob/a78675348aeeb7756014c6d628cd603923f72560/lib/matplotlib/_cm.py
		const scale = [
			[0.32941176470588235,  0.18823529411764706,  0.0196078431372549],
			[0.5490196078431373 ,  0.31764705882352939,  0.0392156862745098],
			[0.74901960784313726,  0.50588235294117645,  0.17647058823529413],
			[0.87450980392156863,  0.76078431372549016,  0.49019607843137253],
			[0.96470588235294119,  0.90980392156862744,  0.76470588235294112],
			[0.96078431372549022,  0.96078431372549022,  0.96078431372549022],
			[0.7803921568627451 ,  0.91764705882352937,  0.89803921568627454],
			[0.50196078431372548,  0.80392156862745101,  0.75686274509803919],
			[0.20784313725490197,  0.59215686274509804,  0.5607843137254902 ],
			[0.00392156862745098,  0.4                ,  0.36862745098039218],
			[0.0                ,  0.23529411764705882,  0.18823529411764706]
		];
		function colorMap(amt) {
			let f = amt * (scale.length-1);
			f = Math.max(0.0, f);
			f = Math.min(scale.length-1 - 1e-6, f);
			let i = Math.floor(f);
			f -= i;
			function cvt(v) {
				return Math.max(0, Math.min(255, Math.round(v * 255)));
			}
			return {
				r:cvt(f * (scale[i+1][0] - scale[i][0]) + scale[i][0]),
				g:cvt(f * (scale[i+1][1] - scale[i][1]) + scale[i][1]),
				b:cvt(f * (scale[i+1][2] - scale[i][2]) + scale[i][2]),
				a: 255
			};
		}

		for (let i = 0; i < spineColors.length; ++i) {
			let v = spineColors[i];
			if (!(v === v)) {
				//"don't care" value
				spineColors[i] = yarn.color;
			} else {
				spineColors[i] = colorMap((v - minVal) / (maxVal - minVal));
			}
		}
		*/
	});

	gl.bindBuffer(gl.ARRAY_BUFFER, perInstanceBuffer);
	gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(RSMEAttribs), gl.STATIC_DRAW);
	instanceCount = RSMEAttribs.length / 10;
	console.log("Have " + instanceCount + " instances.");

	let ParamSideAttribs = [];
	for (let i = 0; i <= 10; ++i) {
		let amt = i / 10.0;
		ParamSideAttribs.push(amt,-1.0);
		ParamSideAttribs.push(amt, 1.0);
	}
	gl.bindBuffer(gl.ARRAY_BUFFER, perVertexBuffer);
	gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(ParamSideAttribs), gl.STATIC_DRAW);
	vertexCount = ParamSideAttribs.length / 2;
	console.log("Have " + vertexCount + " vertices per instance.");

	/*
	//accumulate all the attribs list into one big list:
	let totalAttribs = 0;
	function countBuffer(buffer){
		if (buffer.byteLength == 0) return;
		if (totalAttribs != 0) totalAttribs += 2;
		totalAttribs += buffer.byteLength / BYTES_PER_ATTRIB;
	}
	yarnAttribsList.forEach(countBuffer);
	checkpointAttribsList.forEach(countBuffer);

	let Attribs = new ArrayBuffer(totalAttribs * BYTES_PER_ATTRIB);
	let AttribsArray = new Uint32Array(Attribs);

	let offset = 0;
	function copyBuffer(buffer){
		if (buffer.byteLength == 0) return;
		let array = new Uint32Array(buffer);
		if (offset != 0) {
			AttribsArray.set(AttribsArray.slice(offset - BYTES_PER_ATTRIB/4, offset), offset);
			offset += BYTES_PER_ATTRIB/4;
			AttribsArray.set(array.slice(0, BYTES_PER_ATTRIB/4), offset);
			offset += BYTES_PER_ATTRIB/4;
		}
		AttribsArray.set(array, offset);
		offset += array.length;
	}
	yarnAttribsList.forEach(copyBuffer);
	checkpointAttribsList.forEach(copyBuffer);
	console.assert(offset == Attribs.byteLength/4, "set full array");


	gl.bindBuffer(gl.ARRAY_BUFFER, attribsBuffer);
	gl.bufferData(gl.ARRAY_BUFFER, Attribs, gl.STATIC_DRAW);
	attribsCount = Attribs.byteLength / BYTES_PER_ATTRIB;
	console.log("Have " + Attribs.byteLength + " bytes of attribs == " + attribsCount + " attribs.");
	*/
};


let ext = gl.getExtension("ANGLE_instanced_arrays");

renderer.redraw = function strips_redraw() {

	gl.clearColor(0.1, 0.1, 0.1, 1.0);
	gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

	gl.useProgram(stripProgram.program);

	gl.bindBuffer(gl.ARRAY_BUFFER, perVertexBuffer);

	gl.vertexAttribPointer(stripProgram.attribLocations.ParamSide,
		2, //size
		gl.FLOAT, //type
		false, //normalize
		8, //stride
		0 //offset
	);
	gl.enableVertexAttribArray(stripProgram.attribLocations.ParamSide);


	gl.bindBuffer(gl.ARRAY_BUFFER, perInstanceBuffer);

	gl.vertexAttribPointer(stripProgram.attribLocations.Radius,
		1, //size
		gl.FLOAT, //type
		false, //normalize
		4*10, //stride
		4*0 //offset
	);
	ext.vertexAttribDivisorANGLE(stripProgram.attribLocations.Radius, 1);
	gl.enableVertexAttribArray(stripProgram.attribLocations.Radius);

	gl.vertexAttribPointer(stripProgram.attribLocations.Start,
		3, //size
		gl.FLOAT, //type
		false, //normalize
		4*10, //stride
		4*1 //offset
	);
	ext.vertexAttribDivisorANGLE(stripProgram.attribLocations.Start, 1);
	gl.enableVertexAttribArray(stripProgram.attribLocations.Start);

	gl.vertexAttribPointer(stripProgram.attribLocations.Middle,
		3, //size
		gl.FLOAT, //type
		false, //normalize
		4*10, //stride
		4*4 //offset
	);
	ext.vertexAttribDivisorANGLE(stripProgram.attribLocations.Middle, 1);
	gl.enableVertexAttribArray(stripProgram.attribLocations.Middle);

	gl.vertexAttribPointer(stripProgram.attribLocations.End,
		3, //size
		gl.FLOAT, //type
		false, //normalize
		4*10, //stride
		4*7 //offset
	);
	ext.vertexAttribDivisorANGLE(stripProgram.attribLocations.End, 1);
	gl.enableVertexAttribArray(stripProgram.attribLocations.End);

	gl.uniformMatrix4fv(
		stripProgram.uniformLocations.MVP_mat4,
		false,
		computeMVP()
	);

	
	const eye = camera.computeAt();
	gl.uniform3fv(
		stripProgram.uniformLocations.EYE_vec3,
		new Float32Array([eye.x, eye.y, eye.z])
	);

	const out = camera.computeOut();
	gl.uniform3fv(
		stripProgram.uniformLocations.OUT_vec3,
		new Float32Array([out.x, out.y, eye.z])
	);



	gl.enable(gl.DEPTH_TEST);
	gl.depthFunc(gl.LESS);

	gl.enable(gl.CULL_FACE);
	gl.cullFace(gl.BACK);

	ext.drawArraysInstancedANGLE(gl.TRIANGLE_STRIP, 0, vertexCount, instanceCount);


	//it messes with other renderers to leave these active:
	ext.vertexAttribDivisorANGLE(stripProgram.attribLocations.Radius, 0);
	ext.vertexAttribDivisorANGLE(stripProgram.attribLocations.Start, 0);
	ext.vertexAttribDivisorANGLE(stripProgram.attribLocations.Middle, 0);
	ext.vertexAttribDivisorANGLE(stripProgram.attribLocations.End, 0);

	gl.disable(gl.CULL_FACE);
	gl.disable(gl.DEPTH_TEST);
};


})();
