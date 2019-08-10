"use strict";
(function(){

window.renderers = window.renderers || {};

let renderer = window.renderers.tubes = {};

const colorProgram = initShaderProgram(
`
	uniform mat4 mvp_mat4;

	attribute vec4 Position;
	attribute vec3 Normal;
	attribute vec4 Color;

	varying vec4 color;
	varying vec3 normal;

	void main() {
		gl_Position = mvp_mat4 * Position;
		color = Color;
		normal = Normal;
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

const attribsBuffer = gl.createBuffer();
let attribsCount = 0;

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
		/*
		let inter = 14;
		for (let i = 0; i < inter; ++i) {
			const t = (i + 1.0) / (inter + 1.0);
			const sm = {
				x: t * (m.x - s.x) + s.x,
				y: t * (m.y - s.y) + s.y,
				z: t * (m.z - s.z) + s.z
			};
			const me = {
				x: t * (e.x - m.x) + m.x,
				y: t * (e.y - m.y) + m.y,
				z: t * (e.z - m.z) + m.z
			};
			const sme = {
				x: t * (me.x - sm.x) + sm.x,
				y: t * (me.y - sm.y) + sm.y,
				z: t * (me.z - sm.z) + sm.z
			};
			poly.push(sme.x, sme.y, sme.z);
		}
		*/
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

//build a tubular tristrip around a (one hopes) smooth spine:
//return a ArrayBuffer with attribs inside
function buildTube(spine, radius, color) {
	if (spine.length < 3) return new ArrayBuffer(0);

	//build a cross-section shape:
	const Ring = new Array(5);
	for (let i = 0; i < Ring.length; ++i) {
		const angle = (i / Ring.length) * 2.0 * Math.PI;
		Ring[i] = {x:Math.cos(angle), y:Math.sin(angle)};
	}

	let attribBytes = BYTES_PER_ATTRIB * (
		2 * ( Ring.length + (Ring.length % 2 ? 1 : 0) ) //caps
		+ (spine.length/3 - 1) * 2 * (Ring.length+1) //sweep segments
	);

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

	const ColorA = rgbaToUint32(color);
	const ColorB = rgbaToUint32(color);
	const ColorC = rgbaToUint32(color);
	const ColorD = rgbaToUint32(color);

	/*
	//DEBUG
	attribBytes = (3*4+4*1) * 4;
	tristrip = new ArrayBuffer(attribBytes);
	data = new DataView(tristrip);
	attrib(0,0,0, ColorA);
	attrib(1,0,0, ColorB);
	console.log(new Float32Array(tristrip));
	attrib(0,1,0, ColorC);
	attrib(1,1,0, ColorD);
	console.assert(dataOffset === attribBytes);
	console.log(new Float32Array(tristrip));
	return tristrip;
	*/

	function cap(p1, p2, pt, reverse) {
		//play with index orders so that final cap either
		// (a) if !reverse:
		//   - starts with 0, is oriented along 'along'
		// (b) if reverse
		//   - ends with 0, is oriented along '-along'
		// either way, has *even* attrib count

		let inds = [];
		for (let i = 0; i < Ring.length; ++i) {
			inds.push(i);
		}
		if (inds.length % 2) {
			inds.push(0);
		}
		if (!reverse) {
			inds.push(inds.shift());
		}
		let ord = [];
		for (let i = inds.length/2-1; i >= 0; --i) {
			if (reverse) {
				ord.push(inds[inds.length-1-i]);
				ord.push(inds[i]);
			} else {
				ord.push(inds[i]);
				ord.push(inds[inds.length-1-i]);
			}
		}
		inds = ord;

		if (!reverse) {
			inds.reverse();
		}
		//note: want to always emit even number of verts, start with index 0.
		//console.log(reverse, inds);

		const n = normalize(cross(p1, p2));
		if (reverse) {
			n.x = -n.x;
			n.y = -n.y;
			n.z = -n.z;
		}

		inds.forEach(function(i){
			attrib(
				radius * (Ring[i].x * p1.x + Ring[i].y * p2.x) + pt.x,
				radius * (Ring[i].x * p1.y + Ring[i].y * p2.y) + pt.y,
				radius * (Ring[i].x * p1.z + Ring[i].y * p2.z) + pt.z,
				n.x,
				n.y,
				n.z,
				(reverse ? ColorD : ColorC)
			);
		});
	}

	function sweepTo(nextP1, nextP2, nextPt) {
		//for now, build a ring at each new frame:
		//(eventually, build a connector between frames)
		if (p1.x === p1.x) {
			//build connector:
			for (let i = 0; i < Ring.length; ++i) {
				let j = (Ring.length - i) % Ring.length; //travel cw but start at 0
				attrib(
					radius * (Ring[j].x * p1.x + Ring[j].y * p2.x) + pt.x,
					radius * (Ring[j].x * p1.y + Ring[j].y * p2.y) + pt.y,
					radius * (Ring[j].x * p1.z + Ring[j].y * p2.z) + pt.z,
					Ring[j].x * p1.x + Ring[j].y * p2.x,
					Ring[j].x * p1.y + Ring[j].y * p2.y,
					Ring[j].x * p1.z + Ring[j].y * p2.z,
					ColorA
				);
				attrib(
					radius * (Ring[j].x * nextP1.x + Ring[j].y * nextP2.x) + nextPt.x,
					radius * (Ring[j].x * nextP1.y + Ring[j].y * nextP2.y) + nextPt.y,
					radius * (Ring[j].x * nextP1.z + Ring[j].y * nextP2.z) + nextPt.z,
					Ring[j].x * nextP1.x + Ring[j].y * nextP2.x,
					Ring[j].x * nextP1.y + Ring[j].y * nextP2.y,
					Ring[j].x * nextP1.z + Ring[j].y * nextP2.z,
					ColorB
				);
			}
			prevAttrib(-2 * Ring.length);
			prevAttrib(-2 * Ring.length);
		}

		p1 = {x:nextP1.x, y:nextP1.y, z:nextP1.z};
		p2 = {x:nextP2.x, y:nextP2.y, z:nextP2.z};
		pt = {x:nextPt.x, y:nextPt.y, z:nextPt.z};
	}


	//compute sweep frames by dragging a local coordinate frame along the edges:
	let prevD1 = { x:NaN, y:NaN, z:NaN };
	let prevD2 = { x:NaN, y:NaN, z:NaN };
	let prevAlong = { x:NaN, y:NaN, z:NaN };
	let prevAt = { x:spine[0], y:spine[1], z:spine[2] };

	for (let i = 3; i + 2 < spine.length; i += 3) {
		let at = { x:spine[i+0], y:spine[i+1], z:spine[i+2] };
		//skip duplicated points:
		if (at.x === prevAt.x && at.y === prevAt.y && at.z === prevAt.z) continue;

		let along = {
			x:at.x - prevAt.x,
			y:at.y - prevAt.y,
			z:at.z - prevAt.z
		};
		along = normalize(along);

		let d1;
		if (!(prevD1.x === prevD1.x)) {
			//first segment -- build initial coordinate system:
			if (Math.abs(along.x) <= Math.abs(along.y) && Math.abs(along.x) <= Math.abs(along.z)) {
				d1 = {x:1.0, y:0.0, z:0.0};
			} else if (Math.abs(along.y) <= Math.abs(along.z)) {
				d1 = {x:0.0, y:1.0, z:0.0};
			} else {
				d1 = {x:0.0, y:0.0, z:1.0};
			}
			{ //make d1 orthogonal to along:
				let d = dot(d1, along);
				d1.x -= d * along.x;
				d1.y -= d * along.y;
				d1.z -= d * along.z;
				d1 = normalize(d1);
			}
		} else {
			//not first segment.
			
			//compute next d1 by rotating prevAlong -> along:

			if (Math.abs(dot(prevAlong, along)) < 0.999) {
				//rotate prevD1 to get d1:
				let perp = normalize(cross(prevAlong, along));
				let prevOut = cross(prevAlong, perp);
				let out = cross(along, perp);
				let d1Out = dot(prevD1, prevOut);
				d1 = {
					x:prevD1.x + d1Out * (out.x - prevOut.x),
					y:prevD1.y + d1Out * (out.y - prevOut.y),
					z:prevD1.z + d1Out * (out.z - prevOut.z)
				};
				d1 = normalize(d1);
			} else {
				//prevAlong is pretty close to along, so just transfer and re-normalize:
				d1 = { x:prevD1.x, y:prevD1.y, z:prevD1.z };
				let d = dot(d1, along);
				d1.x -= d * along.x;
				d1.y -= d * along.y;
				d1.z -= d * along.z;
				d1 = normalize(d1);
			}
		}
		let d2 = cross(along, d1);

		if (!(prevD1.x === prevD1.x)) {
			//start sweep:
			cap(d1, d2, prevAt, true);
			sweepTo(d1, d2, prevAt);
		} else {
			//sweep to corner:
			const avgD1 = {
				x:0.5 * (prevD1.x + d1.x),
				y:0.5 * (prevD1.y + d1.y),
				z:0.5 * (prevD1.z + d1.z)
			};
			const avgD2 = {
				x:0.5 * (prevD2.x + d2.x),
				y:0.5 * (prevD2.y + d2.y),
				z:0.5 * (prevD2.z + d2.z)
			};
			sweepTo(avgD1, avgD2, prevAt);
		}

		prevD1 = d1;
		prevD2 = d2;
		prevAlong = along;
		prevAt = at;
	}

	//finish sweep:
	if (prevD1.x === prevD1.x) {
		sweepTo(prevD1, prevD2, prevAt);
		cap(prevD1, prevD2, prevAt);
	}

	console.assert(dataOffset === attribBytes);

	return tristrip;
	
}

//Helper: concatenate ArrayBuffers of attribs:
function concat(a,b) {
	console.assert(a.byteLength % BYTES_PER_ATTRIB === 0, "Expecting buffers of attribs");
	console.assert(b.byteLength % BYTES_PER_ATTRIB === 0, "Expecting buffers of attribs");
	console.assert(BYTES_PER_ATTRIB % 4 === 0, "Expecting attribs to be multiples of 32 bits");

	if (b.byteLength === 0) return a;
	if (a.byteLength === 0) return b;

	const o = new ArrayBuffer(a.byteLength + b.byteLength + 2*BYTES_PER_ATTRIB);

	const av = new Uint32Array(a);
	const bv = new Uint32Array(b);
	const ov = new Uint32Array(o);

	ov.set(av, 0);
	//duplicate last vertex:
	ov.set(av.slice(av.length-BYTES_PER_ATTRIB/4, av.length), av.length);
	//duplicate first vertex:
	ov.set(bv.slice(0, BYTES_PER_ATTRIB/4), av.length + BYTES_PER_ATTRIB/4);
	ov.set(bv, av.length + 2*BYTES_PER_ATTRIB/4);

	return ov;
}

renderer.uploadYarns = function tubes_uploadYarns() {
	//accmulate Position attribs for all yarns:
	let Attribs = new ArrayBuffer(0);

	yarns.yarns.forEach(function(yarn){
		//build checkpoints attribs + record split locations:
		let checkpoints = new ArrayBuffer(0);
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

		//build spine (making sure to split at 'splits' and record lengths):
		const {points:spine, splits:spineSplits, lengths:spineLengths} = quadraticBSpline(yarn.points, 0.5 * yarn.radius, 25.0, splits);

		//console.log(splits, spineSplits, lengths, spineLengths);
		console.assert(spineSplits.length === splits.length, "all splits remapped");

		spineSplits.forEach(function(si, i) {
			const pt = {
				x:spine[3*si+0],
				y:spine[3*si+1],
				z:spine[3*si+2]
			};
			const spineX = [
				pt.x - yarn.radius * 2.0, pt.y, pt.z,
				pt.x + yarn.radius * 2.0, pt.y, pt.z
			];
			const spineY = [
				pt.x, pt.y - yarn.radius * 2.0, pt.z,
				pt.x, pt.y + yarn.radius * 2.0, pt.z
			];
			const spineZ = [
				pt.x, pt.y, pt.z - yarn.radius * 2.0,
				pt.x, pt.y, pt.z + yarn.radius * 2.0
			];
			let color = {r:255, g:255, b:255, a:255};
			if (!(lengths[i] === lengths[i])) {
				//this is a "don't care" point added for vis
				//DEBUG: color = {r:255, g:0, b:0, a:255};
				return;
			}
			const crossX = buildTube(spineX, 0.3 * yarn.radius, color);
			const crossY = buildTube(spineY, 0.3 * yarn.radius, color);
			const crossZ = buildTube(spineZ, 0.3 * yarn.radius, color);
			checkpoints = concat(checkpoints,
				concat(concat(crossX, crossY), crossZ)
			);
		});

		const tube = buildTube(spine, yarn.radius, yarn.color);
	
		//tube is an ArrayBuffer:
		Attribs = concat(Attribs, concat(tube, checkpoints));
	});

	gl.bindBuffer(gl.ARRAY_BUFFER, attribsBuffer);
	gl.bufferData(gl.ARRAY_BUFFER, Attribs, gl.STATIC_DRAW);
	attribsCount = Attribs.byteLength / BYTES_PER_ATTRIB;
	console.log("Have " + Attribs.byteLength + " bytes of attribs == " + attribsCount + " attribs.");
};

renderer.redraw = function tubes_redraw() {

	gl.clearColor(0.1, 0.1, 0.1, 1.0);
	gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

	gl.useProgram(colorProgram.program);
	gl.bindBuffer(gl.ARRAY_BUFFER, attribsBuffer);

	gl.vertexAttribPointer(colorProgram.attribLocations.Color,
		4, //size
		gl.UNSIGNED_BYTE, //type
		true, //normalize
		BYTES_PER_ATTRIB, //stride
		24 //offset
	);
	gl.enableVertexAttribArray(colorProgram.attribLocations.Color);

	gl.vertexAttribPointer(colorProgram.attribLocations.Normal,
		3, //size
		gl.FLOAT, //type
		false, //normalize
		BYTES_PER_ATTRIB, //stride
		12 //offset
	);
	gl.enableVertexAttribArray(colorProgram.attribLocations.Normal);

	gl.vertexAttribPointer(colorProgram.attribLocations.Position,
		3, //size
		gl.FLOAT, //type
		false, //normalize
		BYTES_PER_ATTRIB, //stride
		0 //offset
	);
	gl.enableVertexAttribArray(colorProgram.attribLocations.Position);

	gl.uniformMatrix4fv(
		colorProgram.uniformLocations.mvp_mat4,
		false,
		computeMVP()
	);

	gl.enable(gl.DEPTH_TEST);
	gl.depthFunc(gl.LESS);

	gl.enable(gl.CULL_FACE);
	gl.cullFace(gl.BACK);

	gl.drawArrays(gl.TRIANGLE_STRIP, 0, attribsCount);

	gl.disable(gl.CULL_FACE);
	gl.disable(gl.DEPTH_TEST);
};


})();
