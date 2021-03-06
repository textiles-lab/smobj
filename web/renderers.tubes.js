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
function buildTube(spine, radius, colorOrSplits_, colors_) {
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

	let Color = rgbaToUint32(colors[0]);
	let nextSplit = 1;

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
				Color //DEBUG version: (reverse ? ColorD : ColorC)
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
					Color
				);
				attrib(
					radius * (Ring[j].x * nextP1.x + Ring[j].y * nextP2.x) + nextPt.x,
					radius * (Ring[j].x * nextP1.y + Ring[j].y * nextP2.y) + nextPt.y,
					radius * (Ring[j].x * nextP1.z + Ring[j].y * nextP2.z) + nextPt.z,
					Ring[j].x * nextP1.x + Ring[j].y * nextP2.x,
					Ring[j].x * nextP1.y + Ring[j].y * nextP2.y,
					Ring[j].x * nextP1.z + Ring[j].y * nextP2.z,
					Color
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

		if (nextSplit < splits.length) { //update color:
			const idx = i / 3 - 1;
			while (nextSplit < splits.length && idx >= splits[nextSplit]) { //doing this because may have zero-length segments that got skipped
				Color = rgbaToUint32(colors[nextSplit]);
				++nextSplit;
			}
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

renderer.uploadYarns = function tubes_uploadYarns() {
	//accmulate Position attribs for all yarns:
	let checkpointAttribsList = [];
	let yarnAttribsList = [];

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

		//build spine (making sure to split at 'splits' and record lengths):
		const {points:spine, splits:spineSplits, lengths:spineLengths} = quadraticBSpline(yarn.points, 0.5 * yarn.radius, 25.0, splits);

		//console.log(splits, spineSplits, lengths, spineLengths);
		console.assert(spineSplits.length === splits.length, "all splits remapped");

		console.log((splits.length - 1) + " checkpoint segments.");

		//draw little "+"'s at each checkpoint:
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
			checkpointAttribsList.push( buildTube(spineX, 0.3 * yarn.radius, color) );
			checkpointAttribsList.push( buildTube(spineY, 0.3 * yarn.radius, color) );
			checkpointAttribsList.push( buildTube(spineZ, 0.3 * yarn.radius, color) );
		});

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
		/*
		//don't scale, leave at 0.5-2.0 range.
		spineColors.forEach(function(v){
			if (v === v) {
				minVal = Math.min(minVal,v);
				maxVal = Math.max(maxVal,v);
			}
		});
		*/

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

		yarnAttribsList.push( buildTube(spine, yarn.radius, spineSplits, spineColors) );
	});

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
