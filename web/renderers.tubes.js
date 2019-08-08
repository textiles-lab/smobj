"use strict";
(function(){

window.renderers = window.renderers || {};

let renderer = window.renderers.tubes = {};

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

function quadraticBSpline(points, minDistance, maxAngle) {
	//need at least 3 points (6 coordinates) to make a curve
	if (points.length <= 6) {
		return points;
	}

	//quadratic b-spline through the points (same as bezier through midpoints)

	let poly = [];

	let threshold = Math.cos(maxAngle * Math.PI / 180.0);
	function curveTo(m, e) {
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
			poly.push(e.x, e.y, e.z);
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
			curveTo(sm, sme);
			curveTo(me, e);
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

	//straight segment to first midpoint:
	poly.push(points[0], points[1], points[2]);
	poly.push(0.5 * (points[0] + points[3]), 0.5 * (points[1] + points[4]), 0.5 * (points[2] + points[5]));
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
			z:0.5*(points[i+2]+points[i+5]),
		};
		curveTo(m, e);
	}
	//straight segment from last midpoint to end:
	poly.push(
		points[points.length-3],
		points[points.length-2],
		points[points.length-1]
	);

	return poly;
}

//build a tubular tristrip around a (one hopes) smooth spine:
function buildTube(spine, radius) {
	let tristrip = [];

	if (spine.length < 3) return tristrip;

	//build a cross-section shape:
	const Ring = new Array(5);
	for (let i = 0; i < Ring.length; ++i) {
		const angle = (i / Ring.length) * 2.0 * Math.PI;
		Ring[i] = {x:Math.cos(angle), y:Math.sin(angle)};
	}

	//sweep the cross-section shape along the yarn:
	//..will do this by maintaining a local coordinate system p1,p2 at corners:
	let p1 = { x:NaN, y:NaN, z:NaN};
	let p2 = { x:NaN, y:NaN, z:NaN};
	let at = { x:spine[0], y:spine[1], z:spine[2] };

	function sweepTo(nextP1, nextP2, nextAt) {
		//for now, build a ring at each new frame:
		//(eventually, build a connector between frames)

		p1 = nextP1;
		p2 = nextP2;
		at = nextAt;

		for (let i = 0; i < Ring.length; ++i) {
			tristrip.push(
				radius * (Ring[i].x * p1.x + Ring[i].y * p2.x) + at.x,
				radius * (Ring[i].x * p1.y + Ring[i].y * p2.y) + at.y,
				radius * (Ring[i].x * p1.z + Ring[i].y * p2.z) + at.z
			);
		}
		tristrip.push(tristrip[tristrip.length - Ring.length * 3]);
		tristrip.push(tristrip[tristrip.length - Ring.length * 3]);
		tristrip.push(tristrip[tristrip.length - Ring.length * 3]);
	}


	//compute sweep frames by dragging a local coordinate frame along the edges:
	let prevD1 = { x:NaN, y:NaN, z:NaN };
	let prevD2 = { x:NaN, y:NaN, z:NaN };
	let prevAlong = { x:NaN, y:NaN, z:NaN };

	for (let i = 3; i + 2 < spine.length; i += 3) {
		let next = { x:spine[i+0], y:spine[i+1], z:spine[i+2] };
		//skip duplicated points:
		if (next.x === at.x && next.y === at.y && next.z === at.z) continue;

		let along = {
			x:next.x - at.x,
			y:next.y - at.y,
			z:next.z - at.z
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

		//TODO: compute average?
		sweepTo(d1, d2, at);

		prevD1 = d1;
		prevD2 = d2;
		prevAlong = along;
		at = next;
	}

	//finish sweep:
	if (prevD1.x === prevD1.x) {
		sweepTo(prevD1, prevD2, at);
	}

	//TODO: last frame also

	return tristrip;
	
}

renderer.uploadYarns = function tubes_uploadYarns() {
	//clear yarnStarts array:
	yarnStarts.splice(0,yarnStarts.length);

	//accmulate Position attribs for all yarns:
	let Positions = [];
	yarns.yarns.forEach(function(yarn){
		yarnStarts.push(Positions.length / 3);

		const spine = quadraticBSpline(yarn.points, 0.5 * yarn.radius, 25.0);
		const tube = buildTube(spine, yarn.radius);

		Positions = Positions.concat(tube);
	});
	yarnStarts.push(Positions.length / 3);

	gl.bindBuffer(gl.ARRAY_BUFFER, yarnsBuffer);
	gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(Positions), gl.STATIC_DRAW);
};

renderer.redraw = function tubes_redraw() {

	gl.clearColor(0.1, 0.1, 0.1, 1.0);
	gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

	gl.useProgram(colorProgram.program);
	gl.bindBuffer(gl.ARRAY_BUFFER, yarnsBuffer);

	gl.vertexAttrib4f(colorProgram.attribLocations.Color, 0.0, 1.0, 1.0, 1.0);
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
