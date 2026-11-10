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
 *  - store coordinate system using angle or cos+sin factors? (memory vs compute trade-off)
 */

window.renderers = window.renderers || {};

let renderer = window.renderers.strips = {};

const stripProgram = initShaderProgram(
`
	uniform mat4 MVP_mat4;
	uniform vec3 EYE_vec3;
	uniform vec3 OUT_vec3;

	attribute vec2 PerpNPerpB; //radius + 1st perpendicular direction in terms of normal / binormal
	attribute vec4 Start;  //x,y,z, length
	attribute vec4 Middle; //x,y,z, length
	attribute vec4 End;    //x,y,z, length
	attribute vec3 Normal; //perpendicular to the plane containing start/middle/end

	attribute vec4 PlyMulAddRadSquish; //parameters for current ply -- angle is computed as pi * fract(Mul * len + Add)
	//paramters for current strand:
	attribute vec2 SAMulAdd; //angle is computed as pi * fract(Mul * len + Add)
	attribute vec4 SRMulAddRad1Rad2; //radius is mix(Rad1, Rad2, cos(pi * fract(Mul * len + Add)))

	attribute vec2 ParamSide;

	varying vec4 color;
	varying float pos; //position from -normal to normal across strand
	varying vec3 normal;
	varying vec3 tangent;
	varying vec2 local; //local position in yarn

	#define PI 3.1415926
	#define EPS 0.001 //epsilon for finite difference approximations

	void main() {
		vec4 sm = mix(Start, Middle, ParamSide.x);
		vec4 me = mix(Middle, End, ParamSide.x);
		vec4 sme = mix(sm, me, ParamSide.x);

		vec3 at = sme.xyz;
		float len = sme.w;

		//local frame:
		vec3 Tangent = normalize(mix(Middle.xyz - Start.xyz, End.xyz - Middle.xyz, ParamSide.x)); //notice that there is a factor of two inside the normalize but (shrug)
		vec3 Binormal = cross(Tangent, Normal);
		vec3 p1 = PerpNPerpB.x * Normal + PerpNPerpB.y * Binormal;
		vec3 p2 = -PerpNPerpB.y * Normal + PerpNPerpB.x * Binormal;

		//ply center in local frame:
		float ply_twist = PI * 2.0 * fract(PlyMulAddRadSquish.x * len + PlyMulAddRadSquish.y);
		vec2 ply_dir = vec2(cos(ply_twist), sin(ply_twist));
		vec2 ply_at = PlyMulAddRadSquish.z * ply_dir;

		//strand center in ply frame:
		float strand_twist = PI * 2.0 * fract(SAMulAdd.x * len + SAMulAdd.y);
		float strand_rad = mix(SRMulAddRad1Rad2.z, SRMulAddRad1Rad2.w, 0.5 + 0.5 * cos(PI * 2.0 * fract(SRMulAddRad1Rad2.x * len + SRMulAddRad1Rad2.y)));
		vec2 strand_at = ply_at + strand_rad * (vec2(-ply_dir.y, ply_dir.x) * cos(strand_twist) + ply_dir * PlyMulAddRadSquish.w * sin(strand_twist));

		vec3 pt = strand_at.x * p1 + strand_at.y * p2 + at;


		//ply tangent (lots of chain rule!):
		vec4 d_sme = 2.0 * mix(Middle - Start, End - Middle, ParamSide.x);
		vec3 d_at = d_sme.xyz;
		float d_len = d_sme.w;

		vec3 t = mix(Middle.xyz - Start.xyz, End.xyz - Middle.xyz, ParamSide.x);
		vec3 d_t = Start.xyz - 2.0 * Middle.xyz + End.xyz;
		vec3 d_Tangent =
			(d_t - normalize(t) * dot(d_t, normalize(t))) / length(t) //<-- seems inefficient, can probably save at least one sqrt here
		;
		//given T(t) = A(t) / |A(t)|
		// dT(t)/dt = (dA(t)/dt - T(t)*dot(dA(t)/dt,T(t)) / |A(t)|
		//via: https://www.gamedev.net/forums/topic/553227-vector-calculus-question-derivative-of-a-normalized-vector-from-coutinhos-book/

		vec3 d_Binormal = cross(d_Tangent, Normal) /* constant: + cross(Tangent, d_Normal) */;
		vec3 d_p1 = PerpNPerpB.y * d_Binormal;
		vec3 d_p2 = PerpNPerpB.x * d_Binormal;

		float d_ply_twist = PI * 2.0 * PlyMulAddRadSquish.x * d_len;
		vec2 d_ply_dir = d_ply_twist * vec2(-sin(ply_twist), cos(ply_twist));
		vec2 d_ply_at = PlyMulAddRadSquish.z * d_ply_dir;

		float d_strand_twist = PI * 2.0 * SAMulAdd.x * d_len;
		float d_strand_rad = (SRMulAddRad1Rad2.w - SRMulAddRad1Rad2.z) * 0.5 * -sin(PI * 2.0 * fract(SRMulAddRad1Rad2.x * len + SRMulAddRad1Rad2.y)) * PI * 2.0 * SRMulAddRad1Rad2.x * d_len;

		vec2 d_strand_at = d_ply_at
			+ d_strand_rad * (vec2(-ply_dir.y, ply_dir.x) * cos(strand_twist) + ply_dir * PlyMulAddRadSquish.w * sin(strand_twist))
			+ strand_rad * (
				vec2(-d_ply_dir.y, d_ply_dir.x) * cos(strand_twist) + d_ply_dir * PlyMulAddRadSquish.w * sin(strand_twist)
				+ vec2(-ply_dir.y, ply_dir.x) * d_strand_twist * -sin(strand_twist) + ply_dir * PlyMulAddRadSquish.w * d_strand_twist * cos(strand_twist)
			);

		vec3 d_pt = d_strand_at.x * p1 + strand_at.x * d_p1 + d_strand_at.y * p2 + strand_at.y * d_p2 + d_at;

#ifdef COMPARE_TO_NUMERICAL
		//very lazy way of approximating the derviative:
		vec4 n_sm = mix(Start, Middle, ParamSide.x + EPS);
		vec4 n_me = mix(Middle, End, ParamSide.x + EPS);
		vec4 n_sme = mix(n_sm, n_me, ParamSide.x + EPS);

		vec3 n_at = n_sme.xyz;
		float n_len = n_sme.w;

		vec3 n_Tangent = normalize(mix(Middle.xyz - Start.xyz, End.xyz - Middle.xyz, ParamSide.x + EPS));
		vec3 n_Binormal = cross(n_Tangent, Normal);
		vec3 n_p1 = PerpNPerpB.x * Normal + PerpNPerpB.y * n_Binormal;
		vec3 n_p2 = -PerpNPerpB.y * Normal + PerpNPerpB.x * n_Binormal;

		//ply center in local frame:
		float n_ply_twist = PI * 2.0 * fract(PlyMulAddRadSquish.x * n_len + PlyMulAddRadSquish.y);
		if (n_ply_twist + PI < ply_twist) n_ply_twist += 2.0 * PI;
		if (n_ply_twist - PI > ply_twist) n_ply_twist -= 2.0 * PI;

		vec2 n_ply_dir = vec2(cos(n_ply_twist), sin(n_ply_twist));
		vec2 n_ply_at = PlyMulAddRadSquish.z * n_ply_dir;

		//strand center in ply frame:
		float n_strand_twist = PI * 2.0 * fract(SAMulAdd.x * n_len + SAMulAdd.y);
		if (n_strand_twist + PI < strand_twist) n_strand_twist += 2.0 * PI;
		if (n_strand_twist - PI > strand_twist) n_strand_twist -= 2.0 * PI;

		float n_strand_rad = mix(SRMulAddRad1Rad2.z, SRMulAddRad1Rad2.w, 0.5 + 0.5 * cos(PI * 2.0 * fract(SRMulAddRad1Rad2.x * n_len + SRMulAddRad1Rad2.y)));
		vec2 n_strand_at = n_ply_at + n_strand_rad * (vec2(-n_ply_dir.y, n_ply_dir.x) * cos(n_strand_twist) + n_ply_dir * PlyMulAddRadSquish.w * sin(n_strand_twist));

		vec3 n_pt = n_strand_at.x * n_p1 + n_strand_at.y * n_p2 + n_at;

		//---- set color to check derviatives ---
		//vec3 delta = vec3( (n_sme.xyz - sme.xyz)/EPS - d_sme.xyz);
		//vec3 delta = vec3( (n_ply_twist - ply_twist)/EPS - d_ply_twist);
		//vec3 delta = vec3( (n_ply_dir.xyy - ply_dir.xyy)/EPS - d_ply_dir.xyy);
		//vec3 delta = vec3( (n_ply_at.xyy - ply_at.xyy)/EPS - d_ply_at.xyy);
		//vec3 delta = vec3( (n_strand_rad - strand_rad)/EPS - d_strand_rad);
		//vec3 delta = vec3( (n_strand_twist - strand_twist)/EPS - d_strand_twist);
		//vec3 delta = vec3( (n_strand_at.xyy - strand_at.xyy)/EPS - d_strand_at.xyy);
		//vec3 delta = vec3( (n_Tangent - Tangent)/EPS - d_Tangent);
		//vec3 delta = vec3( (n_Binormal - Binormal)/EPS - d_Binormal);
		//vec3 delta = vec3( (n_p1 - p1)/EPS - d_p1);
		//vec3 delta = vec3( (n_p2 - p2)/EPS - d_p2);
		//vec3 delta = vec3( (n_at - at)/EPS - d_at);
		//vec3 delta = vec3( (n_pt - pt)/EPS - d_pt);
		//color = vec4(10.0*abs(delta), 1); //Color;
#endif

		//shade with direction:
		//color = vec4(normalize(d_pt)*0.5+0.5, 1.0);


		//build strip:
		vec3 along = d_pt;
		vec3 to_eye = EYE_vec3 - pt;
		pt += normalize(cross(along, to_eye)) * ParamSide.y; //now an approx of the perpendicular
		gl_Position = MVP_mat4 * vec4(pt, 1.0);


		color = vec4(1.0, 1.0, 0.0, 1.0); //DEBUG
		normal = normalize(to_eye - dot(to_eye, d_pt) / dot(d_pt, d_pt) * d_pt);
		pos = sign(ParamSide.y);
		tangent = normalize(d_pt);
		local = strand_at;
	}
`,`
	uniform highp vec3 TO_LIGHT_vec3; //distant directional light
	varying lowp vec4 color;
	varying highp vec3 normal;
	varying highp float pos;
	varying highp vec2 at;
	varying highp vec3 tangent;
	varying highp vec2 local;
	void main() {
		highp vec3 n = normalize(normal); //perpendicular to both yarn and direction to eye
		highp vec3 t = normalize(tangent); //tangent to yarn
		highp vec3 b = cross(n, t); //perpendicular to yarn, points toward eye

		highp vec3 s = pos * n + sqrt(1.0 - pos*pos) * b;
		//upper dome light:
		//highp vec3 light = mix(vec3(0.1), vec3(1.0), dot(TO_LIGHT_vec3,s)+0.5);

		highp float e = mix(0.1, 1.0, 1.0 - abs(dot(TO_LIGHT_vec3, t)) );
		//highp float e = mix(0.1, 1.0, 0.5 + dot(TO_LIGHT_vec3, n) );
		e *= mix(0.5, 1.0, smoothstep(0.0, 0.1, length(local)));

		gl_FragColor = vec4(color.rgb * e, color.a);
		//gl_FragColor = vec4(0.5*s+0.5, 1.0);
	}
`);

const perInstanceBuffer = gl.createBuffer();
let instanceCount = 0;
//const perStrandBuffer = gl.createBuffer();
//let strandCount = 0;
const perVertexBuffer = gl.createBuffer();
let vertexCount = 0;

renderer.uploadYarns = function strips_uploadYarns() {
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


		let totalLen = 0.0;
		let perp = null; //perpendicular vector -- will parallel-transport down the segments

		function segment(s, m, e) {
			if (s.x === m.x && s.y === m.y && s.z === m.z && s.x === e.x && s.y === e.y && s.z === e.z) {
				console.log("degenerate curve segment, skipping.");
				return;
			}
			if ((s.x === m.x && s.y === m.y && s.z === m.z)
			 || (m.x === e.x && m.y === e.y && m.z === e.z)) {
				console.log("s == m or m == e curve segment; reparameterizing");
				segment(s, {x:0.5 * (s.x + e.x), y:0.5 * (s.y + e.y), z:0.5 * (s.z + e.z)}, e);
				return;
			}

			let startTangent = normalize({ x:m.x-s.x, y:m.y-s.y, z:m.z-s.z });
			let endTangent = normalize({ x:e.x-m.x, y:e.y-m.y, z:e.z-m.z });
			let normal = cross(startTangent,endTangent);
			if (dot(normal, normal) < 1e-4 * 1e-4) {
				console.log("Flat segment; no natural frame");
				if (Math.abs(startTangent.x) <= Math.abs(startTangent.y) && Math.abs(startTangent.x) <= Math.abs(startTangent.y)) {
					normal = {x:1, y:0, z:0};
				} else if (Math.abs(startTangent.y) <= Math.abs(startTangent.z)) {
					normal = {x:0, y:1, z:0};
				} else {
					normal = {x:0, y:0, z:1};
				}
				let d = dot(normal, startTangent);
				normal.x -= d * startTangent.x;
				normal.y -= d * startTangent.y;
				normal.z -= d * startTangent.z;
			}
			normal = normalize(normal);


			if (perp === null) {
				//just use the normal:
				perp = {x:normal.x, y:normal.y, z:normal.z};
			}

			{ //correct perp in case it isn't exactly perpendicular to the start tangent:
				let d = dot(perp, startTangent);
				perp.x -= d * startTangent.x;
				perp.y -= d * startTangent.y;
				perp.z -= d * startTangent.z;
				perp = normalize(perp);
			}

			let startBinormal = cross(startTangent, normal);

			//record perp in terms of normal and binormal:
			let perpN = dot(normal, perp);
			let perpB = dot(startBinormal, perp);

			let endBinormal = cross(endTangent, normal);

			//parallel transport perp along the path:
			perp = {
				x:perpN * normal.x + perpB * endBinormal.x,
				y:perpN * normal.y + perpB * endBinormal.y,
				z:perpN * normal.z + perpB * endBinormal.z
			};


			let len = Math.sqrt((e.x-s.x)**2 + (e.y-s.y)**2 + (e.z-s.z)**2) * 0.5;

			RSMEAttribs.push( perpN, perpB,
				s.x, s.y, s.z, totalLen,
				m.x, m.y, m.z, totalLen + 0.5 * len,
				e.x, e.y, e.z, totalLen + len,
				normal.x, normal.y, normal.z
			);
			totalLen += len;

		}

		if (yarn.points.length >= 6) {
			let prev = { x:yarn.points[0], y:yarn.points[1], z:yarn.points[2] };
			let next = { x:yarn.points[3], y:yarn.points[4], z:yarn.points[5] };

			segment(prev,
				{x:0.25 * (next.x - prev.x) + prev.x, y:0.25 * (next.y - prev.y) + prev.y, z:0.25 * (next.z - prev.z) + prev.z},
				{x:0.5 * (prev.x + next.x), y:0.5 * (prev.y + next.y), z:0.5 * (prev.z + next.z)}
			);
		}
		for (let i = 8; i < yarn.points.length; i += 3) {
			let prev = { x:yarn.points[i-8], y:yarn.points[i-7], z:yarn.points[i-6] };
			let cur  = { x:yarn.points[i-5], y:yarn.points[i-4], z:yarn.points[i-3] };
			let next = { x:yarn.points[i-2], y:yarn.points[i-1], z:yarn.points[i] };

			segment(
				{x:0.5 * (prev.x + cur.x), y:0.5 * (prev.y + cur.y), z:0.5 * (prev.z + cur.z)},
				cur,
				{x:0.5 * (next.x + cur.x), y:0.5 * (next.y + cur.y), z:0.5 * (next.z + cur.z)}
			);
		}
		if (yarn.points.length > 6) {
			//last segment
			let prev = { x:yarn.points[yarn.points.length-6], y:yarn.points[yarn.points.length-5], z:yarn.points[yarn.points.length-4] };
			let next = { x:yarn.points[yarn.points.length-3], y:yarn.points[yarn.points.length-2], z:yarn.points[yarn.points.length-1] };

			segment(
				{x:0.5 * (prev.x + next.x), y:0.5 * (prev.y + next.y), z:0.5 * (prev.z + next.z)},
				{x:0.75 * (next.x - prev.x) + prev.x, y:0.75 * (next.y - prev.y) + prev.y, z:0.75 * (next.z - prev.z) + prev.z},
				next
			);
		}

	});

	gl.bindBuffer(gl.ARRAY_BUFFER, perInstanceBuffer);
	gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(RSMEAttribs), gl.STATIC_DRAW);
	instanceCount = RSMEAttribs.length / 17;
	console.log("Have " + instanceCount + " instances.");

	let randomIndex = 0;
	function pseudoRandom(min, max) {
		//from: https://stackoverflow.com/questions/4200224/random-noise-functions-for-glsl
		const amt = (Math.sin(randomIndex * 78.233) * 43758.5453) % 1;
		randomIndex += 1;
		return (0.5 * amt + 0.5) * (max - min) + min;
	}

	let ParamSideAttribs = [];
	let radius = yarns.yarns[0].radius; //overall yarn radius
	for (let p = 0; p < 3; ++p) {
		let ply_add = p / 3.0;
		let ply_mul = 1.0 / (5.0 * radius); //one twist every five radii
		let ply_rad = radius * 0.5;
		let ply_squish = 0.8;
		for (let s = 0; s < 20; ++s) {
			let strand_angle_add = s / 20.0;
			let strand_angle_mul = -1.0 / (5.0 * radius); //tighter strand twist than ply twist, opposite direction
			let strand_radius_add = pseudoRandom(0.0,1.0);
			let strand_radius_mul = 1.0 / (pseudoRandom(5.0, 10.0) * radius);
			let strand_rad1 = pseudoRandom(-0.1, 0.5) * ply_rad;
			let strand_rad2 = pseudoRandom(0.5, 1.1) * ply_rad;

			let strand_radius = 0.1 * 0.5 * radius;
			for (let i = 0; i <= 10; ++i) {
				let amt = i / 10.0;
				if (i == 0 && ParamSideAttribs.length > 0) ParamSideAttribs.push( ...ParamSideAttribs.slice(-12) );
				ParamSideAttribs.push(amt,-strand_radius,
					ply_mul, ply_add, ply_rad, ply_squish,
					strand_angle_mul, strand_angle_add,
					strand_radius_mul, strand_radius_add, strand_rad1, strand_rad2);
				if (i == 0 && ParamSideAttribs.length > 12) ParamSideAttribs.push( ...ParamSideAttribs.slice(-12) );
				ParamSideAttribs.push(amt, strand_radius,
					ply_mul, ply_add, ply_rad, ply_squish,
					strand_angle_mul, strand_angle_add,
					strand_radius_mul, strand_radius_add, strand_rad1, strand_rad2);
			}
		}
	}
	gl.bindBuffer(gl.ARRAY_BUFFER, perVertexBuffer);
	gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(ParamSideAttribs), gl.STATIC_DRAW);
	vertexCount = ParamSideAttribs.length / 12;
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
		4*12, //stride
		0 //offset
	);
	gl.enableVertexAttribArray(stripProgram.attribLocations.ParamSide);

	if ('PlyMulAddRadSquish' in stripProgram.attribLocations) {
		gl.vertexAttribPointer(stripProgram.attribLocations.PlyMulAddRadSquish,
			4, //size
			gl.FLOAT, //type
			false, //normalize
			4*12, //stride
			4*2 //offset
		);
		gl.enableVertexAttribArray(stripProgram.attribLocations.PlyMulAddRadSquish);
	}

	if ('SAMulAdd' in stripProgram.attribLocations) {
		gl.vertexAttribPointer(stripProgram.attribLocations.SAMulAdd,
			2, //size
			gl.FLOAT, //type
			false, //normalize
			4*12, //stride
			4*6 //offset
		);
		gl.enableVertexAttribArray(stripProgram.attribLocations.SAMulAdd);
	}

	if ('SRMulAddRad1Rad2' in stripProgram.attribLocations) {
		gl.vertexAttribPointer(stripProgram.attribLocations.SRMulAddRad1Rad2,
			4, //size
			gl.FLOAT, //type
			false, //normalize
			4*12, //stride
			4*8 //offset
		);
		gl.enableVertexAttribArray(stripProgram.attribLocations.SRMulAddRad1Rad2);
	}


	gl.bindBuffer(gl.ARRAY_BUFFER, perInstanceBuffer);

	gl.vertexAttribPointer(stripProgram.attribLocations.PerpNPerpB,
		2, //size
		gl.FLOAT, //type
		false, //normalize
		4*17, //stride
		4*0 //offset
	);
	ext.vertexAttribDivisorANGLE(stripProgram.attribLocations.Radius, 1);
	gl.enableVertexAttribArray(stripProgram.attribLocations.Radius);

	gl.vertexAttribPointer(stripProgram.attribLocations.Start,
		4, //size
		gl.FLOAT, //type
		false, //normalize
		4*17, //stride
		4*2 //offset
	);
	ext.vertexAttribDivisorANGLE(stripProgram.attribLocations.Start, 1);
	gl.enableVertexAttribArray(stripProgram.attribLocations.Start);

	gl.vertexAttribPointer(stripProgram.attribLocations.Middle,
		4, //size
		gl.FLOAT, //type
		false, //normalize
		4*17, //stride
		4*6 //offset
	);
	ext.vertexAttribDivisorANGLE(stripProgram.attribLocations.Middle, 1);
	gl.enableVertexAttribArray(stripProgram.attribLocations.Middle);

	gl.vertexAttribPointer(stripProgram.attribLocations.End,
		4, //size
		gl.FLOAT, //type
		false, //normalize
		4*17, //stride
		4*10 //offset
	);
	ext.vertexAttribDivisorANGLE(stripProgram.attribLocations.End, 1);
	gl.enableVertexAttribArray(stripProgram.attribLocations.End);

	if ('Normal' in stripProgram.attribLocations) {
		gl.vertexAttribPointer(stripProgram.attribLocations.Normal,
			3, //size
			gl.FLOAT, //type
			false, //normalize
			4*17, //stride
			4*14 //offset
		);
		ext.vertexAttribDivisorANGLE(stripProgram.attribLocations.Normal, 1);
		gl.enableVertexAttribArray(stripProgram.attribLocations.Normal);
	}

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

	const to_light = normalize({x:1.0, y:1.0, z:0.5});
	gl.uniform3fv(
		stripProgram.uniformLocations.TO_LIGHT_vec3,
		new Float32Array([to_light.x, to_light.y, to_light.z])
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
	ext.vertexAttribDivisorANGLE(stripProgram.attribLocations.Normal, 0);

	gl.disable(gl.CULL_FACE);
	gl.disable(gl.DEPTH_TEST);
};


})();
