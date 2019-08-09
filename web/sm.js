"use strict";
(function(){

const sm = (typeof(module) === 'undefined' ? (window.sm = {}) : module.exports );

sm.Yarns = function() {
	this.yarns = [];
	this.units = []; //{name:'1', length:1.0} ];
};

/*
 * Return yarns object created from ArrayBuffer,
 * throw on error.
 */

sm.Yarns.fromArrayBuffer = function(buffer) {
	const yarns = new sm.Yarns();

	const view = new DataView(buffer);

	let offset = 0;

	//helpers: read value and update offset
	function getUint8() {
		const got = view.getUint8(offset);
		offset += 1;
		return got;
	}
	function getUint32() {
		const got = view.getUint32(offset, true);
		offset += 4;
		return got;
	}
	function getFloat32() {
		const got = view.getFloat32(offset, true);
		offset += 4;
		return got;
	}

	//helper: read chunk header from view
	function getCount(fourcc, itemSize) {
		console.assert(fourcc.length === 4, "fourcc supplied to getCount should have four characters");
		let gotFourcc = '';
		gotFourcc += String.fromCharCode(getUint8());
		gotFourcc += String.fromCharCode(getUint8());
		gotFourcc += String.fromCharCode(getUint8());
		gotFourcc += String.fromCharCode(getUint8());
		if (gotFourcc !== fourcc) {
			throw new Error("Expecting chunk fourcc '" + fourcc + "', got fourcc '" + gotFourcc + "'.");
		}
		let gotSize = getUint32();
		if (gotSize % itemSize !== 0) {
			throw new Error("Expecting chunk size (" + gotSize + ") to be a multiple of item size (" + itemSize + ").");
		}
		return gotSize / itemSize;
	}

	const pointsCount = getCount('f3..', 12);
	console.log(pointsCount + " points.");
	const points = [];
	for (let i = 0; i < pointsCount; ++i) {
		let x = getFloat32();
		let y = getFloat32();
		let z = getFloat32();
		points.push(x,y,z);
	}
	//helper: get points from range
	function getPoints(begin, end) {
		return points.slice(3*begin, 3*end);
	}

	const sourcesCount = getCount('src.', 4);
	console.log(sourcesCount + " sources.");
	const sources = [];
	for (let i = 0; i < sourcesCount; ++i) {
		const source = getUint32();
		sources.push(source);
	}
	if (3*sources.length !== points.length) {
		throw new Error("Yarns file should have as many sources as points.");
	}

	//helper: get sources from range
	function getSources(begin, end) {
		return sources.slice(begin, end);
	}

	//build map for checkpoint sorting:
	const pointToYarnIndex = new Array(points.length);
	const yarnPointBegin = [];

	const yarnInfosCount = getCount('yarn', 16);
	console.log(yarnInfosCount + " yarnInfos.");
	for (let i = 0; i < yarnInfosCount; ++i) {
		const pointBegin = getUint32();
		const pointEnd = getUint32();
		const radius = getFloat32();
		const r = getUint8();
		const g = getUint8();
		const b = getUint8();
		const a = getUint8();
		if (!(pointBegin <= pointEnd && pointEnd < points.length)) {
			throw new Error("Invalid point indices [" + pointBegin + ", " + pointEnd + ") of " + points.length + " in yarn info.");
		}
		yarns.yarns.push({
			points:getPoints(pointBegin, pointEnd),
			sources:getSources(pointBegin, pointEnd),
			radius:radius,
			color:{ r:r, g:g, b:b, a:a },
			checkpoints:[]
		});

		//update checkpoint lookup map:
		yarnPointBegin.push(pointBegin);
		for (let p = pointBegin; p < pointEnd; ++p) {
			pointToYarnIndex[p] = i;
		}
	}

	const stringsBytes = getCount('strs', 1);
	console.log(stringsBytes + " bytes of strings.");
	const strings = new Uint8Array(buffer, offset, stringsBytes);
	offset += stringsBytes;
	//helper: get string from byte range
	function getString(begin, end) {
		//TODO: consider utf8 decode
		let string = '';
		for (let i = begin; i < end; ++i) {
			string += String.fromCharCode(strings[i]);
		}
		return string;
	}

	const unitsCount = getCount('unit', 12);
	console.log(unitsCount + " units.");
	for (let i = 0; i < unitsCount; ++i) {
		const nameBegin = getUint32();
		const nameEnd = getUint32();
		const length = getFloat32();
		yarns.units.push({
			name:getString(nameBegin, nameEnd),
			length:length
		});
	}

	const checkpointsCount = getCount('chk.', 12);
	console.log(checkpointsCount + " checkpoints.");
	for (let i = 0; i < checkpointsCount; ++i) {
		const point = getUint32();
		const length = getFloat32();
		const unit = getUint32();
		if (unit >= yarns.units.length) {
			throw new Error("Checkpoint with out-of-range unit index " + unit + ".");
		}
		//store checkpoint in proper yarn:
		const yarnIndex = pointToYarnIndex[point];
		if (typeof(yarnIndex) !== 'number') {
			throw new Error("Checkpoint " + point + " doesn't lie in a yarn.");
		}
		yarns.yarns[yarnIndex].checkpoints.push({
			point:point-yarnPointBegin[yarnIndex],
			length:length,
			unit:unit
		});
	}


	yarns.yarns.forEach(function(yarn){
		yarn.sources = sources.slice(yarn.pointsBegin, yarn.pointsEnd);
		delete yarn.pointsBegin;
		delete yarn.pointsEnd;
	});


	return yarns;
};

})();
