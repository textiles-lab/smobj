#include "sm.hpp"

#include <glm/gtx/hash.hpp>
#include <glm/gtx/norm.hpp>
#include <glm/gtx/string_cast.hpp>

#include <iostream>
#include <fstream>
#include <cstring>
#include <unordered_map>
#include <map>


//convert a quadratic bezier through midpoints (i.e., a quadratic b-spline, if I recall correctly) to a Cubic Catmull-Rom Spline, preserving marked points, and keeping the curves within 'tol' of each-other [at all test points]
void midpoint_bezier_to_catmull_rom(std::vector< glm::vec3 > const &yarn, std::vector< uint32_t > yarn_marked, std::vector< glm::vec3 > *rom, std::vector< uint32_t > *rom_marked, const float tol);

int main(int argc, char **argv) {
	if (argc != 3) {
		std::cerr << "Usage:\n\t./yarns-to-bccx <in.yarns> <out.bccx>" << std::endl;
		return 1;
	}
	std::string in_yarns = argv[1];
	std::string out_bccx = argv[2];

	std::cout << "Will convert yarns in '" << in_yarns << "' to a bccx-formatted file '" << out_bccx << "'." << std::endl;

	sm::Yarns yarns = sm::Yarns::load(in_yarns);

	{
		uint32_t total_points = 0;
		uint32_t total_checkpoints = 0;
		for (auto const &yarn : yarns.yarns) {
			total_points += yarn.points.size();
			total_checkpoints += yarn.checkpoints.size();
		}
		std::cout << "Input has " << yarns.yarns.size() << " yarns with a total of " << total_points << " points and " << total_checkpoints << " checkpoints." << std::endl;
		if (total_checkpoints != 0) {
			std::cerr << "WARNING: checkpoint conversion code not yet written." << std::endl;
		}
	}

	//convert yarns-style splines to bccx-style splines:
	std::vector< std::vector< glm::vec3 > > rom_points;
	std::vector< std::vector< uint32_t > > rom_checkpoints;
	rom_points.reserve(yarns.yarns.size());
	rom_checkpoints.reserve(yarns.yarns.size());
	for (auto const &yarn : yarns.yarns) {
		//make a list of checkpoint locations to preserve in the conversion:
		std::vector< uint32_t > checkpoints;
		checkpoints.reserve(yarn.checkpoints.size());
		for (auto const &cp : yarn.checkpoints) {
			if (checkpoints.empty() || cp.point != checkpoints.back()) {
				checkpoints.emplace_back(cp.point);
			}
		}
		//do the conversion:
		std::vector< glm::vec3 > rom;
		std::vector< uint32_t > rom_cp;
		midpoint_bezier_to_catmull_rom(yarn.points, checkpoints, &rom, &rom_cp, 0.1 * yarn.radius);
		assert(rom_cp.size() == checkpoints.size());
		for (auto cp : rom_cp) {
			assert(cp < rom.size());
		}

		rom_points.emplace_back(std::move(rom));
		rom_checkpoints.emplace_back(std::move(rom_cp));

	}


	std::ofstream out(out_bccx, std::ios::binary);

	//BCC(x) format is based on 'sim.h' and 'bcc.cpp' from YarnBallSim
	// as well as http://www.cemyuksel.com/cyCodeBase/soln/using_bcc_files.html
	struct BCCHeader {
		char magic[3] = {'B', 'C', 'C'};
		uint8_t bytes = 0x44; //4-byte integers, 4-byte floating point
		char format[2] = {'C', '0'}; //uniform catmull-rom splines
		uint8_t dimensions = 3; //3d points
		uint8_t up_direction = 2; //1 is y-up, 2 is z-up
		uint64_t num_yarns;
		uint64_t num_nodes;
		char info[40];
	} header;
	static_assert(sizeof(BCCHeader) == 64, "BCC header should be 64 bytes.");

	header.num_yarns = yarns.yarns.size();
	header.num_nodes = 0;
	for (auto const &points : rom_points) {
		header.num_nodes += points.size();
	}
	{ //copy source filename into info:
		memset(header.info, '\0', sizeof(header.info));
		std::string info = "from '" + in_yarns + "'";
		if (info.size() > 40) info = info.substr(0, 40);
		memcpy(header.info, info.data(), info.size());
	}
	out.write(reinterpret_cast< const char * >(&header), sizeof(header));

	//Actually, we fixed that:
	//std::cerr << "WARNING: copying control points directly from yarns to bccx file -- this is technically incorrect in that yarns uses quadratic bezier curves between midpoints of yarns polylines, while bccx files use uniform cubic catmull-rom splines through yarn points." << std::endl;
	for (auto &yarn : yarns.yarns) {
		uint32_t index = &yarn - &yarns.yarns[0];
		int32_t count = rom_points[index].size();
		//NOTE: if count is negative, yarn is a closed loop; yarns format doesn't store closed loops
		//though, if yarns ever were used to store closed loops, it would be by having the first and last point in the yarn occupy exactly the same position:
		if (!yarn.points.empty() && yarn.points[0] == yarn.points.back()) {
			count = -count;
		}
		out.write(reinterpret_cast< const char * >(&count), sizeof(count));
		static_assert(sizeof(rom_points[index][0]) == 12, "yarn points are 3-vectors of floating point values");
		out.write(reinterpret_cast< const char * >(rom_points[index].data()), rom_points[index].size() * sizeof(rom_points[index][0]));
	}

	{ //bccx stuff:
		//yarn radius:
		float radius = 0.1f;
		if (!yarns.yarns.empty()) {
			radius = yarns.yarns[0].radius;
			for (auto const &yarn : yarns.yarns) {
				if (yarn.radius != radius) {
					std::cerr << "WARNING: bccx supports only a single radius, but yarns[0] has radius " << radius << " and yarns[" << (&yarn - &yarns.yarns[0]) << "] has radius " << yarn.radius << ". (Using yarns[0] radius.)" << std::endl;
				}
			}
		}
		out.write(reinterpret_cast< const char * >(&radius), sizeof(radius));
		//yarn point velocities:
		for (auto &yarn : yarns.yarns) {
			uint32_t index = &yarn - &yarns.yarns[0];
			for ([[maybe_unused]] auto const &point : rom_points[index]) {
				glm::vec3 velocity = glm::vec3(0.0f, 0.0f, 0.0f);
				out.write(reinterpret_cast< const char * >(&velocity), sizeof(velocity));
			}
		}
		//yarn point rest lengths:
		// (YarnBallSim's bcc.cpp initializes these to distance-to-next-control-point, so that's what I will use for now; could use checkpoints for this in the future.)
		for (auto &yarn : yarns.yarns) {
			uint32_t index = &yarn - &yarns.yarns[0];
			auto &rp = rom_points[index];
			//auto &cp = rom_checkpoints[index]; //<-- eventually, for checkpoints
			if (yarn.checkpoints.size()) {
				std::cerr << "WARNING: ignoring checkpoints in yarn[" << (&yarn - &yarns.yarns[0]) << "] -- reading lengths from checkpoints is not yet supported by this conversion code." << std::endl;
			}
			for (size_t i = 0; i < rp.size(); ++i) {
				float l = glm::length(rp[i] - rp[(i + 1) % rp.size()]);
				out.write(reinterpret_cast< const char * >(&l), sizeof(l));
			}
		}
	}

	return 0;
}


//based on 'simplify()' function from simplify-yarns.cpp
void midpoint_bezier_to_catmull_rom(std::vector< glm::vec3 > const &yarn, std::vector< uint32_t > yarn_marked, std::vector< glm::vec3 > *rom_, std::vector< uint32_t > *rom_marked_, const float tol) {
	assert(rom_);
	auto &rom = *rom_;
	assert(rom_marked_);
	auto &rom_marked = *rom_marked_;

	rom.clear();
	rom_marked.clear();

	if (yarn.size() == 0) {
		return; //already done!
	} else if (yarn.size() <= 2) {
		//too simple to approximate, just make a straight line:
		rom.emplace_back(yarn[0]);
		rom.insert(rom.end(), yarn.begin(), yarn.end());
		rom.emplace_back(yarn.back());

		rom_marked = yarn_marked;
		for (auto &m : rom_marked) {
			m += 1;
		}
		return;
	}

	const float tol2 = tol*tol;

	//overall idea:
	// (1) establish list of "test points" that are on the yarn and have associated parameter values
	// (2) search through space of possible catmull-rom splines by looking for best spline
	//     - location in search is given by last three control points;
	//     - action in search is picking one new control point;
	//     - cost is number of control points picked.
	//  (could use a heuristic here if something admissible existed...)

	std::vector< glm::vec3 > pts;
	std::vector< uint32_t > pts_marked;
	{ //dice spline into pts:

		float max_step = 10.0f * tol;

		auto line_to = [&pts, &max_step](glm::vec3 const &to) {
			glm::vec3 from = pts.back();
			//NOTE: variable sampling rate is going to do weird things with catmull-rom spline parameter speed (or, hmm, is going to favor nice uniform parameter speed which is maybe fine)
			float length = glm::length(to - from);
			uint32_t steps = std::max(1, int32_t(std::ceil(length / max_step)));
			for (uint32_t i = 1; i < steps; ++i) {
				pts.emplace_back(glm::mix(from, to, float(i) / float(steps)));
			}
			pts.emplace_back(to);
		};

		auto curve_to = [&pts, &max_step](glm::vec3 const &m, glm::vec3 const &e) {
			glm::vec3 s = pts.back();

			//NOTE: using a length approximation *certainly* will mess with catmull-rom spline speed:
			float approx_length = glm::length(m - s) + glm::length(e - m);
			uint32_t steps = std::max(1, int32_t(std::ceil(approx_length / max_step)));

			for (uint32_t i = 1; i < steps; ++i) {
				float t = float(i) / float(steps);
				pts.emplace_back(glm::mix( glm::mix(s,m,t), glm::mix(m,e,t), t));
			}
			pts.emplace_back(e);
		};


		std::vector< uint32_t >::iterator next_marked = yarn_marked.begin();

		//first part (is a straight line):
		if (next_marked != yarn_marked.end() && *next_marked == 0) {
			pts_marked.emplace_back(0);
			++next_marked;
		}
		pts.emplace_back(yarn[0]);
		line_to(glm::mix(yarn[0], yarn[1], 0.5f));

		//next part (is quadratic segments):
		for (uint32_t i = 1; i + 1 < yarn.size(); ++i) {
			glm::vec3 s = pts.back();
			glm::vec3 m = yarn[i];
			glm::vec3 e = glm::mix(yarn[i], yarn[i+1], 0.5f);
			
			glm::vec3 sm = glm::mix(s,m,0.5f);
			glm::vec3 me = glm::mix(m,e,0.5f);
			glm::vec3 sme = glm::mix(sm,me,0.5f);

			//first half of curve:
			curve_to(sm, sme);
			//marked point:
			if (next_marked != yarn_marked.end() && *next_marked == i) {
				pts_marked.emplace_back(pts.size()-1);
				++next_marked;
			}
			//second half of curve:
			curve_to(me, e);
		}

		//final part (is also a straight line):
		line_to(yarn.back());
		if (next_marked != yarn_marked.end() && *next_marked == yarn.size()-1) {
			pts_marked.emplace_back(pts.size()-1);
			++next_marked;
		}
		assert(next_marked == yarn_marked.end());
		assert(pts_marked.size() == yarn_marked.size());
	}
	std::vector< bool > is_marked(pts.size(), false);
	for (auto m : pts_marked) {
		is_marked[m] = true;
	}

	std::cout << "Sampled " << yarn.size() << " (" << yarn_marked.size() << " marked) into " << pts.size() << " test points; fitting..." << std::endl;

	//Now approximate test points with catmull-rom spline:

	//helper to check distance of points in [i1,i2] given control points [i0,i1,i2,i3] have been selected.
	//special cases:
	//  if i0 == -1U: (i1 must be 0) start of spline, point for a is reflection of i2 over i1
	//  if i3 == -1U: (i2 must be pts.size()-1) end of spline, point for i3 is reflection of i1 over i2
	auto check_dis2 = [&pts](uint32_t i0, uint32_t i1, uint32_t i2, uint32_t i3) {
		//look up control points, handling edge conditions:
		assert(i1 < pts.size());
		glm::vec3 p1 = pts[i1];
		assert(i2 < pts.size());
		glm::vec3 p2 = pts[i2];
		glm::vec3 p0;
		if (i0 == -1U) {
			assert(i1 == 0);
			p0 = p1 - (p2 - p1);
		} else {
			assert(i0 < pts.size());
			p0 = pts[i0];
		}
		glm::vec3 p3;
		if (i3 == -1U) {
			assert(i2 == pts.size()-1);
			p3 = p2 - (p1 - p2);
		} else {
			assert(i3 < pts.size());
			p3 = pts[i3];
		}

		float max_dis2 = 0.0f;
		//Note: spline meets i1, i2 exactly, so only testing internal points
		for (uint32_t i = i1 + 1; i < i2; ++i) {
			/* pyramid method:
			//compute point on spline as per https://en.wikipedia.org/wiki/Centripetal_Catmull%E2%80%93Rom_spline
			//time relative to knots at 0 (pa), 1 (pb), 2 (pc), 3 (pd):
			float t = 1.0f + float(i - i1) / float(i2 - i1);
			glm::vec3 a1 = glm::mix(p0,p1, (t-0.0f) / (1.0f - 0.0f));
			glm::vec3 a2 = glm::mix(p1,p2, (t-1.0f) / (2.0f - 1.0f));
			glm::vec3 a3 = glm::mix(p2,p3, (t-2.0f) / (3.0f - 2.0f));
			glm::vec3 b1 = glm::mix(a1,a2, (t-0.0f) / (2.0f - 0.0f));
			glm::vec3 b2 = glm::mix(a2,a3, (t-1.0f) / (3.0f - 1.0f));
			glm::vec3 at = glm::mix(b1,b2, (t-1.0f) / (2.0f - 1.0f));
			*/

			//compute point on spline as per https://www.mvps.org/directx/articles/catmull/
			float s = float(i-i1) / float(i2-i1);
			glm::vec4 coefs = 0.5f * glm::vec4(1.0f, s, s*s, s*s*s);
			/*
			glm::vec3 at = (coefs * glm::transpose(glm::mat4(
				 0.0f, 2.0f, 0.0f, 0.0f,
				-1.0f, 0.0f, 1.0f, 0.0f,
				 2.0f,-5.0f, 4.0f,-1.0f,
				-1.0f, 3.0f,-3.0f, 1.0f
			))) * glm::transpose(glm::mat4x3(p0, p1, p2, p3));
			*/

			//basis polynomial method (without matrices):
			// -- seems fastest in my [very rough] tests
			glm::vec3 at =
				(-coefs[1] + 2.0f * coefs[2] - coefs[3]) * p0
				+ (2.0f * coefs[0] - 5.0f * coefs[2] + 3.0f * coefs[3]) * p1
				+ (coefs[1] + 4.0f * coefs[2] - 3.0f * coefs[3]) * p2
				+ (-coefs[2] + coefs[3]) * p3;

			//assert(glm::length2(at - at2) < 0.001f);

			float dis2 = glm::length2(at - pts[i]);
			max_dis2 = std::max(dis2, max_dis2);
		}

		return max_dis2;
	};

	//Search for best path, where steps go from triples of points to triples of points (extend spline by one segment):
	struct DistanceFrom {
		DistanceFrom(uint32_t distance_, uint32_t from_) : distance(distance_), from(from_) { }
		uint32_t distance;
		uint32_t from;
	};
	std::unordered_map< glm::uvec3, DistanceFrom > distance; //distance == # of points to get to pair
	std::multimap< std::pair< uint32_t, uint32_t >, glm::uvec3 > to_expand; //to expand sorted by heuristic & distance

	auto queue = [&distance, &to_expand](glm::uvec3 at, DistanceFrom df) {
		auto ret = distance.emplace(at, df);
		if (ret.second || ret.first->second.distance > df.distance) {
			ret.first->second = df;
			uint32_t heuristic = df.distance; // + (pts.size() - 1 - at.y); //<-- NOTE: inadmissable heuristic
			to_expand.emplace(std::make_pair(heuristic, df.distance), at);
		}
	};

	constexpr uint32_t MaxStep = 80; //hmm, might want to adjust this?

	//try a bunch of different starting strides:
	for (uint32_t s = 1; s <= MaxStep; ++s) {
		queue(glm::uvec3(-1U, 0, s), DistanceFrom(2, -1U));
	}

	uint32_t best_dis = -1U;
	glm::uvec3 best_end = glm::uvec3(-1U);


	while (!to_expand.empty()) {
		uint32_t dis = to_expand.begin()->first.second;
		glm::uvec3 at = to_expand.begin()->second;
		to_expand.erase(to_expand.begin());
		{
			auto f = distance.find(at);
			assert(f != distance.end());
			assert(f->second.distance <= dis);
			if (f->second.distance < dis) continue;
		}

		if (at.z + 1 == pts.size()) {
			if (dis < best_dis) {
				std::cout << "Reached end in " << dis << " points." << std::endl;
				best_dis = dis;
				best_end = at;
			}
		}


		for (uint32_t s = 1; s < MaxStep; ++s) {
			if (at.z + s >= pts.size()) break;
			//TODO: must also not skip marked points!
			if (check_dis2(at.x, at.y, at.z, at.z + s) <= tol2) {
				//have to also check ending:
				if (at.z + s + 1 != pts.size() || check_dis2(at.y, at.z, at.z + s, -1U) <= tol2) {
					queue(glm::uvec3(at.y, at.z, at.z + s), DistanceFrom(dis + 1, at.x));
				}
			}
			if (is_marked[at.z + s]) break; //can't skip over marked points
		}
	}

	assert(best_end.x < best_end.y && best_end.y < best_end.z && best_end.z < pts.size());
	std::vector< uint32_t > path;
	{ //read back whole path (in reverse):
		path.emplace_back(best_end.z);
		path.emplace_back(best_end.y);
		path.emplace_back(best_end.x);

		uint32_t dis = best_dis;
		
		while (path.back() != 0) {
			auto f = distance.find(glm::uvec3(
				path[path.size()-1],
				path[path.size()-2],
				path[path.size()-3]
			));
			assert(f != distance.end());
			path.emplace_back(f->second.from);
		}
		std::reverse(path.begin(), path.end());
		std::cout << "Read out path of length " << path.size() << "." << std::endl;
		assert(dis == path.size());
	}

	rom.reserve(path.size()+2);
	rom.emplace_back(glm::mix(pts[path[0]], pts[path[1]], -1.0f));
	for (auto &i : path) {
		if (is_marked[i]) rom_marked.emplace_back(rom.size());
		rom.emplace_back(pts[i]);
	}
	rom.emplace_back(glm::mix(pts[path[path.size()-2]], pts[path[path.size()-1]], 2.0f));
	assert(rom_marked.size() == yarn_marked.size());

}
