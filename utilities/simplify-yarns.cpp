#include "sm.hpp"

#include <glm/gtx/hash.hpp>
#include <glm/gtx/norm.hpp>

#include <iostream>
#include <unordered_map>
#include <map>

std::vector< glm::vec3 > simplify(std::vector< glm::vec3 > const &pts, float tolerance) {
	//too short to simplify:
	if (pts.size() <= 2) return pts;

	float tol2 = tolerance * tolerance;

	//where points actually map given the spline-between-midpoints of the yarn:
	std::vector< glm::vec3 > spline_pts;
	spline_pts.reserve(pts.size());
	spline_pts.emplace_back(pts[0]);
	for (uint32_t i = 1; i + 1 < pts.size(); ++i) {
		glm::vec3 s = 0.5f * (pts[i-1] + pts[i]);
		glm::vec3 const &m = pts[i];
		glm::vec3 e = 0.5f * (pts[i] + pts[i+1]);
		glm::vec3 sm = glm::mix(s, m, 0.5f);
		glm::vec3 me = glm::mix(m, e, 0.5f);
		spline_pts.emplace_back(glm::mix(sm, me, 0.5f));
	}
	spline_pts.emplace_back(pts.back());

	//compute the fit given points a,b,c selected:
	auto fit_dis2 = [&pts, &spline_pts](uint32_t a, uint32_t b, uint32_t c) -> float {
		assert(a < b);
		assert(b < c);
		assert(c < pts.size());

		glm::vec3 s = glm::mix(pts[a], pts[b], 0.5f);
		glm::vec3 const &m = pts[b];
		glm::vec3 e = glm::mix(pts[b], pts[c], 0.5f);

		float max_dis2 = 0.0f;
		auto check = [&s, &m, &e, &spline_pts, &max_dis2](uint32_t i, float t) {
			glm::vec3 const &pt = spline_pts[i];
			glm::vec3 test = glm::mix( glm::mix(s,m,t), glm::mix(m,e,t), t);
			float dis2 = glm::length2(test - pt);
			max_dis2 = std::max(max_dis2, dis2);

		};
		//check all interior points:
		for (uint32_t i = a + (b-a+1)/2; i < b; ++i) {
			float t = float(i - a) / float(b - a) - 0.5f;
			assert(t >= 0.0f && t <= 1.0f); //N.b. might slightly miss owning to fp precision
			check(i,t);
		}
		check(b, 0.5f);
		for (uint32_t i = b + 1; i <= b + (c-b)/2; ++i) {
			float t = float(i - b) / float(c - b) + 0.5f;
			assert(t >= 0.0f && t <= 1.0f); //N.b. might slightly miss owning to fp precision
			check(i,t);
		}

		return max_dis2;
	};

	//Search for best path, where steps go from pairs of points to pairs of points (extend spline by one segment):
	std::unordered_map< glm::uvec2, uint32_t > distance; //distance == # of points to get to pair
	std::multimap< std::pair< uint32_t, uint32_t >, glm::uvec2 > to_expand; //to expand sorted by heuristic & distance

	auto queue = [&distance, &to_expand](glm::uvec2 at, uint32_t dis) {
		auto ret = distance.emplace(at, dis);
		if (ret.second || ret.first->second > dis) {
			ret.first->second = dis;
			uint32_t heuristic = dis; // + (pts.size() - 1 - at.y); //<-- NOTE: inadmissable heuristic
			to_expand.emplace(std::make_pair(heuristic, dis), at);
		}
	};

	queue(glm::uvec2(0, 1), 2);

	uint32_t best_dis = -1U;
	glm::uvec2 best_end = glm::uvec2(-1U);

	constexpr uint32_t MaxStep = 80;

	while (!to_expand.empty()) {
		uint32_t dis = to_expand.begin()->first.second;
		glm::uvec2 at = to_expand.begin()->second;
		to_expand.erase(to_expand.begin());
		{
			auto f = distance.find(at);
			assert(f != distance.end());
			assert(f->second <= dis);
			if (f->second < dis) continue;
		}

		//NOTE: could relax condition on at.x if we had a checking function for line segments
		if (at.x + 2 == pts.size() && at.y + 1 == pts.size()) {
			if (dis < best_dis) {
				std::cout << "Reached end in " << dis << " steps." << std::endl;
				best_dis = dis;
				best_end = at;
			}
		}

		for (uint32_t s = 1; s < MaxStep; ++s) {
			if (at.y + s >= pts.size()) break;
			if (fit_dis2(at.x, at.y, at.y + s) <= tol2) {
				queue(glm::uvec2(at.y, at.y + s), dis + 1);
			}
		}
	}

	assert(best_end.x < best_end.y && best_end.y < pts.size());
	std::vector< uint32_t > path;
	{ //read back whole path (in reverse):
		path.emplace_back(best_end.y);
		path.emplace_back(best_end.x);

		uint32_t dis = best_dis;
		
		while (path.back() != 0) {
			uint32_t b = path[path.size()-1];
			bool found = false;
			for (uint32_t s = 1; s < MaxStep; ++s) {
				if (s > b) break;
				uint32_t a = b - s;
				auto f = distance.find(glm::uvec2(a,b));
				if (f == distance.end()) continue;
				if (f->second + 1 > dis) continue;
				assert(f->second + 1 == dis);
				//step back:
				path.emplace_back(a);
				dis -= 1;
				found = true;
				break;
			}
			assert(found);
		}
		std::reverse(path.begin(), path.end());
		std::cout << "Read out path of length " << path.size() << "." << std::endl;
	}

	std::vector< glm::vec3 > new_pts;
	new_pts.reserve(path.size());
	for (auto &i : path) {
		new_pts.emplace_back(pts[i]);
	}

	return new_pts;
}

int main(int argc, char **argv) {
	if (argc != 4) {
		std::cerr << "Usage:\n\t./simplify-yarns <in.yarns> <tolerance|tolerance%> <out.yarns>" << std::endl;
		return 1;
	}
	std::string in_yarns = argv[1];
	float tolerance;
	bool tolerance_is_percent;
	std::string out_yarns = argv[3];
	try {
		std::string tol = argv[2];
		if (tol.size() > 0 && tol[tol.size()-1] == '%') {
			tolerance = std::stof(tol.substr(0,tol.size()-1));
			tolerance_is_percent = true;
		} else {
			tolerance = std::stof(tol);
			tolerance_is_percent = false;
		}
	} catch (std::exception &e) {
		std::cerr << "Failed to parse tolerance '" << argv[2] << "'." << std::endl;
		return 1;
	}

	std::cout << "Will simply yarns in '" << in_yarns << "' by making a path that remains within " << tolerance << (tolerance_is_percent ? " percent" : " units") << " and write to '" << out_yarns << "'." << std::endl;

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
			std::cerr << "WARNING: checkpoint simplification code not yet written." << std::endl;
		}
	}

	for (auto &yarn : yarns.yarns) {
		std::cout << "Simplifying yarn " << (&yarn - &yarns.yarns[0]) << " (radius " << yarn.radius << ", " << yarn.points.size() << " points)..."; std::cout.flush();
		std::vector< uint32_t > subset;
		float yarn_tol = (tolerance_is_percent ? yarn.radius * tolerance / 100.0f : tolerance);
		yarn.points = simplify(yarn.points, yarn_tol);
		yarn.sources.assign(yarn.points.size(), -1U); //NOTE: source info lost :-(
	}

	{
		uint32_t total_points = 0;
		uint32_t total_checkpoints = 0;
		for (auto const &yarn : yarns.yarns) {
			total_points += yarn.points.size();
			total_checkpoints += yarn.checkpoints.size();
		}
		std::cout << "Output has " << yarns.yarns.size() << " yarns with a total of " << total_points << " points and " << total_checkpoints << " checkpoints." << std::endl;
	}

	yarns.save(out_yarns);

	return 0;
}
