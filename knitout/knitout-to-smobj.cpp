
#include <glm/glm.hpp>

#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <vector>
#include <string>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <algorithm>

#include <cassert>

enum Direction : char {
	Left = '+',
	Right = '-',
};

enum Bed : char {
	Front = 'f',
	FrontSliders = 'F',
	BackSliders = 'B',
	Back = 'b',
};

//indices work like this:
// needle   next needle
//     v       v
// -1  0  1 2  3
//  ^-----^-^-- where carriers park (though they travel stitch<->track slightly beside these)

int32_t needle_index(int32_t needle) {
	return needle * 3;
}

int32_t side_index(int32_t needle, Direction dir) {
	return needle_index(needle) + (dir == Left ? -1 : +1);
}

struct FaceEdge {
	enum Flip {
		FlipNo,
		FlipYes
	};
	FaceEdge(uint32_t face_ = -1U, uint32_t edge_ = -1U, uint32_t count_ = 0, Flip flip_ = FlipNo) : face(face_), edge(edge_), count(count_), flip(flip_) { }
	uint32_t face;
	uint32_t edge;
	uint32_t count; //how many items (yarns/loops)
	Flip flip; //should this connection be flipped?
	bool is_valid() const {
		return (face != -1U && edge != -1U);
	}
};

struct Face {
	std::string type; //type edge ... edge
	std::vector< glm::vec3 > vertices;
};

struct Connection {
	Connection(FaceEdge const &a_, FaceEdge const &b_) : a(a_), b(b_) { }
	FaceEdge a,b;
};

//top edges of faces to avoid in the between-bed area:
// (highest-first)
struct TopEdge {
	TopEdge(glm::vec3 const &a_, glm::vec3 const &b_) : a(a_), b(b_) {
		if (a.y < b.y) std::swap(a,b);
	}
	glm::vec3 a,b;
	//kept sorted by highest point:
	bool operator<(TopEdge const &other) const {
		return a.y < other.a.y;
	}
};

struct Gizmo {
	std::vector< uint32_t > faces;
	std::vector< Connection > connections;
	std::vector< TopEdge > top_edges;
	float lift = 0.0f;
};

struct Column {
	FaceEdge top_edge; //by convention, oriented rightwards
	float top_y = 0;
};

struct BedColumns : public std::unordered_map< int32_t, Column > {
	float depth = 0.0f;
};

struct Carrier {
	float depth = 0.0f;
	//info about current position:
	FaceEdge edge; //is a valid edge if carrier is in, otherwise is not
	//by convention, flip on the FaceEdge will be set such that the (flipped) edge points up.
	float x = std::numeric_limits< float >::quiet_NaN();
};

//get ranges of lines within a given threshold of each-other:
bool close_ranges(
	float threshold,
	glm::vec2 const &a1, glm::vec2 const &b1,
	glm::vec2 const &a2, glm::vec2 const &b2,
	std::pair< float, float > *range1_, std::pair< float, float > *range2_ ) {
	assert(range1_);
	auto &range1 = *range1_;
	assert(range2_);
	auto &range2 = *range2_;
	range1.first = 0.0f;
	range1.second = 1.0f;
	range2.first = 0.0f;
	range2.second = 1.0f;

	//clip to range such that all elements of range are below amount when dotted with dir:
	auto clip = [](glm::vec2 const &a, glm::vec2 const &b, glm::vec2 const &dir, float amt, std::pair< float, float > *range) {
		float da = glm::dot(a, dir) - amt;
		float db = glm::dot(b, dir) - amt;
		if (da < 0.0f && db < 0.0f) {
			//entire range is good.
		} else if (da > 0.0f && db > 0.0f) {
			//clip entire range:
			range->first = 1.0f;
			range->second = 0.0f;
		} else if (da <= 0.0f && db > 0.0f) {
			range->second = std::min(range->second, (0.0f - da) / (db-da));
		} else if (da > 0.0f && db <= 0.0f) {
			range->first = std::max(range->first, (0.0f - da) / (db-da));
		} else {
			assert(0 && "unreachable case");
		}
	};

	{ //clip range2 vs (a1-b1):
		glm::vec2 ab1 = glm::normalize(b1-a1);
		glm::vec2 perp1 = glm::vec2(-ab1.y, ab1.x);
		clip(a2, b2, ab1, glm::dot(ab1, b1) + threshold, &range2);
		clip(a2, b2, -ab1, glm::dot(-ab1, a1) + threshold, &range2);
		clip(a2, b2, perp1, glm::dot(perp1, 0.5f*(a1+b1)) + threshold, &range2);
		clip(a2, b2, -perp1, glm::dot(-perp1, 0.5f*(a1+b1)) + threshold, &range2);
	}

	{ //clip range1 vs (a2-b2):
		glm::vec2 ab2 = glm::normalize(b2-a2);
		glm::vec2 perp2 = glm::vec2(-ab2.y, ab2.x);
		clip(a1, b1, ab2, glm::dot(ab2, b2) + threshold, &range1);
		clip(a1, b1, -ab2, glm::dot(-ab2, a2) + threshold, &range1);
		clip(a1, b1, perp2, glm::dot(perp2, 0.5f*(a2+b2)) + threshold, &range1);
		clip(a1, b1, -perp2, glm::dot(-perp2, 0.5f*(a2+b2)) + threshold, &range1);
	}
	return range1.first <= range1.second && range2.first <= range2.second;
}

//depths:
// --- back bed --- @ -1.0
// --- back hooks --- @ -0.75
//  [various crossings here]
// (carriers) @ [0.0 - 0.5]
// --- front hooks --- @ 0.75
// --- front bed ---   @ 1.0

constexpr const float FaceHeight = 1.0f;

struct Translator {
	BedColumns back_bed, back_sliders, front_sliders, front_bed;
	std::unordered_map< std::string, Carrier > carriers;
	float racking = 0.0f;


	Translator(std::vector< std::string > const &carrier_names) {
		back_bed.depth = -1.0f;
		back_sliders.depth = -0.75f;
		front_sliders.depth = 0.75f;
		front_bed.depth = 1.0f;
		for (uint32_t i = 0; i < carrier_names.size(); ++i) {
			Carrier carrier;
			carrier.depth = 0.5f * ((i + 0.5f) / float(carrier_names.size()));
			auto ret = carriers.insert(std::make_pair(carrier_names[i], carrier));
			if (!ret.second) throw std::runtime_error("Carrier '" + carrier_names[i] + "' is named twice.");
		}
	}

	//--- output ---
	std::vector< Face > faces;

	//--- helper functions ---
	std::vector< Carrier * > lookup_carriers(std::vector< std::string > const &cs) {
		std::vector< Carrier * > ret;
		ret.reserve(cs.size());
		std::unordered_set< std::string > seen;
		seen.reserve(cs.size());
		for (auto const &cn : cs) {
			if (!seen.insert(cn).second) {
				throw std::runtime_error("Carrier '" + cn + "' appears twice in carrier set.");
			}
			auto f = carriers.find(cn);
			if (f == carriers.end()) {
				throw std::runtime_error("Carrier '" + cn + "' does not appear in carriers list.");
			}
			ret.emplace_back(&f->second);
		}
		return ret;
	}
	BedColumns &lookup_bed(Bed bed) {
		if (bed == Back) return back_bed;
		if (bed == BackSliders) return back_sliders;
		if (bed == FrontSliders) return front_sliders;
		if (bed == Front) return front_bed;
		assert(0);
	}

	std::multiset< TopEdge > top_edges;
	void check_edge(glm::vec3 const &a, glm::vec3 const &b, float *lift_) {
		assert(lift_);
		auto &lift = *lift_;

		for (auto const &te : top_edges) {
			//okay, all later edges are lower, stop checking:
			if (te.a.y < std::min(a.y,b.y) + lift) break;
			//check for intersection in the xz plane:
			glm::vec2 a_2d(a.x, a.z);
			glm::vec2 b_2d(b.x, b.z);

			glm::vec2 te_a_2d(te.a.x, te.a.z);
			glm::vec2 te_b_2d(te.b.x, te.b.z);

			std::pair< float, float > range, te_range;
			if (!close_ranges(0.1f, a_2d, b_2d, te_a_2d, te_b_2d, &range, &te_range)) continue;

			float y = std::min(
				glm::mix(a.y, b.y, range.first),
				glm::mix(a.y, b.y, range.second)
			);
			float te_y = std::max(
				glm::mix(te.a.y, te.b.y, te_range.first),
				glm::mix(te.a.y, te.b.y, te_range.second)
			);
			lift = std::max(lift, te_y - y);
		}
	}



	//bring carriers together in order to make a stitch in a given direction on a given bed/needle.
	//(so, e.g., if dir is "Right", will bring yarns to the left of the needle to get them ready)
	//adds yarn movement faces and tracks with gizmo; returns the live yarn edge:
	FaceEdge bring_carriers(Direction dir, BedColumns const &bed, int32_t needle, std::vector< Carrier * > const &cs, Gizmo *gizmo_) {
		assert(gizmo_);
		auto &gizmo = *gizmo_;

		//no carriers? nothing to do:
		if (cs.empty()) return FaceEdge();

		//compute x-coordinate to travel along based on bed, racking, and direction:
		float travel_x; //where to travel in depth
		float final_x; //where to arrive on the bed
		final_x = side_index(needle, (dir == Right ? Left : Right));
		if (&bed == &front_bed || &bed == &front_sliders) {
			if (dir == Right) {
				travel_x = side_index(needle, Left) + 0.25f;
			} else { assert(dir == Left);
				travel_x = side_index(needle, Right) - 0.25f;
			}
		} else { assert(&bed == &back_bed || &bed == &back_sliders);
			//account for racking:
			uint32_t r = int32_t(std::floor(racking));
			if (std::floor(racking) == racking) {
				//aligned -- get to the front bed location across from the needle:
				if (dir == Right) {
					travel_x = side_index(needle + r, Left) + 0.25f;
				} else { assert(dir == Left);
					travel_x = side_index(needle + r, Right) - 0.25f;
				}
			} else { assert(std::floor(racking) + 0.25f == racking);
				//quarter pitch -- use the gap inside the indices to the right/left:
				if (dir == Right) {
					travel_x = side_index(needle + r, Right) + 0.25f;
				} else { assert(dir == Left);
					travel_x = side_index(needle + r + 1, Left) - 0.25f;
				}
			}
		}

		//all carriers travel in their tracks to the target index:
		for (auto c : cs) {
			if (c->edge.is_valid()) {
				assert(c->x != travel_x); //carriers are parked on indices, target is always 0.25f away.
			} else {
				c->x = travel_x + 0.25f; //a little "carrier coming in" tail
			}
			//build yarn face to connect things:
			faces.emplace_back();
			Face &face = faces.back();
			gizmo.faces.emplace_back(faces.size()-1);

			if (c->x < travel_x) {
				face.type = "yarn-to-right x +y1 x -y1";
				face.vertices = {
					glm::vec3(c->x, 0.0f, c->depth),
					glm::vec3(travel_x, 0.0f, c->depth),
					glm::vec3(travel_x, FaceHeight, c->depth),
					glm::vec3(c->x, FaceHeight, c->depth),
				};
				check_edge(face.vertices[0], face.vertices[1], &gizmo.lift);
				gizmo.top_edges.emplace_back(face.vertices[2], face.vertices[3]);

				FaceEdge yarn_in = FaceEdge(gizmo.faces.back(), 3, 1, FaceEdge::FlipYes);
				FaceEdge yarn_out = FaceEdge(gizmo.faces.back(), 1, 1, FaceEdge::FlipNo);
				if (c->edge.is_valid()) {
					gizmo.connections.emplace_back(c->edge, yarn_in);
				}
				c->edge = yarn_out;
				c->x = travel_x;

			} else { assert(c->x > travel_x); //can't have exactly equal.
				face.type = "yarn-to-left x -y1 x +y1";
				face.vertices = {
					glm::vec3(travel_x, 0.0f, c->depth),
					glm::vec3(c->x, 0.0f, c->depth),
					glm::vec3(c->x, FaceHeight, c->depth),
					glm::vec3(travel_x, FaceHeight, c->depth),
				};
				check_edge(face.vertices[0], face.vertices[1], &gizmo.lift);
				gizmo.top_edges.emplace_back(face.vertices[2], face.vertices[3]);

				FaceEdge yarn_in = FaceEdge(gizmo.faces.back(), 1, 1, FaceEdge::FlipNo);
				FaceEdge yarn_out = FaceEdge(gizmo.faces.back(), 3, 1, FaceEdge::FlipYes);
				if (c->edge.is_valid()) {
					gizmo.connections.emplace_back(c->edge, yarn_in);
				}
				c->edge = yarn_out;
				c->x = travel_x;
			}
		}

		//Connect all carriers up to the front/back:

		std::vector< Carrier * > cs_by_depth = cs;
		std::sort(cs_by_depth.begin(), cs_by_depth.end(), [](Carrier *a, Carrier *b) {
			return a->depth < b->depth;
		});
		if (cs.size() > 1) {
			std::cerr << "TODO: handle plating when bringing yarns." << std::endl;
		}
		if (&bed == &front_bed || &bed == &front_sliders) {
			FaceEdge prev = cs_by_depth.back()->edge;
			float prev_depth = cs_by_depth.back()->depth;

			//connect forward to bed depth:
			faces.emplace_back();
			Face &face = faces.back();
			gizmo.faces.emplace_back(faces.size()-1);
			std::string Y = std::to_string(prev.count);
			face.type = "yarn-to-right x +y" + Y + " x -y" + Y;
			face.vertices = {
				glm::vec3(travel_x, 0.0f, prev_depth),
				glm::vec3(final_x, 0.0f, bed.depth),
				glm::vec3(final_x, FaceHeight, bed.depth),
				glm::vec3(travel_x, FaceHeight, prev_depth),
			};
			check_edge(face.vertices[0], face.vertices[1], &gizmo.lift);
			gizmo.top_edges.emplace_back(face.vertices[2], face.vertices[3]);

			FaceEdge yarn_in = FaceEdge(gizmo.faces.back(), 3, 1, FaceEdge::FlipYes);
			FaceEdge yarn_out = FaceEdge(gizmo.faces.back(), 1, 1, FaceEdge::FlipNo);

			gizmo.connections.emplace_back(prev, yarn_in);

			return yarn_out;
		} else { assert(&bed == &back_bed || &bed == &back_sliders);
			FaceEdge prev = cs_by_depth.front()->edge;
			float prev_depth = cs_by_depth.front()->depth;

			{ //connect backward to zero depth:
				faces.emplace_back();
				Face &face = faces.back();
				gizmo.faces.emplace_back(faces.size()-1);
				std::string Y = std::to_string(prev.count);
				face.type = "yarn-to-left x -y" + Y + " x +y" + Y;
				face.vertices = {
					glm::vec3(travel_x, 0.0f, 0.0f),
					glm::vec3(travel_x, 0.0f, prev_depth),
					glm::vec3(travel_x, FaceHeight, prev_depth),
					glm::vec3(travel_x, FaceHeight, 0.0f),
				};
				check_edge(face.vertices[0], face.vertices[1], &gizmo.lift);
				gizmo.top_edges.emplace_back(face.vertices[2], face.vertices[3]);

				FaceEdge yarn_in = FaceEdge(gizmo.faces.back(), 1, 1, FaceEdge::FlipNo);
				FaceEdge yarn_out = FaceEdge(gizmo.faces.back(), 3, 1, FaceEdge::FlipYes);

				gizmo.connections.emplace_back(prev, yarn_in);

				prev = yarn_out;
				prev_depth = 0.0f;
			}

			{ //connect further backward to bed:
				faces.emplace_back();
				Face &face = faces.back();
				gizmo.faces.emplace_back(faces.size()-1);
				std::string Y = std::to_string(prev.count);
				face.type = "yarn-to-left x -y" + Y + " x +y" + Y;
				face.vertices = {
					glm::vec3(final_x, 0.0f, bed.depth),
					glm::vec3(travel_x, 0.0f, prev_depth),
					glm::vec3(travel_x, FaceHeight, prev_depth),
					glm::vec3(final_x, FaceHeight, bed.depth),
				};
				check_edge(face.vertices[0], face.vertices[1], &gizmo.lift);
				gizmo.top_edges.emplace_back(face.vertices[2], face.vertices[3]);

				FaceEdge yarn_in = FaceEdge(gizmo.faces.back(), 1, 1, FaceEdge::FlipNo);
				FaceEdge yarn_out = FaceEdge(gizmo.faces.back(), 3, 1, FaceEdge::FlipYes);

				gizmo.connections.emplace_back(prev, yarn_in);

				prev = yarn_out;
				prev_depth = bed.depth;
			}

			return prev;
		}

	};

	void return_carriers(Direction dir, BedColumns const &bed, int32_t needle, std::vector< Carrier * > const &cs, FaceEdge const &stitch_yarn_out, Gizmo *gizmo_) {
		assert(stitch_yarn_out.is_valid());
		assert(gizmo_);
		auto &gizmo = *gizmo_;

		//no carriers? nothing to do:
		if (cs.empty()) return;

		//compute travel location and final resting index based on bed, racking, and direction:
		float start_x = side_index(needle, dir);
		float final_x;
		if (&bed == &front_bed || &bed == &front_sliders) {
			final_x = side_index(needle, dir);
		} else { assert(&bed == &back_bed || &bed == &back_sliders);
			//account for racking:
			uint32_t r = int32_t(std::floor(racking));
			if (std::floor(racking) == racking) {
				//aligned -- get to the front bed location across from the needle:
				final_x = side_index(needle + r, dir);
			} else { assert(std::floor(racking) + 0.25f == racking);
				//quarter pitch -- use the gap inside the indices to the right/left:
				if (dir == Right) {
					final_x = side_index(needle + r + 1, Left);
				} else { assert(dir == Left);
					final_x = side_index(needle + r, Right);
				}
			}
		}
		float travel_x = final_x + (dir == Right ? -0.25f : +0.25f);

		std::vector< Carrier * > cs_by_depth = cs;
		std::sort(cs_by_depth.begin(), cs_by_depth.end(), [](Carrier *a, Carrier *b) {
			return a->depth < b->depth;
		});
		if (cs.size() > 1) {
			std::cerr << "TODO: handle plating when returning yarns." << std::endl;
		}
		if (&bed == &front_bed || &bed == &front_sliders) {
			FaceEdge prev = stitch_yarn_out;
			float prev_depth = bed.depth;
			float prev_x = start_x;

			{ //connect backward to first carrier:
				Carrier *c = cs_by_depth.front();
				faces.emplace_back();
				Face &face = faces.back();
				gizmo.faces.emplace_back(faces.size()-1);
				std::string Y = std::to_string(prev.count);
				face.type = "yarn-to-left x -y" + Y + " x +y" + Y;
				face.vertices = {
					glm::vec3(travel_x, 0.0f, c->depth),
					glm::vec3(prev_x, 0.0f, prev_depth),
					glm::vec3(prev_x, FaceHeight, prev_depth),
					glm::vec3(travel_x, FaceHeight, c->depth),
				};
				check_edge(face.vertices[0], face.vertices[1], &gizmo.lift);
				gizmo.top_edges.emplace_back(face.vertices[2], face.vertices[3]);

				FaceEdge yarn_in = FaceEdge(gizmo.faces.back(), 1, 1, FaceEdge::FlipNo);
				FaceEdge yarn_out = FaceEdge(gizmo.faces.back(), 3, 1, FaceEdge::FlipYes);

				gizmo.connections.emplace_back(prev, yarn_in);
				c->edge = yarn_out;
				c->x = final_x;

				prev = yarn_out;
				prev_depth = c->depth;
				prev_x = travel_x;
			}

			//TODO: actually properly connect these things
			for (uint32_t ci = 1; ci < cs_by_depth.size(); ++ci) {
				cs_by_depth[ci]->edge = FaceEdge();
				cs_by_depth[ci]->x = final_x;
			}

		} else { assert(&bed == &back_bed || &bed == &back_sliders);
			FaceEdge prev = stitch_yarn_out;
			float prev_depth = bed.depth;
			float prev_x = start_x;

			{ //connect forward to zero depth:
				faces.emplace_back();
				Face &face = faces.back();
				gizmo.faces.emplace_back(faces.size()-1);
				std::string Y = std::to_string(prev.count);
				face.type = "yarn-to-right x +y" + Y + " x -y" + Y;
				face.vertices = {
					glm::vec3(prev_x, 0.0f, prev_depth),
					glm::vec3(travel_x, 0.0f, 0.0f),
					glm::vec3(travel_x, FaceHeight, 0.0f),
					glm::vec3(prev_x, FaceHeight, prev_depth),
				};
				check_edge(face.vertices[0], face.vertices[1], &gizmo.lift);
				gizmo.top_edges.emplace_back(face.vertices[2], face.vertices[3]);

				FaceEdge yarn_in = FaceEdge(gizmo.faces.back(), 3, 1, FaceEdge::FlipYes);
				FaceEdge yarn_out = FaceEdge(gizmo.faces.back(), 1, 1, FaceEdge::FlipNo);

				gizmo.connections.emplace_back(prev, yarn_in);

				prev = yarn_out;
				prev_depth = 0.0f;
				prev_x = travel_x;
			}

			{ //connect further forward to first carrier:
				Carrier *c = cs_by_depth.back();
				faces.emplace_back();
				Face &face = faces.back();
				gizmo.faces.emplace_back(faces.size()-1);

				std::string Y = std::to_string(prev.count);
				face.type = "yarn-to-right x +y" + Y + " x -y" + Y;
				face.vertices = {
					glm::vec3(prev_x, 0.0f, prev_depth),
					glm::vec3(travel_x, 0.0f, c->depth),
					glm::vec3(travel_x, FaceHeight, c->depth),
					glm::vec3(prev_x, FaceHeight, prev_depth),
				};
				check_edge(face.vertices[0], face.vertices[1], &gizmo.lift);
				gizmo.top_edges.emplace_back(face.vertices[2], face.vertices[3]);

				FaceEdge yarn_in = FaceEdge(gizmo.faces.back(), 3, 1, FaceEdge::FlipYes);
				FaceEdge yarn_out = FaceEdge(gizmo.faces.back(), 1, 1, FaceEdge::FlipNo);

				gizmo.connections.emplace_back(prev, yarn_in);
				c->edge = yarn_out;
				c->x = final_x;

				prev = yarn_out;
				prev_depth = c->depth;
				prev_x = travel_x;
			}

			//TODO: actually properly connect these things
			for (uint32_t ci = 0; ci + 1 < cs_by_depth.size(); ++ci) {
				cs_by_depth[ci]->edge = FaceEdge();
				cs_by_depth[ci]->x = final_x;
			}
		}

		//TODO: tails from travel_x to final_x for all carriers.


	}


	//--- driver functions ---
	void knit(Direction dir, BedColumns &bed, int32_t needle, std::vector< Carrier * > cs) {
		assert(&bed == &back_bed || &bed == &front_bed);

		Gizmo gizmo;

		FaceEdge yarn_in = bring_carriers(dir, bed, needle, cs, &gizmo);
		FaceEdge loop_in = bed[needle_index(needle)].top_edge;

		std::string L = std::to_string(loop_in.count);
		std::string Y = std::to_string(cs.size());

		FaceEdge stitch_yarn_in;
		FaceEdge stitch_yarn_out;

		gizmo.faces.emplace_back(faces.size());
		faces.emplace_back();
		Face &face = faces.back();
		{ //face type:
			std::string knit_or_purl = ((&bed == &front_bed || &bed == &front_sliders) ? "knit" : "purl");
			if (dir == Right) {
				face.type = knit_or_purl + "-to-right -l" + L + " +y" + Y + " +l" + Y + " -y" + Y;
				stitch_yarn_out = FaceEdge(gizmo.faces.back(), 1, cs.size(), FaceEdge::FlipNo);
				stitch_yarn_in = FaceEdge(gizmo.faces.back(), 3, cs.size(), FaceEdge::FlipYes);
			} else {
				face.type = knit_or_purl + "-to-left -l" + L + " -y" + Y + " +l" + Y + " +y" + Y;
				stitch_yarn_in = FaceEdge(gizmo.faces.back(), 1, cs.size(), FaceEdge::FlipNo);
				stitch_yarn_out = FaceEdge(gizmo.faces.back(), 3, cs.size(), FaceEdge::FlipYes);
			}
		}
		face.vertices = {
				glm::vec3(side_index(needle, Left), 0.0f, bed.depth),
				glm::vec3(side_index(needle, Right), 0.0f, bed.depth),
				glm::vec3(side_index(needle, Right), FaceHeight, bed.depth),
				glm::vec3(side_index(needle, Left), FaceHeight, bed.depth),
		};
		if (loop_in.is_valid()) {
			gizmo.connections.emplace_back(FaceEdge(gizmo.faces.back(), 0), loop_in);
		}
		if (yarn_in.is_valid()) {
			gizmo.connections.emplace_back(stitch_yarn_in, yarn_in);
		}

		//TODO: return carriers! (might need to add to a different gizmo, given potential overlapping travel?)
		return_carriers(dir, bed, needle, cs, stitch_yarn_out, &gizmo);
		
		//increase lift based on edge conflicts in stitch column:
		gizmo.lift = std::max(gizmo.lift, bed[needle_index(needle)].top_y);

		//DEBUG: make sure stitches are separated by a little bit:
		gizmo.lift += 0.05f;

		for (auto fi : gizmo.faces) {
			for (auto &v : faces[fi].vertices) {
				v.y += gizmo.lift;
			}
		}
		//register top edges:
		for (auto &te : gizmo.top_edges) {
			te.a.y += gizmo.lift;
			te.b.y += gizmo.lift;
			top_edges.insert(te);
		}

		//TODO: proper glue faces for yarn and loop connections, if lift is big(?)

		//register top edge with column:
		bed[needle_index(needle)].top_edge = FaceEdge(gizmo.faces.back(), 2, cs.size());
		bed[needle_index(needle)].top_y = faces[gizmo.faces.back()].vertices[2].y;

		assert(bed[needle_index(needle)].top_edge.is_valid()); //DEBUG
	}
};

int main(int argc, char **argv) {

	if (argc != 3) {
		std::cerr << "Usage:\n\t./knitout-to-smobj <in.knitout> <out.smobj>" << std::endl;
		return 1;
	}

	std::string in_knitout = argv[1];
	std::string out_smobj = argv[2];
	std::cout << "Will interpret knitout in '" << in_knitout << "' and write smobj to '" << out_smobj << "'" << std::endl;

	//parse knitout:

	std::ifstream knitout(in_knitout, std::ios::binary);

	enum {
		MagicValue,
		CommentHeaders,
		Body
	} section = MagicValue;

	std::unique_ptr< Translator > translator;

	std::string line;
	while (std::getline(knitout, line, '\n')) {
		if (section == MagicValue) {
			if (line == ";!knitout-2") {
				section = CommentHeaders;
				continue;
			} else {
				std::cerr << "ERROR: expecting magic value line, got '" << line << "'" << std::endl;
				return 1;
			}
		}
		if (section == CommentHeaders) {
			if (line.size() >= 2 && line.substr(0,2) == ";;") {
				auto idx = line.find(": ");
				if (idx == std::string::npos) {
					std::cerr << "ERROR: no ': ' in comment header '" << line << "'" << std::endl;
					return 1;
				}
				std::string header = line.substr(2,idx-2);
				std::string value = line.substr(idx+2);
				if (header == "Carriers") {
					std::vector< std::string > tokens;
					std::string tok;
					std::istringstream str(value);
					while (str >> tok) {
						tokens.emplace_back(tok);
					}
					translator.reset(new Translator(tokens));
				} else {
					std::cerr << "NOTE: ignoring header '" << line << "'" << std::endl;
				}
				continue;
			} else {
				if (!translator) {
					std::cerr << "ERROR: ';;Carriers: ...' header is required." << std::endl;
					return 1;
				}
				section = Body;
			}
		}
		assert(section == Body);
		assert(translator);

		{ //remove comments:
			auto idx = line.find(';');
			if (idx != std::string::npos) {
				if (idx + 1 < line.size() && line[idx+1] == ';') {
					std::cerr << "WARNING: comment-header-like comment outside header section '" << line << "'" << std::endl;
				}
				line = line.substr(0,idx);
			}
		}

		//split into tokens:
		std::vector< std::string > tokens;
		{
			std::string tok;
			std::istringstream str(line);
			while (str >> tok) {
				tokens.emplace_back(tok);
			}
		}

		auto parse_direction = [](std::string const &tok, Direction *dir) {
			if (tok == "+") {
				*dir = Right;
			} else if (tok == "-") {
				*dir = Left;
			} else {
				throw std::runtime_error("invalid direction '" + tok + "'");
			}
		};

		auto parse_bedneedle = [](std::string const &tok, Bed *bed, int32_t *needle) {
			std::string number;
			if        (tok.size() >= 2 && tok.substr(0,2) == "fs") {
				*bed = FrontSliders;
				number = tok.substr(2);
			} else if (tok.size() >= 2 && tok.substr(0,2) == "bs") {
				*bed = BackSliders;
				number = tok.substr(2);
			} else if (tok.size() >= 1 && tok.substr(0,1) == "b") {
				*bed = Back;
				number = tok.substr(1);
			} else if (tok.size() >= 1 && tok.substr(0,1) == "f") {
				*bed = Front;
				number = tok.substr(1);
			} else {
				throw std::runtime_error("invalid bed '" + tok + "'");
			}
			*needle = std::stoi(number);
		};

		if (tokens.size() == 0) {
			//ignore empty line
		} else if (tokens[0] == "knit") {
			//knit D N CS
			if (tokens.size() < 3) throw std::runtime_error("knit must have at least two parameters.");
			Direction dir;
			parse_direction(tokens[1], &dir);
			Bed bed;
			int32_t needle;
			parse_bedneedle(tokens[2], &bed, &needle);
			std::vector< std::string > carriers(tokens.begin() + 3, tokens.end());

			translator->knit(
				dir,
				translator->lookup_bed(bed),
				needle,
				translator->lookup_carriers(carriers)
			);

		} else if (tokens[0] == "tuck") {
		} else if (tokens[0] == "split") {
		} else if (tokens[0] == "split") {
		} else {
			std::cout << "WARNING: unsupported operation '" << line << "'" << std::endl;
		}
	}

	if (!translator) {
		std::cerr << "ERROR: ';;Carriers: ...' header is required." << std::endl;
		return 1;
	}

	//write smobj file:
	std::cout << "Made " << translator->faces.size() << " faces." << std::endl;

	std::ofstream out(out_smobj, std::ios::binary);

	//face type library:
	std::map< std::string, uint32_t > type_index;
	for (auto const &f : translator->faces) {
		if (type_index.insert(std::make_pair(f.type, type_index.size()+1)).second) {
			out << "L " << f.type << "\n";
		}
	}

	//faces:
	uint32_t vertex_index = 0;
	for (auto const &f : translator->faces) {
		std::string f_line = "f";
		for (auto const &v : f.vertices) {
			out << "v " << v.x << " " << v.y << " " << v.z << "\n";
			f_line += " " + std::to_string(++vertex_index);
		}
		out << f_line << "\n";
		out << "T " << type_index[f.type] << "\n"; //DEBUG: << " # " << f.type << "\n";
	}

	return 0;
}
