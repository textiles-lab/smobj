#ifndef GLM_ENABLE_EXPERIMENTAL
#define GLM_ENABLE_EXPERIMENTAL
#endif
#include <glm/glm.hpp>
#include <glm/gtx/hash.hpp>

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
#include <functional>
#include <list>

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
// needle     next needle
//     v         v
// -1  0  1 2 3  4
//  ^-----^---^---- where carriers travel vertically
//     ^----^----^-- where carriers travel front-to-back (4*n+2 are for 1/4 pitch racking)

int32_t needle_index(int32_t needle) {
	return needle * 4;
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
	Face(uint32_t source_) : source(source_) { }
	std::string type; //type edge ... edge
	std::vector< glm::vec3 > vertices;
	uint32_t source = 0; //1-based source line number
};

struct Connection {
	Connection(FaceEdge const &a_, FaceEdge const &b_, std::string const &why_) : a(a_), b(b_), why(why_) {
		//DEBUG:
		if (a.count != b.count) {
			std::cerr << "Count mis-match " << a.count << " vs " << b.count << " for \"" << why << "\"" << std::endl;
			assert(a.count == b.count);
		}
	}
	FaceEdge a,b;
	std::string why; //DEBUG
};

struct YarnVertical {
	YarnVertical(int32_t index_, FaceEdge const &a_, Direction da_, FaceEdge const &b_, Direction db_) : index(index_), a(a_), b(b_), da(da_), db(db_) { }
	int32_t index;
	FaceEdge a,b;
	Direction da, db;
};

struct LoopVertical {
	LoopVertical(int32_t index_, FaceEdge const &a_, FaceEdge const &b_) : index(index_), a(a_), b(b_) { }
	int32_t index;
	FaceEdge a,b;
};

struct Gizmo {
	std::vector< uint32_t > faces;
	std::vector< Connection > connections;
	std::vector< YarnVertical > yarn_verticals;
	std::vector< LoopVertical > loop_verticals;
	float lift = 0.0f;
	std::vector< std::function< void(float) > > set_lift; //called in order after lift is finalized
};

struct Column {
	FaceEdge top_edge; //by convention, oriented rightwards
	float top_y = 0;
};

struct BedColumns : public std::unordered_map< int32_t, Column > {
	float depth = 0.0f;
};

struct Horizon {
	//returns maximum value in range (begin,end):
	float get_value(int32_t begin, int32_t end) const {
		float value = -std::numeric_limits< float >::infinity();
		for (int32_t i = begin; i < end; ++i) {
			auto f = values.find(i);
			if (f != values.end()) {
				value = std::max(value, f->second);
			}
		}
		return value;
	}
	//sets range (begin,end) to a value:
	void raise_value(int32_t begin, int32_t end, float value) {
		for (int32_t i = begin; i < end; ++i) {
			auto ret = values.insert(std::make_pair(i, value));
			assert(ret.first->second <= value);
			ret.first->second = value;
		}
	}
	//simple/slow version: track as list of values at points:
	std::unordered_map< int32_t, float > values;
	//TODO: could be tracked as a piecewise-constant list of open intervals at depths:
	//std::map< int32_t, float > endpoints;
};


struct Unit {
	Unit(std::string const &_name = "", float _length = 1.0f) : name(_name), length(_length) { }
	std::string name;
	float length;
};

struct Checkpoint {
	Checkpoint(uint32_t face_ = -1U, uint32_t edge_ = -1U, uint32_t crossing_ = -1U, float _length = 0.0f, Unit const *_unit = nullptr) : face(face_), edge(edge_), crossing(crossing_), length(_length), unit(_unit) { }
	uint32_t face;
	uint32_t edge;
	uint32_t crossing;

	float length;
	Unit const *unit;

	bool is_valid() const {
		return (face != -1U && edge != -1U && crossing != -1U);
	}
};

struct Carrier {
	std::string name; //for debugging/error reporting purposes.

	//information about carrier track:
	float depth = 0.0f; //depth (z) for track
	Horizon horizon; //height of tallest face along track at every point.

	//track is divided into "travel zones" near stitches and "lift zones" between them
	//carriers are always parked in lift zones.

	//info about current parked position:
	bool ready = false; //has 'in' been called? I.e. is the carrier ready to be in?

	FaceEdge parked_edge; //is a valid edge if carrier is in, otherwise is not
	//by convention, edge points upward
	int32_t parked_index = 0;
	Direction parked_direction = Right; //was it moving left or right when parked?
	//height of parked position is not tracked (or queried?)

	//last stitch is used to kick carriers as needed:
	BedColumns const *last_bed = nullptr;
	int32_t last_needle = 0;
	Direction last_dir = Right;

	//used when adding lengths:
	Checkpoint last_checkpoint;

	//the 'anchor' is the connected stitch:
	//(this gets a bit fuzzy when splits happen, but estimating inter-stitch distances is generally going to be a but fuzzy)
	BedColumns const *anchor_bed = nullptr;
	int32_t anchor_needle = 0;
	Direction anchor_side = Right;
	std::map< Unit const *, float > anchor_depth;
};

struct Crossings {
	//crossing indicies work with 3n being loop crossings and 3n+/-1 being yarn crossings
	float check_crossing(int32_t front, int32_t back) const {
		float value = -std::numeric_limits< float >::infinity();
		//TODO: could use range lookups to do better, or use height-sorting to do better:
		for (auto const &c : crossings) {
			if ( (c.first.first == front || c.first.second == back)
			  || (c.first.first <= front && c.first.second >= back)
			  || (c.first.first >= front && c.first.second <= back) ) {
				value = std::max(value, c.second);
			}
		}
		return value;
	}
	void add_crossing(int32_t front, int32_t back, float value) {
		auto ret = crossings.insert(std::make_pair(std::make_pair(front, back), value));
		assert(ret.first->second <= value);
		ret.first->second = value;
	}
	std::map< std::pair< int32_t, int32_t >, float > crossings;
};

//positions:

//        v--------v--------v-v----------- where yarns/loops travel between stitches
// travel | stitch | travel | | travel |
//   ^---------^--------^----------^------ where yarns/loops travel vertically

constexpr float StitchWidth = 1.0f;
constexpr float TravelWidth = 0.25f;
constexpr float GapWidth = 0.125f;

float stitch_x(int32_t index, Direction side) {
	int32_t needle = index/4;
	assert(needle * 4 == index && "stitches are always at 4n");
	return needle * (StitchWidth + TravelWidth + GapWidth + TravelWidth)
		+ (side == Left ? -0.5f * StitchWidth : 0.5f * StitchWidth);
}

float travel_x(int32_t index, Direction side) {
	int32_t needle = index/4;
	int32_t ofs = index - needle*4;
	if (ofs < 0) {
		ofs += 4;
		needle -= 1;
	}
	assert(ofs >= 0 && ofs < 4);
	assert(needle * 4 + ofs == index);

	if (ofs == 1) {
		return stitch_x(needle*4, Right) + (side == Left ? 0.0f : TravelWidth);
	} else if (ofs == 3) {
		return stitch_x((needle+1)*4, Left) + (side == Left ?-TravelWidth : 0.0f);
	} else {
		assert(0 && "Travels are always at 4n +/- 1");
	}
}

//depths:
// --- back bed --- @ -1.0
// --- back hooks --- @ -0.75
//  [various crossings here @ [-0.5 - 0.0]
// (carriers) @ [0.0 - 0.5]
// --- front hooks --- @ 0.75
// --- front bed ---   @ 1.0

constexpr const float ShearBack = -0.5f;
constexpr const float ShearFront = 0.0f;
constexpr const float FaceHeight = 1.0f;

struct Translator {
	BedColumns back_bed, back_sliders, front_sliders, front_bed;
	std::unordered_map< std::string, Carrier > carriers;
	Crossings crossings;
	float racking = 0.0f;
	uint32_t stitch_index = 5; //current stitch table setting (machine units)

	uint32_t source = 0; //current line number (applied to created faces)

	struct {
		uint32_t empty_drops = 0; //stitches that did nothing!
	} stats;


	Translator(std::vector< std::string > const &carrier_names) {
		back_bed.depth = -1.0f;
		back_sliders.depth = -0.75f;
		front_sliders.depth = 0.75f;
		front_bed.depth = 1.0f;
		for (uint32_t i = 0; i < carrier_names.size(); ++i) {
			Carrier carrier;
			carrier.name = carrier_names[i];
			//carrier_names is front-to-back order:
			carrier.depth = 0.5f * ((carrier_names.size() - 1 - i + 0.5f) / float(carrier_names.size()));
			auto ret = carriers.insert(std::make_pair(carrier.name, carrier));
			if (!ret.second) throw std::runtime_error("Carrier '" + carrier_names[i] + "' is named twice.");
		}
		units.emplace_back("1", 1.0);
	}

	//--- output ---
	std::vector< Face > faces;
	std::vector< Connection > connections;
	std::list< Unit > units;
	std::vector< Checkpoint > checkpoints;

	//--- helper functions ---

	Unit const *get_unit(std::string name, float default_length = 1.0f) {
		for (auto const &u : units) {
			if (u.name == name) return &u;
		}
		units.emplace_back(name, default_length);
		return &units.back();
	}
	//some specific units:

	//loop length of current stitch:
	Unit const *get_loop_length_unit() {
		return get_unit("l" + std::to_string(stitch_index), 1.0f);
	}
	//height of stitch for anchor_depth:
	Unit const *get_stitch_height_unit() {
		return get_unit("h" + std::to_string(stitch_index), 0.7f);
	}
	//needle units:
	// at quarter-pitch, let's say it looks like this:
	//  back:       |---0---|   |---1---|
	// front: |---0---|   |---1---|
	// units:     |a|b|a|a|b|a|
	Unit const *get_needle_width_a_unit() {
		return get_unit("na", 0.1f);
	}
	Unit const *get_needle_width_b_unit() {
		return get_unit("nb", 0.3f);
	}
	//
	Unit const *get_bed_gap_unit() {
		return get_unit("g", 0.6f);
	}

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

	enum SetupSpecial {
		SpecialNone,
		SpecialStopAtShearFront
	};

	//bring carriers together in order to make a stitch in a given direction on a given bed/needle.
	//adds both the bring-to and the come-from for the stitch to the gizmo,
	//sets yarn_to_stitch and yarn_from_stitch properly for connection to the stitch.
	void setup_carriers(Direction dir, BedColumns const &bed, int32_t needle, std::vector< Carrier * > const &cs,
		Gizmo *gizmo_,
		FaceEdge *yarn_to_stitch_, FaceEdge *yarn_from_stitch_,
		SetupSpecial special = SpecialNone) {
		assert(gizmo_);
		auto &gizmo = *gizmo_;
		assert(yarn_to_stitch_);
		auto &yarn_to_stitch = *yarn_to_stitch_;
		assert(yarn_from_stitch_);
		auto &yarn_from_stitch = *yarn_from_stitch_;

		if (special == SpecialStopAtShearFront) {
			assert(&bed == &back_bed || &bed == &back_sliders);
		}

		yarn_to_stitch = FaceEdge();
		yarn_from_stitch = FaceEdge();

		//no carriers? nothing to do:
		if (cs.empty()) return;

		for (auto c : cs) {
			if (!c->ready) {
				throw std::runtime_error("Making a stitch with carrier '" + c->name + "', which is not in.");
			}
		}

		//yarn is going to have two parts:
		// - the "bring", which gets carriers to the lift lanes at the proper side of the stitch (and moves other carriers out of the way)
		// - the "surround", which gets carriers to/from the stitch along the carrier travel lanes
		//The "bring" is going to be assembled and committed first,
		//  and the "surround" will be added to the same gizmo as the stitch.

		//where does bring need to bring the carriers?
		int32_t bring_index;
		//where does surround transit to/from carriers?
		int32_t surround_index;
		//where does carrier end up:
		int32_t leave_index;

		carrier_indices(dir, bed, needle, &bring_index, &surround_index, &leave_index);

		//First, build 'bring':
		Gizmo bring_gizmo;

		//all non-active carriers get kicked if the surround index (i.e., where the yarn is moved) is between that carrier's parked position and its logical position.
		{
			std::unordered_set< Carrier * > used(cs.begin(), cs.end());
			for (auto &cn : carriers) {
				Carrier *c = &cn.second;
				if (used.count(c)) continue;
				if (!c->last_bed) continue;
				int32_t logical_index;
				carrier_indices(c->last_dir, *c->last_bed, c->last_needle, nullptr, nullptr, &logical_index);
				/*std::cout << "Carrier " << c->name << ":" << std::endl;
				std::cout << "    parked index: " << c->parked_index << std::endl;
				std::cout << "    logical index: " << logical_index << std::endl;
				std::cout << " vs surround index: " << surround_index << std::endl;*/
				assert(c->parked_index != surround_index);
				assert(logical_index != surround_index);
				if ((c->parked_index < surround_index) == (logical_index < surround_index)) continue;
				//std::cout << "    !!KICKING!!" << std::endl;
				bring_carrier(c, logical_index, &bring_gizmo);
			}
		}

		//all carriers travel in their tracks to the target index:
		for (auto c : cs) {
			//If carrier is out, it appears at bring_index:
			if (!c->parked_edge.is_valid()) {

				c->parked_index = bring_index;
				c->parked_direction = dir;

				//build a parked face coming from the right:
				faces.emplace_back(source);
				bring_gizmo.faces.emplace_back(faces.size()-1);
				Face &face = faces.back();
				face.vertices = {
					glm::vec3(travel_x(c->parked_index, Left), 0.0f, c->depth),
					glm::vec3(travel_x(c->parked_index, Right), 0.0f, c->depth),
					glm::vec3(travel_x(c->parked_index, Right), FaceHeight, c->depth),
					glm::vec3(travel_x(c->parked_index, Left), FaceHeight, c->depth),
				};

				c->last_checkpoint.face = faces.size()-1;
				c->last_checkpoint.crossing = 0;

				if (dir == Left) {
					face.type = "yarn-to-left x -y1 x +y1";
				
					c->parked_edge = FaceEdge(faces.size()-1, 3, 1, FaceEdge::FlipYes);
					c->last_checkpoint.edge = 1;

				} else {
					face.type = "yarn-to-right x +y1 x -y1";

					c->parked_edge = FaceEdge(faces.size()-1, 1, 1, FaceEdge::FlipNo);
					c->last_checkpoint.edge = 3;
				}

				c->last_checkpoint.unit = get_unit("in");
				c->last_checkpoint.length = 1.0f;
				checkpoints.emplace_back(c->last_checkpoint);

				bring_gizmo.lift = std::max(bring_gizmo.lift, c->horizon.get_value(bring_index, bring_index+1));
				bring_gizmo.set_lift.emplace_back([c,bring_index](float lift){
					c->horizon.raise_value(bring_index, bring_index+1, lift); //<-- lift set to bottom of stitch at travel enter
				});

			}
			bring_carrier(c, bring_index, &bring_gizmo);

			//also (preemptively) update last_* info:
			c->last_dir = dir;
			c->last_bed = &bed;
			c->last_needle = needle;
		}


		//Carriers are all now parked at the correct edge, so commit the "bring" gizmo / resolve lift:

		commit(bring_gizmo);

		//Now build the surround:
		gizmo.lift = std::max(gizmo.lift, bring_gizmo.lift); //surround must be at least as high as the bring gizmo

		//Connect all carriers up to the front/back:

		std::vector< Carrier * > cs_by_depth = cs;
		std::sort(cs_by_depth.begin(), cs_by_depth.end(), [](Carrier *a, Carrier *b) {
			return a->depth < b->depth;
		});

		std::map< Carrier *, uint32_t > cs_to_index;
		for (uint32_t i = 0; i < cs.size(); ++i) {
			auto ret = cs_to_index.insert(std::make_pair(cs[i], i));
			assert(ret.second && "c can only appear once in cs");
		}

		if (&bed == &front_bed || &bed == &front_sliders) {
			//font bed:

			//check and mark for update all of the relevant surround-travel-lifts:
			for (auto &nc : carriers) {
				Carrier *c2 = &nc.second;
				//carriers' horizons own the depth *behind* them, so this is a strict inequality:
				if (c2->depth > cs_by_depth.front()->depth) {
					gizmo.lift = std::max(gizmo.lift, c2->horizon.get_value(surround_index, surround_index+1));
					gizmo.set_lift.emplace_back([c2,surround_index](float lift){
						c2->horizon.raise_value(surround_index, surround_index+1, lift + FaceHeight);
					});
				}
			}
			//NOTE: I think it's okay not to check/update heights on front_bed / front_sliders here because stitches share those values and the drive code will deal with it.

			//track previously assembled edge:
			FaceEdge live_to_surround, live_to_travel;
			float live_depth = std::numeric_limits< float >::quiet_NaN();

			//walk back-to-front:
			for (uint32_t i = 0; i < cs_by_depth.size(); ++i) {
				Carrier *c = cs_by_depth[i];
				assert(c->parked_index == bring_index);
				std::vector< Carrier * > live_above;
				std::vector< Carrier * > live_below;
				for (uint32_t j = 0; j < i; ++j) {
					if (cs_to_index[cs_by_depth[j]] < cs_to_index[cs_by_depth[i]]) {
						live_below.emplace_back(cs_by_depth[j]);
					} else {
						assert(cs_to_index[cs_by_depth[j]] > cs_to_index[cs_by_depth[i]]);
						live_above.emplace_back(cs_by_depth[j]);
					}
				}

				//build face to connect from/to vertical travel:
				FaceEdge travel_to_surround;
				FaceEdge surround_to_travel;
				{
					MCTParams params;
					params.c = c;
					params.bring_index = bring_index;
					params.leave_index = leave_index;
					make_carrier_tabs(params, &gizmo, &travel_to_surround, &surround_to_travel);
				}

				if (live_above.empty() && live_below.empty()) {
					//first/only yarn in plating:
					assert(!live_to_surround.is_valid());
					assert(!live_to_travel.is_valid());

					live_to_surround = travel_to_surround;
					live_to_travel = surround_to_travel;
					live_depth = c->depth;

				} else {
					assert(live_to_surround.count == live_above.size() + live_below.size());
					assert(live_to_travel.count == live_above.size() + live_below.size());

					PSMParams psm;

					psm.start_z = live_depth;
					psm.add_z = c->depth;
					//NOTE: when moving to "to-back" version, will need to edit this:
					psm.end_z = glm::mix(psm.add_z, (i + 1 < cs_by_depth.size() ? cs_by_depth[i+1]->depth : bed.depth), 0.5f);
					psm.split_z = glm::mix(psm.start_z, psm.add_z, 0.5f);
					psm.merge_z = glm::mix(psm.add_z, psm.end_z, 0.5f);
					psm.live_above = live_above.size();
					psm.live_below = live_below.size();

					//std::cout << "(" << psm.start_z << ", " << psm.split_z << ", " << psm.add_z << ", " << psm.merge_z << ", " << psm.end_z << ")" << std::endl; //DEBUG

					psm.x = stitch_x(surround_index, (dir == Right ? Left : Right));
					psm.inward = true;
					psm.live = live_to_surround;
					psm.travel = travel_to_surround;
					live_to_surround = plating_split_merge(psm, &gizmo);

					psm.x = stitch_x(surround_index, dir);
					psm.inward = false;
					psm.live = live_to_travel;
					psm.travel = surround_to_travel;
					live_to_travel = plating_split_merge(psm, &gizmo);

					live_depth = psm.end_z;
				}
			}

			//connect live edge to stitch:
			{
				assert(live_to_surround.is_valid());

				faces.emplace_back(source);
				gizmo.faces.emplace_back(faces.size()-1);
				Face &face = faces.back();
				std::string Y = std::to_string(live_to_surround.count);
				face.type = "yarn-to-right x +y" + Y + " x -y" + Y;
				float x = stitch_x(surround_index, (dir == Right ? Left : Right));
				face.vertices = {
					glm::vec3(x, 0.0f, live_depth),
					glm::vec3(x, 0.0f, bed.depth),
					glm::vec3(x, FaceHeight, bed.depth),
					glm::vec3(x, FaceHeight, live_depth),
				};
				FaceEdge yarn_in = FaceEdge(faces.size()-1, 3, live_to_surround.count, FaceEdge::FlipYes);
				FaceEdge yarn_out = FaceEdge(faces.size()-1, 1, live_to_surround.count, FaceEdge::FlipNo);

				gizmo.connections.emplace_back(live_to_surround, yarn_in, "A");

				//NOTE: I don't think this needs lift management

				yarn_to_stitch = yarn_out;
			}

			//connect live edge *from* stitch:
			{
				assert(live_to_travel.is_valid());

				faces.emplace_back(source);
				gizmo.faces.emplace_back(faces.size()-1);
				Face &face = faces.back();
				std::string Y = std::to_string(live_to_travel.count);
				face.type = "yarn-to-left x -y" + Y + " x +y" + Y;
				float x = stitch_x(surround_index, dir);
				face.vertices = {
					glm::vec3(x, 0.0f, live_depth),
					glm::vec3(x, 0.0f, bed.depth),
					glm::vec3(x, FaceHeight, bed.depth),
					glm::vec3(x, FaceHeight, live_depth),
				};
				FaceEdge yarn_in = FaceEdge(faces.size()-1, 1, live_to_travel.count, FaceEdge::FlipNo);
				FaceEdge yarn_out = FaceEdge(faces.size()-1, 3, live_to_travel.count, FaceEdge::FlipYes);

				gizmo.connections.emplace_back(yarn_out, live_to_travel, "B");

				//NOTE: I don't think this needs lift management

				yarn_from_stitch = yarn_in;
				if (yarn_from_stitch.count != yarn_to_stitch.count) {
					std::cerr << "NOTE: DEBUG -- not connecting yarn_from_stitch" << std::endl;
					yarn_from_stitch = FaceEdge(); //DEBUG!!
				}
			}

		} else { assert(&bed == &back_bed || &bed == &back_sliders);
			//back bed code walks backward from front bed:

			//check and mark for update all of the relevant surround-travel-lifts:
			for (auto &nc : carriers) {
				Carrier *c2 = &nc.second;
				//carriers' horizons own the depth *behind* them, so this is a non-strict inequality:
				if (c2->depth <= cs_by_depth.back()->depth) {
					gizmo.lift = std::max(gizmo.lift, c2->horizon.get_value(surround_index, surround_index+1));
					gizmo.set_lift.emplace_back([c2,surround_index](float lift){
						c2->horizon.raise_value(surround_index, surround_index+1, lift + FaceHeight);
					});
				}
			}
			//NOTE: I think it's okay not to check/update heights on back_bed / back_sliders here because stitches share those values and the drive code will deal with it.

			//track previously assembled edge:
			FaceEdge live_to_surround, live_to_travel;
			float live_depth = std::numeric_limits< float >::quiet_NaN();

			//walk front-to-back:
			for (uint32_t i = cs_by_depth.size()-1; i < cs_by_depth.size(); --i) {
				Carrier *c = cs_by_depth[i];
				assert(c->parked_index == bring_index);
				std::vector< Carrier * > live_above;
				std::vector< Carrier * > live_below;
				for (uint32_t j = i+1; j < cs_by_depth.size(); ++j) {
					if (cs_to_index[cs_by_depth[j]] < cs_to_index[cs_by_depth[i]]) {
						live_below.emplace_back(cs_by_depth[j]);
					} else {
						assert(cs_to_index[cs_by_depth[j]] > cs_to_index[cs_by_depth[i]]);
						live_above.emplace_back(cs_by_depth[j]);
					}
				}

				//build face to connect from/to vertical travel:
				FaceEdge travel_to_surround;
				FaceEdge surround_to_travel;
				{
					MCTParams params;
					params.c = c;
					params.bring_index = bring_index;
					params.leave_index = leave_index;
					make_carrier_tabs(params, &gizmo, &travel_to_surround, &surround_to_travel);
				}

				if (live_above.empty() && live_below.empty()) {
					//first/only yarn in plating:
					assert(!live_to_surround.is_valid());
					assert(!live_to_travel.is_valid());

					live_to_surround = travel_to_surround;
					live_to_travel = surround_to_travel;
					live_depth = c->depth;

				} else {
					assert(live_to_surround.count == live_above.size() + live_below.size());
					assert(live_to_travel.count == live_above.size() + live_below.size());

					PSMParams psm;

					psm.start_z = live_depth;
					psm.add_z = c->depth;
					psm.end_z = glm::mix(psm.add_z, (i > 0 ? cs_by_depth[i-1]->depth : 0.0f), 0.5f);
					psm.split_z = glm::mix(psm.start_z, psm.add_z, 0.5f);
					psm.merge_z = glm::mix(psm.add_z, psm.end_z, 0.5f);
					psm.live_above = live_above.size();
					psm.live_below = live_below.size();

					//std::cout << "(" << psm.start_z << ", " << psm.split_z << ", " << psm.add_z << ", " << psm.merge_z << ", " << psm.end_z << ")" << std::endl; //DEBUG

					psm.x = stitch_x(surround_index, (dir == Right ? Left : Right));
					psm.inward = true;
					psm.live = live_to_surround;
					psm.travel = travel_to_surround;
					live_to_surround = plating_split_merge(psm, &gizmo);

					psm.x = stitch_x(surround_index, dir);
					psm.inward = false;
					psm.live = live_to_travel;
					psm.travel = surround_to_travel;
					live_to_travel = plating_split_merge(psm, &gizmo);

					live_depth = psm.end_z;

				}
			}

			//connect live edges to shear front:
			{
				assert(live_to_surround.is_valid());

				faces.emplace_back(source);
				gizmo.faces.emplace_back(faces.size()-1);
				Face &face = faces.back();
				std::string Y = std::to_string(live_to_surround.count);
				face.type = "yarn-to-left x -y" + Y + " x +y" + Y;
				float x = travel_x(bring_index, (dir == Right ? Right : Left));
				face.vertices = {
					glm::vec3(x, 0.0f, ShearFront),
					glm::vec3(x, 0.0f, live_depth),
					glm::vec3(x, FaceHeight, live_depth),
					glm::vec3(x, FaceHeight, ShearFront),
				};
				FaceEdge yarn_in = FaceEdge(faces.size()-1, 1, live_to_surround.count, FaceEdge::FlipNo);
				FaceEdge yarn_out = FaceEdge(faces.size()-1, 3, live_to_surround.count, FaceEdge::FlipYes);

				gizmo.connections.emplace_back(live_to_surround, yarn_in, "live to shear front");

				live_to_surround = yarn_out;

				//NOTE: no lift management -- owned by area behind last yarn
			}

			{
				assert(live_to_travel.is_valid());

				faces.emplace_back(source);
				gizmo.faces.emplace_back(faces.size()-1);
				Face &face = faces.back();
				std::string Y = std::to_string(live_to_travel.count);
				face.type = "yarn-to-right x +y" + Y + " x -y" + Y;
				//float x = stitch_x(surround_index, dir);
				float x = travel_x(leave_index, (dir == Right ? Left : Right));
				face.vertices = {
					glm::vec3(x, 0.0f, ShearFront),
					glm::vec3(x, 0.0f, live_depth),
					glm::vec3(x, FaceHeight, live_depth),
					glm::vec3(x, FaceHeight, ShearFront),
				};
				FaceEdge yarn_in = FaceEdge(faces.size()-1, 3, live_to_travel.count, FaceEdge::FlipYes);
				FaceEdge yarn_out = FaceEdge(faces.size()-1, 1, live_to_travel.count, FaceEdge::FlipNo);

				gizmo.connections.emplace_back(yarn_out, live_to_travel, "shear front to live");

				//NOTE: I don't think this needs lift management

				live_to_travel = yarn_in;
			}

			if (special == SpecialStopAtShearFront) {
				yarn_to_stitch = live_to_surround;
				yarn_from_stitch = live_to_travel;
				return;
			}

			//connect shear front to shear back:
			{
				assert(live_to_surround.is_valid());

				faces.emplace_back(source);
				gizmo.faces.emplace_back(faces.size()-1);
				Face &face = faces.back();
				std::string Y = std::to_string(live_to_surround.count);
				face.type = "yarn-to-left x -y" + Y + " x +y" + Y;
				float back_x = stitch_x(needle_index(needle), (dir == Right ? Left : Right));
				//float x = stitch_x(surround_index, (dir == Right ? Left : Right));
				float x = travel_x(bring_index, (dir == Right ? Right : Left));
				face.vertices = {
					glm::vec3(back_x, 0.0f, ShearBack),
					glm::vec3(x, 0.0f, ShearFront),
					glm::vec3(x, FaceHeight, ShearFront),
					glm::vec3(back_x, FaceHeight, ShearBack),
				};
				FaceEdge yarn_in = FaceEdge(faces.size()-1, 1, live_to_surround.count, FaceEdge::FlipNo);
				FaceEdge yarn_out = FaceEdge(faces.size()-1, 3, live_to_surround.count, FaceEdge::FlipYes);

				gizmo.connections.emplace_back(live_to_surround, yarn_in, "shear front to back");

				live_to_surround = yarn_out;

				//lift management for bed crossing:
				int32_t back_bring_index = side_index(needle, (dir == Right ? Left : Right));
				gizmo.lift = std::max(gizmo.lift, crossings.check_crossing(bring_index, back_bring_index));
				gizmo.set_lift.emplace_back([this,bring_index,back_bring_index](float lift){
					crossings.add_crossing(bring_index, back_bring_index, lift + FaceHeight); //<-- set to top of stitch for crossing
				});

			}

			{

				assert(live_to_travel.is_valid());

				faces.emplace_back(source);
				gizmo.faces.emplace_back(faces.size()-1);
				Face &face = faces.back();
				std::string Y = std::to_string(live_to_travel.count);
				face.type = "yarn-to-right x +y" + Y + " x -y" + Y;
				float back_x = stitch_x(needle_index(needle), dir);
				//float x = stitch_x(surround_index, dir);
				float x = travel_x(leave_index, (dir == Right ? Left : Right));
				face.vertices = {
					glm::vec3(back_x, 0.0f, ShearBack),
					glm::vec3(x, 0.0f, ShearFront),
					glm::vec3(x, FaceHeight, ShearFront),
					glm::vec3(back_x, FaceHeight, ShearBack),
				};
				FaceEdge yarn_in = FaceEdge(faces.size()-1, 3, live_to_travel.count, FaceEdge::FlipYes);
				FaceEdge yarn_out = FaceEdge(faces.size()-1, 1, live_to_travel.count, FaceEdge::FlipNo);

				gizmo.connections.emplace_back(yarn_out, live_to_travel, "shear front to live");

				//NOTE: I don't think this needs lift management

				live_to_travel = yarn_in;

				//lift management for bed crossing:
				int32_t back_leave_index = side_index(needle, dir);
				gizmo.lift = std::max(gizmo.lift, crossings.check_crossing(leave_index, back_leave_index));
				gizmo.set_lift.emplace_back([this,leave_index,back_leave_index](float lift){
					crossings.add_crossing(leave_index, back_leave_index, lift + FaceHeight); //<-- set to top of stitch for crossing
				});

			}


			//connect live edges to stitch:
			{
				assert(live_to_surround.is_valid());

				faces.emplace_back(source);
				gizmo.faces.emplace_back(faces.size()-1);
				Face &face = faces.back();
				std::string Y = std::to_string(live_to_surround.count);
				face.type = "yarn-to-left x -y" + Y + " x +y" + Y;
				float back_x = stitch_x(needle_index(needle), (dir == Right ? Left : Right));
				face.vertices = {
					glm::vec3(back_x, 0.0f, bed.depth),
					glm::vec3(back_x, 0.0f, ShearBack),
					glm::vec3(back_x, FaceHeight, ShearBack),
					glm::vec3(back_x, FaceHeight, bed.depth),
				};
				FaceEdge yarn_in = FaceEdge(faces.size()-1, 1, live_to_surround.count, FaceEdge::FlipNo);
				FaceEdge yarn_out = FaceEdge(faces.size()-1, 3, live_to_surround.count, FaceEdge::FlipYes);

				gizmo.connections.emplace_back(live_to_surround, yarn_in, "shear back to bed");

				//NOTE: I don't think this needs lift management

				yarn_to_stitch = yarn_out;
			}

			//connect live edge *from* stitch:
			{
				assert(live_to_travel.is_valid());

				faces.emplace_back(source);
				gizmo.faces.emplace_back(faces.size()-1);
				Face &face = faces.back();
				std::string Y = std::to_string(live_to_travel.count);
				face.type = "yarn-to-right x +y" + Y + " x -y" + Y;
				float back_x = stitch_x(needle_index(needle), dir);
				face.vertices = {
					glm::vec3(back_x, 0.0f, bed.depth),
					glm::vec3(back_x, 0.0f, ShearBack),
					glm::vec3(back_x, FaceHeight, ShearBack),
					glm::vec3(back_x, FaceHeight, bed.depth),
				};
				FaceEdge yarn_in = FaceEdge(faces.size()-1, 3, live_to_travel.count, FaceEdge::FlipYes);
				FaceEdge yarn_out = FaceEdge(faces.size()-1, 1, live_to_travel.count, FaceEdge::FlipNo);

				gizmo.connections.emplace_back(yarn_out, live_to_travel, "bed to shear back");

				//NOTE: I don't think this needs lift management

				yarn_from_stitch = yarn_in;
			}

		}

		//add yarn checkpoints:
		for (auto &c : cs) {
			uint32_t ci = &c - &cs[0];
			Checkpoint before_stitch(yarn_to_stitch.face, yarn_to_stitch.edge, (yarn_to_stitch.flip == FaceEdge::FlipYes ? cs.size() - 1 - ci : ci));
			Checkpoint after_stitch(yarn_from_stitch.face, yarn_from_stitch.edge, (yarn_from_stitch.flip == FaceEdge::FlipYes ? cs.size() - 1 - ci : ci));
			if (c->last_checkpoint.is_valid()) {
				//previous stitch -> this stitch:

				//compute distance to previous attachment point:

				//start with depth accumulated into anchor:
				std::map< Unit const *, float > length = c->anchor_depth;
				auto add = [&length](float amt, Unit const *unit) {
					length.insert(std::make_pair(unit, 0.0f)).first->second += amt;
				};

				//add in the bed gap, if needed:
				if ((&bed == &back_bed || &bed == &back_sliders)
				 != (c->anchor_bed == &back_bed || c->anchor_bed == &back_sliders)) {
					add(1.0f, get_bed_gap_unit());
				}

				//estimate steps along bed using a modified index system:
				//
				//  n*6  -- needle n
				//     (a)
				//       +1  -- left side (back bed at quarter pitch)
				//     (b)
				//   +2  -- right side (n)
				//     (a)
				//       n*6+3 -- middle == quarter-pitch aligned racking
				//     (a)
				//   +4  -- left side n+1
				//     (b)
				//       +5  -- right side (back bed at quarter pitch)
				//     (a)
				// (n+1)*6 -- needle n+1
				//

				auto make_index = [this](BedColumns const &bed, int32_t needle, Direction side) {
					int32_t index = 6 * needle;
					if (&bed == &back_bed || &bed == &back_sliders) {
						int32_t r = std::floor(racking);
						if (r == racking) {
							index += 6*r;
						} else { assert(racking - r == 0.25f);
							index += 6*r+3;
						}
					}
					index += (side == Left ? -2 : 2);
					return index;
				};

				int32_t min_index = make_index(*c->anchor_bed, c->anchor_needle, c->anchor_side);
				int32_t max_index = make_index(bed, needle, (dir == Right ? Left : Right));
				if (min_index > max_index) std::swap(min_index, max_index);
				for (int32_t i = min_index; i < max_index; ++i) {
					int32_t step = i % 6;
					if (step < 0) step += 6;
					if (step == 1 || step == 4) {
						add(1.0f, get_needle_width_b_unit());
					} else {
						add(1.0f, get_needle_width_a_unit());
					}
				}

				for (auto const &ul : length) {
					c->last_checkpoint.unit = ul.first;
					c->last_checkpoint.length = ul.second;
					checkpoints.emplace_back(c->last_checkpoint);
				}
			}
			//this stitch:
			before_stitch.unit = get_loop_length_unit();
			before_stitch.length = 1.0f;
			checkpoints.emplace_back(before_stitch);

			//last checkpoint:
			after_stitch.unit = get_unit("1");
			after_stitch.length = 0.0f;
			checkpoints.emplace_back(after_stitch);

			c->last_checkpoint = after_stitch; //store this (blank) checkpoint

			//re-anchor yarn:
			c->anchor_bed = &bed;
			c->anchor_needle = needle;
			c->anchor_depth.clear();
			c->anchor_side = dir;
		}

	}

	void carrier_indices(Direction dir, BedColumns const &bed, int32_t needle,
		int32_t *bring_index, int32_t *surround_index, int32_t *leave_index) const {

		if (&bed == &front_bed || &bed == &front_sliders) {
			//front bed: just align with needle
			if (bring_index) *bring_index = side_index(needle, (dir == Right ? Left : Right));
			if (surround_index) *surround_index = needle_index(needle);
			if (leave_index) *leave_index = side_index(needle, dir);
		} else { assert(&bed == &back_bed || &bed == &back_sliders);
			//back bed: account for racking
			uint32_t r = int32_t(std::floor(racking));
			if (std::floor(racking) == racking) { //aligned racking
				if (bring_index) *bring_index = side_index(needle+r, (dir == Right ? Left : Right));
				if (surround_index) *surround_index = needle_index(needle+r);
				if (leave_index) *leave_index = side_index(needle+r, dir);
			} else { assert(std::floor(racking) + 0.25f == racking); //quarter pitch -- use 4*n+2 index
				if (dir == Right) {
					if (bring_index) *bring_index = side_index(needle + r, Right);
					if (leave_index) *leave_index = side_index(needle + r + 1, Left);
				} else { assert(dir == Left);
					if (bring_index) *bring_index = side_index(needle + r + 1, Left);
					if (leave_index) *leave_index = side_index(needle + r, Right);
				}
				if (surround_index) *surround_index = needle_index(needle+r)+2;
			}
		}
	}


	//bring a carrier to a given index.
	//NOTE: carrier *must already* be in(!!)
	//NOTE2: DOES NOT SET c->last_* (though does set c->parked_*)
	void bring_carrier(Carrier *c, int32_t bring_index, Gizmo *bring_gizmo_) {
		assert(bring_gizmo_);
		auto &bring_gizmo = *bring_gizmo_;

		assert(c->parked_edge.is_valid());

		//build bring face:
		if (c->parked_index < bring_index) {
			//need to lift + move to right:

			faces.emplace_back(source);
			bring_gizmo.faces.emplace_back(faces.size()-1);
			Face &face = faces.back();
			face.type = "yarn-to-right x +y1 x -y1";
			face.vertices = {
				glm::vec3(travel_x(c->parked_index, Left), 0.0f, c->depth),
				glm::vec3(travel_x(bring_index, Right), 0.0f, c->depth),
				glm::vec3(travel_x(bring_index, Right), FaceHeight, c->depth),
				glm::vec3(travel_x(c->parked_index, Left), FaceHeight, c->depth),
			};

			FaceEdge yarn_in = FaceEdge(faces.size()-1, 3, 1, FaceEdge::FlipYes);
			FaceEdge yarn_out = FaceEdge(faces.size()-1, 1, 1, FaceEdge::FlipNo);

			//mark for merge w/ optional yarn-*-up-* face:
			//  c->parked_edge (from c->parked_direction)
			// up to
			//  yarn_in (to Right)
			bring_gizmo.yarn_verticals.emplace_back(c->parked_index,
				c->parked_edge, c->parked_direction,
				yarn_in, Right);

			int32_t begin = c->parked_index;
			int32_t end = bring_index + 1;

			bring_gizmo.lift = std::max(bring_gizmo.lift, c->horizon.get_value(begin, end));
			bring_gizmo.set_lift.emplace_back([c,begin,end](float lift){
				c->horizon.raise_value(begin, end-1, lift + FaceHeight);
				c->horizon.raise_value(end-1, end, lift); //<-- set to bottom of stitch at travel enter
			});

			c->parked_index = bring_index;
			c->parked_edge = yarn_out;
			c->parked_direction = Right;

		} else if (c->parked_index > bring_index) {
			//need to lift + move to left:

			faces.emplace_back(source);
			bring_gizmo.faces.emplace_back(faces.size()-1);
			Face &face = faces.back();
			face.type = "yarn-to-left x -y1 x +y1";
			face.vertices = {
				glm::vec3(travel_x(bring_index, Left), 0.0f, c->depth),
				glm::vec3(travel_x(c->parked_index, Right), 0.0f, c->depth),
				glm::vec3(travel_x(c->parked_index, Right), FaceHeight, c->depth),
				glm::vec3(travel_x(bring_index, Left), FaceHeight, c->depth),
			};

			FaceEdge yarn_in = FaceEdge(faces.size()-1, 1, 1, FaceEdge::FlipNo);
			FaceEdge yarn_out = FaceEdge(faces.size()-1, 3, 1, FaceEdge::FlipYes);

			//mark for merge w/ optional yarn-*-up-* face:
			//  c->parked_edge (going c->parked_direction)
			// up to
			//  yarn_in (going Left)
			bring_gizmo.yarn_verticals.emplace_back(c->parked_index,
				c->parked_edge, c->parked_direction,
				yarn_in, Left);


			int32_t begin = bring_index;
			int32_t end = c->parked_index + 1;

			bring_gizmo.lift = std::max(bring_gizmo.lift, c->horizon.get_value(begin, end));
			bring_gizmo.set_lift.emplace_back([c,begin,end](float lift){
				c->horizon.raise_value(begin+1,end, lift + FaceHeight);
				c->horizon.raise_value(begin,begin+1, lift); //<-- set to bottom of stitch at travel enter
			});

			c->parked_index = bring_index;
			c->parked_edge = yarn_out;
			c->parked_direction = Left;
		}

	}

	//make to/from carrier tabs for a given carrier at a given index:
	struct MCTParams {
		Carrier *c; //the carrier to make tabs for -- parked_edge *will* be modified.
		int32_t bring_index;
		int32_t leave_index;
	};
	void make_carrier_tabs(MCTParams const &p, Gizmo *gizmo_, FaceEdge *travel_to_surround_, FaceEdge *surround_to_travel_) {
		assert(gizmo_);
		auto &gizmo = *gizmo_;
		assert(travel_to_surround_);
		auto &travel_to_surround = *travel_to_surround_;
		assert(surround_to_travel_);
		auto &surround_to_travel = *surround_to_travel_;

		assert(p.c);
		assert(p.bring_index != p.leave_index);

		Carrier *c = p.c;
		int32_t bring_index = p.bring_index;
		int32_t leave_index = p.leave_index;

		FaceEdge left_yarn_in, left_yarn_out, right_yarn_in, right_yarn_out;

		{ //left tab:
			faces.emplace_back(source);
			gizmo.faces.emplace_back(faces.size()-1);
			Face &face = faces.back();
			face.type = (bring_index < leave_index ? "yarn-to-right x +y1 x -y1" : "yarn-to-left x -y1 x +y1");
			int32_t left_index = std::min(bring_index, leave_index);
			face.vertices = {
				glm::vec3(travel_x(left_index, Left), 0.0f, p.c->depth),
				glm::vec3(travel_x(left_index, Right), 0.0f, p.c->depth),
				glm::vec3(travel_x(left_index, Right), FaceHeight, p.c->depth),
				glm::vec3(travel_x(left_index, Left), FaceHeight, p.c->depth),
			};
			if (bring_index < leave_index) {
				left_yarn_in = FaceEdge(faces.size()-1, 3, 1, FaceEdge::FlipYes);
				left_yarn_out = FaceEdge(faces.size()-1, 1, 1, FaceEdge::FlipNo);
			} else {
				left_yarn_in = FaceEdge(faces.size()-1, 1, 1, FaceEdge::FlipNo);
				left_yarn_out = FaceEdge(faces.size()-1, 3, 1, FaceEdge::FlipYes);
			}
		}

		
		{ //right tab:
			faces.emplace_back(source);
			gizmo.faces.emplace_back(faces.size()-1);
			Face &face = faces.back();
			face.type = (bring_index < leave_index ? "yarn-to-right x +y1 x -y1" : "yarn-to-left x -y1 x +y1");
			int32_t right_index = std::max(bring_index, leave_index);
			face.vertices = {
				glm::vec3(travel_x(right_index, Left), 0.0f, c->depth),
				glm::vec3(travel_x(right_index, Right), 0.0f, c->depth),
				glm::vec3(travel_x(right_index, Right), FaceHeight, c->depth),
				glm::vec3(travel_x(right_index, Left), FaceHeight, c->depth),
			};

			if (bring_index < leave_index) {
				right_yarn_in = FaceEdge(faces.size()-1, 3, 1, FaceEdge::FlipYes);
				right_yarn_out = FaceEdge(faces.size()-1, 1, 1, FaceEdge::FlipNo);
			} else {
				right_yarn_in = FaceEdge(faces.size()-1, 1, 1, FaceEdge::FlipNo);
				right_yarn_out = FaceEdge(faces.size()-1, 3, 1, FaceEdge::FlipYes);
			}
		}

		if (bring_index < leave_index) {
			//mark for merge w/ optional yarn-*-up-* face:
			//  c->parked_edge (from c->parked_direction)
			// up to
			//  yarn_in (to Right)
			assert(c->parked_index == bring_index); //c should be here
			gizmo.yarn_verticals.emplace_back(bring_index,
				c->parked_edge, c->parked_direction,
				left_yarn_in, Right);

			travel_to_surround = left_yarn_out;

			surround_to_travel = right_yarn_in;

			c->parked_edge = right_yarn_out;
			c->parked_direction = Right;
			c->parked_index = leave_index;
		} else {
			//mark for merge w/ optional yarn-*-up-* face:
			//  c->parked_edge (from c->parked_direction)
			// up to
			//  yarn_in (to Left)
			assert(c->parked_index == bring_index); //c should be here
			gizmo.yarn_verticals.emplace_back(bring_index,
				c->parked_edge, c->parked_direction,
				right_yarn_in, Left);
	
			travel_to_surround = right_yarn_out;

			surround_to_travel = left_yarn_in;

			c->parked_edge = left_yarn_out;
			c->parked_direction = Left;
			c->parked_index = leave_index;
		}


		gizmo.lift = std::max(gizmo.lift, c->horizon.get_value(bring_index, bring_index+1));
		gizmo.lift = std::max(gizmo.lift, c->horizon.get_value(leave_index, leave_index+1));
		gizmo.set_lift.emplace_back([c,leave_index,bring_index](float lift){
			c->horizon.raise_value(bring_index, bring_index+1, lift + FaceHeight); //<-- set to top of stitch at travel exit
			c->horizon.raise_value(leave_index, leave_index+1, lift); //<-- set to bottom of stitch at travel enter
		});
	}

	// Split works like this (viewed from the left):
	//    ,|,
	// ,|' |-|-|
	//|-|--|--'
	// ^---------- split into above/below
	//    ^------- route to above/below new edge
	//       ^-^-- two-step merge
	//^-- start_z
	//  ^-- split_z
	//     ^-- add_z (== carrier depth)
	//       ^-- merge_z
	//         ^-- end_z

	struct PSMParams {
		float x = std::numeric_limits< float >::quiet_NaN();
		bool inward = false; //is yarn heading toward surround or away from surround?

		float start_z = std::numeric_limits< float >::quiet_NaN();
		float split_z = std::numeric_limits< float >::quiet_NaN();
		float add_z = std::numeric_limits< float >::quiet_NaN();
		float merge_z = std::numeric_limits< float >::quiet_NaN();
		float end_z = std::numeric_limits< float >::quiet_NaN();

		FaceEdge live; //edge with what came before
		uint32_t live_above = -1U;
		uint32_t live_below = -1U;
		FaceEdge travel; //edge with travel
	};


	FaceEdge plating_split_merge(PSMParams const &p, Gizmo *gizmo_) {
		assert(p.x == p.x);
		assert(p.start_z == p.start_z);
		assert(p.split_z == p.split_z);
		assert(p.add_z == p.add_z);
		assert(p.merge_z == p.merge_z);
		assert(p.end_z == p.end_z);

		float const &x = p.x;
		float const &start_z = p.start_z;
		float const &split_z = p.split_z;
		float const &add_z = p.add_z;
		float const &merge_z = p.merge_z;
		float const &end_z = p.end_z;

/*
		assert(start_z < split_z);
		assert(split_z < add_z);
		assert(add_z < merge_z);
		assert(merge_z < end_z);
*/

		assert(p.live.count == p.live_above + p.live_below);

		assert(p.travel.count); //no much point in doing otherwise.

		assert(gizmo_);
		auto &gizmo = *gizmo_;

		std::string travel = (p.inward ? "yarn-to-right" : "yarn-to-left");
		std::string split = (p.inward ? "yarn-unplate-to-right" : "yarn-plate-to-left");
		std::string merge = (p.inward ? "yarn-plate-to-right" : "yarn-unplate-to-left");
		std::string in_sign = (p.inward ? "-" : "+");
		std::string out_sign = (p.inward ? "+" : "-");

		FaceEdge top_edge, bottom_edge; //<-- used to track live edges
		{ //split face
			faces.emplace_back(source);
			gizmo.faces.emplace_back(faces.size()-1);
			Face &face = faces.back();
			FaceEdge yarn_in;
			if (p.live_above && p.live_below) {
				std::string Y = std::to_string(p.live_above + p.live_below);
				std::string YT = std::to_string(p.live_above);
				std::string YB = std::to_string(p.live_below);
				face.type = split + " x " + out_sign + "y" + YB + " " + out_sign + "y" + YT + " x " + in_sign + "y" + Y;
				face.vertices = {
					glm::vec3(x, 0.0f, start_z),
					glm::vec3(x, 0.0f, split_z),
					glm::vec3(x, 0.5f * FaceHeight, split_z),
					glm::vec3(x, FaceHeight, split_z),
					glm::vec3(x, FaceHeight, start_z),
				};
				yarn_in = FaceEdge(faces.size()-1, 4, p.live_above + p.live_below, FaceEdge::FlipYes);
				bottom_edge = FaceEdge(faces.size()-1, 1, p.live_below, FaceEdge::FlipNo);
				top_edge = FaceEdge(faces.size()-1, 2, p.live_above, FaceEdge::FlipNo);
			} else if (p.live_above) {
				std::string Y = std::to_string(p.live_above);
				face.type = travel + " x " + out_sign + "y" + Y + " x " + in_sign + "y" + Y;
				face.vertices = {
					glm::vec3(x, 0.0f, start_z),
					glm::vec3(x, 0.5f * FaceHeight, split_z),
					glm::vec3(x, FaceHeight, split_z),
					glm::vec3(x, FaceHeight, start_z),
				};
				yarn_in = FaceEdge(faces.size()-1, 3, p.live_above, FaceEdge::FlipYes);
				bottom_edge = FaceEdge();
				top_edge = FaceEdge(faces.size()-1, 1, p.live_above, FaceEdge::FlipNo);
			} else { assert(p.live_below);
				std::string Y = std::to_string(p.live_below);
				face.type = travel + " x " + out_sign + "y" + Y + " x " + in_sign + "y" + Y;
				face.vertices = {
					glm::vec3(x, 0.0f, start_z),
					glm::vec3(x, 0.0f, split_z),
					glm::vec3(x, 0.5f * FaceHeight, split_z),
					glm::vec3(x, FaceHeight, start_z),
				};
				yarn_in = FaceEdge(faces.size()-1, 3, p.live_below, FaceEdge::FlipYes);
				bottom_edge = FaceEdge(faces.size()-1, 1, p.live_below, FaceEdge::FlipNo);
				top_edge = FaceEdge();
			}

			gizmo.connections.emplace_back(p.live, yarn_in, "plate split");
		}

		//routing faces:
		if (top_edge.count) { //top travel face
			faces.emplace_back(source);
			gizmo.faces.emplace_back(faces.size()-1);
			Face &face = faces.back();
			std::string Y = std::to_string(top_edge.count);
			face.type = travel + " x " + out_sign + "y" + Y + " x " + in_sign + "y" + Y;
			face.vertices = {
				glm::vec3(x, 0.5f*FaceHeight, split_z),
				glm::vec3(x, (2.0f/3.0f)*FaceHeight, add_z),
				glm::vec3(x, 1.0f*FaceHeight, add_z),
				glm::vec3(x, 1.0f*FaceHeight, split_z),
			};
			//in:
			FaceEdge yarn_in = FaceEdge(faces.size()-1, 3, top_edge.count, FaceEdge::FlipYes);
			gizmo.connections.emplace_back(top_edge, yarn_in, "split top");

			//out:
			FaceEdge yarn_out = FaceEdge(faces.size()-1, 1, top_edge.count, FaceEdge::FlipNo);
			top_edge = yarn_out;
		}
		if (bottom_edge.count) { //bottom travel face
			faces.emplace_back(source);
			gizmo.faces.emplace_back(faces.size()-1);
			Face &face = faces.back();
			std::string Y = std::to_string(bottom_edge.count);
			face.type = travel + " x " + out_sign + "y" + Y + " x " + in_sign + "y" + Y;
			face.vertices = {
				glm::vec3(x, 0.0f*FaceHeight, split_z),
				glm::vec3(x, 0.0f*FaceHeight, add_z),
				glm::vec3(x, (1.0f/3.0f)*FaceHeight, add_z),
				glm::vec3(x, 0.5f*FaceHeight, split_z),
			};
			//in:
			FaceEdge yarn_in = FaceEdge(faces.size()-1, 3, bottom_edge.count, FaceEdge::FlipYes);
			gizmo.connections.emplace_back(bottom_edge, yarn_in, "split bottom");

			//out:
			FaceEdge yarn_out = FaceEdge(faces.size()-1, 1, bottom_edge.count, FaceEdge::FlipNo);
			bottom_edge = yarn_out;
		}

		{ //fix middle edge to be shorter:
			auto get_vertex = [this](FaceEdge const &fe, uint8_t vertex) -> glm::vec3 & {
				if ((fe.flip == FaceEdge::FlipNo) == (vertex == 0)) {
					return faces[fe.face].vertices[fe.edge];
				} else {
					return faces[fe.face].vertices[(fe.edge+1)%faces[fe.face].vertices.size()];
				}
			};
			glm::vec3 &bottom = get_vertex(p.travel, 0);
			glm::vec3 &top = get_vertex(p.travel, 1);
			assert(bottom.y == 0.0f);
			bottom.y = (1.0f/3.0f)*FaceHeight;
			assert(top.y == FaceHeight);
			top.y = (2.0f/3.0f)*FaceHeight;
		}

		{ //merge: top step
			faces.emplace_back(source);
			gizmo.faces.emplace_back(faces.size()-1);
			Face &face = faces.back();

			FaceEdge yarn_bottom_in, yarn_out;
			if (top_edge.count) {
				std::string YT = std::to_string(top_edge.count);
				std::string YB = std::to_string(p.travel.count);
				std::string Y = std::to_string(top_edge.count+p.travel.count);
				face.type = merge + " x " + out_sign + "y" + Y + " x " + in_sign + "y" + YT + " " + in_sign + "y" + YB;
				face.vertices = {
					glm::vec3(x, (1.0f/3.0f)*FaceHeight, add_z),
					glm::vec3(x, 0.5f*FaceHeight, merge_z),
					glm::vec3(x, 1.0f*FaceHeight, merge_z),
					glm::vec3(x, 1.0f*FaceHeight, add_z),
					glm::vec3(x, (2.0f/3.0f)*FaceHeight, add_z),
				};
				//in:
				yarn_bottom_in = FaceEdge(faces.size()-1, 4, p.travel.count, FaceEdge::FlipYes);
				FaceEdge yarn_top_in = FaceEdge(faces.size()-1, 3, top_edge.count, FaceEdge::FlipYes);
				gizmo.connections.emplace_back(top_edge, yarn_top_in, "plate merge top");

				yarn_out = FaceEdge(faces.size()-1, 1, p.travel.count+top_edge.count, FaceEdge::FlipNo);
			} else {
				std::string Y = std::to_string(p.travel.count);
				face.type = travel + " x " + out_sign + "y" + Y + " x " + in_sign + "y" + Y;
				face.vertices = {
					glm::vec3(x, (1.0f/3.0f)*FaceHeight, add_z),
					glm::vec3(x, 0.5f*FaceHeight, merge_z),
					glm::vec3(x, 1.0f*FaceHeight, merge_z),
					glm::vec3(x, (2.0f/3.0f)*FaceHeight, add_z),
				};
				//in:
				yarn_bottom_in = FaceEdge(faces.size()-1, 3, p.travel.count, FaceEdge::FlipYes);
				yarn_out = FaceEdge(faces.size()-1, 1, p.travel.count, FaceEdge::FlipNo);
			}
			gizmo.connections.emplace_back(p.travel, yarn_bottom_in, "plate merge main");

			//out:
			top_edge = yarn_out;
		}
		if (bottom_edge.count) { //merge: bottom travel face
			faces.emplace_back(source);
			gizmo.faces.emplace_back(faces.size()-1);
			Face &face = faces.back();
			std::string Y = std::to_string(bottom_edge.count);
			face.type = travel + " x " + out_sign + "y" + Y + " x " + in_sign + "y" + Y;
			face.vertices = {
				glm::vec3(x, 0.0f*FaceHeight, add_z),
				glm::vec3(x, 0.0f*FaceHeight, merge_z),
				glm::vec3(x, 0.5f*FaceHeight, merge_z),
				glm::vec3(x, (1.0f/3.0f)*FaceHeight, add_z),
			};
			//in:
			FaceEdge yarn_in = FaceEdge(faces.size()-1, 3, bottom_edge.count, FaceEdge::FlipYes);
			gizmo.connections.emplace_back(bottom_edge, yarn_in, "merge bottom travel");

			//out:
			FaceEdge yarn_out = FaceEdge(faces.size()-1, 1, bottom_edge.count, FaceEdge::FlipNo);
			bottom_edge = yarn_out;
		}

		{ //merge: second step
			faces.emplace_back(source);
			gizmo.faces.emplace_back(faces.size()-1);
			Face &face = faces.back();
			FaceEdge yarn_top_in, yarn_out;
			if (bottom_edge.count) {
				std::string YT = std::to_string(top_edge.count);
				std::string YB = std::to_string(bottom_edge.count);
				std::string Y = std::to_string(top_edge.count+bottom_edge.count);
				face.type = merge + " x " + out_sign + "y" + Y + " x " + in_sign + "y" + YT + " " + in_sign + "y" + YB;
				face.vertices = {
					glm::vec3(x, 0.0f*FaceHeight, merge_z),
					glm::vec3(x, 0.0f*FaceHeight, end_z),
					glm::vec3(x, 1.0f*FaceHeight, end_z),
					glm::vec3(x, 1.0f*FaceHeight, merge_z),
					glm::vec3(x, 0.5f*FaceHeight, merge_z),
				};
				FaceEdge yarn_bottom_in = FaceEdge(faces.size()-1, 4, bottom_edge.count, FaceEdge::FlipYes);
				yarn_top_in = FaceEdge(faces.size()-1, 3, top_edge.count, FaceEdge::FlipYes);
				yarn_out = FaceEdge(faces.size()-1, 1, bottom_edge.count+top_edge.count, FaceEdge::FlipNo);

				gizmo.connections.emplace_back(bottom_edge, yarn_bottom_in, "plate merge bottom");
			} else {
				std::string Y = std::to_string(top_edge.count);
				face.type = travel + " x " + out_sign + "y" + Y + " x " + in_sign + "y" + Y;
				face.vertices = {
					glm::vec3(x, 0.5f*FaceHeight, merge_z),
					glm::vec3(x, 0.0f*FaceHeight, end_z),
					glm::vec3(x, 1.0f*FaceHeight, end_z),
					glm::vec3(x, 1.0f*FaceHeight, merge_z),
				};
				yarn_top_in = FaceEdge(faces.size()-1, 3, top_edge.count, FaceEdge::FlipYes);
				yarn_out = FaceEdge(faces.size()-1, 1, top_edge.count, FaceEdge::FlipNo);
			}

			assert(top_edge.count >= 1);
			gizmo.connections.emplace_back(top_edge, yarn_top_in, "plate merge top2");

			return yarn_out;
		}

	}

	//apply gizmo lift to all faces:
	void commit(Gizmo &gizmo) {
		//apply lift value:
		for (auto fi : gizmo.faces) {
			for (auto &v : faces[fi].vertices) {
				v.y += gizmo.lift;
			}
		}

		//build loop verticals:
		for (auto const &lv : gizmo.loop_verticals) {
			glm::vec3 top_left;
			glm::vec3 top_right;
			glm::vec3 bottom_left;
			glm::vec3 bottom_right;

			
			{ //low edge:
				FaceEdge const &edge = lv.a;
				assert(edge.is_valid());

				Face &face = faces[edge.face];

				bottom_left = face.vertices[(edge.edge + (edge.flip == FaceEdge::FlipNo ? 0 : 1)) % face.vertices.size()];
				bottom_right = face.vertices[(edge.edge + (edge.flip == FaceEdge::FlipNo ? 1 : 0)) % face.vertices.size()];
				assert(bottom_left.x < bottom_right.x);
			}
			{ //high edge:
				FaceEdge const &edge = lv.b;
				assert(edge.is_valid());

				Face &face = faces[edge.face];

				top_left = face.vertices[(edge.edge + (edge.flip == FaceEdge::FlipNo ? 0 : 1)) % face.vertices.size()];
				top_right = face.vertices[(edge.edge + (edge.flip == FaceEdge::FlipNo ? 1 : 0)) % face.vertices.size()];
				assert(top_left.x < top_right.x);
			}

			if (top_left == bottom_left && top_right == bottom_right) {
				gizmo.connections.emplace_back(lv.a, lv.b, "merged loop");
			} else {
				faces.emplace_back(source);
				Face &face = faces.back();
				assert(lv.a.count == lv.b.count);
				std::string L = std::to_string(lv.a.count);
				face.type = "loop -l" + L + " x +l" + L + " x";
				face.vertices = {
					bottom_left,
					bottom_right,
					top_right,
					top_left
				};
				//DEBUG:
				//std::cout << face.type << std::endl;
				//for (auto v : face.vertices) {
				//	std::cout << v.x << " " << v.y << " " << v.z << std::endl; //DEBUG
				//}

				FaceEdge loop_in(faces.size()-1, 0, lv.a.count, FaceEdge::FlipNo);
				FaceEdge loop_out(faces.size()-1, 2, lv.a.count, FaceEdge::FlipYes);

				gizmo.connections.emplace_back(lv.a, loop_in, "lv in");
				gizmo.connections.emplace_back(loop_out, lv.b, "lv out");
			}


		}


		//build yarn verticals:
		for (auto const &yv : gizmo.yarn_verticals) {
			float left = travel_x(yv.index, Left);
			float right = travel_x(yv.index, Right);

			std::string type = "yarn-";

			//oriented bottom-to-top:
			FaceEdge to_exit;
			FaceEdge from_entrance;

			glm::vec3 top_left;
			glm::vec3 top_right;
			glm::vec3 bottom_left;
			glm::vec3 bottom_right;

			{ //mitre low edge:
				FaceEdge const &edge = yv.a;
				Direction const &dir = yv.da;
				assert(edge.is_valid());

				Face &face = faces[edge.face];

				glm::vec3 &low = face.vertices[(edge.edge + (edge.flip == FaceEdge::FlipNo ? 0 : 1)) % face.vertices.size()];
				glm::vec3 &high = face.vertices[(edge.edge + (edge.flip == FaceEdge::FlipNo ? 1 : 0)) % face.vertices.size()];
				assert(low.y < high.y);
				if (dir == Right) {
					assert(low.x == right);
					assert(high.x == right);
					high.x = glm::mix(left, right, 0.125f);
					low.x = glm::mix(left, right, 0.875f);
					type += "right";
					from_entrance = edge;
					bottom_left = high;
					bottom_right = low;
					assert(bottom_left.x < bottom_right.x);
				} else {
					assert(low.x == left);
					assert(high.x == left);
					high.x = glm::mix(left, right, 0.875f);
					low.x = glm::mix(left, right, 0.125f);
					type += "left";
					from_entrance = edge;
					bottom_left = low;
					bottom_right = high;
					assert(bottom_left.x < bottom_right.x);
				}
			}
			type += "-up-";

			{ //mitre high edge:
				FaceEdge const &edge = yv.b;
				Direction const &dir = yv.db;
				assert(edge.is_valid());

				Face &face = faces[edge.face];

				glm::vec3 &low = face.vertices[(edge.edge + (edge.flip == FaceEdge::FlipNo ? 0 : 1)) % face.vertices.size()];
				glm::vec3 &high = face.vertices[(edge.edge + (edge.flip == FaceEdge::FlipNo ? 1 : 0)) % face.vertices.size()];
				assert(low.y < high.y);
				if (dir == Right) {
					assert(low.x == left);
					assert(high.x == left);
					high.x = glm::mix(left, right, 0.125f);
					low.x = glm::mix(left, right, 0.875f);
					type += "right";
					to_exit = edge;
					top_left = high;
					top_right = low;
					assert(top_left.x < top_right.x);
				} else {
					assert(low.x == right);
					assert(high.x == right);
					high.x = glm::mix(left, right, 0.875f);
					low.x = glm::mix(left, right, 0.125f);
					type += "left";
					to_exit = edge;
					top_left = low;
					top_right = high;
					assert(top_left.x < top_right.x);
				}
			}
			type += " -y1 x +y1 x";

			if (top_left == bottom_left && top_right == bottom_right) {
				//std::cout << "NOTE: merged case." << std::endl;
				gizmo.connections.emplace_back(from_entrance, to_exit, "merged");
			} else {
				faces.emplace_back(source);
				Face &face = faces.back();
				face.type = type;
				face.vertices = {
					bottom_left,
					bottom_right,
					top_right,
					top_left
				};
				//DEBUG:
				//std::cout << face.type << std::endl;
				//for (auto v : face.vertices) {
				//	std::cout << v.x << " " << v.y << " " << v.z << std::endl; //DEBUG
				//}

				FaceEdge yarn_in(faces.size()-1, 0, 1, (yv.da == Right ? FaceEdge::FlipYes : FaceEdge::FlipNo));
				if (yarn_in.flip == FaceEdge::FlipNo) {
					assert(face.vertices[yarn_in.edge].y < face.vertices[(yarn_in.edge+1)%4].y);
				} else {
					assert(face.vertices[yarn_in.edge].y > face.vertices[(yarn_in.edge+1)%4].y);
				}
				FaceEdge yarn_out(faces.size()-1, 2, 1, (yv.db == Right ? FaceEdge::FlipNo : FaceEdge::FlipYes));
				if (yarn_out.flip == FaceEdge::FlipNo) {
					assert(face.vertices[yarn_out.edge].y < face.vertices[(yarn_out.edge+1)%4].y);
				} else {
					assert(face.vertices[yarn_out.edge].y > face.vertices[(yarn_out.edge+1)%4].y);
				}


				gizmo.connections.emplace_back(from_entrance, yarn_in, "yv in");
				gizmo.connections.emplace_back(yarn_out, to_exit, "yv out");
			}
		}

		//add connections: (NOTE: why are these in gizmo anyway?)
		connections.insert(connections.end(), gizmo.connections.begin(), gizmo.connections.end());

		//report on lift values to update horizons:
		float lift = gizmo.lift;
		for (auto const &fn : gizmo.set_lift) {
			fn(lift);
		}
	}

	//--- driver functions ---
	void in(std::vector< Carrier * > cs) {
		for (auto c : cs) {
			if (c->ready) {
				throw std::runtime_error("Calling in/inhook on carrier " + c->name + " which is already in.");
			}
			c->ready = true;
		}
	}

	void out(std::vector< Carrier * > cs) {
		for (auto c : cs) {
			if (!c->ready) {
				throw std::runtime_error("Calling out/outhook on carrier " + c->name + " which is not in.");
			}

			//---
			//build checkpoints for end of yarn:
			if (c->parked_edge.is_valid()) {
				assert(c->parked_edge.count == 1);

				c->last_checkpoint.unit = get_unit("out");
				c->last_checkpoint.length = 1.0f;
				checkpoints.emplace_back(c->last_checkpoint);

				//TODO: could use anchor depth here, potentially

				Checkpoint end;
				end.face = c->parked_edge.face;
				end.edge = c->parked_edge.edge;
				end.crossing = 0;
				end.unit = get_unit("1");
				end.length = 0.0f;
				checkpoints.emplace_back(end);
			}
			//---
			//reset carrier data:

			c->ready = false;
			c->parked_edge = FaceEdge();
			c->parked_index = 0;
			c->parked_direction = Right;

			c->last_dir = Right;
			c->last_needle = 0;
			c->last_bed = nullptr;

			c->last_checkpoint = Checkpoint();

			c->anchor_bed = nullptr;
			c->anchor_needle = 0;
			c->anchor_side = Right;
			c->anchor_depth.clear();
		}
	}


	void tuck(Direction dir, BedColumns &bed, int32_t needle, std::vector< Carrier * > cs) {
		assert(&bed == &back_bed || &bed == &front_bed);

		//a-miss doesn't show up in output:
		if (cs.empty()) return;

		Gizmo gizmo;

		//set up yarn to/from stitch:
		FaceEdge yarn_to_stitch, yarn_from_stitch;
		setup_carriers(dir, bed, needle, cs, &gizmo, &yarn_to_stitch, &yarn_from_stitch);

		FaceEdge loop_in = bed[needle_index(needle)].top_edge;


		{ //build stitch face:
			FaceEdge stitch_yarn_in, stitch_yarn_out;

			faces.emplace_back(source);
			gizmo.faces.emplace_back(faces.size()-1);
			Face &face = faces.back();
			{ //face type:
				std::string L = std::to_string(loop_in.count);
				std::string Y = std::to_string(yarn_to_stitch.count);
				std::string S = std::to_string(yarn_to_stitch.count + loop_in.count);

				std::string behind_or_infront = (&bed == &front_bed ? "behind" : "infront");
				if (dir == Right) {
					face.type = "tuck-" + behind_or_infront + "-to-right -l" + L + " +y" + Y + " +l" + S + " -y" + Y;
					stitch_yarn_out = FaceEdge(gizmo.faces.back(), 1, cs.size(), FaceEdge::FlipNo);
					stitch_yarn_in = FaceEdge(gizmo.faces.back(), 3, cs.size(), FaceEdge::FlipYes);
				} else {
					face.type = "tuck-" + behind_or_infront + "-to-left -l" + L + " -y" + Y + " +l" + S + " +y" + Y;
					stitch_yarn_in = FaceEdge(gizmo.faces.back(), 1, cs.size(), FaceEdge::FlipNo);
					stitch_yarn_out = FaceEdge(gizmo.faces.back(), 3, cs.size(), FaceEdge::FlipYes);
				}
			}
			int32_t stitch_index = needle_index(needle);
			face.vertices = {
					glm::vec3(stitch_x(stitch_index, Left), 0.0f, bed.depth),
					glm::vec3(stitch_x(stitch_index, Right), 0.0f, bed.depth),
					glm::vec3(stitch_x(stitch_index, Right), FaceHeight, bed.depth),
					glm::vec3(stitch_x(stitch_index, Left), FaceHeight, bed.depth),
			};
			
			if (loop_in.is_valid()) {
				gizmo.loop_verticals.emplace_back(stitch_index, loop_in, FaceEdge(gizmo.faces.back(), 0, loop_in.count));

				//mark for loop connection:
				// loop_in
				// up to
				// FaceEdge(gizmo.faces.back(), 0)
			}

			if (yarn_to_stitch.is_valid()) {
				gizmo.connections.emplace_back(yarn_to_stitch, stitch_yarn_in, "tuck in");
			}
			if (yarn_from_stitch.is_valid()) {
				gizmo.connections.emplace_back(stitch_yarn_out, yarn_from_stitch, "tuck out");
			}

			//increase lift based on edge conflicts in stitch column:
			gizmo.lift = std::max(gizmo.lift, bed[needle_index(needle)].top_y);

			//apply lift and create merge faces:
			commit(gizmo);

			//register loop with column:
			bed[needle_index(needle)].top_edge = FaceEdge(gizmo.faces.back(), 2, yarn_to_stitch.count + loop_in.count, FaceEdge::FlipYes);
			bed[needle_index(needle)].top_y = faces[gizmo.faces.back()].vertices[2].y;
		}

	}

	void miss(Direction dir, BedColumns &bed, int32_t needle, std::vector< Carrier * > cs) {
		if (cs.empty()) {
			std::cout << "NOTE: skipping a 'miss' that has no yarns." << std::endl;
			return;
		}

		for (auto c : cs) {
			if (!c->ready) {
				throw std::runtime_error("Making a miss with carrier '" + c->name + "', which is not in.");
			}
		}

		int32_t leave_index;

		carrier_indices(dir, bed, needle, nullptr, nullptr, &leave_index);

		Gizmo gizmo;

		//move each carrier and reset its logical position:
		for (auto c : cs) {
			if (!c->parked_edge.is_valid()) {
				std::cout << "NOTE: not generating yarn for miss with '" << c->name << "' before it makes a stitch." << std::endl;
				continue;
			}

			bring_carrier(c, leave_index, &gizmo);

			c->last_dir = dir;
			c->last_bed = &bed;
			c->last_needle = needle;
		}

		commit(gizmo);
	}

	void knit(Direction dir, BedColumns &bed, int32_t needle, std::vector< Carrier * > cs) {
		assert(&bed == &back_bed || &bed == &front_bed);

		Gizmo gizmo;

		//update anchor depths (all carriers):
		for (auto &cn : carriers) {
			Carrier *c = &cn.second;
			if (!c->last_bed) continue;
			if (c->anchor_bed == &bed && c->anchor_needle == needle) {
				c->anchor_depth.insert(std::make_pair(get_stitch_height_unit(), 0.0f)).first->second += 1.0f;
			}
		}

		//set up yarn to/from stitch:
		FaceEdge yarn_to_stitch, yarn_from_stitch;
		setup_carriers(dir, bed, needle, cs, &gizmo, &yarn_to_stitch, &yarn_from_stitch);

		FaceEdge loop_in = bed[needle_index(needle)].top_edge;

		if (loop_in.count == 0 && cs.size() == 0) {
			//skip a no-op stitch.
			stats.empty_drops += 1;
			return;
		}


		{ //build stitch face:
			FaceEdge stitch_yarn_in, stitch_yarn_out;

			faces.emplace_back(source);
			gizmo.faces.emplace_back(faces.size()-1);
			Face &face = faces.back();
			{ //face type:
				std::string L = std::to_string(loop_in.count);
				std::string Y = std::to_string(cs.size());
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
			int32_t stitch_index = needle_index(needle);
			face.vertices = {
					glm::vec3(stitch_x(stitch_index, Left), 0.0f, bed.depth),
					glm::vec3(stitch_x(stitch_index, Right), 0.0f, bed.depth),
					glm::vec3(stitch_x(stitch_index, Right), FaceHeight, bed.depth),
					glm::vec3(stitch_x(stitch_index, Left), FaceHeight, bed.depth),
			};
			
			if (loop_in.is_valid()) {
				gizmo.loop_verticals.emplace_back(stitch_index, loop_in, FaceEdge(gizmo.faces.back(), 0, loop_in.count));

				//mark for loop connection:
				// loop_in
				// up to
				// FaceEdge(gizmo.faces.back(), 0)
			}

			if (yarn_to_stitch.is_valid()) {
				gizmo.connections.emplace_back(yarn_to_stitch, stitch_yarn_in, "knit in");
			}
			if (yarn_from_stitch.is_valid()) {
				gizmo.connections.emplace_back(stitch_yarn_out, yarn_from_stitch, "knit out");
			}

			//increase lift based on edge conflicts in stitch column:
			gizmo.lift = std::max(gizmo.lift, bed[needle_index(needle)].top_y);

			//apply lift and create merge faces:
			commit(gizmo);

			//register loop with column:
			bed[needle_index(needle)].top_edge = FaceEdge(gizmo.faces.back(), 2, cs.size(), FaceEdge::FlipYes);
			bed[needle_index(needle)].top_y = faces[gizmo.faces.back()].vertices[2].y;
		}
	}

	void split(Direction dir, BedColumns &from_bed, int32_t from_needle, BedColumns &to_bed, int32_t to_needle, std::vector< Carrier * > cs) {
		assert(&from_bed == &back_bed || &from_bed == &back_sliders || &from_bed == &front_sliders || &from_bed == &front_bed);
		assert(&to_bed == &back_bed || &to_bed == &back_sliders || &to_bed == &front_sliders || &to_bed == &front_bed);

		bool is_front_to_back;
		if ( (&from_bed == &back_bed || &from_bed == &back_sliders)
		  && (&to_bed == &front_bed || &to_bed == &front_sliders) ) {
			is_front_to_back = false;
		} else if ( (&from_bed == &front_bed || &from_bed == &front_sliders)
		         && (&to_bed == &back_bed || &to_bed == &back_sliders) ) {
			is_front_to_back = true;
		} else {
			throw std::runtime_error("xfer must be f <-> b");
		}
		if ( (&from_bed == &front_sliders || &from_bed == &back_sliders)
		  && (&to_bed == &front_sliders || &to_bed == &back_sliders) ) {
			throw std::runtime_error("xfer can't be fs <-> bs");
		}

		BedColumns &front_bed = (is_front_to_back ? from_bed : to_bed);
		BedColumns &back_bed = (is_front_to_back ? to_bed : from_bed);
		int32_t front_needle = (is_front_to_back ? from_needle : to_needle);
		int32_t back_needle = (is_front_to_back ? to_needle : from_needle);
		int32_t front_index = needle_index(front_needle);
		int32_t back_index = needle_index(back_needle);

		if (back_needle + racking != front_needle) {
			throw std::runtime_error("xfer with improper racking.");
		}

		Gizmo gizmo;

		//if this is a split,
		//update anchor depths (all carriers):
		if (!cs.empty()) {
			for (auto &cn : carriers) {
				Carrier *c = &cn.second;
				if (!c->last_bed) continue;
				if (c->anchor_bed == &from_bed && c->anchor_needle == from_needle) {
					c->anchor_depth.insert(std::make_pair(get_stitch_height_unit(), 0.0f)).first->second += 1.0f;
				}
			}
		}

		uint32_t from_front_count = front_bed[needle_index(front_needle)].top_edge.count;
		uint32_t from_back_count = back_bed[needle_index(back_needle)].top_edge.count;

		//set up yarn to/from stitch @ shear plane:
		FaceEdge yarn_to_stitch, yarn_from_stitch;
		setup_carriers(dir, back_bed, back_needle, cs, &gizmo, &yarn_to_stitch, &yarn_from_stitch, SpecialStopAtShearFront);

		{ //determine lift based on a clear lane to front and back:
			//(this is a bit heavy-handed, as it could clear based on front/back loop existence)
			float R_from_back = 0.125f * from_back_count;
			float R_from_front = 0.125f * from_front_count; //lift a bit more to avoid yarns overlapping the top of other travelling loops
			gizmo.lift = std::max(gizmo.lift, this->front_bed[front_index].top_y + R_from_front);
			gizmo.lift = std::max(gizmo.lift, this->front_sliders[front_index].top_y + R_from_front);
			gizmo.lift = std::max(gizmo.lift, this->back_sliders[back_index].top_y + R_from_back);
			gizmo.lift = std::max(gizmo.lift, this->back_bed[back_index].top_y + R_from_back);
			for (auto const &nc : carriers) {
				gizmo.lift = std::max(gizmo.lift, nc.second.horizon.get_value(front_index, front_index+1) + R_from_front);
			}
			gizmo.lift = std::max(gizmo.lift, crossings.check_crossing(front_index, back_index) + R_from_back);

			gizmo.set_lift.emplace_back([this,front_index,back_index](float lift){
				float top = lift + FaceHeight;
				//TODO: could set some lower if no to-front / to-back output loop

				this->front_bed[front_index].top_y = top;
				this->front_sliders[front_index].top_y = top;
				this->back_bed[back_index].top_y = top;
				this->back_sliders[back_index].top_y = top;
				for (auto &nc : carriers) {
					nc.second.horizon.raise_value(front_index, front_index+1, top);
				}
				crossings.add_crossing(front_index, back_index, top);
			});
		}

		//set up connections to front/back beds:
		FaceEdge front_loop_to_stitch, back_loop_to_stitch;
		FaceEdge front_loop_from_stitch, back_loop_from_stitch;
		{
			FaceEdge front_loop_in = front_bed[needle_index(front_needle)].top_edge;
			FaceEdge back_loop_in = back_bed[needle_index(back_needle)].top_edge;
			if (front_loop_in.is_valid() && front_loop_in.count > 0) {
				faces.emplace_back(source);
				gizmo.faces.emplace_back(faces.size()-1);
				Face &face = faces.back();

				std::string L = std::to_string(front_loop_in.count);
				face.type = "loop -l" + L +" x +l" + L + " x";
				float left = stitch_x(front_index, Left);
				float right = stitch_x(front_index, Right);
				face.vertices = {
					glm::vec3(left, 0.0f, front_bed.depth),
					glm::vec3(right, 0.0f, front_bed.depth),
					glm::vec3(glm::mix(left, right, 0.5f), 0.0f, ShearFront),
					glm::vec3(left, 0.0f, ShearFront),
				};

				FaceEdge loop_in(gizmo.faces.back(), 0, front_loop_in.count, FaceEdge::FlipNo);
				FaceEdge loop_out(gizmo.faces.back(), 2, front_loop_in.count, FaceEdge::FlipYes);

				gizmo.loop_verticals.emplace_back(front_index, front_loop_in, loop_in);

				front_loop_to_stitch = loop_out;
			}

			if (back_loop_in.is_valid() && back_loop_in.count > 0) {
				faces.emplace_back(source);
				gizmo.faces.emplace_back(faces.size()-1);
				Face &face = faces.back();

				std::string L = std::to_string(back_loop_in.count);
				face.type = "loop -l" + L +" x +l" + L + " x";
				float back_left = stitch_x(back_index, Left);
				float back_right = stitch_x(back_index, Right);
				float left = stitch_x(front_index, Left);
				float right = stitch_x(front_index, Right);
				face.vertices = {
					glm::vec3(back_left, 0.0f, back_bed.depth),
					glm::vec3(back_right, 0.0f, back_bed.depth),
					glm::vec3(right, 0.0f, ShearFront),
					glm::vec3(glm::mix(left, right, 0.5f), 0.0f, ShearFront),
				};

				FaceEdge loop_in(gizmo.faces.back(), 0, back_loop_in.count, FaceEdge::FlipNo);
				FaceEdge loop_out(gizmo.faces.back(), 2, back_loop_in.count, FaceEdge::FlipYes);

				gizmo.loop_verticals.emplace_back(back_index, back_loop_in, loop_in);

				back_loop_to_stitch = loop_out;
			}

			uint32_t front_loop_out_count = (is_front_to_back ? cs.size() : front_loop_in.count + back_loop_in.count);
			if (front_loop_out_count > 0) {
				faces.emplace_back(source);
				gizmo.faces.emplace_back(faces.size()-1);
				Face &face = faces.back();

				std::string L = std::to_string(front_loop_out_count);
				face.type = "loop -l" + L +" x +l" + L + " x";
				float left = stitch_x(front_index, Left);
				float right = stitch_x(front_index, Right);
				face.vertices = {
					glm::vec3(left, FaceHeight, ShearFront),
					glm::vec3(glm::mix(left, right, 0.5f), FaceHeight, ShearFront),
					glm::vec3(right, FaceHeight, front_bed.depth),
					glm::vec3(left, FaceHeight, front_bed.depth),
				};

				FaceEdge loop_in(gizmo.faces.back(), 0, front_loop_out_count, FaceEdge::FlipNo);
				FaceEdge loop_out(gizmo.faces.back(), 2, front_loop_out_count, FaceEdge::FlipYes);

				front_loop_from_stitch = loop_in;

				front_bed[front_index].top_edge = loop_out;
			} else {
				front_bed[front_index].top_edge = FaceEdge();
			}

			uint32_t back_loop_out_count = (is_front_to_back ? front_loop_in.count + back_loop_in.count : cs.size());
			if (back_loop_out_count > 0) {
				faces.emplace_back(source);
				gizmo.faces.emplace_back(faces.size()-1);
				Face &face = faces.back();

				std::string L = std::to_string(back_loop_out_count);
				face.type = "loop -l" + L +" x +l" + L + " x";
				float left = stitch_x(front_index, Left);
				float right = stitch_x(front_index, Right);
				float back_left = stitch_x(back_index, Left);
				float back_right = stitch_x(back_index, Right);
				face.vertices = {
					glm::vec3(glm::mix(left, right, 0.5f), FaceHeight, ShearFront),
					glm::vec3(right, FaceHeight, ShearFront),
					glm::vec3(back_right, FaceHeight, back_bed.depth),
					glm::vec3(back_left, FaceHeight, back_bed.depth),
				};

				FaceEdge loop_in(gizmo.faces.back(), 0, back_loop_out_count, FaceEdge::FlipNo);
				FaceEdge loop_out(gizmo.faces.back(), 2, back_loop_out_count, FaceEdge::FlipYes);

				back_loop_from_stitch = loop_in;

				back_bed[back_index].top_edge = loop_out;
			} else {
				back_bed[back_index].top_edge = FaceEdge();
			}


		}

		
		//build stitch @ shear plane:
		{
			faces.emplace_back(source);
			gizmo.faces.emplace_back(faces.size()-1);
			Face &face = faces.back();

			face.type = "split-";
			face.type += (is_front_to_back ? "front" : "back");
			face.type += "-to-";
			face.type += (dir == Right ? "right" : "left");
			face.type += " -l" + std::to_string(front_loop_to_stitch.count);
			face.type += " -l" + std::to_string(back_loop_to_stitch.count);
			face.type += (dir == Right ? " +y" : " -y") + std::to_string(cs.size());
			face.type += " +l" + std::to_string(back_loop_from_stitch.count);
			face.type += " +l" + std::to_string(front_loop_from_stitch.count);
			face.type += (dir == Right ? " -y" : " +y") + std::to_string(cs.size());

			float left = stitch_x(front_index, Left);
			float right = stitch_x(front_index, Right);

			face.vertices = {
					glm::vec3(left, 0.0f, ShearFront),
					glm::vec3(glm::mix(left,right,0.5f), 0.0f, ShearFront),
					glm::vec3(right, 0.0f, ShearFront),
					glm::vec3(right, FaceHeight, ShearFront),
					glm::vec3(glm::mix(left,right,0.5f), FaceHeight, ShearFront),
					glm::vec3(left, FaceHeight, ShearFront),
			};

			//connect to loops/yarns:
			if (dir == Right) {
				if (yarn_to_stitch.is_valid()) {
					gizmo.connections.emplace_back(yarn_to_stitch, FaceEdge(gizmo.faces.back(), 5, cs.size(), FaceEdge::FlipYes), "yarn to split-right");
				}
				if (yarn_from_stitch.is_valid()) {
					gizmo.connections.emplace_back(yarn_from_stitch, FaceEdge(gizmo.faces.back(), 2, cs.size(), FaceEdge::FlipNo), "yarn to split-right");
				}
			} else { assert(dir == Left);
				if (yarn_to_stitch.is_valid()) {
					gizmo.connections.emplace_back(yarn_to_stitch, FaceEdge(gizmo.faces.back(), 2, cs.size(), FaceEdge::FlipNo), "yarn to split-left");
				}
				if (yarn_from_stitch.is_valid()) {
					gizmo.connections.emplace_back(yarn_from_stitch, FaceEdge(gizmo.faces.back(), 5, cs.size(), FaceEdge::FlipYes), "yarn to split-left");
				}
			}

			if (front_loop_to_stitch.is_valid()) {
				gizmo.connections.emplace_back(front_loop_to_stitch, FaceEdge(gizmo.faces.back(), 0, front_loop_to_stitch.count, FaceEdge::FlipNo), "front loop to split");
			}
			if (back_loop_to_stitch.is_valid()) {
				gizmo.connections.emplace_back(back_loop_to_stitch, FaceEdge(gizmo.faces.back(), 1, back_loop_to_stitch.count, FaceEdge::FlipNo), "back loop to split");
			}
			if (front_loop_from_stitch.is_valid()) {
				gizmo.connections.emplace_back(front_loop_from_stitch, FaceEdge(gizmo.faces.back(), 4, front_loop_from_stitch.count, FaceEdge::FlipYes), "front loop from split");
			}
			if (back_loop_from_stitch.is_valid()) {
				gizmo.connections.emplace_back(back_loop_from_stitch, FaceEdge(gizmo.faces.back(), 3, back_loop_from_stitch.count, FaceEdge::FlipYes), "back loop from split");
			}

		}

		//apply lift and create merge faces:
		commit(gizmo);

		//if this is an xfer, anchor positions follow the moved loops:
		//update anchor positions (all carriers):
		if (cs.empty()) {
			for (auto &cn : carriers) {
				Carrier *c = &cn.second;
				if (!c->last_bed) continue;
				if (c->anchor_bed == &from_bed && c->anchor_needle == from_needle) {
					c->anchor_bed = &to_bed;
					c->anchor_needle = to_needle;
				}
			}
		}

	}


	//----- DEBUG helpers -----
	void DEBUG_add_horizons() {
		auto add_horizon = [this](float left, float right, float y, float depth) {
			faces.emplace_back(source);
			Face &face = faces.back();
			face.type = "horizon x x x";
			face.vertices = {
				glm::vec3(left, y, depth),
				glm::vec3(right, y, depth),
				glm::vec3(0.5f * (left + right), y + 0.25f * FaceHeight, depth),
			};
		};
		for (auto const &nc : carriers) {
			Carrier const *c = &nc.second;
			for (auto const &iv : c->horizon.values) {
				int32_t needle = iv.first / 4;
				int32_t ofs = iv.first - 4*needle;
				if (ofs < 0) {
					needle -= 1;
					ofs += 4;
				}
				assert(needle*4 + ofs == iv.first);
				float left, right;
				if (ofs == 0) {
					left = stitch_x(iv.first, Left);
					right = stitch_x(iv.first, Right);
				} else if (ofs == 1 || ofs == 3) {
					left = travel_x(iv.first, Left);
					right = travel_x(iv.first, Right);
				} else if (ofs == 2) {
					left = travel_x(4*needle+1, Right);
					right = travel_x(4*(needle+1)-1, Left);
				} else {
					assert(0 && "unreachable case");
				}
				add_horizon(left, right, iv.second, c->depth);
			}
		}
		for (auto bed : {&back_bed, &back_sliders, &front_sliders, &front_bed}) {
			for (auto const &ic : *bed) {
				int32_t needle = ic.first / 4;
				int32_t ofs = ic.first - 4*needle;
				if (ofs < 0) {
					needle -= 1;
					ofs += 4;
				}
				assert(needle*4 + ofs == ic.first);
				float left, right;
				if (ofs == 0) {
					left = stitch_x(ic.first, Left);
					right = stitch_x(ic.first, Right);
				} else if (ofs == 1 || ofs == 3) {
					assert(0 && "I don't think we use non-needle inds on bedcolumns");
				} else if (ofs == 2) {
					assert(0 && "I don't think we use non-needle inds on bedcolumns");
				} else {
					assert(0 && "unreachable case");
				}
				add_horizon(left, right, ic.second.top_y, bed->depth);
			}
		}
	}
};

int main(int argc, char **argv) {

	struct {
		std::string in_knitout = "";
		std::string out_smobj = "";
		bool drop_all = true;
		bool compact = false;
	} args;

	{ //parse arguments:
		bool end_of_options = false;
		std::vector< std::string > files;

		bool usage = false;
		for (int argi = 1; argi < argc; ++argi) {
			std::string arg = argv[argi];
			if (!end_of_options && arg.substr(0,2) == "--") {
				std::string option = arg.substr(2);
				if (option == "") {
					end_of_options = true;
				} else if (option == "compact") {
					args.compact = true;
				} else if (option == "no-drop") {
					args.drop_all = false;
				} else {
					std::cerr << "Unrecognized option '" << arg << "'" << std::endl;
					usage = true;
				}
			} else {
				files.push_back(arg);
			}
		}

		if (files.size() == 2) {
			args.in_knitout = files[0];
			args.out_smobj = files[1];
		} else {
			usage = true;
		}

		if (usage) {
			std::cerr << "Usage:\n\t./knitout-to-smobj [--compact] [--no-drop] [--] <in.knitout> <out.smobj>" << std::endl;
			return 1;
		}
	}

	std::cout << "Will interpret knitout in '" << args.in_knitout << "' and write smobj to '" << args.out_smobj << "'.\n";
	if (args.compact) std::cout << " The output file will be compacted (no helpful comments) (--compact).\n";
	if (args.drop_all) std::cout << " Any loops remaining at the end of the file will be dropped.\n";
	else std::cout << " Any loops remaining at the end of the file will be gaps in the yarn (--no-drop).\n";
	std::cout.flush();

	//parse knitout:

	std::ifstream knitout(args.in_knitout, std::ios::binary);

	enum {
		MagicValue,
		CommentHeaders,
		Body
	} section = MagicValue;

	std::unique_ptr< Translator > translator;

	std::string line;
	uint32_t line_number = 0;
	std::vector< std::string > lines;
	while (std::getline(knitout, line, '\n')) {
		line_number += 1;
		lines.emplace_back(line);
		if (translator) {
			translator->source = line_number;
		}

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
		} else if (tokens[0] == "in" || tokens[0] == "inhook") {
			std::vector< std::string > carriers(tokens.begin() + 1, tokens.end());
			translator->in(translator->lookup_carriers(carriers));
		} else if (tokens[0] == "releasehook") {
			//TODO: could provide some debugging output if releasehook not called when using inhook.
		} else if (tokens[0] == "out" || tokens[0] == "outhook") {
			std::vector< std::string > carriers(tokens.begin() + 1, tokens.end());
			translator->out(translator->lookup_carriers(carriers));
		} else if (tokens[0] == "miss") {
			//miss D N CS
			if (tokens.size() < 3) throw std::runtime_error("miss must have at least two parameters.");
			Direction dir;
			parse_direction(tokens[1], &dir);
			Bed bed;
			int32_t needle;
			parse_bedneedle(tokens[2], &bed, &needle);
			std::vector< std::string > carriers(tokens.begin() + 3, tokens.end());

			translator->miss(
				dir,
				translator->lookup_bed(bed),
				needle,
				translator->lookup_carriers(carriers)
			);

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
		} else if (tokens[0] == "drop") {
			//drop N
			if (tokens.size() != 2) throw std::runtime_error("drop should just have a needle parameter.");
			Bed bed;
			int32_t needle;
			parse_bedneedle(tokens[1], &bed, &needle);

			translator->knit(
				Right,
				translator->lookup_bed(bed),
				needle,
				std::vector< Carrier * >()
			);

		} else if (tokens[0] == "tuck") {
			//tuck D N CS
			if (tokens.size() < 3) throw std::runtime_error("tuck must have at least two parameters.");
			Direction dir;
			parse_direction(tokens[1], &dir);
			Bed bed;
			int32_t needle;
			parse_bedneedle(tokens[2], &bed, &needle);
			std::vector< std::string > carriers(tokens.begin() + 3, tokens.end());

			translator->tuck(
				dir,
				translator->lookup_bed(bed),
				needle,
				translator->lookup_carriers(carriers)
			);
		} else if (tokens[0] == "split") {
			//split D N N2 CS
			if (tokens.size() < 4) throw std::runtime_error("split must have at least three parameters.");
			Direction dir;
			parse_direction(tokens[1], &dir);
			Bed bed;
			int32_t needle;
			parse_bedneedle(tokens[2], &bed, &needle);
			Bed bed2;
			int32_t needle2;
			parse_bedneedle(tokens[3], &bed2, &needle2);
			std::vector< std::string > carriers(tokens.begin() + 4, tokens.end());

			translator->split(
				dir,
				translator->lookup_bed(bed), needle,
				translator->lookup_bed(bed2), needle2,
				translator->lookup_carriers(carriers)
			);
		} else if (tokens[0] == "xfer") {
			//xfer N N2
			if (tokens.size() != 3) throw std::runtime_error("xfer must have exactly three parameters.");
			Bed bed;
			int32_t needle;
			parse_bedneedle(tokens[1], &bed, &needle);
			Bed bed2;
			int32_t needle2;
			parse_bedneedle(tokens[2], &bed2, &needle2);

			translator->split(
				Right,
				translator->lookup_bed(bed), needle,
				translator->lookup_bed(bed2), needle2,
				std::vector< Carrier * >()
			);
		} else if (tokens[0] == "rack") {
			//rack R
			if (tokens.size() != 2) throw std::runtime_error("rack should have exactly one parameter.");
			std::istringstream str(tokens[1]);
			float rack;
			if (!(str >> rack)) throw std::runtime_error("rack value should be a number.");
			std::string temp;
			if (str >> temp) throw std::runtime_error("rack value had trailing junk: '" + temp + "...'.");
			if (!(std::floor(rack) == rack || std::floor(rack) + 0.25f == rack)) throw std::runtime_error("expecting either aligned or quarter-pitch racking.");
			translator->racking = rack;
		} else if (tokens[0] == "x-stitch-number") {
			if (tokens.size() != 2) throw std::runtime_error("x-stitch-number should have exactly one parameter.");
			std::istringstream str(tokens[1]);
			int32_t stitch_number;
			if (!(str >> stitch_number)) throw std::runtime_error("x-stitch-number should be an integer");
			std::string temp;
			if (str >> temp) throw std::runtime_error("x-stitch-number value had trailing junk: '" + temp + "...'.");
			if (!(0 <= stitch_number && stitch_number <= 120)) throw std::runtime_error("expecting stitch number in the interval [0,120]");
			translator->stitch_index = stitch_number;
		} else {
			std::cout << "WARNING: unsupported operation '" << line << "'" << std::endl;
		}
	}

	if (!translator) {
		std::cerr << "ERROR: ';;Carriers: ...' header is required." << std::endl;
		return 1;
	}

	if (translator->stats.empty_drops) {
		std::cout << "NOTE: pattern contained " << translator->stats.empty_drops << " instances where an empty needle was knit with an empty carrier set." << std::endl;
	}

	//reset source before auto-drop:
	line_number = 0;
	translator->source = 0;

	//out any remaining carriers:
	for (auto &nc : translator->carriers) {
		Carrier *c = &nc.second;
		if (c->ready) {
			std::cerr << "NOTE: taking out carrier '" << c->name << "', which was left in at the end of pattern." << std::endl;
			translator->out(std::vector< Carrier * >(1, c));
		}
	}

	if (args.drop_all) {
		//drop any remaining loops:
		int32_t min_index = std::numeric_limits< int32_t >::max();
		int32_t max_index = std::numeric_limits< int32_t >::min();
		auto check_bed = [&min_index,&max_index](BedColumns const &bed, int32_t offset) {
			for (auto const &pc : bed) {
				min_index = std::min(min_index, pc.first + offset);
				max_index = std::max(max_index, pc.first + offset);
			}
		};
		int32_t back_offset;
		{ //determine back_offset given racking:
			uint32_t r = int32_t(std::floor(translator->racking));
			if (std::floor(translator->racking) == translator->racking) { //aligned racking
				back_offset = r*4;
			} else { assert(std::floor(translator->racking) + 0.25f == translator->racking); //quarter pitch
				back_offset = r*4+2;
			}
		}

		check_bed(translator->back_bed, back_offset);
		check_bed(translator->back_sliders, back_offset);
		check_bed(translator->front_sliders, 0);
		check_bed(translator->front_bed, 0);
		
		std::vector< std::string > dropped;

		auto check_drop = [&translator,&dropped](BedColumns &bed, int32_t index, std::string const &bed_name) {
			if ((index % 4) != 0) return; //not a needle
			int32_t needle = index / 4;
			auto f = bed.find(index);
			if (f == bed.end()) return; //nothing here
			if (f->second.top_edge.count == 0) return; //no loop to drop
			//okay, need to drop:
			translator->knit(
				Right,
				bed,
				needle,
				std::vector< Carrier * >()
			);
			dropped.emplace_back(bed_name + std::to_string(needle));
		};

		for (int32_t idx = min_index; idx <= max_index; ++idx) {
			check_drop(translator->back_sliders, idx - back_offset, "bs");
			check_drop(translator->back_bed, idx - back_offset, "b");
			check_drop(translator->front_sliders, idx, "fs");
			check_drop(translator->front_bed, idx, "f");
		}
		if (!dropped.empty()) {
			std::cerr << "NOTE: dropped " << dropped.size() << " loops that remained on the bed at end of pattern:\n";
			for (auto const &d : dropped) {
				std::cerr << " " << d;
			}
			std::cerr << "\n";
		}

	}

	//write smobj file:
	std::cout << "Made " << translator->faces.size() << " faces." << std::endl;

	//This will add some blank "horizon" faces to make it clear where the tops of stacks are.
	//translator->DEBUG_add_horizons();

	std::ofstream out(args.out_smobj, std::ios::binary);

	//face type library:
	std::map< std::string, uint32_t > type_index;
	for (auto const &f : translator->faces) {
		if (type_index.insert(std::make_pair(f.type, type_index.size()+1)).second) {
			out << "L " << f.type << "\n";
		}
	}

	//faces:
	std::unordered_set< uint32_t > did_output; //did source code for a given line get dumpped yet?
	uint32_t vertex_index = 0;
	for (auto const &f : translator->faces) {
		std::string f_line = "f";
		for (auto const &v : f.vertices) {
			out << "v " << v.x << " " << v.y << " " << v.z << "\n";
			f_line += " " + std::to_string(++vertex_index);
		}
		out << f_line << "\n";
		out << "T " << type_index[f.type];
		if (!args.compact) out << " # " << f.type;
		out << "\n";
		out << "N " << f.source;
		if (!args.compact) {
			if (f.source > 0 && did_output.insert(f.source).second) {
				assert(f.source <= lines.size());
				out << " # \"" << lines[f.source-1] << "\"";
			}
		}
		out << "\n";
	}
	
	{ //PARANOIA: build list of edge types and check against connections!
		std::map< std::pair< uint32_t, uint32_t >, std::string > face_edge_type;
		for (auto const &f : translator->faces) {
			std::vector< std::string > edge_types;
			{
				std::istringstream str(f.type);
				std::string tok;
				str >> tok; //discard name
				while (str >> tok) {
					edge_types.emplace_back(tok);
				}
			}
			assert(edge_types.size() == f.vertices.size());
			uint32_t fi = &f - &(translator->faces[0]);
			for (uint32_t ei = 0; ei < edge_types.size(); ++ei) {
				auto ret = face_edge_type.insert(std::make_pair(std::make_pair(fi, ei), edge_types[ei]));
				assert(ret.second);
			}
		}
		for (auto const &c : translator->connections) {
			auto fa = face_edge_type.find(std::make_pair(c.a.face, c.a.edge));
			assert(fa != face_edge_type.end());
			auto fb = face_edge_type.find(std::make_pair(c.b.face, c.b.edge));
			assert(fb != face_edge_type.end());
			std::string ta = fa->second;
			std::string tb = fb->second;
			//std::cout << "Found " << ta << " to " << tb << " (" << c.why << ")" << std::endl;
			assert(ta.size() == tb.size());
			assert(ta.substr(1) == tb.substr(1));
			assert((ta[0] == '+' && tb[0] == '-') || (ta[0] == '-' && tb[0] == '+'));
		}

	}

	//edges:
	for (auto const &c : translator->connections) {
		out << "e " << (c.a.face+1) << "/" << (c.a.edge+1)
		    << " " << (c.b.face+1) << "/" << (c.a.flip == c.b.flip ? "-" : "") << (c.b.edge+1);
		if (!args.compact) out << " # " << c.why;
		out << "\n";
	}

	//units:
	for (auto const &u : translator->units) {
		out << "U " << u.name << " " << u.length << "\n";
	}

	{ //checkpoints:
		std::unordered_map< Unit const *, uint32_t > unit_index;
		for (auto const &u : translator->units) {
			unit_index[&u] = unit_index.size() + 1; //1-based
		}

		std::unordered_map< glm::uvec3, std::map< uint32_t, float > > checkpoints;
		//aggregate checkpoints:
		for (auto const &c : translator->checkpoints) {
			glm::uvec3 key = glm::uvec3(c.face, c.edge, c.crossing);
			checkpoints[key].insert( std::make_pair( unit_index[c.unit], 0.0f )).first->second += c.length;
		}

		//output aggregated checkpoints:
		for (auto const &c : translator->checkpoints) {
			glm::uvec3 key = glm::uvec3(c.face, c.edge, c.crossing);
			auto f = checkpoints.find(key);
			if (f == checkpoints.end()) continue; //I guess we already got this one

			out << "c " << (c.face+1) << "/" << (c.edge+1) << "/" << (c.crossing+1);
			for (auto const &l : f->second) {
				if (l.first == 0 && l.second == 0.0f) {
					continue; //remove 0.0*'1' -- any other zero-length checkpoints shouldn't exist but we'll leave 'em in anyway.
				}
				assert(l.first != 0 && "this process should *never* use an absolute unit.");
				out << ' ' << l.second << ' ' << l.first;
			}
			out << '\n';
			checkpoints.erase(f);
		}
		
	}

	return 0;
}
