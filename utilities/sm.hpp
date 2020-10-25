#pragma once

#include <glm/glm.hpp>
#include <set>
#include <map>
#include <string>
#include <vector>
#include <optional>
#include <variant>

/*
 * Data structures for dealing with common smobj-related tasks:
 * sm::Mesh represents the contents of '.smobj' files
 * sm::Library represents a face library
 * sm::Code represents a face program library
 */

namespace sm {

	struct BedNeedle{
		char bed = 'x'; // todo enum 'f','b'(hooks) 'F','B'(sliders) 'x' (dontcare)
		int needle = 0;
		// todo have a convention that nudge is always to the right of the lower needle index
		// to make it easier to compare without bothering with floats
		int8_t nudge = 0; // location = needle + 0.5 * nudge
		uint32_t depth = 0; // maybe a way to hold more needles at the same location without overlap, 0 is the actual needle
		BedNeedle(){}
		BedNeedle (char b, int n) {
			bed = b; needle = n;
		}
		BedNeedle (char b, float n) {
			bed = b; needle = n;
			nudge = (n - needle)*2;
		}
		float location() const{
			return 0.5f*nudge + needle;
		}
		bool operator< (const BedNeedle& rhs) const {
			float loc = location();
			float rloc = rhs.location();
			return std::tie(bed, loc) < std::tie(rhs.bed, rloc);
		}
		bool operator== (const BedNeedle& rhs) const {
			float loc = location();
			float rloc = rhs.location();
			return std::tie(bed, loc) == std::tie(rhs.bed, rloc);
		}
		std::string to_string() const{
			int n = location();
			return bed + std::to_string(n);
		}
		bool dontcare() const{
			return bed == 'x';
		}
	};

	// awkward that sm.hpp seems to have more and more knitout specific information, but well...
	struct Instr{
		enum Operation : char {
			In = 'i', Out = 'o', Knit = 'k', Split = 's', Xfer = 'x', Miss = 'm',  Tuck = 't', Drop = 'd', Unknown = 'u'
		} op = Unknown;
		enum Direction : char {
			Left = '-', Right = '+', None = '*'
		} direction = None;
		BedNeedle src; // location consumed (n/a: tuck)
		BedNeedle tgt;   // location produced (n/a: drop)
		BedNeedle tgt2; // aux produced (only for split)
		std::string yarns;
		int rack() const {
			//TODO enum this
			
			if(src.bed == 'b' && tgt.bed == 'f'){
				return src.needle - tgt.needle;
			}
			else if(tgt.bed == 'b' && src.bed == 'f'){
				return tgt.needle - src.needle;
			}
			return 0;
		}
		uint32_t face = -1U; // which smobj face (not code face) does this instruction come from, (-1U==> xfers for scheduling)
		
		std::string to_string(bool include_racking = false) const{
			std::string ret = "";
			switch(op){
				case Instr::In :
					ret += "in ";
					ret += yarns;
					break;
				case Instr::Out :
					ret += "out ";
					ret += yarns;
					break;
				case Instr::Knit:
					ret += "knit ";
					ret += (char)direction; ret += " ";
					ret += tgt.to_string() +  " ";
					ret += yarns;
					break;
				case Instr::Tuck:
					ret += "tuck ";
					ret += (char)direction; ret += " ";
					ret += tgt.to_string() + " ";
					ret += yarns;
					break;
				case Instr::Miss:
					ret += "miss ";
					ret += (char)direction; ret += " ";
					ret += tgt.to_string() + " ";
					ret += yarns;
					break;
				case Instr::Xfer:
					// insert rack, here
					if(include_racking){
						ret += "rack " + std::to_string(rack()) + "\n";
					}
					ret += "xfer ";
					ret += src.to_string() + " ";
					ret += tgt.to_string() + " ";
					break;
				case Instr::Split:
					// insert rack here
					if(include_racking){
						ret += "rack " + std::to_string(rack()) + "\n";
					}
					ret += "split ";
					ret += (char)direction; ret += " ";
					ret += tgt.to_string() + " "; // src, tgt or tgt,tgt2
					ret += tgt2.to_string() + " ";
					ret += yarns;
					break;
				case Instr::Drop:
					ret += "drop ";
					ret += src.to_string();
					break;
				default:
					assert(false && "unknown instruction");
			}
			return ret;
		}

		void translate(int offset){
			if(!src.dontcare()) src.needle += offset;
			if(!tgt.dontcare()) tgt.needle += offset;
			if(!tgt2.dontcare()) tgt2.needle += offset;
		}
		bool is_move() const{
			return (op == Knit || op == Tuck || op == Split || op == Drop);
		}
		bool is_loop() const{
			return (op == Xfer);
		}
	};


	

	struct MachineState{
		int racking = 0;
		int tension  = 0; // stitch value maybe
		// need to know if loop exists or not
		// need to know yarn location
		// TODO maintain connectivity
		std::map<BedNeedle, uint32_t> active_needles; // number of loops on a particular bed-needle
		std::vector< std::vector<Instr> > passes;
		bool make(Instr instr);
		bool empty();
	};


	struct Library;
	//The contents of a .smobj can be loaded as a "Mesh":
struct Mesh {
	//Library signatures ('L' lines):
	std::vector< std::string > library;

	//Vertex locations ('v' lines):
	std::vector< glm::vec3 > vertices;

	//Faces ('f' lines):
	//NOTE: vertex and type indices are 0-based, unlike in .smobj file:
	struct Face : public std::vector< uint32_t > {
		uint32_t type = -1U;
		uint32_t source = 0; //1-based source line number; 0 means 'unknown'
	};
	std::vector< Face > faces;

	//Connections ('e' lines):
	struct FaceEdge {
		uint32_t face = -1U;
		uint32_t edge = -1U;
		bool operator==(FaceEdge const &o) const { return (face == o.face && edge == o.edge); }
		bool operator< (const FaceEdge& rhs) const {
			return std::tie(face, edge) < std::tie(rhs.face, rhs.edge);
		}
	};
	struct Connection {
		FaceEdge a,b;
		bool flip = false;
		//!flip:  a[0] glued to b[0], a[1] glued to b[1]
		// flip: a[0] glued to b[1], a[1] glued to b[0]
		bool operator==(Connection const &o) const { return (a == o.a && b == o.b && flip == o.flip); }
	};
	std::vector< Connection > connections;

	// Structures for maintianing a xfer-only instruction stream + associating it with connections:
	
	std::vector< Instr > move_instructions; //the rest of the "stream" Is maintaining this within smobj the best idea?
	
	// is maintaining this inside smobj the best idea?
	struct MoveConnection {
		uint32_t c_idx = -1U;
		uint32_t i_idx = -1U;
		Connection connection; // keep a copy for sanity checking, for now
	};
	// slightly awkward but need to maintain 2 connections-> 1 instruction 
	std::vector< MoveConnection > move_connections; // a mapping of connections associated with move_instructions, the same instruction can have multiple connections

	std::vector< std::pair<uint32_t, uint32_t> > total_order; // face-id/instruction-id face_id = -1U, order from move_instructions
	struct Hint {
		enum HintType : char{
			Resource = 'r', // needle
			Order = 'o', // Instruction order
			Variant = 'v', // code face variant
		} type;
		enum HintSource : char{
			User = 'u',
			Inferred = 'i',
			Heuristic = 'h',
		} src;
		FaceEdge lhs;;
		std::variant<BedNeedle,  FaceEdge, std::string> rhs;
		uint32_t inferred_from = -1U; // optional index to maintain what constraint this was inferred from if at all
		std::string to_string() const{
			std::string str = std::to_string(lhs.face) + "/" + std::to_string(lhs.edge) + " ";
			if(type == Resource){
				BedNeedle r = std::get<BedNeedle>(rhs);
				str += std::to_string(r.bed) + "/" + std::to_string(r.needle);
			}
			else if(type == Order){
				FaceEdge r = std::get<FaceEdge>(rhs);
				str += std::to_string(r.face) + "/" + std::to_string(r.edge);
			}
			else if(type == Variant){
				std::string r = std::get<std::string>(rhs);
				str += r;
			}
			else{
				str += " NONE.";
			}
			return str;
		}
		bool operator== (const Hint& other) const {
		return std::tie(type,src,lhs,rhs) == std::tie(other.type, other.src, other.lhs, other.rhs);
	}

	};
	std::vector< Hint > hints; // can be an unordered_map, but hints could be multiple?

	//Units ('U' lines):
	struct Unit {
		std::string name;
		float length = 1.0f;
	};
	std::vector< Unit > units;

	//Checkpoints ('c' lines):
	//NOTE: ideally, *last* checkpoint on yarn will always have zero following length, but will exist. This code does not properly load a zero-length checkpoints. So it will, effectively, vanish.
	struct Checkpoint {
		uint32_t face = -1U;
		uint32_t edge = -1U;
		uint32_t crossing = -1U;
		float length = 0.0f;
		uint32_t unit = -1U;
	};
	std::vector< Checkpoint > checkpoints;

	//read a mesh from a ".smobj" file, throws on error:
	static Mesh load(std::string const &filename);

	//write an ".smobj" (precondition: mesh is valid)
	void save(std::string const &filename) const;

	//write a text version of the instructions, also needs library to keep track of edge types..
	void save_instructions(std::string const &filename, sm::Library const &library) const;

	//rip mesh : associate every face with its own vertex
	void rip(uint32_t from = 0);


	void remove_inferred_hints(uint32_t from = -1U);
	void remove_heuristic_hints(uint32_t from = -1U);
};


//The contents of an '.sf' face library are stored as a 'Library':
struct Library {
	//----- db management: load/save to a file -----
	static Library load(std::string const &filename); //throws on error
	void save(std::string const &filename) const;

	//----- library stores faces -----
	struct Face {
		std::string name; //descriptive name for type

		struct EdgePoint {
			EdgePoint(uint32_t edge_ = 0, float along_ = 0.5f, float z_ = 0.0f) : edge(edge_), along(along_), z(z_) { }
			uint32_t edge; //edge for point
			float along; //offset along edge (0 or 1 means yarn end => don't use this point when drawing)
			float z; //offset in normal direction
		};
		struct Yarn {
			EdgePoint begin, end;
			std::vector< glm::vec3 > middle;
		};
		struct Edge {
			glm::vec2 vertex;
			enum Direction : char {
				In = '-', Out = '+', Any = '*'
			} direction = Any;
			std::string type;
		};

		//face, with vertices in CCW order on the xy plane:
		std::vector< Edge > edges;

		//yarns passing through the face:
		std::vector< Yarn > yarns;

		//as a convenience, faces may be defined as "derived" from other faces:
		struct Derive {
			uint8_t by = 0; //an 'or' of ByBits
			std::string from = ""; //key() for the source face, or "" if not derived from any other face
			std::string expect_key = ""; //key() for the finished face or "" if not specified
			enum ByBit : uint8_t {
				MirrorXBit = (1 << 0),
				MirrorZBit = (1 << 1),
				ReverseYarnBit = (1 << 2),
			};
		} derive;

		//key() produces the corresponding .smobj library ('L' line) signature for the face:
		std::string key() const {
			std::string ret = name;
			for (auto const &e : edges) {
				ret += ' ';
				if (e.direction == Edge::In) ret += '-';
				else if (e.direction == Edge::Out) ret += '+';
				ret += e.type;
			}
			return ret;
		}
	};
	std::vector< Face > faces;
};

//The contents of an '.code' face library are stored as a 'Code':
// TODO: This makes sm very knitting specific, should this be structured differently?
struct Code {
	//----- db management: load/save to a file -----
	static Code load(std::string const &filename); //throws on error
	void save(std::string const &filename) const;


	//----- code per faces -----
	struct Face {
		std::string name; //descriptive name for type
		std::string variant="";
		bool check_bn = true; // if false, the template bedneedles are overwritten by hints, not checked
		std::vector<std::string> carriers;
		struct Edge {
			enum Direction : char {
				In = '-', Out = '+', Any = '*'
			} direction = Any;
			std::string type;
			BedNeedle bn;
			std::string yarns=""; //string (z-index?)
		};

		//face, with vertices in CCW order:
		std::vector< Edge > edges;
		//Instructions
		std::vector< Instr > instrs;

		//TODO: derive face code for opposite bed variant, opp direction variant
		//key() produces the corresponding .smobj library ('L' line) signature for the face:
		std::string key() const {
			std::string ret = name;
			for (auto const &e : edges) {
				ret += ' ';
				if (e.direction == Edge::In) ret += '-';
				else if (e.direction == Edge::Out) ret += '+';
				ret += e.type;
			}
			ret += ' ' + variant; //add variant to the key
			return ret;
		}
		std::string key_library() const {
			std::string ret = name;
			for (auto const &e : edges) {
				ret += ' ';
				if (e.direction == Edge::In) ret += '-';
				else if (e.direction == Edge::Out) ret += '+';
				ret += e.type;
			}
			return ret;
		}
		std::string knitout_string(int translate_to = 0, int index = -1, bool verbose = false) const {
			// verbose generates explicit instructions for rack + additional comments
			std::string ret ="";
			if(verbose){
				ret += ";from: " +  key() + ":" + "instr  " + std::to_string(index)+  "\n";
			}
			// concatenate all instructions
			for(auto const &i_ : instrs){
				int32_t idx = &i_ - &instrs[0];
				auto i = i_;
				if(index >= 0 && index != idx) continue;  // get instruction by index
				i.src.needle += translate_to;
				i.tgt.needle += translate_to;
				i.tgt2.needle += translate_to;
				ret += i.to_string(verbose);
				
				if(index < 0){
					ret += "\n";
				}
			}

			return ret;
		}
	};
	std::vector< Face > faces;
};




//The contents of a '.yarns' file are represented as a 'Yarns' structure:
struct Yarns {

	//load from a '.yarns' file, throw on error:
	static Yarns load(std::string const &filename);

	//save to a '.yarns' file:
	void save(std::string const &filename) const;

	//yarns:
	struct Yarn {
		std::vector< glm::vec3 > points;
		std::vector< uint32_t > sources; //1-based source line number; 0 means 'unknown'
		struct Checkpoint {
			uint32_t point;
			float length;
			uint32_t unit;
		};
		std::vector< Checkpoint > checkpoints;
		float radius = 0.1f; //yarns are radius-0.1f in canonical faces, but this can get scaled if faces are shrunk
		glm::u8vec4 color = glm::u8vec4(0xff, 0xff, 0xff, 0xff);
	};
	std::vector< Yarn > yarns;

	//units info:
	struct Unit {
		std::string name;
		float length;
	};
	std::vector< Unit > units;
};

//------ helper functions ------

// order face instructions to generate a total order (base order, may be reordered for efficiency)

bool compute_total_order(sm::Mesh &mesh, sm::Code const &code);

// order the mesh faces in a partial order that follows dependencies using the face library
bool can_order_faces(sm::Mesh const  &mesh,  sm::Library const  &library, std::vector<uint32_t> *_order);

// generate knitout code from an ordered set of faces using the code library
std::string knitout(sm::Mesh const &mesh, sm::Code const &code);

// can Hint h be added to mesh m without offending any existing hints?
bool add_hint(sm::Mesh::Hint h, sm::Mesh *mesh, sm::Library const &library, sm::Code const &code, std::vector<sm::Mesh::Hint> *offenders);

// verify existing hints
bool verify(sm::Mesh const &mesh,  Code const &code, std::vector<Mesh::Hint> *_offenders, bool strict=false);


// transfer and slack helpers
// insert a custom face to satisfy some slack and/or xfer constraint
// given the entire structure of the mesh.
bool create_out_slack_xfer(uint32_t face_id,  sm::Mesh &mesh, sm::Library &library, sm::Code &code);
bool create_in_slack_xfer(uint32_t face_id,  sm::Mesh &mesh, sm::Library &library, sm::Code &code);
//bool create_out_xfer();	// TODO
//bool create_in_xfer();	// TODO
//bool create_out_slack(); // this is not allowed, right (unless explicit)?
//bool create_in_slack();  // this is not allowed, right (unless explicit)?



//build from a mesh and library, connecting yarns over edges:
// will warn on inconsistent edge types/directions, and missing face types
// will erase output yarns structure
void mesh_and_library_to_yarns(Mesh const &mesh, Library const &library, Yarns *yarns);

//visualization helper: build tubes around a yarn spine:
// will append to attribs array
struct YarnAttribs {
	YarnAttribs(glm::vec3 const &Position_, glm::vec3 const &Normal_, glm::u8vec4 const &Color_) : Position(Position_), Normal(Normal_), Color(Color_) { }
	YarnAttribs(YarnAttribs const &) = default;
	glm::vec3 Position;
	glm::vec3 Normal;
	glm::u8vec4 Color;
};
enum Quality {
	QualityLow, //uncapped, low-poly-count
	QualityMedium, //capped, low-poly-count
	QualityHigh //capped, high-poly-count
};
//end_attrib, if supplied, is filled with the size of the attribute vector after finishing each yarn. So the attribs in range [ (*end_attrib)[i-1], (*end_attrib)[i] ) correspond to yarn i. (yarn 0 starts at 0, of course)
void yarns_to_tristrip(Yarns const &yarns, std::vector< YarnAttribs > *attribs, Quality quality = QualityLow, std::vector< size_t > *end_attrib = nullptr);


//face construction helper: derive one face from another by mirroring/reversing:
// note: contents of 'target' (except for name) will be over-written
void derive_face(sm::Library::Face const &source, uint8_t by_bits, sm::Library::Face *target);

// helper
// todo maybe this should be templated
bool  partial_order_to_sequence(std::set<std::pair<uint32_t, uint32_t>> partial, std::vector<uint32_t> *_sequence);

} //namespace sm
