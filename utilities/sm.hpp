#pragma once

#include <glm/glm.hpp>

#include <map>
#include <string>
#include <vector>

/*
 * Data structures for dealing with common smobj-related tasks:
 * sm::Mesh represents the contents of '.smobj' files
 * sm::Library represents a face library
 */

namespace sm {

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
	};
	struct Connection {
		FaceEdge a,b;
		bool flip = false;
		//!flip:  a[0] glued to b[0], a[1] glued to b[1]
		// flip: a[0] glued to b[1], a[1] glued to b[0]
	};
	std::vector< Connection > connections;

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

	//rip mesh : associate every face with its own vertex
	void rip(uint32_t from = 0);

};


//The contents of an '.sf' face library are stored as a 'Library':
struct Library {
	//----- db management: load/save to a file -----
	static Library load(std::string const &filename); //throws on error
	void save(std::string const &filename);

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

} //namespace sm
