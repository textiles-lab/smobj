#include "sm.hpp"

#include <glm/gtx/hash.hpp>

#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <deque>
#include <random>


//---------------------------------
//Mesh

sm::Mesh sm::Mesh::load(std::string const &filename) {
	sm::Mesh mesh;

	struct Signature {
		std::string name;
		std::vector< std::string > edges;
		std::string key() const {
			std::string ret = name;
			for (auto const &e : edges) {
				ret += ' ';
				ret += e;
			}
			return ret;
		}
	};

	std::vector< Signature > library; //name edge0 .. edgeN
	std::vector< uint32_t > types;
	std::vector< uint32_t > lines;

	std::ifstream in(filename, std::ios::binary);
	std::string line;
	while (std::getline(in, line)) {
		{ //trim comments
			auto i = line.find('#');
			if (i != std::string::npos) line = line.substr(0, i);
		}
		std::istringstream str(line);
		std::string cmd;
		if (!(str >> cmd)) continue; //ignore blank lines
		if (cmd == "L") {
			Signature signature;
			if (!(str >> signature.name)) throw std::runtime_error("no name in L row");
			std::string tok;
			while (str >> tok) {
				signature.edges.emplace_back(tok);
			}
			if (signature.edges.size() < 3) throw std::runtime_error("degenerate face type");
			library.emplace_back(signature);
			mesh.library.emplace_back(signature.key());
		} else if (cmd == "v") {
			glm::vec3 v;
			if (!(str >> v.x >> v.y >> v.z)) throw std::runtime_error("failed to read x y z after v");
			std::string temp;
			if (str >> temp) throw std::runtime_error("trailing junk (" + temp + "...) in v line");
			mesh.vertices.emplace_back(v);
		} else if (cmd == "f") {
			std::vector< uint32_t > inds;
			std::string tok;
			while (str >> tok) {
				std::istringstream tokstr(tok);
				int32_t vi;
				if (!(tokstr >> vi)) throw std::runtime_error("No vertex index in f token.");
				if (vi < 0) {
					vi = int32_t(mesh.vertices.size()) + vi;
				} else {
					vi = vi - 1;
				}
				if (!(vi >= 0 && vi < int32_t(mesh.vertices.size()))) throw std::runtime_error("Out-of-range vertex index in f token.");
				char slash;
				if (tokstr >> slash) {
					if (slash != '/') throw std::runtime_error("trailing (non-slash) junk in f token");
					int32_t ti;
					if (tokstr >> ti) {
						//TODO: deal with texture coordinates
					}
				}
				inds.emplace_back(vi);
			}
			if (inds.size() < 3) throw std::runtime_error("face with only " + std::to_string(inds.size()) + " vertices.");
			mesh.faces.emplace_back();
			static_cast< std::vector< uint32_t > & >(mesh.faces.back()) = inds;
		} else if (cmd == "T") {
			int32_t idx;
			if (!(str >> idx)) throw std::runtime_error("expecting 1-based type index after T");
			if (idx < 1 || idx > (int32_t)library.size()) throw std::runtime_error("invalid type index '" + std::to_string(idx) + "'");
			types.emplace_back(idx-1);
			std::string temp;
			if (str >> temp) throw std::runtime_error("trailing junk (" + temp + "...) in T line");
		} else if (cmd == "ln") {
			int32_t idx;
			if (!(str >> idx)) throw std::runtime_error("expecting line number after ln");
			if (idx < 0) throw std::runtime_error("invalid line number index '" + std::to_string(idx) + "'");
			lines.emplace_back(idx);
			std::string temp;
			if (str >> temp) throw std::runtime_error("trailing junk (" + temp + "...) in ln line");

		} else if (cmd == "e") {
			char slash1, slash2;
			int32_t a_face, a_edge;
			int32_t b_face, b_edge;
			if (!(str >> a_face >> slash1 >> a_edge >> b_face >> slash2 >> b_edge) || slash1 != '/' || slash2 != '/') throw std::runtime_error("e line that doesn't look like 1/4 3/1");
			bool flip = false;
			if (a_edge < 0) {
				a_edge = -a_edge;
				flip = !flip;
			}
			if (b_edge < 0) {
				b_edge = -b_edge;
				flip = !flip;
			}
			Connection con;
			con.flip = flip;

			if (a_face < 1 || a_face > int32_t(mesh.faces.size())) throw std::runtime_error("face out of range in e line.");
			con.a.face = a_face - 1;
			if (a_edge < 1 || a_edge > int32_t(mesh.faces[con.a.face].size())) throw std::runtime_error("edge out of range in e line.");
			con.a.edge = a_edge - 1;

			if (b_face < 1 || b_face > int32_t(mesh.faces.size())) throw std::runtime_error("face out of range in e line.");
			con.b.face = b_face - 1;
			if (b_edge < 1 || b_edge > int32_t(mesh.faces[con.b.face].size())) throw std::runtime_error("edge out of range in e line.");
			con.b.edge = b_edge - 1;

			mesh.connections.emplace_back(con);

			std::string temp;
			if (str >> temp) throw std::runtime_error("trailing junk (" + temp + "...) in e line");
		} else if (cmd == "U") { //unit definition
			std::string name;
			float length;
			if (!(str >> name >> length)) throw std::runtime_error("expecting name and length after U");
			if (length < 0.0f) throw std::runtime_error("invalid length '" + std::to_string(length) + "'");
			for (auto const &u : mesh.units) {
				if (u.name == name) throw std::runtime_error("duplicated unit name '" + u.name + "'");
			}
			if (mesh.units.empty() && (name != "1" || length != 1.0f)) throw std::runtime_error("First unit must be named '1' and have length 1.0");
			mesh.units.emplace_back();
			mesh.units.back().name = name;
			mesh.units.back().length = length;
			std::string temp;
			if (str >> temp) throw std::runtime_error("trailing junk (" + temp + "...) in U line");
		} else if (cmd == "c") { //checkpoint definition
			char slash1, slash2;
			int32_t face, edge, crossing;
			if (!(str >> face >> slash1 >> edge >> slash2 >> crossing) || slash1 != '/' || slash2 != '/') {
				throw std::runtime_error("c line that doesn't start like 2/3/1");
			}

			//validate face/edge (crossing hard to validate without library)
			if (face < 1 || face > int32_t(mesh.faces.size())) throw std::runtime_error("c line with invalid face '" + std::to_string(face) + "'");
			if (edge < 1 || edge > int32_t(mesh.faces[face].size())) throw std::runtime_error("c line with invalid edge '" + std::to_string(edge) + "'");
			if (crossing < 1) throw std::runtime_error("c line with invalid crossing '" + std::to_string(crossing) + "'");


			//read pairs of lengths and unit indices:
			float length;
			int32_t idx;
			while (str >> length >> idx) {
				if (length < 0.0f) throw std::runtime_error("c line has invalid length '" + std::to_string(length) + "'");
				if (idx < 1 || idx > int32_t(mesh.units.size())) throw std::runtime_error("c line with invalid unit index '" + std::to_string(idx) + "'");
				mesh.checkpoints.emplace_back();
				mesh.checkpoints.back().face = face-1;
				mesh.checkpoints.back().edge = edge-1;
				mesh.checkpoints.back().crossing = crossing-1;
				mesh.checkpoints.back().length = length;
				mesh.checkpoints.back().unit = idx-1;
			}

			std::string temp;
			if (str >> temp) throw std::runtime_error("trailing junk (" + temp + "...) in c line");
		} else {
			throw std::runtime_error("Unrecognized command '" + cmd + "'");
		}
	}

	std::unordered_set< glm::uvec2 > connected;
	for (auto const &c : mesh.connections) {
		if (!connected.insert(glm::uvec2(c.a.face, c.a.edge)).second) {
			throw std::runtime_error("Multiple connection to edge " + std::to_string(c.a.face+1) + "/" + std::to_string(c.a.edge+1) + ".");
		}
		if (!connected.insert(glm::uvec2(c.b.face, c.b.edge)).second) {
			throw std::runtime_error("Multiple connection to edge " + std::to_string(c.b.face+1) + "/" + std::to_string(c.b.edge+1) + ".");
		}
	}

	std::unordered_set< std::string > keys;
	for (auto const &s : library) {
		if (!keys.insert(s.key()).second) {
			throw std::runtime_error("Library signature '" + s.key() + "' is duplicated");
		}
	}

	if (types.size() != mesh.faces.size()) throw std::runtime_error("should be a 'T' for every face.");
	if (lines.empty()) lines.resize(mesh.faces.size(), 0);
	if (lines.size() != mesh.faces.size()) throw std::runtime_error("should be a 'ln' for every face or for no faces.");

	for (uint32_t i = 0; i < types.size(); ++i) {
		if (mesh.faces[i].size() != library[types[i]].edges.size()) {
			throw std::runtime_error("face/library edge count mismatch.");
		}
		mesh.faces[i].type = types[i];
		mesh.faces[i].line_number = lines[i];
	}


	std::vector< glm::vec3 > normals(mesh.vertices.size(), glm::vec3(0.0f));

	for (auto const &face : mesh.faces) {
		for (uint32_t vi = 0; vi < face.size(); ++vi) {
			glm::vec3 const &v = mesh.vertices[face[vi]];
			for (uint32_t i = 1; i + 1 < face.size(); ++i) {
				glm::vec3 const &a = mesh.vertices[face[(vi+i)%face.size()]];
				glm::vec3 const &b = mesh.vertices[face[(vi+i+1)%face.size()]];
				//lazy area-weighted normals:
				normals[face[vi]] += glm::cross(b-a, v-a);
			}
		}
	}
	for (auto &n : normals) {
		if (n == glm::vec3(0.0f)) {
			throw std::runtime_error("Vertex " + std::to_string(&n-&normals[0]+1) + " has sum-zero surrounding face area.");
		}
	}


	return mesh;
}

sm::Library sm::Library::load(std::string const &filename) {
	sm::Library library;
	sm::Library::Face *current = nullptr;

	std::ifstream in(filename, std::ios::binary);
	uint32_t line_number = 0;
	auto line_info = [&]() -> std::string {
		return "[" + filename + ":" + std::to_string(line_number) + "] ";
	};
	std::string line;
	while (std::getline(in, line)) {
		++line_number;
		{ //trim comments
			auto i = line.find('#');
			if (i != std::string::npos) line = line.substr(0, i);
		}
		std::istringstream str(line);
		auto expect_char = [&](char c) {
			char got;
			if (!(str >> got)) {
				throw std::runtime_error(line_info() + "Expecting '" + c + "' but failed to read");
			} else if (got != c) {
				throw std::runtime_error(line_info() + "Expecting '" + c + "' but got '" + got + "'");
			}
		};
		std::string tok;
		if (!(str >> tok)) continue; //ignore blank lines
		if (tok == "face") {
			std::string name;
			if (!(str >> name)) throw std::runtime_error(line_info() + "Failed to read name field from face line");
			std::string temp;
			if (str >> temp) throw std::runtime_error(line_info() + "Trailing junk (" + temp + "...) in face line");

			library.faces.emplace_back();
			library.faces.back().name = name;
			current = &library.faces.back();
		} else if (tok == "edge") {
			if (!current) throw std::runtime_error(line_info() + "edge line without face line");
			Library::Face::Edge edge;
			expect_char('(');
			if (!(str >> edge.vertex.x)) throw std::runtime_error(line_info() + "Failed to read x-coodinate");
			expect_char(',');
			if (!(str >> edge.vertex.y)) throw std::runtime_error(line_info() + "Failed to read x-coodinate");
			expect_char(')');
			std::string type;
			if (!(str >> type)) throw std::runtime_error(line_info() + "Failed to read [+-]?type in edge line");
			if (type[0] == '+') {
				edge.direction = sm::Library::Face::Edge::Out;
				edge.type = type.substr(1);
			} else if (type[0] == '-') {
				edge.direction = sm::Library::Face::Edge::In;
				edge.type = type.substr(1);
			} else {
				edge.direction = sm::Library::Face::Edge::Any;
				edge.type = type;
			}
			std::string temp;
			if (str >> temp) throw std::runtime_error(line_info() + "Trailing junk (" + temp + "...) in edge line");

			current->edges.emplace_back(edge);
		} else if (tok == "yarn") {
			if (!current) throw std::runtime_error(line_info() + "yarn line without face line");
			Library::Face::Yarn yarn;
			{ //begin
				expect_char('[');
				if (!(str >> yarn.begin.edge)) throw std::runtime_error(line_info() + "Failed to read edge index in yarn begin");
				if (!(yarn.begin.edge < current->edges.size())) throw std::runtime_error(line_info() + "Yarn begin edge is out-of-range.");
				expect_char(':');
				if (!(str >> yarn.begin.along)) throw std::runtime_error(line_info() + "Failed to read edge along");
				expect_char(',');
				if (!(str >> yarn.begin.z)) throw std::runtime_error(line_info() + "Failed to read edge z-value");
				expect_char(']');
			}
			//middle:
			while (true) {
				while (str.peek() != str.eof() && isspace(str.peek())) str.get();
				if (str.peek() == '[') break;
				expect_char('(');
				glm::vec3 pt;
				if (!(str >> pt.x)) throw std::runtime_error(line_info() + "Failed to middle point x");
				expect_char(',');
				if (!(str >> pt.y)) throw std::runtime_error(line_info() + "Failed to middle point x");
				expect_char(',');
				if (!(str >> pt.z)) throw std::runtime_error(line_info() + "Failed to middle point x");
				expect_char(')');

				yarn.middle.emplace_back(pt);
			}
			{ //end
				expect_char('[');
				if (!(str >> yarn.end.edge)) throw std::runtime_error(line_info() + "Failed to read edge index in yarn end");
				if (!(yarn.end.edge < current->edges.size())) throw std::runtime_error(line_info() + "Yarn end edge is out-of-range.");
				expect_char(':');
				if (!(str >> yarn.end.along)) throw std::runtime_error(line_info() + "Failed to read edge along");
				expect_char(',');
				if (!(str >> yarn.end.z)) throw std::runtime_error(line_info() + "Failed to read edge z-value");
				expect_char(']');
			}
			current->yarns.emplace_back(yarn);
		} else {
			throw std::runtime_error(line_info() + "Unrecognized line-start token '" + tok + "'");
		}
	}

	std::set< std::string > keys;

	for (auto const &face : library.faces) {
		
		auto ret = keys.insert(face.key());
		if (!ret.second) throw std::runtime_error("Duplicate face signature: '" + *ret.first + "'");
	}

	return library;
}

#if 0

void StitchFaces::save(std::string const &filename) {
	std::ofstream out(filename, std::ios::binary);
	for (auto const &name_face : faces) {
		std::string const &name = name_face.first;
		StitchFace const &face = name_face.second;
		out << "face " << name << '\n';
		for (auto const &edge : face.edges) {
			out << "\tedge ";
			out << '(' << edge.vertex.x << ',' << edge.vertex.y << ')';
			out << ' ';
			if      (edge.direction == StitchFace::Edge::In) out << "-";
			else if (edge.direction == StitchFace::Edge::Out) out << "+";
			out << edge.type;
			out << '\n';
		}
		for (auto const &yarn : face.yarns) {
			out << "\tyarn ";
			out << '[' << yarn.begin.edge << ':' << yarn.begin.along << ',' << yarn.begin.z << ']';
			for (auto const &v : yarn.middle) {
				out << ' ' << '(' << v.x << ',' << v.y << ',' << v.z << ')';
			}
			out << ' ' << '[' << yarn.end.edge << ':' << yarn.end.along << ',' << yarn.end.z << ']';
			out << '\n';
		}
	}
}
#endif

//------------------------------------------------

//helper for writing vectors of data:
template< typename T >
static void write(std::ostream &out, std::string magic, std::vector< T > const &data) {
	assert(magic.size() == 4);
	uint32_t size = sizeof(T) * data.size();
	out.write(magic.c_str(), 4);
	out.write(reinterpret_cast< const char * >(&size), sizeof(uint32_t));
	out.write(reinterpret_cast< const char * >(data.data()), sizeof(T)*data.size());
}

void sm::Yarns::save(std::string const &filename) const {
	std::ofstream out(filename, std::ios::binary);

	//arrays-of-structures that will be written:
	static_assert(sizeof(glm::vec3) == 12, "vec3 is packed");
	std::vector< glm::vec3 > out_points;

	struct YarnInfo {
		uint32_t point_begin;
		uint32_t point_end;
		float radius;
		glm::u8vec4 color;
	};
	static_assert(sizeof(YarnInfo) == 16, "YarnInfo is packed");
	std::vector< YarnInfo > out_yarns;

	std::vector< char > out_strings;

	struct UnitInfo {
		uint32_t name_begin;
		uint32_t name_end;
		float length;
	};
	static_assert(sizeof(UnitInfo) == 12, "UnitInfo is packed");
	std::vector< UnitInfo > out_units;

	struct CheckpointInfo {
		uint32_t point;
		float length;
		uint32_t unit;
	};
	static_assert(sizeof(CheckpointInfo) == 12, "CheckpointInfo is packed");
	std::vector< CheckpointInfo > out_checkpoints;

	std::vector< uint32_t > out_line_numbers;

	//fill the arrays:

	for (auto const &unit : units) {
		//PERHAPS: if (&unit == &units[0]) assert(unit.name == "1" && unit.length == 1.0f);
		out_units.emplace_back();
		out_units.back().name_begin = out_strings.size();
		out_strings.insert(out_strings.end(), unit.name.begin(), unit.name.end());
		out_units.back().name_end = out_strings.size();
		out_units.back().length = unit.length;
	}

	for (auto const &yarn : yarns) {
		assert(yarn.points.size() == yarn.line_numbers.size());

		for (auto const &cp : yarn.checkpoints) {
			assert(cp.unit < units.size());

			out_checkpoints.emplace_back();
			out_checkpoints.back().point = cp.point + out_points.size();
			out_checkpoints.back().length = cp.length;
			out_checkpoints.back().unit = cp.unit;
		}
		out_yarns.emplace_back();
		out_yarns.back().point_begin = out_points.size();
		out_points.insert(out_points.end(), yarn.points.begin(), yarn.points.end());
		out_line_numbers.insert(out_line_numbers.end(), yarn.line_numbers.begin(), yarn.line_numbers.end());
		out_yarns.back().point_end = out_points.size();
		out_yarns.back().radius = yarn.radius;
		out_yarns.back().color = yarn.color;
	}

	write(out, "f3..", out_points);
	write(out, "yarn", out_yarns);
	write(out, "strs", out_strings);
	write(out, "unit", out_units);
	write(out, "chk.", out_checkpoints);
	write(out, "src.", out_line_numbers);

}


//------------------------------------------------

float twice_area(glm::vec2 const &a, glm::vec2 const &b, glm::vec2 const &c) {
	return glm::dot( c - a, glm::vec2( -(b.y - a.y), b.x - a.x ) );
}

void sm::mesh_and_library_to_yarns(sm::Mesh const &mesh, sm::Library const &library, sm::Yarns *yarns_) {
	assert(yarns_);
	auto &yarns = *yarns_;
	yarns = Yarns();

	auto conv_along = [](float along) {
		return std::round(along * 1000);
	};

	auto inv_along = [](uint32_t along) {
		return 1000 - along;
	};

	auto conv_z = [](float z) {
		return std::round(z * 1000);
	};

	//look up all mesh face types in library:
	std::vector< sm::Library::Face const * > mesh_library;
	{
		std::map< std::string, sm::Library::Face const * > signature_to_face;
		for (auto const &face : library.faces) {
			signature_to_face.insert(std::make_pair(face.key(), &face));
		}
		mesh_library.reserve(mesh.library.size());
		for (auto const &sig : mesh.library) {
			auto f = signature_to_face.find(sig);
			if (f == signature_to_face.end()) {
				std::cerr << "WARNING: library is missing face with signature '" << sig << "'" << std::endl;
				mesh_library.emplace_back(nullptr);
			} else {
				mesh_library.emplace_back(f->second);
			}
		}
		assert(mesh_library.size() == mesh.library.size());
	}

	//generate long connected chains of yarn segments by looking at edge labels:
	struct ChainSegment {
		ChainSegment(Mesh::Face const *face_ = nullptr, uint32_t yarn_ = -1U, bool reverse_ = false) : face(face_), yarn(yarn_), reverse(reverse_) { }
		Mesh::Face const *face = nullptr;
		uint32_t yarn = -1U;
		bool reverse = false;
	};
	std::vector< std::vector< ChainSegment > > chains;
	{
		//build map connecting v0,v1,along,z pairs to other v0,v1,along,z via ChainSegment:
		std::vector< ChainSegment > segs;
		std::unordered_map< glm::ivec4, std::pair< glm::ivec4, uint32_t > > forward, reverse;

		for (auto const &face : mesh.faces) {
			//no yarns in face without type:
			assert(face.type < mesh_library.size());
			if (!mesh_library[face.type]) continue;
			auto const &sf = *mesh_library[face.type];
			for (uint32_t yi = 0; yi < sf.yarns.size(); ++yi) {
				segs.emplace_back(&face, yi, false);
				glm::ivec4 from;
				from.x = face[sf.yarns[yi].begin.edge];
				from.y = face[(sf.yarns[yi].begin.edge+1)%face.size()];
				from.z = conv_along(sf.yarns[yi].begin.along);
				from.w = conv_z(sf.yarns[yi].begin.z);
				glm::ivec4 to;
				to.x = face[sf.yarns[yi].end.edge];
				to.y = face[(sf.yarns[yi].end.edge+1)%face.size()];
				to.z = conv_along(sf.yarns[yi].end.along);
				to.w = conv_z(sf.yarns[yi].end.z);

				if (from.x > from.y) {
					std::swap(from.x,from.y);
					from.z = inv_along(from.z);
				}
				if (to.x > to.y) {
					std::swap(to.x,to.y);
					to.z = inv_along(to.z);
				}

				auto ret = forward.insert(std::make_pair(from, std::make_pair(to, segs.size()-1)));
				assert(ret.second);
				ret = reverse.insert(std::make_pair(to, std::make_pair(from, segs.size()-1)));
				assert(ret.second);
			}
		}

	
		std::unordered_map< glm::uvec2, glm::uvec3 > connections;
		for (auto const &c : mesh.connections) {
			glm::uvec2 from = glm::uvec2(
				mesh.faces[c.a.face][c.a.edge],
				mesh.faces[c.a.face][(c.a.edge+1)%mesh.faces[c.a.face].size()]
			);
			glm::uvec2 to = glm::uvec2(
				mesh.faces[c.b.face][c.b.edge],
				mesh.faces[c.b.face][(c.b.edge+1)%mesh.faces[c.b.face].size()]
			);

			if (c.flip) std::swap(to.x, to.y);

			if (from.x > from.y) {
				std::swap(from.x, from.y);
				std::swap(to.x, to.y);
			}
			{
				auto ret = connections.insert(std::make_pair(from, glm::ivec3(to, c.flip ? 1 : 0)));
				if (!ret.second) {
					std::cerr << "WARNING: multi-connected edge." << std::endl;
				}
			}

			std::swap(from, to);
			if (from.x > from.y) {
				std::swap(from.x, from.y);
				std::swap(to.x, to.y);
			}
			{
				auto ret = connections.insert(std::make_pair(from, glm::ivec3(to, c.flip ? 1 : 0)));
				if (!ret.second) {
					std::cerr << "WARNING: multi-connected edge." << std::endl;
				}
			}
		}
		auto paired_edge = [&](glm::ivec4 *at_) {
			assert(at_);
			auto &at = *at_;
			auto f = connections.find(glm::uvec2(at.x, at.y));
			if (f == connections.end()) {
				at.x = at.y = -1;
			} else {
				at.x = f->second.x;
				at.y = f->second.y;
				at.z = inv_along(at.z);

				if (at.x > at.y) {
					std::swap(at.x, at.y);
					at.z = inv_along(at.z);
				}
			}
		};

		/*
		//DEBUG:dump connections:
		for (auto arr : {&forward, &reverse})
		for (auto const &s : *arr) {
			std::cout << "[" << s.first.x << "," << s.first.y << " @ " << s.first.z << "]"
			      << " -> [" << s.second.first.x << "," << s.second.first.y << " @ " << s.second.first.z << "]";

			glm::ivec3 p = s.second.first;
			paired_edge(&p);
			if (p.x >= 0 && p.y >= 0) {
				std::cout << " -> [" << p.x << "," << p.y << " @ " << p.z << "]";
			}
			std::cout << "\n";
		}
		std::cout.flush();
		*/


		//now read + delete from the map:
		while (!forward.empty()) {
			assert(forward.size() == reverse.size());

			std::deque< ChainSegment > chain;

			glm::ivec4 seed = forward.begin()->first;
			//pass 0: read chain in one direction
			//pass 1: read chain in the opposite direction
			for (uint32_t pass = 0; pass < 2; ++pass) {
				glm::ivec4 at = seed;
				if (pass == 1) {
					paired_edge(&at);
				}
				//std::cout << "p" << pass << " [" << at.x << "," << at.y << " @ " << at.z << "]" << std::endl; //DEBUG
				while (true) {
					glm::ivec4 to;
					uint32_t via;
					bool via_reverse;
					{ //find a next segment:
						auto f = forward.find(at);
						if (f == forward.end()) {
							auto r = reverse.find(at);
							if (r == reverse.end()) {
								break;
							} else {
								//record:
								to = r->second.first;
								via = r->second.second;
								via_reverse = true;
								//consume:
								reverse.erase(r);
								f = forward.find(to);
								assert(f != forward.end());
								assert(f->second.first == at);
								assert(f->second.second == via);
								forward.erase(f);
							}
						} else {
							//record:
							to = f->second.first;
							via = f->second.second;
							via_reverse = false;
							//consume:
							forward.erase(f);
							auto r = reverse.find(to);
							assert(r != reverse.end());
							assert(r->second.first == at);
							assert(r->second.second == via);
							reverse.erase(r);
						}
					}
					at = to;
					//std::cout << "   -> [" << at.x << "," << at.y << " @ " << at.z << "]"; //DEBUG
					paired_edge(&at);
					//std::cout << " -> [" << at.x << "," << at.y << " @ " << at.z << "]" << std::endl; //DEBUG
					if (pass == 0) {
						chain.emplace_back(segs[via]);
						chain.back().reverse = via_reverse;
					} else {
						chain.emplace_front(segs[via]);
						chain.front().reverse = !via_reverse;
					}
				}
			}

			chains.emplace_back(chain.begin(), chain.end());
		}

		{ //show chain sizes breakdown:
			std::map< uint32_t, uint32_t > hist;
			for (auto const &chain : chains) {
				hist.insert(std::make_pair(chain.size(), 0)).first->second += 1;
			}
			std::cout << "Have " << chains.size() << " chains:\n";
			for (auto const &sc : hist) {
				std::cout << "  " << sc.second << " of size " << sc.first << "\n";
			}
			std::cout.flush();
		}
	}

	//Idea: express xy coords of yarn points in each face in terms of generalized barycentric coordinates.
	//particularly, use the formulation of [Wachpress 1975], as given in equation (6) of:
	//  http://geometry.caltech.edu/pubs/MHBD02.pdf
	//
	// for polygon  .. q_i .. with point p,
	//  weight for vertex j is the product of the area of triangle triangle q_{j-1}, q_{j}, q_{j+1} and the areas
	//  of all the triangles between p and edge segments that don't include j.

	struct BaryFace {
		BaryFace() = default;
		BaryFace(sm::Library::Face const &face) {
			auto get_coords = [&face](glm::vec2 const &p) {
				std::vector< float > areas;
				areas.reserve(face.edges.size());
				for (uint32_t i = 0; i < face.edges.size(); ++i) {
					glm::vec2 const &a = face.edges[i].vertex;
					glm::vec2 const &b = face.edges[(i + 1 < face.edges.size() ? i + 1 : 0)].vertex;
					float area = twice_area(p,a,b);
					areas.emplace_back(area);
				}
				std::vector< float > coords;
				float total_coord = 0.0f;
				coords.reserve(face.edges.size());
				for (uint32_t i = 0; i < face.edges.size(); ++i) {
					uint32_t pi = (i > 0 ? i - 1 : face.edges.size() - 1);
					glm::vec2 const &prev = face.edges[pi].vertex;
					glm::vec2 const &cur = face.edges[i].vertex;
					glm::vec2 const &next = face.edges[(i + 1 < face.edges.size() ? i + 1 : 0)].vertex;
					float coord = twice_area(prev, cur, next);
					for (uint32_t a = 0; a < areas.size(); ++a) {
						if (a != i && a != pi) coord *= areas[a];
					}
					total_coord += coord;
					coords.emplace_back(coord);
				}
				assert(total_coord > 0.0f);
				float inv_total_coord = 1.0f / total_coord;
				for (auto &coord : coords) {
					coord *= inv_total_coord;
				}
				return coords;
			};
			for (auto const &yarn : face.yarns) {
				yarns.emplace_back();
				auto do_edge = [&](sm::Library::Face::EdgePoint const &ep) {
					if (ep.along == 0.0f || ep.along == 1.0f) {
						//don't add edge points for yarn ends
						return;
					}
					yarns.back().emplace_back();
					YarnPoint &yp = yarns.back().back();
					yp.resize(face.edges.size(), 0.0f);
					yp[ep.edge] = 1.0f - ep.along;
					yp[(ep.edge + 1)%yp.size()] = ep.along;
					yp.z = ep.z;
				};
				do_edge(yarn.begin);
				for (auto const &m : yarn.middle) {
					yarns.back().emplace_back();
					YarnPoint &yp = yarns.back().back();
					static_cast< std::vector< float > & >(yp) = get_coords(glm::vec2(m));
					yp.z = m.z;
				}
				do_edge(yarn.end);
			}
		}
		struct YarnPoint : public std::vector< float > { float z = 0.0f; };

		std::vector< std::vector< YarnPoint > > yarns;
	};

	std::vector< BaryFace > mesh_baryfaces;
	mesh_baryfaces.reserve(mesh_library.size());

	for (auto fp : mesh_library) {
		if (!fp) {
			mesh_baryfaces.emplace_back();
			continue;
		}
		mesh_baryfaces.emplace_back(BaryFace(*fp));
	}
	assert(mesh_baryfaces.size() == mesh_library.size());


	//---------------
	//compute scaling factor for yarn radius and z offset based on average face areas:
	float radius_scale = 1.0f;
	{
		std::vector< float > mesh_proto_area;
		mesh_proto_area.reserve(mesh_library.size());
		for (auto fp : mesh_library) {
			if (!fp) {
				mesh_proto_area.emplace_back(0.0f);
				continue;
			}
			auto const &face = *fp;
			float area = 0.0f;
			for (uint32_t i = 2; i < face.edges.size(); ++i) {
				area += twice_area(face.edges[0].vertex, face.edges[i-1].vertex, face.edges[i].vertex);
			}
			mesh_proto_area.emplace_back(area);
		}
		assert(mesh_proto_area.size() == mesh_library.size());

		float avg_ratio = 0.0f;
		uint32_t avg_ratio_count = 0;

		for (auto const &face : mesh.faces) {
			if (!mesh_library[face.type]) continue;
			glm::vec3 middle = glm::vec3(0.0f);
			for (uint32_t vi = 0; vi < face.size(); ++vi) {
				middle += mesh.vertices[face[vi]];
			}
			middle /= face.size();

			float area = 0.0f;
			for (uint32_t vi = 0; vi < face.size(); ++vi) {
				glm::vec3 const &a = mesh.vertices[face[vi]];
				glm::vec3 const &b = mesh.vertices[face[(vi+1)%face.size()]];
				area += glm::length(glm::cross(a - middle, b - middle));
			}

			float ratio = area / mesh_proto_area[face.type];

			avg_ratio += ratio;
			avg_ratio_count += 1;
		}
		if (avg_ratio_count == 0) {
			avg_ratio = 1.0f;
		} else {
			avg_ratio /= avg_ratio_count;
		}
		radius_scale = std::sqrt(avg_ratio);
		std::cout << "Idea radius scaling, based on " << avg_ratio_count << " area ratios: " << radius_scale << std::endl;
	}

	//---------------
	//compute vertex normals:

	std::vector< glm::vec3 > normals(mesh.vertices.size(), glm::vec3(0.0f));

	for (auto const &face : mesh.faces) {
		for (uint32_t vi = 0; vi < face.size(); ++vi) {
			glm::vec3 const &v = mesh.vertices[face[vi]];
			for (uint32_t i = 1; i + 1 < face.size(); ++i) {
				glm::vec3 const &a = mesh.vertices[face[(vi+i)%face.size()]];
				glm::vec3 const &b = mesh.vertices[face[(vi+i+1)%face.size()]];

				//area-weighted is okay with degenerate triangles: assert(v != a && v != b);

				//lazy area-weighted normals:
				normals[face[vi]] += glm::cross(b-a, v-a);

				/*
				//this appears to be broken

				//triangle a, b, v contributes to normal at v by:
				// cot( v-a, b-a ) * (v-b)
				// cot( v-b, a-b ) * (v-a)

				// n.b. cot = cos / sin

				auto vec_cot = [](glm::vec3 const &A, glm::vec3 const &B) {
					float sin = glm::dot(A,B); //sin * |A| * |B|
					std::cout << A.x << " " << A.y << " " << A.z << " / " << B.x << " " << B.y << " " << B.z << std::endl;
					assert(sin != 0.0f);
					float cos = glm::length(glm::cross(A,B)); //cos * |A| * |B|
					float ret = cos / sin;

					return ret;
				};
				assert(face[vi] < normals.size());

				normals[face[vi]] += vec_cot(v-a, b-a) * (b-v);
				normals[face[vi]] += vec_cot(v-b, a-b) * (a-v);
				*/
			}
		}
	}
	for (auto &n : normals) {
		n = glm::normalize(n);
	}

	/*
	//DEBUG: just show normals.
	for (uint32_t i = 0; i < mesh.verts.size(); ++i) {
		yarns.emplace_back();
		std::vector< glm::vec3 > &yarn = yarns.back();
		yarn.emplace_back(mesh.verts[i]);
		yarn.emplace_back(mesh.verts[i] + 0.3f * normals[i]);
	}
	{
		yarns.emplace_back();
		yarns.back().emplace_back(-1.0f, 0.0f, 0.0f);
		yarns.back().emplace_back( 1.0f, 0.0f, 0.0f);
		yarns.emplace_back();
		yarns.back().emplace_back(0.0f,-1.0f, 0.0f);
		yarns.back().emplace_back(0.0f, 1.0f, 0.0f);
	}*/

	//---------------
	// face yarns -> mesh

	auto yarn_color = []() -> glm::u8vec4 {
		static std::mt19937 mt(0x15469519);
		glm::vec3 col;
		float h = mt() / float(mt.max());

		h = h - std::floor(h);
		h *= 6.0f;
		if (h < 1.0f) {
			col.r = 1.0f;
			col.g = h - 0.0f;
			col.b = 0.0f;
		} else if (h < 2.0f) {
			col.r = 2.0f - h;
			col.g = 1.0f;
			col.b = 0.0f;
		} else if (h < 3.0f) {
			col.r = 0.0f;
			col.g = 1.0f;
			col.b = h - 2.0f;
		} else if (h < 4.0f) {
			col.r = 0.0f;
			col.g = 4.0f - h;
			col.b = 1.0f;
		} else if (h < 5.0f) {
			col.r = h - 4.0f;
			col.g = 0.0f;
			col.b = 1.0f;
		} else { //(h < 6.0f)
			col.r = 1.0f;
			col.g = 0.0f;
			col.b = 6.0f - h;
		}

		return glm::u8vec4(
			std::max(0, std::min(255, int32_t(col.r * 255))),
			std::max(0, std::min(255, int32_t(col.g * 255))),
			std::max(0, std::min(255, int32_t(col.b * 255))),
			255
		);
	};

	float mismatch = 0.0f;
	for (auto const &chain : chains) {
		yarns.yarns.emplace_back();
		yarns.yarns.back().radius *= radius_scale;
		yarns.yarns.back().color = yarn_color();

		for (auto const &seg : chain) {
			assert(seg.face);
			auto const &face = *seg.face;

			auto &bf = mesh_baryfaces[face.type];

			auto yarn = bf.yarns[seg.yarn];
			if (seg.reverse) {
				std::reverse(yarn.begin(), yarn.end());
			}

			for (auto const &yp : yarn) {
				assert(yp.size() == face.size());
				glm::vec3 acc = glm::vec3(0.0f);
				glm::vec3 normal_acc = glm::vec3(0.0f);
				for (uint32_t i = 0; i < yp.size(); ++i) {
					acc += yp[i] * mesh.vertices[face[i]];
					normal_acc += yp[i] * normals[face[i]];
				}
				normal_acc = glm::normalize(normal_acc); //maybe
				glm::vec3 pt = acc + normal_acc * (yp.z * radius_scale);
				if (&yp == &yarn[0] && !yarns.yarns.back().points.empty()) {
					//don't duplicate boundary points.
					mismatch += glm::length(pt - yarns.yarns.back().points.back());
					//average just in case
					yarns.yarns.back().points.back() = 0.5f * (yarns.yarns.back().points.back() + pt);
					yarns.yarns.back().line_numbers.back() = face.line_number;
				} else {
					yarns.yarns.back().points.emplace_back(pt);
					yarns.yarns.back().line_numbers.emplace_back(face.line_number);
				}
			}
		}
		//TODO: for circular yarns, make sure first/last points match exactly.
	}
	std::cout << "Total boundary mis-match from chains (should be very small): " << mismatch << std::endl;

	/*
	for (auto const &face : mesh.faces) {
		auto f = baryfaces.find(face.type);
		if (f == baryfaces.end()) continue;
		auto &bf = f->second;
		for (auto const &yarn : bf.yarns) {
			yarns.emplace_back();
			yarns.back().radius *= radius_scale;
			for (auto const &yp : yarn) {
				assert(yp.size() == face.size());
				glm::vec3 acc = glm::vec3(0.0f);
				glm::vec3 normal_acc = glm::vec3(0.0f);
				for (uint32_t i = 0; i < yp.size(); ++i) {
					acc += yp[i] * mesh.verts[face[i]];
					normal_acc += yp[i] * normals[face[i]];
				}
				normal_acc = glm::normalize(normal_acc); //maybe
				yarns.back().emplace_back(acc + normal_acc * (yp.z * radius_scale));
			}
		}
	}*/

	//---------------

	glm::vec3 min = glm::vec3(std::numeric_limits< float >::infinity());
	glm::vec3 max = glm::vec3(-std::numeric_limits< float >::infinity());
	for (auto const &yarn : yarns.yarns) {
		for (auto const &v : yarn.points) {
			min = glm::min(min, v);
			max = glm::max(max, v);
		}
	}

	/*
	{	//DEBUG: "bounding box"
		yarns.emplace_back();
		yarns.back().emplace_back(min.x, min.y, min.z);
		yarns.back().emplace_back(min.x, min.y, min.z);
		yarns.back().emplace_back(max.x, min.y, min.z);
		yarns.back().emplace_back(max.x, min.y, min.z);
		yarns.back().emplace_back(max.x, max.y, min.z);
		yarns.back().emplace_back(max.x, max.y, min.z);
		yarns.back().emplace_back(min.x, max.y, min.z);
		yarns.back().emplace_back(min.x, max.y, min.z);
		yarns.back().emplace_back(min.x, min.y, max.z);
		yarns.back().emplace_back(min.x, min.y, max.z);
		yarns.back().emplace_back(max.x, min.y, max.z);
		yarns.back().emplace_back(max.x, min.y, max.z);
		yarns.back().emplace_back(max.x, max.y, max.z);
		yarns.back().emplace_back(max.x, max.y, max.z);
		yarns.back().emplace_back(min.x, max.y, max.z);
		yarns.back().emplace_back(min.x, max.y, max.z);
		yarns.back().emplace_back(min.x, min.y, max.z);
		yarns.back().emplace_back(min.x, min.y, max.z);
	}
	*/


	std::cout << "Generated yarns lie in [" << min.x << "," << max.x << "]x[" << min.y << ", " << max.y << "]x[" << min.z << "," << max.z << "]." << std::endl;

}

#if 0

//------------------------------------------------

void StitchFaces::yarns_to_tristrip(std::vector< Yarn > const &yarns, std::vector< YarnAttribs > *attribs_, Quality quality) {
	assert(attribs_);
	auto &attribs = *attribs_;

	uint32_t Angles = (quality != QualityLow ? 16 : 8);
	std::vector< glm::vec2 > Circle;
	Circle.reserve(Angles);
	for (uint32_t a = 0; a < Angles; ++a) {
		double ang = M_PI * 2.0 * double(a) / double(Angles);
		Circle[a].x = std::cos(ang);
		Circle[a].y = std::sin(ang);
	}

	for (auto const &yarn : yarns) {
		if (yarn.size() < 2) continue;

		/*//DEBUG: no smoothing:
		for (uint32_t i = 1; i < yarn.size(); ++i) {
			add_tube(yarn[i-1], yarn[i], yarn.radius, yarn.color, 1.0f, 1.0f);
		}
		continue;
		*/
		//what if, instead, we drew with quadratic b-splines?
		//basic function:
		// x^2/2  0..1
		// (-2x^2 + 6x - 3)/2 1..2
		// (3-x)^2/2  2..3

		//this is the same as a quadratic bezier through the midpoints.

		std::vector< glm::vec3 > poly;
		poly.emplace_back(yarn[0]);
		poly.emplace_back(0.5f * (yarn[1] + yarn[0]));

		float thresh = std::cos((quality == QualityHigh ? 5.0f : 15.0f) / 180.0f * float(M_PI));

		std::function< void(glm::vec3 const &, glm::vec3 const &) > curve_to = [&](glm::vec3 const &m, glm::vec3 const &e) {
			{ //is curve already "flat enough"?
				glm::vec3 a = (m - poly.back());
				glm::vec3 b = (e - m);
				if (glm::dot(a,b) > std::sqrt(glm::dot(a,a) * glm::dot(b,b)) * thresh
				|| glm::dot(b-a,b-a) < 0.5 * 0.5 * yarn.radius * yarn.radius) {
					poly.emplace_back(m);
					poly.emplace_back(e);
					return;
				}
			}
			{ //subdivide:
				glm::vec3 sm = 0.5f * (poly.back() + m);
				glm::vec3 me = 0.5f * (m + e);
				glm::vec3 sme = 0.5f * (sm + me);
				curve_to(sm, sme);
				curve_to(me, e);
			}
			/*
			
			glm::vec3 s = at;
			for (uint32_t i = 1; i <= 5; ++i) {
				float t = (i / 5.0f);
				glm::vec3 a = glm::mix(s,m,t);
				glm::vec3 b = glm::mix(m,e,t);
				line_to(glm::mix(a,b,t));
			}*/
		};
		for (uint32_t i = 2; i < yarn.size(); ++i) {
			curve_to(yarn[i-1], 0.5f * (yarn[i-1] + yarn[i]));
		}
		poly.emplace_back(yarn.back());

		for (auto const &v : poly) {
			assert(v.x == v.x); //PARANOIA: are we getting NaNs?
		}

		//Now draw a cylinder around the polyline:

		{
			glm::vec3 a = poly[0];
			glm::vec3 along, p1, p2;
			{ //initial frame for verts:
				along = glm::normalize(poly[1]-poly[0]);
				if (std::abs(along.x) <= std::abs(along.y) && std::abs(along.x) <= std::abs(along.z)) {
					p1 = glm::vec3(1.0f, 0.0f, 0.0f);
				} else if (std::abs(along.y) <= std::abs(along.z)) {
					p1 = glm::vec3(0.0f, 1.0f, 0.0f);
				} else {
					p1 = glm::vec3(0.0f, 0.0f, 1.0f);
				}
				p1 = glm::normalize(p1 - glm::dot(p1, along) * along);
				p2 = glm::cross(along, p1);
			}

			//starting cap:
			for (uint32_t r = 1; r <= Angles / 4; ++r) {
				if (!attribs.empty()) attribs.emplace_back(attribs.back());
				for (uint32_t i = 0; i < Angles; ++i) {
					glm::vec3 n0 = -Circle[r-1].y * along
					             + Circle[r-1].x * ( p1 * Circle[i].x + p2 * Circle[i].y );
					glm::vec3 n1 = -Circle[r].y * along
					             + Circle[r].x * ( p1 * Circle[i].x + p2 * Circle[i].y );	
					attribs.emplace_back(a + yarn.radius * n0, n0, yarn.color);
					if (i == 0 && attribs.size() > 1) attribs.emplace_back(attribs.back());
					attribs.emplace_back(a + yarn.radius * n1, n1, yarn.color);
				}
				attribs.emplace_back(attribs[attribs.size() - 2*Angles]);
				attribs.emplace_back(attribs[attribs.size() - 2*Angles]);
			}

			glm::vec3 d1 = p1;
			glm::vec3 d2 = p2;
			glm::vec3 at = a;
			auto connect_to = [&](glm::vec3 const &next_d1, glm::vec3 const &next_d2, glm::vec3 const &next_at) {
				attribs.emplace_back(attribs.back());
				for (uint32_t i = 0; i < Angles; ++i) {
					glm::vec3 n = d1 * Circle[i].x + d2 * Circle[i].y;
					glm::vec3 next_n = next_d1 * Circle[i].x + next_d2 * Circle[i].y;
					attribs.emplace_back(next_at + yarn.radius * next_n, next_n, yarn.color);
					if (i == 0) attribs.emplace_back(attribs.back());
					attribs.emplace_back(at + yarn.radius * n, n, yarn.color);
				}
				attribs.emplace_back(attribs[attribs.size() - 2*Angles]);
				attribs.emplace_back(attribs[attribs.size() - 2*Angles]);
				
				d1 = next_d1;
				d2 = next_d2;
				at = next_at;
			};

			for (uint32_t bi = 1; bi < poly.size(); ++bi) {
				glm::vec3 const &b = poly[bi];
				if (b == a) continue; //<-- shouldn't be needed(?)

				glm::vec3 next_along = glm::normalize(b - a);

				glm::vec3 next_p1, next_p2;
				{ //rotate p1 to get next_p1, compute next_p2:
					if (std::abs(glm::dot(along, next_along)) < 0.999f) {
						glm::vec3 perp = glm::normalize(glm::cross(along, next_along));
						assert(perp.x == perp.x); //No NaNs!
						glm::vec3 out = glm::cross(along, perp);
						glm::vec3 next_out = glm::cross(next_along, perp);
						next_p1 = p1 + glm::dot(p1, out) * (next_out - out);
					} else {
						//NOTE: if along/next_along are 180 degrees apart, this may not work:
						next_p1 = p1;
					}

					next_p1 = glm::normalize(next_p1 - glm::dot(next_p1, next_along) * next_along);
					next_p2 = glm::cross(next_along, next_p1);

				}

				if (bi + 1 == poly.size()) {
					connect_to(next_p1, next_p2, b);
				} else {
					//compute half-way frame:
					//NOTE: could do some scaling of half_p[12] to account for cross-section loss, but the effect seems small enough to ignore.
					glm::vec3 half_along = glm::normalize(along + next_along);
					glm::vec3 half_p1 = glm::normalize(p1 + next_p1);
					glm::vec3 half_p2 = glm::cross(half_along, half_p1);
					connect_to(half_p1, half_p2, b);
				}

				along = next_along;
				p1 = next_p1;
				p2 = next_p2;
				a = b;

			}

			for (uint32_t r = 1; r <= Angles / 4; ++r) {
				attribs.emplace_back(attribs.back());
				for (uint32_t i = 0; i < Angles; ++i) {
					glm::vec3 n0 = Circle[r-1].y * along
					             + Circle[r-1].x * ( p1 * Circle[i].x + p2 * Circle[i].y );
					glm::vec3 n1 = Circle[r].y * along
					             + Circle[r].x * ( p1 * Circle[i].x + p2 * Circle[i].y );
	
					attribs.emplace_back(a + yarn.radius * n1, n1, yarn.color);
					if (i == 0) attribs.emplace_back(attribs.back());
					attribs.emplace_back(a + yarn.radius * n0, n0, yarn.color);
				}
				attribs.emplace_back(attribs[attribs.size() - 2*Angles]);
				attribs.emplace_back(attribs[attribs.size() - 2*Angles]);
			}

		}


	}
}
#endif
