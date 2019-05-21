
#include <glm/glm.hpp>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <map>

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

int32_t needle_index(int32_t needle) {
	return needle * 3;
}

int32_t side_index(int32_t needle, Direction dir) {
	return needle_index(needle) + (dir == Left ? -1 : +1);
}

struct FaceEdge {
	FaceEdge(uint32_t face_ = -1U, uint32_t edge_ = -1U, uint32_t count_ = 0) : face(face_), edge(edge_), count(count_) { }
	uint32_t face;
	uint32_t edge;
	uint32_t count; //how many items (yarns/loops)
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

struct Gizmo {
	std::vector< uint32_t > faces;
	std::vector< Connection > connections;
};

struct Column {
	FaceEdge top_edge;
	float top_y = 0;
};

struct BedColumns : public std::unordered_map< int32_t, Column > {
	float depth = 0.0f;
};

struct Carrier {
	FaceEdge edge; //is a valid edge if carrier is in, otherwise is not
	int32_t index = 0;
	float depth = 0.0f;
};

//depths:
// --- back bed --- @ -1.0
// --- back hooks --- @ -0.75
//  [various crossings here]
// (carriers) @ [0.0 - 0.5]
// --- front hooks --- @ 0.75
// --- front bed ---   @ 1.0

struct Translator {
	BedColumns back_bed, back_sliders, front_sliders, front_bed;
	std::unordered_map< std::string, Carrier > carriers;

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
		std::unordered_set< std::string > seen;
		for (auto const &cn : cs) {
			if (!seen.insert(cn).second) {
				throw std::runtime_error("Carrier '" + cn + "' appears twice in carrier set.");
			}
			auto f = carriers.find(cn);
			if (f == carriers.end()) {
				throw std::runtime_error("Carrier '" + cn + "' does not appear in carriers list.");
			}
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

	//--- driver functions ---
	void knit(Direction dir, BedColumns &bed, int32_t needle, std::vector< Carrier * > cs) {
		assert(&bed == &back_bed || &bed == &front_bed);

		Gizmo gizmo;

		//FaceEdge edge = bring_carriers(cs, side_index(needle, (dir === Right ? Left : Right)), &gizmo);

		FaceEdge loop_in = bed[needle_index(needle)].top_edge;

		std::string L = std::to_string(loop_in.count);
		std::string Y = std::to_string(cs.size());

		const float H = 1.0f;

		if (dir == Right) {
			faces.emplace_back();
			Face &face = faces.back();
			face.type = "knit-to-right -l" + L + " +y" + Y + " +l" + Y + " -y" + Y;
			face.vertices = {
				glm::vec3(side_index(needle, Left), 0.0f, bed.depth),
				glm::vec3(side_index(needle, Right), 0.0f, bed.depth),
				glm::vec3(side_index(needle, Right), H, bed.depth),
				glm::vec3(side_index(needle, Left), H, bed.depth),
			};
			faces.emplace_back(face);

			gizmo.faces.emplace_back(faces.size() - 1);
			if (loop_in.is_valid()) {
				gizmo.connections.emplace_back(FaceEdge(gizmo.faces.back(), 0), loop_in);
			}

			float lift = 0.0f;
			{
				lift = std::max(lift, bed[needle_index(needle)].top_y);
			}
			//TODO: increase lift based on bottom/top edge conflicts in crossing area...

			//DEBUG: make sure stitches are separated by a little bit:
			lift += 0.05f;

			for (auto fi : gizmo.faces) {
				for (auto &v : faces[fi].vertices) {
					v.y += lift;
				}
			}

			//TODO: proper glue faces for yarn.

			//register top edge with column:
			bed[needle_index(needle)].top_edge = FaceEdge(gizmo.faces.back(), 2);
			bed[needle_index(needle)].top_y = faces[gizmo.faces.back()].vertices[2].y;
		} else { assert(dir == Left);
			
		}
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
		type_index.insert(std::make_pair(f.type, type_index.size()+1));
	}
	for (auto const &ti : type_index) {
		out << "L " << ti.first << "\n";
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
		out << "T " << type_index[f.type] << "\n";
	}

	return 0;
}
