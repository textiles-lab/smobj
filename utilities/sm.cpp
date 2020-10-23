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
#include <iterator>


//---------------------------------
//Machine

bool sm::MachineState::make(sm::Instr instr){

	if(instr.op == sm::Instr::Xfer || instr.op == sm::Instr::Split){
		if(racking != instr.rack()){
			//rack changed, might need to keep track of pass?
			racking = instr.rack();
			passes.emplace_back();
		}
	}
	if(passes.empty()) passes.emplace_back();
	auto &curr_pass = passes.back();

	// remove the loop(s) from the src
	// TODO maintain loops and their stacking order
	if(instr.op == sm::Instr::Xfer || instr.op == sm::Instr::Split){
		// move src to tgt
		assert(!instr.src.dontcare());
		assert(!instr.tgt.dontcare());
		if(active_needles[instr.src] == 0){
			std::cerr << "[xfer/split] from empty location. " << std::endl;
			std::cerr << instr.to_string() << std::endl;
			return false;
		}
		if(instr.face == -1U && active_needles[instr.tgt]){
			std::cerr << "[xfer*] to non-empty location. " << std::endl;
			std::cerr << instr.to_string() << std::endl;
			return false;
		}
		active_needles[instr.tgt] += active_needles[instr.src];
	}
	else if(!instr.tgt.dontcare()) {
		std::istringstream stream(instr.yarns);
		std::vector<std::string> tokens{std::istream_iterator<std::string>{stream},
                      std::istream_iterator<std::string>{}};
		int nl = tokens.size();
		active_needles[instr.tgt] += nl;
	}
	if(!instr.src.dontcare()){
		active_needles[instr.src] = 0;
	}
	if(!instr.tgt2.dontcare()){
		std::istringstream stream(instr.yarns);
		std::vector<std::string> tokens{std::istream_iterator<std::string>{stream},
                      std::istream_iterator<std::string>{}};
		int nl = tokens.size();
		active_needles[instr.tgt2] += nl;

	}

	if(instr.face == -1U){
		// Does xfer* instruction create yarn tangle  (not local)
		{
		}

		// Does xfer* instruction create loop tangle (not local)
		{
		}

		// Does xfer* violate slack on the bed
		{
		}
	} else {
		// Does instruction use locations not available to it ? (needs more info smobj+code)
	}

	curr_pass.emplace_back(instr);

	return true;
}

// are there any loops left on the machine?
bool sm::MachineState::empty(){
	for(auto const &pr : active_needles){
		if(pr.second != 0) return false;
	}
	return true;
}
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
	std::vector< uint32_t > sources;

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
		} else if (cmd == "N") {
			int32_t idx;
			if (!(str >> idx)) throw std::runtime_error("expecting line number after ln");
			if (idx < 0) throw std::runtime_error("invalid line number index '" + std::to_string(idx) + "'");
			sources.emplace_back(idx);
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
		} else if (cmd == "h"){
			char slash;
			int32_t face, edge;
			if( !(str >> face >> slash >> edge) || slash != '/') throw std::runtime_error("h line doesn't begin like 1/2.");
			Hint h;
			h.type = sm::Mesh::Hint::Resource;
			h.src = sm::Mesh::Hint::User;
			if( face < 1 || face > int32_t(mesh.faces.size())) throw std::runtime_error("face out of range in h line.");
			if( edge < 1 || edge > int32_t(mesh.faces[face-1].size())) throw std::runtime_error("edge out of range in h line.");

			h.lhs.face = face - 1;
			h.lhs.edge = edge -1;

			sm::BedNeedle bn;
			std::string next;
			if (str >> next){
				std::istringstream str2(next);
				if(std::isalpha(next[0])){
					char bed;
					str2 >> bed;
					bn.bed = bed;
				}
				float needle;
				if(str2 >> needle){
					bn.needle = needle;
					if(needle - bn.needle > 0.25) bn.nudge = 1;
					if(bn.needle - needle > 0.25) bn.nudge = -1;
				}
				h.rhs = bn;
			}

			char src;
			if(str >> src){
				if(src == 'u') h.src =  sm::Mesh::Hint::User;
				else if(src == 'i') h.src =  sm::Mesh::Hint::Inferred;
				else if(src == 'h') h.src =  sm::Mesh::Hint::Heuristic;
			}

			mesh.hints.emplace_back(h);

		} else if(cmd == "o"){

			Hint h;
			h.type = sm::Mesh::Hint::Order;
			h.src  = sm::Mesh::Hint::User;
			char slash;
			int32_t face1, face2, edge1, edge2;
			if(!(str >> face1 >> slash >> edge1) || slash != '/'){
				throw std::runtime_error("order hint not as expected f/i");
			}
			if(!(str >> face2 >> slash >> edge2) || slash != '/'){
				throw std::runtime_error("order hint not as expected f/i");
			}
			h.lhs.face = face1-1; h.lhs.edge = edge1-1;
			FaceEdge rhs; rhs.face = face2-1; rhs.edge  = edge2-1;
			h.rhs = rhs;

			char src;
			if(str >> src){
				if(src == 'u') h.src =  sm::Mesh::Hint::User;
				else if(src == 'i') h.src =  sm::Mesh::Hint::Inferred;
				else if(src == 'h') h.src =  sm::Mesh::Hint::Heuristic;
			}
			std::cout << "Hint: " << h.to_string() << std::endl;
			mesh.hints.emplace_back(h);

		} else if(cmd == "t"){

			Hint h;
			h.type = sm::Mesh::Hint::Variant;
			h.src  = sm::Mesh::Hint::User;
			int32_t face;
			if(!(str >> face)){
				throw std::runtime_error("variant hint does not have a face");
			}
			h.lhs.face = face-1;
			h.lhs.edge = -1;
			std::string var;
			if(!(str >> var)) throw std::runtime_error("variant hint does not have a variant name");
			h.rhs = var;
			char src;
			if(str >> src){
				if(src == 'u') h.src =  sm::Mesh::Hint::User;
				else if(src == 'i') h.src =  sm::Mesh::Hint::Inferred;
				else if(src == 'h') h.src =  sm::Mesh::Hint::Heuristic;
			}
			mesh.hints.emplace_back(h);
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
			if (edge < 1 || edge > int32_t(mesh.faces[face-1].size())) throw std::runtime_error("c line with invalid edge '" + std::to_string(edge) + "' (face '" + std::to_string(face) + "' has " + std::to_string(mesh.faces[face-1].size()) + " edges)");
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
		}
		else if(cmd == "I"){
			uint32_t f, i;
			char slash;
			if(!(str >> f >> slash >> i) || slash != '/'){
				throw std::runtime_error("Instruction order is not in expected format face/instruction");
			}
			mesh.total_order.emplace_back(std::make_pair(f-1, i-1));
		}
		else if(cmd == "vn"){
			//ignore
		}
		else if(cmd == "vt"){
			//ignore
		}
		else if(cmd == "usemtl"){
			//ignore
		}
		else if(cmd == "mtllib"){
			//ignore
		}
		else {
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


	if (types.empty()) {
		types.resize(mesh.faces.size(),0); //empty
		std::unordered_map<std::string, uint32_t> key_type;
		// add an empty string to the library
		assert(library.empty() && "Library must be fully specified or empty" );
		for(auto &f : mesh.faces){
			Signature empty_signature;
			empty_signature.name = "empty"+std::to_string(f.size());
			for(uint32_t i = 0; i < f.size(); ++i){
				empty_signature.edges.emplace_back("-");
			}
			if(!key_type.count(empty_signature.key())){

				library.emplace_back(empty_signature);
				mesh.library.emplace_back(empty_signature.key());
				key_type[empty_signature.key()] = key_type.size();
				std::cout << "Added empty string to library (" << library.size() << ")" << std::endl;
			}
			types[&f-&mesh.faces[0]] = key_type[empty_signature.key()];
		}
	}
	if (types.size() != mesh.faces.size()) throw std::runtime_error("should be a 'T' for every face.");
	if (sources.empty()) sources.resize(mesh.faces.size(), 0);
	if (sources.size() != mesh.faces.size()) throw std::runtime_error("should be a 'ln' for every face or for no faces.");

	{
		// do all the smobj checks
		for (uint32_t i = 0; i < types.size(); ++i) {
			if (mesh.faces[i].size() != library[types[i]].edges.size()) {
				throw std::runtime_error("face/library edge count mismatch.");
			}
			mesh.faces[i].type = types[i];
			mesh.faces[i].source = sources[i];
		}
	}
	std::vector< glm::vec3 > normals(mesh.vertices.size(), glm::vec3(std::numeric_limits< float >::quiet_NaN()));

	for (auto const &face : mesh.faces) {
		for (uint32_t vi = 0; vi < face.size(); ++vi) {
			glm::vec3 const &v = mesh.vertices[face[vi]];
			for (uint32_t i = 1; i + 1 < face.size(); ++i) {
				glm::vec3 const &a = mesh.vertices[face[(vi+i)%face.size()]];
				glm::vec3 const &b = mesh.vertices[face[(vi+i+1)%face.size()]];
				//lazy area-weighted normals:
				//note: doing this "init-to-zero" thing here so that stray vertices can be recognized and assigned arbitrary normals:
				if (!(normals[face[vi]].x == normals[face[vi]].x)) {
					normals[face[vi]] = glm::vec3(0.0f);
				}
				normals[face[vi]] += glm::cross(b-a, v-a);
			}
		}
	}

	uint32_t stray_vertices = 0;
	for (auto &n : normals) {
		if (!(n.x == n.x)) {
			++stray_vertices;
			n = glm::vec3(0.0f, 0.0f, 1.0f); //arbitrary
		} else if (n == glm::vec3(0.0f)) {
			throw std::runtime_error("Vertex " + std::to_string(&n-&normals[0]+1) + " has sum-zero surrounding face area.");
		}
		n = normalize(n);
	}

	if (stray_vertices > 0) {
		std::cout << "NOTE: there are " << stray_vertices << " stray (unreferenced) vertices; they have been assigned +z normals." << std::endl;
	}

	return mesh;
}

void sm::Mesh::save_instructions(std::string const &filename, sm::Library const &face_library) const{

	std::ofstream out(filename, std::ios::binary);
	//slow but for now
	auto out_loops = [&](sm::Mesh::Face f)->uint32_t{
		for(auto l : face_library.faces){
			if(l.key() == library[f.type]){
				uint32_t c = 0;
				for(auto e : l.edges){
					// yarn edges are getting over counted, but they will get cancelled:
					if(e.direction == sm::Library::Face::Edge::Out ){
						c++;
					}
				}
				return c;
			}
		}
		return 1;
	};
	auto in_loops = [&](sm::Mesh::Face f)->uint32_t{
		for(auto l : face_library.faces){
			if(l.key() == library[f.type]){
				uint32_t c = 0;
				for(auto e : l.edges){
					if(e.direction == sm::Library::Face::Edge::In){
						c++;
					}
				}
				return c;
			}
		}
		return 1;
	};
	std::string prev_face = "--- pattern start ---- (" + std::to_string(faces.size()) + " stitch-faces )";
	uint32_t count = 1;
	int active_loops = 0;
	// requires : faces appear in construction sequence
	uint32_t STEP  = 10;
	for (auto const &f : faces) {
		if(f.type < this->library.size()){
			// go through all connections that have this face, adding the +edges, removing the -edges
			// slow, but okay
			active_loops += -in_loops(f) + out_loops(f);
			if(this->library[f.type] == prev_face && (&f-&faces[0])%STEP != 0){
				count++;
			}
			else{
				out << prev_face << " ... " << count << " times. \n";
				prev_face = this->library[f.type];
				count = 1;
				if( f.size() != 4 || (&f - &faces[0])%STEP == 0 ){
					out << "\t\tactive loops: " << active_loops << "\n";
				}
			}
		}
	}
	out << prev_face << " ...  " << count << " times.\n";
	out << "---- pattern end ----\n";
	std::cout << "active_loops " << active_loops << std::endl;
}

void sm::Mesh::save(std::string const &filename) const {
	std::ofstream out(filename, std::ios::binary);

	for (auto const &l : library) {
		out << "L " << l << "\n";
	}

	for (auto const &v : vertices) {
		out << "v " << v.x << " " << v.y << " " << v.z << "\n";
	}

	for (auto const &f : faces) {
		out << "f";
		for (auto const i : f) {
			out << " " << (i+1);
		}
		out << "\n";
		out << "T " << (f.type+1) << "\n";
		out << "N " << (f.source) << "\n";
	}

	for (auto const &c : connections) {
		out << "e " << (c.a.face+1) << "/" << (c.a.edge+1) << " " << (c.b.face+1) << "/" << (c.flip ? "-" : "") << (c.b.edge+1) << "\n";
	}
	for (auto const &h : hints) {
		if(h.type == sm::Mesh::Hint::Resource){
		out << "h " << (h.lhs.face + 1) << "/" << (h.lhs.edge+1) << " ";
		auto bn = std::get<sm::BedNeedle>(h.rhs);
		out  << bn.bed;
		out  << bn.needle + 0.5f*bn.nudge;
		out  << ' ' << (char)h.src;
		out <<"\n";
		}
		else if(h.type == sm::Mesh::Hint::Order){
		out << "o " << (h.lhs.face + 1) << "/" << (h.lhs.edge+1) << " ";
		auto fe = std::get<sm::Mesh::FaceEdge>(h.rhs);
		out  << (fe.face+1) << "/" << (fe.edge +1) << " ";
		out  << (char)h.src;
		out <<"\n";
		}
		else if(h.type == sm::Mesh::Hint::Variant){
		out << "t " << (h.lhs.face + 1) << " " << std::get<std::string>(h.rhs) << " ";
		out  << (char)h.src;
		out <<"\n";
		}
	}
	for (auto const &u : units) {
		out << "U " << u.name << " " << u.length << "\n";
	}
	std::unordered_map< std::string, std::string > following_length;
	for (auto const &c : checkpoints) {
		auto &l = following_length[std::to_string(c.face+1) + "/" + std::to_string(c.edge+1) + "/" + std::to_string(c.crossing+1)];
		if (l != "") l += " ";
		l += std::to_string(c.length) + " " + std::to_string(c.unit+1);
	}
	for (auto const &fl : following_length) {
		out << "c " << fl.first << " " << fl.second << "\n";
	}
}

void sm::Mesh::rip(uint32_t fixed_id){

	if(this->faces.empty()) return;
	bool do_rotating = true;
	if(fixed_id > this->faces.size()){
		do_rotating = false;
	}
	std::vector<glm::vec3> updated_vertices;
	float avg_len = 0;
	uint32_t count  = 0;
	std::vector<glm::vec3> normals;
	for(auto &f : this->faces){
		glm::vec3 n = glm::vec3(0,0,0);
		for(uint32_t i = 0; i < f.size(); ++i){
			avg_len += glm::length(this->vertices[f[i]] - this->vertices[f[(i+1)%f.size()]]);
			n += glm::cross(this->vertices[f[i]] - this->vertices[f[(i+1)%f.size()]], -this->vertices[f[(f.size()+i-1)%f.size()]] + this->vertices[f[i]]);
		}
		normals.emplace_back(glm::normalize(n));
		count+= f.size();
	}

	avg_len /= count;

	std::map<std::pair<uint32_t, uint32_t>, uint32_t> edge_to_face;
	std::unordered_set<uint32_t>  fixed_set;

	for(auto &f : this->faces){
		for(uint32_t i = 0; i < f.size(); ++i){
			edge_to_face[std::make_pair(f[i], f[(i+1)%f.size()])] = &f - &this->faces[0];
		}
	}
	// figure out what the appropriate face for rotation should be given a face with known rotation
	// face 0 for default otherwise pass from interface
	// need some rules for basic face types, this is for quads..
	if(do_rotating){
		auto index_in_face  = [](sm::Mesh::Face f, uint32_t v)->uint32_t{
			for(uint32_t i = 0; i < f.size(); ++i){
				if(f[i] == v) return i;
			}
			return -1U;
		};

		auto set_first = [](sm::Mesh::Face &f, uint32_t v){
			uint32_t c = 0;
			while(f[0] != v){
				std::rotate(f.begin(), f.begin()+1, f.end());
				c++;
				if( c > f.size()+1) break; // infinite loop
			}
			if(c > f.size()) assert(false);

		};

		(void)index_in_face;
		fixed_set.insert(fixed_id);

		while(fixed_set.size() != this->faces.size()){
			auto old_size = fixed_set.size();
			for(auto id : fixed_set){
				auto ff = this->faces[id];
				if(ff.size() == 4){ // some basic rules for other shapes?
					auto left = edge_to_face[ std::make_pair(ff[0], ff[3])];
					auto down = edge_to_face[ std::make_pair(ff[1], ff[0])];
					auto right = edge_to_face[ std::make_pair(ff[2], ff[1])];
					auto up = edge_to_face[ std::make_pair(ff[3], ff[2])];
					auto lf = this->faces[left];
					auto df = this->faces[down];
					auto rf = this->faces[right];
					auto uf = this->faces[up];

					auto l = index_in_face(this->faces[left],ff[0]);
					if(l != -1U){
						set_first(this->faces[left], lf[(lf.size()+l-1)%lf.size()]);

						if(!fixed_set.count(left))
							fixed_set.insert(left);
					}
					auto d = index_in_face(this->faces[down],ff[0]);
					if(d != -1U){
						set_first(this->faces[down], df[(df.size()+d+1)%df.size()]);
						if(!fixed_set.count(right))
							fixed_set.insert(right);
					}

					auto r = index_in_face(this->faces[right],ff[1]);
					if(r != -1U){
						set_first(this->faces[right], rf[(rf.size()+r)%rf.size()]);
						if(!fixed_set.count(down))
							fixed_set.insert(down);
					}

					auto u = index_in_face(this->faces[up],ff[3]);
					if(u != -1U){
						set_first(this->faces[up], uf[(uf.size()+u)%uf.size()]);
						if(!fixed_set.count(up))
							fixed_set.insert(up);
					}

				}
			}
			if(old_size == fixed_set.size()) break;
		}


	}
	std::map<std::pair<uint32_t, uint32_t>, uint32_t> old_to_new_facevertex;
	for(auto &f : this->faces){
		glm::vec3 n = normals[&f - &this->faces[0]];
		auto ff = f;
		for(auto &v: f){
			float offset_eps = 0.01f;
			glm::vec3 e = -this->vertices[v] + this->vertices[ ff[(&v-&f[0]+1)%f.size()]];
			glm::vec3 offset_a = glm::normalize(glm::cross(n,e))*offset_eps*avg_len;
			e = -this->vertices[v] + this->vertices[ ff[(f.size()+&v-&f[0]-1)%f.size()]];
			glm::vec3 offset_b = glm::normalize(glm::cross(n,e))*offset_eps*avg_len;
			assert(offset_a == offset_a);
			assert(offset_b == offset_b);
			updated_vertices.emplace_back(this->vertices[v]+offset_a - offset_b);
			old_to_new_facevertex[std::make_pair(&f- &this->faces[0],v)] = (updated_vertices.size()-1);
			v = updated_vertices.size()-1;


		}
	}
	this->connections.clear(); // clear old connections?
	this->vertices = updated_vertices;

	// could have just used find but..
	auto edge_index = [](sm::Mesh::Face f, uint32_t v)->uint32_t{
		for(uint32_t i = 0; i < f.size(); ++i){
			if(f[i] == v) return i;
		}
		assert(false);
		return -1U;
	};
	std::set<std::pair<uint32_t, uint32_t>> done;
	// make new connections
	for(auto pr : edge_to_face){
		if(done.count(pr.first)) continue;
		auto opp = std::make_pair(pr.first.second, pr.first.first);
		if(done.count(opp)) continue;
		if(edge_to_face.count(opp)){

			Connection c;
			c.a.face = pr.second;
			if(!old_to_new_facevertex.count(std::make_pair(c.a.face, pr.first.first))) continue;
			uint32_t v1 = old_to_new_facevertex[std::make_pair(c.a.face, pr.first.first)];
			c.a.edge = edge_index(this->faces[c.a.face], v1);
			c.b.face = edge_to_face[opp];
			if(!old_to_new_facevertex.count(std::make_pair(c.b.face, pr.first.second))) continue;
			uint32_t v2 = old_to_new_facevertex[std::make_pair(c.b.face, pr.first.second)];
			c.b.edge = edge_index(this->faces[c.b.face], v2);
			c.flip = true;

			done.insert(pr.first);
			done.insert(opp);

			this->connections.emplace_back(c);
		}

	}

}


void sm::Mesh::remove_inferred_hints(uint32_t from){ // default-1U
	hints.erase(std::remove_if(hints.begin(), hints.end(),[from](sm::Mesh::Hint const &h)->bool{
				return (h.src == sm::Mesh::Hint::Inferred && (from == -1U || h.inferred_from == from));
				}), hints.end());
}

void sm::Mesh::remove_heuristic_hints(uint32_t from){ // default-1U

	hints.erase(std::remove_if(hints.begin(), hints.end(),[from](sm::Mesh::Hint const &h)->bool{
				return (h.src == sm::Mesh::Hint::Heuristic && (from == -1U || h.inferred_from == from));
				}), hints.end());
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
		} else if (tok == "derive") {
			current = nullptr; //clear current face

			std::string name;
			if (!(str >> name)) throw std::runtime_error(line_info() + "Failed to read name field from face line");
			sm::Library::Face::Derive derive;

			std::string by_or_from;
			{
				std::string type;
				while (str >> type) {
					if (type == "by" || type == "from") {
						by_or_from = type;
						break;
					}
					derive.expect_key += ' ';
					derive.expect_key += type;
				}
			}

			if (derive.expect_key != "") {
				derive.expect_key = name + derive.expect_key;
			} else {
				std::cerr << line_info() + " WARNING: derive line using just face name may be brittle; consider adding edge types." << std::endl;
			}

			if (!(by_or_from == "by" || by_or_from == "from")) throw std::runtime_error(line_info() + "derive line should have face name and edge labels followed by 'by' or 'from'");
			if (by_or_from == "by") {
				std::string tok;
				while (str >> tok) {
					if (tok == "from") {
						break;
					} else if (tok == "mirror-x") {
						if (derive.by & sm::Library::Face::Derive::MirrorXBit) throw std::runtime_error(line_info() + "derive shouldn't mention mirror-x twice.");
						derive.by |= sm::Library::Face::Derive::MirrorXBit;
					} else if (tok == "mirror-z") {
						if (derive.by & sm::Library::Face::Derive::MirrorZBit) throw std::runtime_error(line_info() + "derive shouldn't mention mirror-z twice.");
						derive.by |= sm::Library::Face::Derive::MirrorZBit;
					} else if (tok == "reverse-yarn") {
						if (derive.by & sm::Library::Face::Derive::ReverseYarnBit) throw std::runtime_error(line_info() + "derive shouldn't mention reverse-yarn twice.");
						derive.by |= sm::Library::Face::Derive::ReverseYarnBit;
					} else {
						throw std::runtime_error(line_info() + "unknown derivation operation '" + tok + "'");
					}
				}
				if (tok != "from") throw std::runtime_error(line_info() + "derive line should have 'from' after operations");
			} else {
				assert(by_or_from == "from");
			}
			{ //read source:
				std::string tok;
				while (str >> tok) {
					if (derive.from != "") derive.from += ' ';
					derive.from += tok;
				}
			}
			library.faces.emplace_back();
			library.faces.back().name = name;
			library.faces.back().derive = derive;

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
				while (str.peek() != EOF && isspace(str.peek())) str.get();
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

	std::map< std::string, Library::Face const * > keys;

	std::vector< Library::Face * > pending;
	pending.reserve(library.faces.size());

	for (auto &face : library.faces) {

		//defer running derivation logic:
		if (face.derive.from != "") {
			pending.emplace_back(&face);
			continue;
		}

		auto ret = keys.emplace(face.key(), &face);
		if (!ret.second) throw std::runtime_error("Duplicate face signature: '" + ret.first->first + "'");
	}

	while (!pending.empty()) {
		std::vector< Library::Face * > next_pending;
		for (auto fp : pending) {

			auto source = keys.find(fp->derive.from);
			if (source == keys.end()) {
				next_pending.emplace_back(fp);
				continue;
			}

			derive_face(*source->second, fp->derive.by, fp);

			if (fp->derive.expect_key != "" && fp->derive.expect_key != fp->key()) {
				throw std::runtime_error("Deriving face with expected key '" + fp->derive.expect_key + "', but got key '" + fp->key() + "'.");
			}

			//std::cout << "Derived '" << fp->key() << "' from '" << fp->derive.from << "'" << std::endl; //DEBUG

			auto ret = keys.emplace(fp->key(), fp);
			if (!ret.second) throw std::runtime_error("Duplicate face signature: '" + ret.first->first + "'");
		}

		if (next_pending.size() == pending.size()) {
			//failed to make progress:
			std::string message = "";
			for (auto fp : pending) {
				if (message != "") message += "; ";
				message += "derived face '" + fp->name + "' missing source '" + fp->derive.from + "'";
			}
			throw std::runtime_error(message);
		}

		pending = std::move(next_pending);
	}

	return library;
}

void sm::Library::save(std::string const &filename) const{
	std::ofstream out(filename, std::ios::binary);
	for (auto const &face : faces) {
		if (face.derive.from != "") {
			out << "derive ";
			if (face.derive.expect_key != "") out << face.derive.expect_key;
			else out << face.name;
			if (face.derive.by) {
				out << " by";
				if (face.derive.by & sm::Library::Face::Derive::MirrorXBit) out << " mirror-x";
				if (face.derive.by & sm::Library::Face::Derive::MirrorZBit) out << " mirror-z";
				if (face.derive.by & sm::Library::Face::Derive::ReverseYarnBit) out << " reverse-yarn";
			}
			out << " from " << face.derive.from << '\n';
			continue;
		}
		out << "face " << face.name << '\n';
		for (auto const &edge : face.edges) {
			out << "\tedge ";
			out << '(' << edge.vertex.x << ',' << edge.vertex.y << ')';
			out << ' ';
			if      (edge.direction == sm::Library::Face::Edge::In) out << "-";
			else if (edge.direction == sm::Library::Face::Edge::Out) out << "+";
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

//------------------------------------------------

// Code

sm::Code sm::Code::load(std::string const &filename) {
	sm::Code code_library;
	std::ifstream in(filename, std::ios::binary);
	sm::Code::Face *current = nullptr;
	std::vector< sm::Instr > *current_instrs = nullptr;
	{

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
			std::string tok;
			if (!(str >> tok)) continue; //ignore blank lines
			if (tok == "face") {
				std::string name;
				if (!(str >> name)) throw std::runtime_error(line_info() + "Failed to read name field from face line");
				std::string temp;
				if (str >> temp) throw std::runtime_error(line_info() + "Trailing junk (" + temp + "...) in face line");

				code_library.faces.emplace_back();
				code_library.faces.back().name = name;
				current = &code_library.faces.back();
				current_instrs = nullptr;
			} else if (tok == "variant") {
				if (!current) throw std::runtime_error(line_info() + "variant line without face line");
				std::string ver;
				if(!(str >> ver)) throw std::runtime_error(line_info() + "Failed to read name in variant line");
				current->variant = ver;
			} else if (tok == "edge") {
				if (!current) throw std::runtime_error(line_info() + "edge line without face line");
				Code::Face::Edge edge;
				std::string type;
				if (!(str >> type)) throw std::runtime_error(line_info() + "Failed to read [+-]?type in edge line");
				if (type[0] == '+') {
					edge.direction = sm::Code::Face::Edge::Out;
					edge.type = type.substr(1);
				} else if (type[0] == '-') {
					edge.direction = sm::Code::Face::Edge::In;
					edge.type = type.substr(1);
				} else {
					edge.direction = sm::Code::Face::Edge::Any;
					edge.type = type;
				}
				if(edge.type != "x"){
					char bed;
					if (!(str >> bed)) throw std::runtime_error(line_info() + "edge without bed");
					edge.bn.bed = bed;
					float needle;
					if (!(str >> needle) && !edge.bn.dontcare()) throw std::runtime_error(line_info() + "edge without needle");
					if(!edge.bn.dontcare()){
						edge.bn.needle = needle;
						edge.bn.nudge = (needle - edge.bn.needle)*2;
					}
					std::string y = "";
					while (str >> y ){
						if (!edge.yarns.empty()) edge.yarns += " ";
						edge.yarns += y;
					}
				} else {
					std::string temp;
					if( str >> temp ) throw std::runtime_error(line_info() + "Trailing junk (" + temp + "...) with dontcare edge line");
				}
				// any trailing junk will be treated as yarn...
				current->edges.emplace_back(edge);
			} else if (tok == "code") {
				if (!current) throw std::runtime_error(line_info() + "code line without face line");
				current_instrs = &(current->instrs);
			} else if (tok == ";;Carriers:") {
				if (!current_instrs) throw std::runtime_error(line_info() + "Carrier header line without code line");
				//todo save local carrier order
				std::string c;
				while(str >> c){
					current->carriers.emplace_back(c);
				}
			} else if (tok == "knit") {
				if (!current_instrs) throw std::runtime_error(line_info() + "instr line without code line");
				current_instrs->emplace_back();
				auto &ins = current_instrs->back();
				ins.op = sm::Instr::Knit;
				char direction;
				if(!(str >> direction)){
					throw std::runtime_error(line_info() + "knit line without direction");
				}
				if( direction == '+'){
					ins.direction = sm::Instr::Right;
				} else if( direction == '-'){
					ins.direction = sm::Instr::Left;
				} else{
					throw std::runtime_error(line_info() + "knit without +/- direction");
				}
				char bed;
				int needle;
				if(!(str >> bed)){
					throw std::runtime_error(line_info() + "knit line without bed");
				}
				if(!(str >> needle)){
					throw std::runtime_error(line_info() + "knit line without needle");
				}
				ins.src.bed = bed; ins.src.needle = needle;
				ins.tgt = ins.src;
				std::string y;
				while(str >> y){
					if(!ins.yarns.empty()) ins.yarns += " ";
					ins.yarns += y;
				}
				//ins.rack = rack_from(ins.src, ins.tgt);
			} else if (tok == "tuck") {
				if (!current_instrs) throw std::runtime_error(line_info() + "instr line without code line");
				current_instrs->emplace_back();
				auto &ins = current_instrs->back();
				ins.op = sm::Instr::Tuck;
				char direction;
				if(!(str >> direction)){
					throw std::runtime_error(line_info() + "tuck line without direction");
				}
				if( direction == '+'){
					ins.direction = sm::Instr::Right;
				} else if( direction == '-'){
					ins.direction = sm::Instr::Left;
				} else{
					throw std::runtime_error(line_info() + "tuck without +/- direction");
				}
				char bed;
				int needle;
				if(!(str >> bed)){
					throw std::runtime_error(line_info() + "tuck line without bed");
				}
				if(!(str >> needle)){
					throw std::runtime_error(line_info() + "tuck line without needle");
				}
				ins.tgt.bed = bed; ins.tgt.needle = needle;
				std::string y;
				while(str >> y){
					if(!ins.yarns.empty()) ins.yarns += " ";
					ins.yarns += y;
				}
				//ins.rack = rack_from(ins.src, ins.tgt);
			} else if (tok == "miss") {
				if (!current_instrs) throw std::runtime_error(line_info() + "instr line without code line");
				current_instrs->emplace_back();
				auto &ins = current_instrs->back();
				ins.op = sm::Instr::Miss;
				char direction;
				if(!(str >> direction)){
					throw std::runtime_error(line_info() + "miss line without direction");
				}
				if( direction == '+'){
					ins.direction = sm::Instr::Right;
				} else if( direction == '-'){
					ins.direction = sm::Instr::Left;
				} else{
					throw std::runtime_error(line_info() + "miss without +/- direction");
				}
				char bed;
				int needle;
				if(!(str >> bed)){
					throw std::runtime_error(line_info() + "miss line without bed");
				}
				if(!(str >> needle)){
					throw std::runtime_error(line_info() + "miss line without needle");
				}
				ins.src.bed = bed; ins.src.needle = needle;
				ins.tgt = ins.src;
				std::string y;
				while(str >> y){
					if(!ins.yarns.empty()) ins.yarns += " ";
					ins.yarns += y;
				}
				//ins.rack = rack_from(ins.src, ins.tgt);
			} else if (tok == "xfer") {
				if (!current_instrs) throw std::runtime_error(line_info() + "instr line without code line");
				current_instrs->emplace_back();
				auto &ins = current_instrs->back();
				ins.op = sm::Instr::Xfer;
				char from_bed;
				int from_needle;
				char direction;
				if(!(str >> direction)) throw std::runtime_error(line_info() + "xfer line without direction or bed");
				if (direction == '+' || direction == '-'){
					ins.direction = (direction == '+' ? sm::Instr::Right : sm::Instr::Left);
					if(!(str >> from_bed)){
						throw std::runtime_error(line_info() + "xfer line without (from)bed");
					}
				}
				else{
					from_bed = direction; // no direction was specified so that must have been the from bed
				}

				if(!(str >> from_needle)){
					throw std::runtime_error(line_info() + "xfer line without (from)needle");
				}
				ins.src.bed = from_bed; ins.src.needle = from_needle;
				char to_bed;
				int to_needle;
				if(!(str >> to_bed)){
					throw std::runtime_error(line_info() + "xfer line without (to)bed");
				}
				if(!(str >> to_needle)){
					throw std::runtime_error(line_info() + "xfer line without (to)needle");
				}
				ins.tgt.bed = to_bed;
				ins.tgt.needle = to_needle;
				//ins.rack = rack_from(ins.src, ins.tgt);
				std::string temp;
				if (str >> temp) throw std::runtime_error(line_info() + "Trailing junk (" + temp + "...) in xfer line");
			} else if (tok == "split") {
				if (!current_instrs) throw std::runtime_error(line_info() + "instr line without code line");
				current_instrs->emplace_back();
				auto &ins = current_instrs->back();
				ins.op = sm::Instr::Split;
				char direction;
				if(!(str >> direction)){
					throw std::runtime_error(line_info() + "split line without direction");
				}
				if( direction == '+'){
					ins.direction = sm::Instr::Right;
				} else if( direction == '-'){
					ins.direction = sm::Instr::Left;
				} else{
					throw std::runtime_error(line_info() + "split without +/- direction");
				}
				char bed;
				int needle;
				if(!(str >> bed)){
					throw std::runtime_error(line_info() + "split line without bed");
				}
				if(!(str >> needle)){
					throw std::runtime_error(line_info() + "split line without needle");
				}
				ins.src.bed = bed; ins.src.needle = needle;

				char to_bed;
				int to_needle;
				if(!(str >> to_bed)) throw std::runtime_error(line_info() + "split without target bed");
				if(!(str >> to_needle))  throw std::runtime_error(line_info() + "split without target needle");
				ins.tgt.bed = to_bed; ins.tgt.needle = to_needle;
				ins.tgt2 = ins.src;
				std::string y;
				while(str >> y){
					if(!ins.yarns.empty()) ins.yarns += " ";
					ins.yarns += y;
				}
				//ins.rack = rack_from(ins.src, ins.tgt);
			} else if (tok == "drop"){
				if (!current_instrs) throw std::runtime_error(line_info() + "instr line without code line");
				current_instrs->emplace_back();
				auto &ins = current_instrs->back();
				ins.op = sm::Instr::Drop;
				char bed; int needle;

				if(!(str >> bed)) throw std::runtime_error(line_info() + "drop without bed");
				if(!(str >> needle))  throw std::runtime_error(line_info() + "drop without needle");
				ins.src.bed = bed;
				ins.src.needle = needle;
				std::string temp;
				if (str >> temp) throw std::runtime_error(line_info() + "Trailing junk (" + temp + "...) in drop line");
			} else if (tok == "in") {
				if (!current_instrs) throw std::runtime_error(line_info() + "instr line without code line");
				current_instrs->emplace_back();
				auto &ins = current_instrs->back();
				ins.op = sm::Instr::In;

				std::string y;
				while(str >> y){
					if(!ins.yarns.empty()) ins.yarns += " ";
					ins.yarns += y;
				}
			}  else if (tok == "out") {
				if (!current_instrs) throw std::runtime_error(line_info() + "instr line without code line");
				current_instrs->emplace_back();
				auto &ins = current_instrs->back();
				ins.op = sm::Instr::Out;

				std::string y;
				while(str >> y){
					if(!ins.yarns.empty()) ins.yarns += " ";
					ins.yarns += y;
				}
			}
			else {
				throw std::runtime_error(line_info() + "Unrecognized line-start token '" + tok + "'");
			}
		}

		std::map< std::string, Code::Face const * > keys;

		for (auto &face : code_library.faces) {
			auto ret = keys.emplace(face.key(), &face);
			if (!ret.second) throw std::runtime_error("Duplicate face signature in code library: '" + ret.first->first + "'");
		}

		// verify that instructions are valid (does this need to ssa?)
		// todo: also verify yarn positions
		for (auto const &face : code_library.faces){
			std::set<sm::BedNeedle> temporaries, incoming, outgoing;
			std::set<std::string> carriers; // todo split into incoming and outgoing carriers
			for (auto const &edge : face.edges){
				if(edge.direction == Code::Face::Edge::In && edge.type[0] == 'l' && !edge.bn.dontcare()){
					incoming.insert(edge.bn);
				}
				else if(edge.direction == Code::Face::Edge::Out &&  edge.type[0] == 'l' && !edge.bn.dontcare()){
					outgoing.insert(edge.bn);
				}
				if(!edge.yarns.empty()) {
					carriers.insert(edge.yarns);
					// also include split up carriers
					std::istringstream str(edge.yarns);
					std::string yarn;
					while(str >> yarn) {
						carriers.insert(yarn);
					}
				}
			}

			for(auto const &ins : face.instrs){
				if(!ins.src.dontcare() && !(incoming.count(ins.src) || temporaries.count(ins.src))){
					throw std::runtime_error("src bed-needle " + ins.src.to_string() + " is not an input resource:" + face.key());
				}
				temporaries.erase(ins.src);
				if(!ins.tgt.dontcare()) {
					temporaries.insert(ins.tgt);
				}
				if(!ins.tgt2.dontcare()) {
					temporaries.insert(ins.tgt2);
				}
				if(!ins.yarns.empty()){
					std::istringstream str(ins.yarns);
					std::string yarn;
					while (str >> yarn){
						if (!carriers.count(yarn)) throw std::runtime_error("yarn [" + yarn + " of " + ins.yarns + "] not specified on edge label by code:" + face.key());
					}
				}
			}

			for(auto bn : outgoing){
				if(!temporaries.count(bn))
					throw std::runtime_error("tgt bed-needle [" + bn.to_string() + "] is not produced by code:" + face.key());
				temporaries.erase(bn);
			}


			if(!temporaries.empty()){
				std::string err = "";
				for(auto bn : temporaries) err += " " + bn.to_string();
				throw std::runtime_error("Temporary locations [ " + err + " ] are not cleared in " + face.key());
			}


		}
	}
	return code_library;
}

void sm::Code::save(std::string const &filename) const {
	std::ofstream out(filename);
	for(auto const &face : faces){
		out << "face " << face.name << '\n';
		if(face.variant != ""){
			out << "\tvariant " << face.variant << '\n';
		}
		for(auto const &edge : face.edges) {
			out << "\tedge ";
			if (edge.direction == sm::Code::Face::Edge::In) out << "-";
			else if (edge.direction == sm::Code::Face::Edge::Out) out << "+";
			out << edge.type << ' ';
			if(edge.type != "x"){ // convention none type ==> no bed-needle specified
				if(edge.bn.bed != 'x')
					out << edge.bn.bed <<  (edge.bn.needle + edge.bn.nudge * 0.5) << ' ';
				else
					out << 'x' << ' ';
			}
			out << edge.yarns;
			out << '\n';
		}
		// don't use instr.knitout_string() here, it uses translation and  also inserts racking
		out << "\tcode \n";
		out << "\t\t;;Carriers: "; for(auto c : face.carriers) out << c << " "; out << '\n';
		for(auto const &instr : face.instrs) {
			out <<"\t\t";
			switch(instr.op){
				case sm::Instr::Knit:
					out <<"knit ";
					out << (char)instr.direction << ' ';
					out << instr.tgt.bed << instr.tgt.needle << ' ';
					out << instr.yarns << '\n';
					break;
				case sm::Instr::Tuck:
					out <<"tuck ";
					out << (char)instr.direction << ' ';
					out << instr.tgt.bed << instr.tgt.needle << ' ';
					out << instr.yarns << '\n';
					break;
				case sm::Instr::Miss:
					out <<"miss ";
					out << (char)instr.direction << ' ';
					out << instr.src.bed << instr.src.needle << ' ';
					out << instr.yarns << '\n';
					break;
				case sm::Instr::Xfer:
					out <<"xfer ";
					if(instr.direction != sm::Instr::None){
						out << (char)instr.direction << ' ';
					}
					out << instr.src.bed << instr.src.needle << ' ';
					out << instr.tgt.bed << instr.tgt.needle << ' ';
					out << '\n';
					break;
				case sm::Instr::Split:
					out <<"split ";
					out <<(char)instr.direction << ' ';
					out << instr.src.bed << instr.src.needle << ' ';
					out << instr.tgt.bed << instr.tgt.needle << ' ';
					out << instr.yarns << '\n';
					break;
				case sm::Instr::Drop:
					out <<"drop ";
					out << instr.src.bed << instr.src.needle << '\n';
					break;
				case sm::Instr::In:
					out <<"in ";
					out << instr.yarns << '\n';
					break;
				case sm::Instr::Out:
					out <<"out ";
					out << instr.yarns << '\n';
					break;
				default:
					assert(false && "unknown knitout operation");
			}
		}
		out << '\n';
	}
}

//------------------------------------------------


namespace {
	//internal structures used for yarns file format:
	struct YarnInfo {
		uint32_t point_begin;
		uint32_t point_end;
		float radius;
		glm::u8vec4 color;
	};
	static_assert(sizeof(YarnInfo) == 16, "YarnInfo is packed");
	struct UnitInfo {
		uint32_t name_begin;
		uint32_t name_end;
		float length;
	};
	static_assert(sizeof(UnitInfo) == 12, "UnitInfo is packed");
	struct CheckpointInfo {
		uint32_t point;
		float length;
		uint32_t unit;
	};
	static_assert(sizeof(CheckpointInfo) == 12, "CheckpointInfo is packed");
}


//helper for reading vectors of data:
template< typename T >
static void read(std::istream &in, std::string magic, std::vector< T > *data_) {
	assert(magic.size() == 4);
	assert(data_);
	auto &data = *data_;

	struct {
		char magic[4];
		uint32_t size;
	} header;
	static_assert(sizeof(header) == 8, "header is packed");

	if (!in.read(reinterpret_cast< char * >(&header), sizeof(header))) {
		throw std::runtime_error("Failed to read header for '" + magic + "' chunk.");
	}
	if (std::string(header.magic, 4) != magic) {
		throw std::runtime_error("Expected '" + magic + "' chunk, got '" + std::string(header.magic, 4) + "'.");
	}
	if (header.size % sizeof(T) != 0) {
		throw std::runtime_error("Size of '" + magic + "' chunk not divisible by " + std::to_string(sizeof(T)) +".");
	}

	data.resize(header.size / sizeof(T));
	if (!in.read(reinterpret_cast< char * >(data.data()), data.size()*sizeof(T))) {
		throw std::runtime_error("Failed to read " + std::to_string(data.size()) + " elements (" + std::to_string(header.size) + " bytes) from '" + magic + "' chunk.");
	}
}

sm::Yarns sm::Yarns::load(std::string const &filename) {

	//arrays-of-structures that will be read:
	static_assert(sizeof(glm::vec3) == 12, "vec3 is packed");
	std::vector< glm::vec3 > in_points;
	std::vector< uint32_t > in_sources;
	std::vector< YarnInfo > in_yarns;
	std::vector< char > in_strings;
	std::vector< UnitInfo > in_units;
	std::vector< CheckpointInfo > in_checkpoints;

	{ //read from file:
		std::ifstream in(filename, std::ios::binary);
		read(in, "f3..", &in_points);
		read(in, "src.", &in_sources);
		read(in, "yarn", &in_yarns);
		read(in, "strs", &in_strings);
		read(in, "unit", &in_units);
		read(in, "chk.", &in_checkpoints);
	}

	if (in_sources.size() != in_points.size()) {
		throw std::runtime_error("Points and sources size mismatch in '" + filename + "'.");
	}

	sm::Yarns ret;

	ret.units.reserve(in_units.size());
	for (auto const &unit : in_units) {
		if (!(unit.name_begin <= unit.name_end && unit.name_end <= in_strings.size())) {
			throw std::runtime_error("Incorrect unit name indices in '" + filename + "'.");
		}
		ret.units.emplace_back();
		ret.units.back().name = std::string(in_strings.begin() + unit.name_begin, in_strings.begin() + unit.name_end);
		ret.units.back().length = unit.length;
	}
	/*
	struct YarnInfo {
		uint32_t point_begin;
		uint32_t point_end;
		float radius;
		glm::u8vec4 color;
	};
	*/

	auto cpi = in_checkpoints.begin();
	ret.yarns.reserve(in_yarns.size());
	for (auto const &yarn : in_yarns) {
		if (!(yarn.point_begin <= yarn.point_end && yarn.point_end <= in_points.size())) {
			throw std::runtime_error("Incorrect yarn indices in '" + filename + "'.");
		}
		ret.yarns.emplace_back();
		ret.yarns.back().points.assign(in_points.begin() + yarn.point_begin, in_points.begin() + yarn.point_end);
		ret.yarns.back().sources.assign(in_sources.begin() + yarn.point_begin, in_sources.begin() + yarn.point_end);
		ret.yarns.back().radius = yarn.radius;
		ret.yarns.back().color = yarn.color;

		//figure out the range of checkpoints for the yarn:
		auto checkpoint_begin = cpi;
		while (cpi != in_checkpoints.end() && cpi->point < yarn.point_end) {
			if (cpi->point < yarn.point_begin) {
				throw std::runtime_error("Out-of-order checkpoint in '" + filename + "'.");
			}
			if (!(cpi->unit < ret.units.size())) {
				throw std::runtime_error("Invalid unit in checkpoint in '" + filename + "'.");
			}
			++cpi;
		}
		auto checkpoint_end = cpi;

		//copy checkpoints into the yarn:
		ret.yarns.back().checkpoints.reserve(checkpoint_end - checkpoint_begin);
		for (auto cp = checkpoint_begin; cp != checkpoint_end; ++cp) {
			ret.yarns.back().checkpoints.emplace_back();
			ret.yarns.back().checkpoints.back().point = cp->point;
			ret.yarns.back().checkpoints.back().unit = cp->unit;
			ret.yarns.back().checkpoints.back().length = cp->length;
		}
	}

	if (cpi != in_checkpoints.end()) {
		throw std::runtime_error("Unused checkpoints in '" + filename + "'.");
	}

	return ret;
}

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
	std::vector< uint32_t > out_sources;
	std::vector< YarnInfo > out_yarns;
	std::vector< char > out_strings;
	std::vector< UnitInfo > out_units;
	std::vector< CheckpointInfo > out_checkpoints;

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
		assert(yarn.points.size() == yarn.sources.size());

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
		out_sources.insert(out_sources.end(), yarn.sources.begin(), yarn.sources.end());
		out_yarns.back().point_end = out_points.size();
		out_yarns.back().radius = yarn.radius;
		out_yarns.back().color = yarn.color;
	}

	write(out, "f3..", out_points);
	write(out, "src.", out_sources);
	write(out, "yarn", out_yarns);
	write(out, "strs", out_strings);
	write(out, "unit", out_units);
	write(out, "chk.", out_checkpoints);

}


//------------------------------------------------

static float twice_area(glm::vec2 const &a, glm::vec2 const &b, glm::vec2 const &c) {
	return glm::dot( c - a, glm::vec2( -(b.y - a.y), b.x - a.x ) );
}

void sm::mesh_and_library_to_yarns(sm::Mesh const &mesh, sm::Library const &library, sm::Yarns *yarns_) {
	assert(yarns_);
	auto &yarns = *yarns_;
	yarns = Yarns();

	auto conv_along = [](float along) -> int32_t {
		return int32_t(std::round(along * 1000));
	};

	auto inv_along = [](uint32_t along) {
		return 1000 - along;
	};

	auto conv_z = [](float z) -> int32_t {
		return int32_t(std::round(z * 1000));
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
		int8_t direction_vote = 0; //negative => yarn should go against this segment; 0 => no preference; positive => yarn should go along this segment
	};

	std::vector< std::vector< ChainSegment > > chains;
	{
		//build map connecting v0,v1,along,z pairs to other v0,v1,along,z via ChainSegment:
		std::vector< ChainSegment > segs;
		std::unordered_map< glm::ivec4, std::pair< glm::ivec4, uint32_t > > forward, reverse;

		uint32_t invalid_assignment = 0;

		for (auto const &face : mesh.faces) {
			//no yarns in face without type:
			assert(face.type < mesh_library.size());
			if (!mesh_library[face.type]) continue;
			auto const &sf = *mesh_library[face.type];
			if (sf.edges.size() != face.size()) {
				++invalid_assignment;
				continue;
			}
			for (uint32_t yi = 0; yi < sf.yarns.size(); ++yi) {
				segs.emplace_back(&face, yi, false);

				//construct direction vote from edge labels:
				auto is_yarn_type = [](std::string const &type) {
					return (!type.empty() && type[0] == 'y');
				};

				if (is_yarn_type(sf.edges[sf.yarns[yi].begin.edge].type)) {
					if (sf.edges[sf.yarns[yi].begin.edge].direction == Library::Face::Edge::In) {
						segs.back().direction_vote += 1;
					} else if (sf.edges[sf.yarns[yi].begin.edge].direction == Library::Face::Edge::Out) {
						segs.back().direction_vote -= 1;
					}
				}

				if (is_yarn_type(sf.edges[sf.yarns[yi].end.edge].type)) {
					if (sf.edges[sf.yarns[yi].end.edge].direction == Library::Face::Edge::In) {
						segs.back().direction_vote -= 1;
					} else if (sf.edges[sf.yarns[yi].end.edge].direction == Library::Face::Edge::Out) {
						segs.back().direction_vote += 1;
					}
				}


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

		if (invalid_assignment) {
			std::cerr << "WARNING: have " << invalid_assignment << " faces whose types have a different number of edges than the face." << std::endl;
		}


		std::unordered_map< glm::uvec2, glm::uvec2 > connections;
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
				auto ret = connections.insert(std::make_pair(from,to));
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
				auto ret = connections.insert(std::make_pair(from,to));
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


		uint32_t flipped_chains = 0;
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

			//look at orientation votes for chain:
			int32_t direction_vote = 0;
			for (auto const &seg : chain) {
				direction_vote += (seg.reverse ? -seg.direction_vote : seg.direction_vote);
			}
			uint32_t concurring = 0;
			uint32_t dissenting = 0;
			for (auto const &seg : chain) {
				int8_t vote = (seg.reverse ? -seg.direction_vote : seg.direction_vote);
				if (vote < 0) {
					if (direction_vote < 0) {
						++concurring;
					} else {
						++dissenting;
					}
				}
				if (vote > 0) {
					if (direction_vote > 0) {
						++concurring;
					} else {
						++dissenting;
					}
				}
			}
			if (dissenting > 0) {
				std::cerr << "WARNING: had chain with direction vote " << direction_vote << " with " << concurring << " concurring and " << dissenting << " dissenting (would generally expect *no* dissenting votes)." << std::endl;
			}

			if (direction_vote < 0) {
				flipped_chains += 1;
				std::reverse(chain.begin(), chain.end());
				for (auto &seg : chain) {
					seg.reverse = !seg.reverse;
				}
			}

			chains.emplace_back(chain.begin(), chain.end());
		}
		if (flipped_chains) {
			std::cerr << "NOTE: flipped " << flipped_chains << " chains owing to direction vote." << std::endl;
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

	//for every (non-null) face in library, compute edge + crossing for each yarn start/end:
	// (will be used for checkpoint lookup)
	struct EdgeCrossing {
		uint32_t edge = -1U;
		uint32_t crossing = -1U;
	};
	std::vector< std::vector< std::pair< EdgeCrossing, EdgeCrossing > > > mesh_library_yarn_crossings(mesh_library.size());
	{
		for (auto const &sfp : mesh_library) {
			if (sfp == nullptr) continue;

			auto const &sf = *sfp;
			auto &yarn_crossings = mesh_library_yarn_crossings[&sfp - &mesh_library[0]];
			yarn_crossings.resize(sf.yarns.size());

			//list edge points by edges:
			std::vector< std::vector< std::pair< Library::Face::EdgePoint, EdgeCrossing * > > > edge_points(sf.edges.size());

			for (uint32_t yi = 0; yi < sf.yarns.size(); ++yi) {
				assert(sf.yarns[yi].begin.edge < edge_points.size());
				edge_points[sf.yarns[yi].begin.edge].emplace_back(sf.yarns[yi].begin, &yarn_crossings[yi].first);

				assert(sf.yarns[yi].end.edge < edge_points.size());
				edge_points[sf.yarns[yi].end.edge].emplace_back(sf.yarns[yi].end, &yarn_crossings[yi].second);
			}

			//sort crossings at each edge, assign to proper crossings:
			for (auto &ep : edge_points) {
				std::stable_sort(ep.begin(), ep.end(), [](
					std::pair< Library::Face::EdgePoint, EdgeCrossing * > const &a,
					std::pair< Library::Face::EdgePoint, EdgeCrossing * > const &b) -> bool {
					if (a.first.along != b.first.along) return a.first.along < b.first.along;
					return a.first.z < b.first.z;
				});
				for (auto &p : ep) {
					p.second->edge = uint32_t(&ep - &edge_points[0]);
					p.second->crossing = uint32_t(&p - &ep[0]);
				}
			}

			//quick paranoid check:
			for (uint32_t yi = 0; yi < sf.yarns.size(); ++yi) {
				assert(sf.yarns[yi].begin.edge == yarn_crossings[yi].first.edge);
				assert(sf.yarns[yi].end.edge == yarn_crossings[yi].second.edge);
			}

		}
	}

	//lookup table for checkpoints:
	std::unordered_multimap< glm::uvec3, Mesh::Checkpoint const * > fec_to_checkpoints;
	for (auto const &c : mesh.checkpoints) {
		fec_to_checkpoints.insert(std::make_pair(
			glm::uvec3(c.face, c.edge, c.crossing),
			&c
		));
	}
	assert(fec_to_checkpoints.size() == mesh.checkpoints.size());

	//---------------
	// face yarns -> mesh

	std::mt19937 mt(0x15469519);
	auto yarn_color = [&mt]() -> glm::u8vec4 {
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
	float max_mismatch = 0.0f;
	for (auto const &chain : chains) {
		yarns.yarns.emplace_back();
		yarns.yarns.back().radius *= radius_scale;
		yarns.yarns.back().color = yarn_color();

		std::vector< std::pair< uint32_t, Mesh::Checkpoint const * > > checkpoints;

		for (auto const &seg : chain) {
			assert(seg.face);
			auto const &face = *seg.face;

			glm::uvec3 begin_fec = glm::uvec3(
				seg.face - &mesh.faces[0],
				mesh_library_yarn_crossings[face.type][seg.yarn].first.edge,
				mesh_library_yarn_crossings[face.type][seg.yarn].first.crossing
			);

			glm::uvec3 end_fec = glm::uvec3(
				seg.face - &mesh.faces[0],
				mesh_library_yarn_crossings[face.type][seg.yarn].second.edge,
				mesh_library_yarn_crossings[face.type][seg.yarn].second.crossing
			);

			auto &bf = mesh_baryfaces[face.type];

			auto yarn = bf.yarns[seg.yarn];
			if (seg.reverse) {
				std::reverse(yarn.begin(), yarn.end());
				std::swap(begin_fec, end_fec);
			}

			{ //add checkpoints for beginning of chain segment
				uint32_t point = ( yarns.yarns.back().points.empty() ? 0 : yarns.yarns.back().points.size() - 1);
				auto r = fec_to_checkpoints.equal_range(begin_fec);
				for (auto i = r.first; i != r.second; /* later */ ) {
					//assign checkpoint based on yarn point:
					checkpoints.emplace_back(std::make_pair(point, i->second));

					//remove checkpoint from lookup structure:
					auto old = i;
					++i;
					fec_to_checkpoints.erase(old);
				}
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
					float m = glm::length(pt - yarns.yarns.back().points.back());
					mismatch += m;
					max_mismatch = std::max(max_mismatch, m);
					//average just in case
					yarns.yarns.back().points.back() = 0.5f * (yarns.yarns.back().points.back() + pt);
					yarns.yarns.back().sources.back() = face.source;
				} else {
					yarns.yarns.back().points.emplace_back(pt);
					yarns.yarns.back().sources.emplace_back(face.source);
				}
			}

			{ //add checkpoints for end of chain segment
				assert(!yarns.yarns.back().points.empty());
				uint32_t point = yarns.yarns.back().points.size() - 1;
				auto r = fec_to_checkpoints.equal_range(end_fec);
				for (auto i = r.first; i != r.second; /* later */ ) {
					//assign checkpoint based on yarn point:
					checkpoints.emplace_back(std::make_pair(point, i->second));

					//remove checkpoint from lookup structure:
					auto old = i;
					++i;
					fec_to_checkpoints.erase(old);
				}
			}
		}
		//(old) TODO: for circular yarns, make sure first/last points match exactly.

		//Make sure yarn is properly oriented for checkpoints:
		if (!checkpoints.empty()) {
			assert(checkpoints.size() >= 2 && "yarns with any checkpoints *must* be covered with checkpoints -- thus, should have at least two checkpoints");
			assert(checkpoints[0].first == 0 && "yarns with any checkpoints *must* be covered with checkpoints -- thus, should always start at 0");
			assert(checkpoints.back().first == yarns.yarns.back().points.size()-1 && "yarns with any checkpoints *must* be covered with checkpoints -- thus, should always end at last point");
			//check if checkpoints are increasing or decreasing:
			bool is_increasing = false;
			bool is_decreasing = false;

			for (uint32_t i = 1; i < checkpoints.size(); ++i) {
				if (checkpoints[i-1].first != checkpoints[i].first) { //ignore order at the same point
					if (checkpoints[i-1].second < checkpoints[i].second) {
						is_increasing = true;
						assert(!is_decreasing && "checkpoints should be ordered along yarns");
					} else {
						is_decreasing = true;
						assert(!is_increasing && "checkpoints should be ordered along yarns");
					}
				}
			}
			assert(is_decreasing || is_increasing);

			if (is_decreasing) {
				std::cerr << "WARNING: checkpoint-enforced yarn order does not match with edge-label-suggested yarn order. Flipping to match checkpoints." << std::endl;

				//need to flip yarn orientation so that checkpoints line up properly:
				std::reverse(yarns.yarns.back().points.begin(), yarns.yarns.back().points.end());
				std::reverse(yarns.yarns.back().sources.begin(), yarns.yarns.back().sources.end());
				//move source labels to the start of their segments:
				yarns.yarns.back().sources.erase(yarns.yarns.back().sources.begin());
				yarns.yarns.back().sources.emplace_back(0);

				std::reverse(checkpoints.begin(), checkpoints.end());
				for (auto &pc : checkpoints) {
					pc.first = yarns.yarns.back().points.size() - 1 - pc.first;
				}

				//PARANOIA: check that flipping worked
				for (uint32_t i = 1; i < checkpoints.size(); ++i) {
					if (checkpoints[i-1].first != checkpoints[i].first) { //ignore order at the same point
						assert(checkpoints[i-1].second < checkpoints[i].second);
					}
				}
			}

			for (auto const &pc : checkpoints) {
				yarns.yarns.back().checkpoints.emplace_back();
				yarns.yarns.back().checkpoints.back().point = pc.first;
				yarns.yarns.back().checkpoints.back().length = pc.second->length;
				yarns.yarns.back().checkpoints.back().unit = pc.second->unit;
			}
		}

	}
	std::cout << "Total/max boundary mis-match from chains (should be very small): " << mismatch << "/" << max_mismatch << std::endl;

	//transfer units:
	//TODO: could consider eliminating any unused units
	for (auto const &u : mesh.units) {
		yarns.units.emplace_back();
		yarns.units.back().name = u.name;
		yarns.units.back().length = u.length;
	}

	std::cout << "Of " << mesh.checkpoints.size() << " checkpoints, " << fec_to_checkpoints.size() << " remain unassigned." << std::endl;

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

//------------------------------------------------

void sm::yarns_to_tristrip(sm::Yarns const &yarns, std::vector< sm::YarnAttribs > *attribs_, sm::Quality quality, std::vector< size_t > *end_attrib_) {
	assert(attribs_);
	auto &attribs = *attribs_;

	if (end_attrib_) end_attrib_->clear();

	uint32_t Angles = (quality != QualityLow ? 16 : 8);
	std::vector< glm::vec2 > Circle;
	Circle.reserve(Angles);
	for (uint32_t a = 0; a < Angles; ++a) {
		double ang = M_PI * 2.0 * double(a) / double(Angles);
		Circle[a].x = std::cos(ang);
		Circle[a].y = std::sin(ang);
	}

	for (auto const &yarn_struct : yarns.yarns) {
		std::vector< glm::vec3 > const &yarn = yarn_struct.points;
		float yarn_radius = yarn_struct.radius;
		glm::u8vec4 yarn_color = yarn_struct.color;
		if (yarn.size() < 2) {
			if (end_attrib_) end_attrib_->emplace_back(attribs.size());
			continue;
		}

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
				|| glm::dot(b-a,b-a) < 0.5 * 0.5 * yarn_radius * yarn_radius) {
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
					attribs.emplace_back(a + yarn_radius * n0, n0, yarn_color);
					if (i == 0 && attribs.size() > 1) attribs.emplace_back(attribs.back());
					attribs.emplace_back(a + yarn_radius * n1, n1, yarn_color);
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
					attribs.emplace_back(next_at + yarn_radius * next_n, next_n, yarn_color);
					if (i == 0) attribs.emplace_back(attribs.back());
					attribs.emplace_back(at + yarn_radius * n, n, yarn_color);
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

					attribs.emplace_back(a + yarn_radius * n1, n1, yarn_color);
					if (i == 0) attribs.emplace_back(attribs.back());
					attribs.emplace_back(a + yarn_radius * n0, n0, yarn_color);
				}
				attribs.emplace_back(attribs[attribs.size() - 2*Angles]);
				attribs.emplace_back(attribs[attribs.size() - 2*Angles]);
			}

		}

		if (end_attrib_) end_attrib_->emplace_back(attribs.size());
	}
	if (end_attrib_) assert(end_attrib_->size() == yarns.yarns.size());
}

void sm::derive_face(sm::Library::Face const &face, uint8_t by_bits, sm::Library::Face *face2_) {

	assert(face2_);
	auto &face2 = *face2_;

	face2.derive.from = face.key();
	face2.derive.by = by_bits;

	face2.edges.resize(face.edges.size());

	uint32_t r = 0;
	if (by_bits & sm::Library::Face::Derive::MirrorXBit) {
		//new first vertex should be after the run of loop-in edges on the side:
		auto remove_number_suffix = [&](std::string str) {
			while (!str.empty() && str[str.size()-1] >= '0' && str[str.size()-1] <= '9') {
				str.erase(str.size()-1);
			}
			return str;
		};
		while (r < face.edges.size() && face.edges[r].direction == face.edges[0].direction && remove_number_suffix(face.edges[r].type) == remove_number_suffix(face.edges[0].type)) ++r;

		//should not call on edges with all in edges:
		if (r >= face.edges.size()) {
			throw std::runtime_error("Cannot derive face '" + face2.name + "' from face '" + face.key() + "' via mirror-x because source has no out edges.");
		}

		for (uint32_t i = 0; i < face2.edges.size(); ++i) {
			//vertices should appear in reversed order, and should start with 'r':
			uint32_t v = (r + face.edges.size() - i) % face.edges.size();
			face2.edges[i].vertex = glm::vec2(-face.edges[v].vertex.x, face.edges[v].vertex.y);
			//edge info gets copied from previous edge:
			uint32_t e = (v + face.edges.size() - 1) % face.edges.size();
			face2.edges[i].direction = face.edges[e].direction;
			face2.edges[i].type = face.edges[e].type;
		}
	} else {
		face2.edges = face.edges;
	}

	if (by_bits & sm::Library::Face::Derive::ReverseYarnBit) {
		for (auto &e : face2.edges) {
			if (e.type.size() >= 1 && e.type[0] == 'y') {
				if (e.direction == sm::Library::Face::Edge::In) {
					e.direction = sm::Library::Face::Edge::Out;
				} else if (e.direction == sm::Library::Face::Edge::Out) {
					e.direction = sm::Library::Face::Edge::In;
				}
			}
		}
	}

	face2.yarns.clear();
	for (auto const &yarn : face.yarns) {
		face2.yarns.emplace_back();
		auto copy_end = [&](sm::Library::Face::EdgePoint const &from, sm::Library::Face::EdgePoint &to) {
			if (by_bits & sm::Library::Face::Derive::MirrorXBit) {
				// e = (r - i - 1) -> r - e - 1 = i
				to.edge = (r + face.edges.size() - from.edge - 1) % face.edges.size();
				to.along = 1.0f - from.along;
			} else {
				to.edge = from.edge;
				to.along = from.along;
			}
			if (by_bits & sm::Library::Face::Derive::MirrorZBit) {
				to.z = -from.z;
			} else {
				to.z = from.z;
			}
		};
		copy_end(yarn.begin, face2.yarns.back().begin);
		copy_end(yarn.end, face2.yarns.back().end);
		for (auto const &m : yarn.middle) {
			face2.yarns.back().middle.emplace_back(
				(by_bits & sm::Library::Face::Derive::MirrorXBit ? -m.x : m.x),
				m.y,
				(by_bits & sm::Library::Face::Derive::MirrorZBit ? -m.z : m.z)
			);
		}
	}

}

//-----------------------------------
// partial order of faces
// returns true if successfully ordered
// maybe should return a reordered Mesh..
// else false

sm::Mesh sm::order_faces(sm::Mesh const &mesh, sm::Library const &library){
	sm::Mesh out = mesh;
	std::map<std::string, uint32_t> name_to_lib_idx;
	std::vector<uint32_t> order;
	for(auto const &l : library.faces){
		name_to_lib_idx[l.key()] = &l-&library.faces[0];
	}
	uint32_t iterations = 0;
	std::unordered_set<uint32_t> completed_faces;
	// Build up candidate faces instead of this ridiculousness...but works
	while(true){
		for(auto const &f : mesh.faces){
			uint32_t fid = &f - &mesh.faces[0];
			if(completed_faces.count(fid)) continue;
			if(name_to_lib_idx.count(mesh.library[f.type]) == 0) {
				// probably not a throw?
				std::cout << "Face " << fid << " does not have a valid library name " << mesh.library[f.type] << std::endl;
				break;
			}
			auto l = library.faces[name_to_lib_idx[mesh.library[f.type]]]; // library face
			// does it have an "in", is the "in" connection done?
			bool ins_available = true;
			for(auto &e : l.edges){
				auto eid = &e - &l.edges[0];
				if(e.direction == sm::Library::Face::Edge::In){
					sm::Mesh::FaceEdge fe; fe.face = fid; fe.edge = eid;
					for(auto c : mesh.connections){
						if(c.a == fe){
							if(!completed_faces.count(c.b.face)) ins_available = false;
						}
						else if(c.b == fe){
							if(!completed_faces.count(c.a.face)) ins_available = false;
						}
					}
				}
			}
			if (ins_available){
				completed_faces.insert(&f - &mesh.faces[0]);
				order.emplace_back(&f - &mesh.faces[0]);
			}
		}
		if(completed_faces.size() == mesh.faces.size()) break;
		++iterations;
		if(iterations > mesh.faces.size()) break;
	}
	//DEBUG
	if(0){
		std::cout << "ORDER:"; for(auto o : order) std::cout << o << " "; std::cout << std::endl; // DEBUG just print this out
		if(mesh.faces.size() == completed_faces.size()){
		}
		else{
			std::cout <<"WARNING: faces could NOT be ordered!" << std::endl;
		}
	}
	if(completed_faces.size() == mesh.faces.size()){
		for(uint32_t i = 0; i < mesh.faces.size(); ++i){
			out.faces[i] = mesh.faces[order[i]];
		}
		for(auto &c : out.connections){
			c.a.face = std::distance(order.begin(), std::find(order.begin(), order.end(), c.a.face));
			c.b.face = std::distance(order.begin(), std::find(order.begin(), order.end(), c.b.face));
		}
		for(auto &h : out.hints){
			h.lhs.face = std::distance(order.begin(), std::find(order.begin(), order.end(), h.lhs.face));
		}
		for(auto &c : out.checkpoints){
			if(c.face != -1U)
				c.face = std::distance(order.begin(), std::find(order.begin(), order.end(), c.face));
		}
	}
	// A mesh where faces are correctly ordered if ordering was possible, else original..
	return out;
}
bool sm::add_hint(sm::Mesh::Hint const h, sm::Mesh *_mesh, sm::Library const &library, sm::Code const &code, std::vector<sm::Mesh::Hint> * _offenders){
	assert(_mesh);
	assert(_offenders);
	sm::Mesh &mesh = *_mesh;
	sm::Mesh temp = *_mesh;
	//std::vector<sm::Mesh::Hint> &offenders = *_offenders;
	temp.hints.emplace_back(h);
	if(verify(temp, code, _offenders)){
		mesh.hints.emplace_back(h);
		return true;
	}

	return false;
}

// if is_fully_hinted is true, checks also if mesh is fully hinted
// else only checks if hints are consistent
bool sm::verify(sm::Mesh const &mesh, sm::Code const &code, std::vector<sm::Mesh::Hint> *_offenders,  bool is_fully_hinted){

	assert(_offenders);
	std::vector<sm::Mesh::Hint> &offenders = *_offenders;

	offenders.clear();
	// TODO fill this up with offending hints

	// process all the hints into easy to access formats
	std::map<std::string, uint32_t> name_to_code_idx;
	std::map<uint32_t, std::string> face_variant;
	std::map<int32_t, int32_t> face_translation;
	for(auto const &c : code.faces){
		name_to_code_idx[c.key()] = &c - &code.faces[0];
	}



	// ----------------- hint verification -------------
	// Variant hints:
	for(auto h : mesh.hints){
		// TODO check that variant hints are consistent
		if(h.type == sm::Mesh::Hint::Variant){
			if(face_variant.count(h.lhs.face) && face_variant[h.lhs.face] != std::get<std::string>(h.rhs)){
				std::cerr << "[Variant] Hints for face " << h.lhs.face << " are inconsistent: " << face_variant[h.lhs.face] << ", " << std::get<std::string>(h.rhs) << std::endl;
				offenders.emplace_back(h);
				return false;
			}
			face_variant[h.lhs.face] = std::get<std::string>(h.rhs);
		}
	}

	{
		// hints don't currently specify yarn carriers, but should
	}
	auto face_to_code_key = [&](const sm::Mesh::Face &f)->std::string{
		std::string variant = "";
		if(face_variant.count(&f- &mesh.faces[0])){
			variant = face_variant[&f - &mesh.faces[0]];
		}
		std::string signature = mesh.library[f.type] + ' ' + variant;
		return signature;
	};

	// faces are associated with code signatures that exist in the code library
	if(is_fully_hinted){
		for(auto const &f : mesh.faces){
			std::string signature = face_to_code_key(f);
			if(!name_to_code_idx.count(signature)){
			// nothing offends, just missing hint -- can create an empty hint signature to indicate which face is missing?
			std::cerr << "Code does not exist for face signature " << mesh.library[f.type] << " code signature " << signature << std::endl;
			return false;
			}
		}
	}

	// Types of offences:
	// 1. Variant: Face has a variant name that does not exists (done)
	// 2. Variant: Face has a variant name that is inconsistent with resource hints. (done)
	// 3. Resource: Edge has a resource hint that is inconsistent with other edge resource hints (done)
	// 4. Resource: Edge has a resource hint that is inconsistent with connections.

	{
		// 1.
		for(auto const &h : mesh.hints){
			if(h.type != sm::Mesh::Hint::Variant) continue;
			std::string variant_name = std::get<std::string>(h.rhs);
			auto f = mesh.faces[h.lhs.face];
			std::string signature = mesh.library[f.type] + " " + variant_name;
			if(!name_to_code_idx.count(signature)){
				offenders.emplace_back(h); // variant name
			}
		}
		// 4.
		for(auto c : mesh.connections){
			// if there exists a resource hint for c.a and a resource hint for c.b and they are inconsistent
			// TODO only for loop edges

			auto ha_it = std::find_if(mesh.hints.begin(), mesh.hints.end(),[&](const sm::Mesh::Hint &h)->bool{return (h.type == sm::Mesh::Hint::Resource && h.lhs == c.a);} );
			auto hb_it = std::find_if(mesh.hints.begin(), mesh.hints.end(),[&](const sm::Mesh::Hint &h)->bool{return (h.type == sm::Mesh::Hint::Resource && h.lhs == c.b);} );
			if(ha_it != mesh.hints.end() && hb_it != mesh.hints.end()){
				sm::BedNeedle rhs_a = std::get<sm::BedNeedle>(ha_it->rhs);
				sm::BedNeedle rhs_b = std::get<sm::BedNeedle>(hb_it->rhs);
				if(!(rhs_a == rhs_b)){
					offenders.emplace_back(*ha_it);
					offenders.emplace_back(*hb_it);
				}
			}
		}

	}

	// Resource hints
	for(auto h: mesh.hints){
		if(h.type == sm::Mesh::Hint::Resource){
			const sm::Mesh::Face &f = mesh.faces[h.lhs.face];
			std::string signature = face_to_code_key(f);
			auto const &l = code.faces[name_to_code_idx[signature]];
			if(l.edges.size() != f.size()){
				// resource hint is not offending, code and face don't match -- must be caught much earlier than this...
				std::cerr <<"Code face and mesh face have different number of edges." << std::endl;
				return false;
			}
			const sm::BedNeedle bn = std::get<sm::BedNeedle>(h.rhs);
			const sm::BedNeedle bn_template = l.edges[h.lhs.edge].bn;
			if(bn_template.dontcare()) continue;
			if(bn.bed != bn_template.bed){
				offenders.emplace_back(h); // edge resource conflicts with variant.
				// also find variant hint for this face and add it to offence.
				for(auto hh : mesh.hints){
					if(hh.type == sm::Mesh::Hint::Variant && hh.lhs.face == h.lhs.face){
						offenders.emplace_back(hh);
					}
				}
			}
			int offset = bn.location() - bn_template.location();
			if( face_translation.count(h.lhs.face) && face_translation[h.lhs.face] != offset){
				std::cerr <<"Resource hints are not consistent!" << h.lhs.face << "/" << h.lhs.edge<< std::endl;
				std::cerr <<"\tOld translation: " << face_translation[h.lhs.face] << " new translation: " << offset << std::endl;
				std::cerr <<"\tLoc Template: " << bn_template.location() << " Hint: " << bn.location() << std::endl;

				offenders.emplace_back(h);

				return false;
			}
			face_translation[h.lhs.face] = offset;
		}
	}
	// all faces have a resource hint
	for(auto &f : mesh.faces){
		if(!face_translation.count(&f - &mesh.faces[0]) && is_fully_hinted){
			// no offenders, missing hints
			std::cerr << "All faces have not been hinted. "<< &f - &mesh.faces[0] << std::endl;
			return false;
		}
	}


	std::vector<std::string> carriers;
	{
		// Go over all the yarn carriers that appear in each face code
		// to generate a carrier use string
		for(auto const &f : mesh.faces){
			std::string signature = face_to_code_key(f);
			auto const &l = code.faces[name_to_code_idx[signature]];
			for(auto c : l.carriers){
				if(std::find(carriers.begin(), carriers.end(), c) == carriers.end()){
					carriers.emplace_back(c);
				}
			}
		}
		// Check that carriers are globally consistent
		for(auto const &f : mesh.faces){
			std::string signature = face_to_code_key(f);
			auto const &l = code.faces[name_to_code_idx[signature]];
			for(uint32_t k = 1; k < l.carriers.size(); ++k){
				auto d1 = std::distance(carriers.begin(), std::find(carriers.begin(), carriers.end(),l.carriers[k-1]));
				auto d2 = std::distance(carriers.begin(), std::find(carriers.begin(), carriers.end(), l.carriers[k]));
				if(d2 < d1){
					// maybe this should be a wrong variant offence TODO
					std::cerr << "Local face carrier order not consistent with global order. " << std::endl;
					return false;
				}
			}
		}
	}

	// Order hints
	{

		//check for consistency independent of total order for verification purposes
		{
			uint32_t N = 0;
			std::map<sm::Mesh::FaceEdge, std::set<sm::Mesh::FaceEdge>> depends_on, inverse_depends_on;
			for(auto const &f : mesh.faces){
				std::string signature = face_to_code_key(f);
				auto const &l = code.faces[name_to_code_idx[signature]];
				N += l.instrs.size();
				for(uint32_t i = 0; i < l.instrs.size(); ++i){
					sm::Mesh::FaceEdge fe; fe.face = &f - &mesh.faces[0];
					fe.edge = i;
					depends_on[fe] = std::set<sm::Mesh::FaceEdge>();
					inverse_depends_on[fe] = std::set<sm::Mesh::FaceEdge>();
				}
			}
			std::vector<sm::Mesh::FaceEdge> total_order;
			std::set<sm::Mesh::FaceEdge> done;
			for(auto const &h: mesh.hints){
				if(h.type == sm::Mesh::Hint::Order){
					depends_on[h.lhs].insert(std::get<sm::Mesh::FaceEdge>(h.rhs));
					inverse_depends_on[std::get<sm::Mesh::FaceEdge>(h.rhs)].insert(h.lhs);
				}
			}
			bool did_update = true;
			while(did_update){ //lazy-version
				did_update = false;
				for(auto ba : depends_on){
					if(ba.second.empty() && !done.count(ba.first)){
						total_order.emplace_back(ba.first);
						done.insert(ba.first);
						for(auto x: inverse_depends_on[ba.first]){
							assert(depends_on[x].count(ba.first));
							depends_on[x].erase(ba.first);
						}
						did_update = true;
					}
				}
			}
			if(total_order.size() < N){
				std::cerr << "Order hints are inconsistent. " << std::endl;
				//TODO order hints are inconsistent, need to figure out which one and add to offending list
				return false;
			}
		}
		// respect total order, as long as total order is consistent, this should be consistent
		auto order = mesh.total_order;
		for(auto const &h: mesh.hints){
			if(h.type == sm::Mesh::Hint::Order){
				std::pair<uint32_t, uint32_t> before, after;
				before.first = h.lhs.face;
				before.second = h.lhs.edge;
				after.first = (std::get<sm::Mesh::FaceEdge>(h.rhs)).face;
				after.second = (std::get<sm::Mesh::FaceEdge>(h.rhs)).edge;
				auto a_it = std::find(order.begin(), order.end(), after);
				auto b_it = std::find(order.begin(), order.end(), before);
				if(std::distance(order.begin(), a_it) < std::distance(order.begin(), b_it)){
					std::cerr << "Order hint inconsistent with total order " << std::endl;
					offenders.emplace_back(h);
					return false;
				}
			}
		}
	}

	// Slack is respected
	for(auto const &c : mesh.connections){
		int32_t a_offset  = face_translation[c.a.face];
		int32_t b_offset  = face_translation[c.b.face];
		auto const &ca = code.faces[name_to_code_idx[face_to_code_key(mesh.faces[c.a.face])]];
		auto const &cb = code.faces[name_to_code_idx[face_to_code_key(mesh.faces[c.b.face])]];
		auto const &ea = ca.edges[c.a.edge];
		auto const &eb = cb.edges[c.b.edge];
		float a = ea.bn.location() + a_offset;
		float b = eb.bn.location() + b_offset;
		if(std::abs(a-b) > 0.5){
			// connection does not respect slack/loop constraint
			std::cerr << "Slack is not respected at connection between " << c.a.face << "/" << c.a.edge << " and " << c.b.face << "/" << c.b.edge << " offsets: " << a << ", " << b  << std::endl;
			//TODO find the resource hints for these connections and add them to offending
			return false;
		}
	}
	// sanity of total order
	{

		// no repetition
		{
			auto order = mesh.total_order;
			auto it = std::unique(order.begin(), order.end());
			if(it != order.end()){
				std::cerr << "Total order has duplicate entries!"  << std::endl;
				return false;
			}
		}
		// all face instructions covered
		for(auto const &f : mesh.faces){
			std::string signature = face_to_code_key(f);
			auto const &l = code.faces[name_to_code_idx[signature]];
			for(uint32_t i = 0; i < l.instrs.size(); ++i){
				const std::pair<uint32_t, uint32_t> fi = std::make_pair(&f - &mesh.faces[0], i);
				auto it = std::find(mesh.total_order.begin(), mesh.total_order.end(), fi);
				if(it == mesh.total_order.end() && is_fully_hinted){
					std::cerr << "Total order is incomplete, does not feature face/instr : "<< &f - &mesh.faces[0] << "/" << i << " ." << std::endl;
					return false;
				}
			}
		}
		// if total order exists and covers all face/instructions without repetition then it must be acyclic
	}
	// --------------End of hint verification -------------

	// if !is_fully_hinted, then true only indicates that no inconsistency exists
	if(!offenders.empty()) return false;
	return true;
}

bool sm::compute_total_order(sm::Mesh &mesh, sm::Code const &code){

	mesh.total_order.clear();

	// go through all the faces and just add them as-is
	// process all the hints into easy to access formats
	std::map<std::string, uint32_t> name_to_code_idx;
	std::map<uint32_t, std::string> face_variant;

	for(auto const &c : code.faces){
		name_to_code_idx[c.key()] = &c - &code.faces[0];
	}
	for(auto h : mesh.hints){
		// TODO check that variant hints are consistent
		if(h.type == sm::Mesh::Hint::Variant){
			if(face_variant.count(h.lhs.face) && face_variant[h.lhs.face] != std::get<std::string>(h.rhs)){
				std::cerr << "[Variant] Hints for face " << h.lhs.face << " are inconsistent: " << face_variant[h.lhs.face] << ", " << std::get<std::string>(h.rhs) << std::endl;
				return false;
			}
			face_variant[h.lhs.face] = std::get<std::string>(h.rhs);
		}
	}


	for(uint32_t fid = 0; fid < mesh.faces.size(); ++fid){
		if(face_variant.count(fid) == 0){
			std::cout << "Variant missing for face " << fid << std::endl;
			return false;
		}
	}
	for(uint32_t fid = 0; fid < mesh.faces.size(); ++fid){
		std::string name = mesh.library[mesh.faces[fid].type] + " " + face_variant[fid];
		if(!name_to_code_idx.count(name)){
			std::cout << "Code/Variant missing for face " << fid << std::endl;
			return false;
		}
		uint32_t cid = name_to_code_idx[name];
		auto const &c = code.faces[cid];
		for(uint32_t i = 0; i < c.instrs.size(); ++i){
			std::pair<uint32_t, uint32_t> fi;
			fi.first = fid;
			fi.second = i;
			mesh.total_order.emplace_back(fi);
		}
	}
	return true;
}

// knitout from faces
// assumes hinting is complete and valid
std::string sm::knitout(sm::Mesh const &mesh, sm::Code const &code){


	std::vector<Mesh::Hint> offenders;
	bool strict = true;
	bool verified = verify(mesh, code, &offenders, strict);

	if(!verified){
		std::cout << "Constraints are not complete; #offending constraints = " << offenders.size() << std::endl;
		return "";
	}
	assert(offenders.empty());

	// constraints check out, now "run" the machine

	// process all the hints into easy to access formats
	std::map<std::string, uint32_t> name_to_code_idx;
	std::map<uint32_t, std::string> face_variant;
	std::map<int32_t, int32_t> face_translation;
	std::vector<std::string> carriers;

	auto face_to_code_key = [&](const sm::Mesh::Face &f)->std::string{
		std::string variant = "";
		if(face_variant.count(&f- &mesh.faces[0])){
			variant = face_variant[&f - &mesh.faces[0]];
		}
		std::string signature = mesh.library[f.type] + ' ' + variant;
		return signature;
	};


	for(auto const &c : code.faces){
		name_to_code_idx[c.key()] = &c - &code.faces[0];
	}
	for(auto h : mesh.hints){
		// TODO check that variant hints are consistent
		if(h.type == sm::Mesh::Hint::Variant){
			if(face_variant.count(h.lhs.face) && face_variant[h.lhs.face] != std::get<std::string>(h.rhs)){
				std::cerr << "[Variant] Hints for face " << h.lhs.face << " are inconsistent: " << face_variant[h.lhs.face] << ", " << std::get<std::string>(h.rhs) << std::endl;
				return "";
			}
			face_variant[h.lhs.face] = std::get<std::string>(h.rhs);
		}
	}

	// Resource hints
	for(auto h: mesh.hints){
		if(h.type == sm::Mesh::Hint::Resource){
			const sm::Mesh::Face &f = mesh.faces[h.lhs.face];
			std::string signature = face_to_code_key(f);
			auto const &l = code.faces[name_to_code_idx[signature]];
			const sm::BedNeedle bn = std::get<sm::BedNeedle>(h.rhs);
			const sm::BedNeedle bn_template = l.edges[h.lhs.edge].bn;
			if(bn_template.dontcare()) continue;
			int offset = bn.location() - bn_template.location();
			// verify checks that these are indeed valid.
			face_translation[h.lhs.face] = offset;
		}
	}
	{
		// Go over all the yarn carriers that appear in each face code
		// to generate a carrier header string
		// TODO this needs to be order consistent.
		// a face with c:1,3 c:2,3 c:1,2 should generate carriers:1,2,3
		// i.e., sort with a function that goes over all the code faces
		for(auto const &f : mesh.faces){
			std::string signature = face_to_code_key(f);
			auto const &l = code.faces[name_to_code_idx[signature]];
			for(auto c : l.carriers){
				if(std::find(carriers.begin(), carriers.end(), c) == carriers.end()){
					carriers.emplace_back(c);
				}
			}
		}
	}
	// --------------Code Generation------------------------
	std::string knitout_string = "";
	// knitout header
	knitout_string += ";!knitout-2\n";


	{
		// append all the carriers
		knitout_string += ";;Carriers: ";
		for(auto c : carriers) knitout_string += c + " ";
		knitout_string += "\n";
	}
	// what to blame:
	knitout_string += ";auto-generated from smobj, see sm.cpp,  sm::knitout()\n";


	// go by order, plug translate by and hint index to get knitout instruction
	for(auto const &fi : mesh.total_order){
		if(!face_translation.count(fi.first)) {
		std::cerr << "Mesh is not sufficiently hinted, total order missing " << std::endl; break;
		}
		std::string signature = mesh.library[mesh.faces[fi.first].type] +  ' ' + (face_variant.count(fi.first) ? face_variant[fi.first] : "");
		auto const &l = code.faces[name_to_code_idx[signature]];
		int t = face_translation[fi.first];
		knitout_string += l.knitout_string(t, fi.second, true);
		knitout_string += "\n";

		//check if instruction can be made, make
		//make_instruction_on_machine()

	}

	return knitout_string;
}


bool sm::create_out_slack_xfer(uint32_t face_id,  sm::Mesh &mesh, sm::Library &library, sm::Code &code){

	std::cout << "Creating out-slack-xfer face " << std::endl;
	assert(face_id < mesh.faces.size());

	auto f = mesh.faces[face_id];

	uint32_t yarn_out_edge = -1U;
	uint32_t loop_out_edge = -1U;

	uint32_t active_loop_edges = 0;

	// this might be > rack-limit, assumes infinite-racking machine
	int required_shift = 0;

	std::map<std::string, uint32_t> name_to_idx;
	for(auto const &l : library.faces){
		name_to_idx[l.key()] = &l - &library.faces[0];
	}

	//count active loops after face_id is constructed (not including loop_out)
	//expects faces to be globally partially ordered.
	//find edge connections for last yarn-out and last loop-out
	auto const &fl = library.faces[name_to_idx[mesh.library[mesh.faces[face_id].type]]];
	for(uint32_t i = 0; i < fl.edges.size(); ++i){
		auto e1 = fl.edges[i];
		auto e2 = fl.edges[(i+1)%fl.edges.size()];
		if(e1.direction != sm::Library::Face::Edge::Out) continue;
		if(e2.direction != sm::Library::Face::Edge::Out) continue;
		if(e1.type[0] == 'l' && e2.type[0] == 'y'){
			yarn_out_edge = (i+1)%fl.edges.size();
			loop_out_edge = i;
		}
		else if(e1.type[0] == 'y' && e2.type[0] == 'l'){
			loop_out_edge = (i+1)%fl.edges.size();
			yarn_out_edge = i;

		}
	}

	assert(loop_out_edge != -1U);
	assert(yarn_out_edge != -1U);
	std::cout << "Found active loop-out and yarn-out edges. " << std::endl;
	std::set<sm::Mesh::FaceEdge> active_set;

	{
		sm::Mesh::FaceEdge yfe, lfe;
		yfe.face = lfe.face = face_id;
		yfe.edge = yarn_out_edge;
		lfe.edge = loop_out_edge;
		active_set.insert(yfe);
		active_set.insert(lfe);
	}

	auto find_needle_from_constraint = [&](sm::Mesh::FaceEdge lhs)->int{
		for(auto const &h : mesh.hints){
			if(h.type == sm::Mesh::Hint::Resource && h.lhs == lhs){
				return std::get<sm::BedNeedle>(h.rhs).needle;
			}
		}
		std::cerr << "Hint does not exist for yarn-out edges!" << std::endl;
		return 0;
	};
	//find connections of type 'loop' that have one face < face_id && one > face_id
	for(auto const &c : mesh.connections){
		auto la = library.faces[name_to_idx[mesh.library[mesh.faces[c.a.face].type]]];
		auto lb = library.faces[name_to_idx[mesh.library[mesh.faces[c.b.face].type]]];
		auto a = c.a;
		auto b = c.b;
		if(a.face == face_id && a.edge == yarn_out_edge){
			// get required slack from c.b
			required_shift = find_needle_from_constraint(c.b) - find_needle_from_constraint(c.a);
		}
		else if(b.face == face_id && b.edge == yarn_out_edge){
			// get required slack from c.a
			required_shift = find_needle_from_constraint(c.a) - find_needle_from_constraint(c.b);
		}

		//if not loop-loop continue
		{
			if(la.edges[c.a.edge].type[0] != 'l'){
				assert(lb.edges[c.b.edge].type[0] != 'l'); // types must match, right?
				continue;
			}
		}

		if(c.b.face <= face_id && c.a.face > face_id){
			std::swap(a,b);
			std::swap(la,lb);
		}
		// only case we care about
		if(a.face <= face_id && b.face > face_id){
			if(a.face == face_id && a.edge == loop_out_edge){
				// this is already recorded.
			}
			else{
				assert(la.edges[a.edge].direction == sm::Library::Face::Edge::Out);
				++active_loop_edges;
				active_set.insert(a);
			}
		}
	}
	std::cout << "Required shift:  " << required_shift << std::endl;
	std::cout << "Identified active edges: "; for(auto fe : active_set) std::cout << fe.face << "/" << fe.edge << " "; std::cout << std::endl;
	std::vector<sm::Mesh::FaceEdge> edge_list;
	// order edges to form the appropriate connections
	if(!active_set.empty()){
		std::cout << "Ordering " << active_set.size() << " active edges. " << std::endl;
		//TODO prune active_set, remove pieces not connected in the connected components of face_id.
		{
		}

		//order edges based in some reasonable manner
		//are edges always consistently ordered ?
		edge_list.emplace_back(*active_set.begin());
		active_set.erase(edge_list[0]);

		while(!active_set.empty()){
		//std::cout << "\t(b)Ordered edge-list "; for(auto e : edge_list) std::cout << e.face << "/" << e.edge << " "; std::cout << std::endl;
			sm::Mesh::FaceEdge next;
			next.face = -1U;
			auto prev = edge_list.back();
			float d = std::numeric_limits<float>::infinity();
			bool begin = false;
			for(auto fe : active_set){
				uint32_t v0 = mesh.faces[prev.face][(prev.edge + 1)%mesh.faces[prev.face].size()];
				uint32_t v1 = mesh.faces[fe.face][fe.edge];
				uint32_t v2 = mesh.faces[prev.face][prev.edge];
				uint32_t v3 = mesh.faces[fe.face][(fe.edge + 1)%mesh.faces[fe.face].size()];
				if(glm::distance(mesh.vertices[v0], mesh.vertices[v1]) < d){
					d = glm::distance(mesh.vertices[v0], mesh.vertices[v1]);
					next = fe;
					begin  = false;
				}
				if(glm::distance(mesh.vertices[v2], mesh.vertices[v3]) < d){
					d = glm::distance(mesh.vertices[v2], mesh.vertices[v3]);
					next = fe;
					begin = true;
				}
			}
			assert(next.face != -1U);
			if(begin) edge_list.insert(edge_list.begin(), next);
			else edge_list.emplace_back(next);
			active_set.erase(next);
		}

		//std::cout << "\tOrdered edge-list "; for(auto e : edge_list) std::cout << e.face << "/" << e.edge << " "; std::cout << std::endl;
	}

	std::cout << "Ordered edge-list "; for(auto e : edge_list) std::cout << e.face << "/" << e.edge << " "; std::cout << std::endl;

	// face has 2*active_out_edges + 2*(1 loop out edge) + 2*(1 yarn out edge) + 2 boundary no-op edges
	sm::Mesh::Face mesh_face;
	sm::Library::Face lib_face;
	sm::Code::Face code_face;
	std::string name = "slack_xfer_after_"+std::to_string(face_id);
	//mesh_face
	{
		mesh_face.type = mesh.library.size();
		mesh.library.emplace_back(name);
		std::cout << "Added " << name << " to mesh library" << std::endl;
		// ... add mesh vertices...
		std::vector<glm::vec3> vertices;
		for(uint32_t i = 0; i < edge_list.size(); ++i){
			auto fe = edge_list[i];
			uint32_t v = mesh.faces[fe.face][fe.edge];
			vertices.emplace_back(mesh.vertices[v] + glm::vec3(0.f,0.f, 0.5f));
			if(i+1 ==edge_list.size()){
				uint32_t w = mesh.faces[fe.face][(fe.edge+1)%mesh.faces[fe.face].size()];
				vertices.emplace_back(mesh.vertices[w] + glm::vec3(0.f,0.f,0.5f));
			}
		}
		std::cout << "added " << vertices.size() << " vertices.. " << std::endl;
		int n = vertices.size()-1;
		if(n > 0){
			for(int i = n ; i >= 0; i--){
				//TODO offset these a bit based on face/edge normal
				//arbitrary offset to avoid NaNs in yarn interpolation
				vertices.emplace_back(vertices[i]+glm::vec3(0.0f,0.25f,0.25f));

			}
		}
		std::cout << "adding " << vertices.size() << " vertices." << std::endl;
		for(uint32_t i = 0; i < vertices.size(); ++i){
			mesh_face.emplace_back(i + mesh.vertices.size());
		}
		std::cout << "adding indices to face " << mesh_face.size()  << std::endl;
		for(auto v : vertices) mesh.vertices.emplace_back(v);
		mesh.faces.emplace_back(mesh_face);

		std::cout << "Added mesh face. " << std::endl;
	}
	//lib_face
	{
		lib_face.name = name;
		// ...edges...
		glm::vec2 pos(0,0);
		for(auto fe : edge_list){
			if(fe.face == face_id && fe.edge == yarn_out_edge){
				lib_face.edges.emplace_back();
				lib_face.edges.back().direction = sm::Library::Face::Edge::In;
				lib_face.edges.back().type = "y";
				lib_face.edges.back().vertex = pos;
			}
			else{
				lib_face.edges.emplace_back();
				lib_face.edges.back().direction = sm::Library::Face::Edge::In;
				lib_face.edges.back().type = "l";
				lib_face.edges.back().vertex = pos;
			}
			pos.x += 1.f;
		}
		{
			lib_face.edges.emplace_back();
			lib_face.edges.back().direction = sm::Library::Face::Edge::Any;
			lib_face.edges.back().type = "x";
			lib_face.edges.back().vertex = pos;
		}
		pos.y = 1.f;
		for(auto it = edge_list.rbegin(); it != edge_list.rend(); ++it){
			auto fe = *it;
			if(fe.face == face_id && fe.edge == yarn_out_edge){
				lib_face.edges.emplace_back();
				lib_face.edges.back().direction = sm::Library::Face::Edge::Out;
				lib_face.edges.back().type = "y";
				lib_face.edges.back().vertex = pos;
			}
			else{
				lib_face.edges.emplace_back();
				lib_face.edges.back().direction = sm::Library::Face::Edge::Out;
				lib_face.edges.back().type = "l";
				lib_face.edges.back().vertex = pos;
			}
			pos.x -= 1.f;
		}
		{
			lib_face.edges.emplace_back();
			lib_face.edges.back().direction = sm::Library::Face::Edge::Any;
			lib_face.edges.back().type = "x";
			lib_face.edges.back().vertex = pos;
		}


		for(uint32_t i = 0; i < edge_list.size(); ++i){
			//fixed values for connections to work sensibly:
			auto fe = edge_list[i];

			if(fe.face == face_id && fe.edge == yarn_out_edge){
				sm::Library::Face::Yarn yarn;
				yarn.begin.edge = i;
				yarn.begin.along = 0.25f;
				yarn.end.edge = (edge_list.size() - i +edge_list.size());
				yarn.end.along = 0.75f;

				lib_face.yarns.emplace_back(yarn);
			}
			else{
				sm::Library::Face::Yarn yarn;
				yarn.begin.edge = i;
				yarn.begin.along = 0.275f;
				yarn.end.edge = (edge_list.size() - i +edge_list.size());
				yarn.end.along = 0.725f;
				lib_face.yarns.emplace_back(yarn);
				yarn.begin.along = 0.725f;
				yarn.end.along = 0.275f;
				lib_face.yarns.emplace_back(yarn);
			}
			// for each edge, add yarn segments
		}

		library.faces.emplace_back(lib_face);

		std::cout << "Added lib face. " << std::endl;
	}
	//code_face
	{
		code_face.name = name;
		code_face.check_bn = false;
		code_face.variant = "v";
		for(auto e : lib_face.edges){
			code_face.edges.emplace_back();
			if(e.direction == sm::Library::Face::Edge::In) code_face.edges.back().direction = sm::Code::Face::Edge::Out;
			else if(e.direction == sm::Library::Face::Edge::Out) code_face.edges.back().direction = sm::Code::Face::Edge::Out;
			code_face.edges.back().type = e.type;
			//code_face.edges.back().bn = ?? what template value

			// add instrs for as many edges
		}
		code.faces.emplace_back(code_face);
		std::cout << "Added code face. " << std::endl;
	}
	{
		assert(lib_face.edges.size() == mesh_face.size());
	}

	{
		// update all the connections between involving active_set
		uint32_t nc = mesh.connections.size();
		uint32_t fid = mesh.faces.size()-1;
		//todo flip loop connections?
		for(uint32_t i = 0; i < nc; ++i){
			auto &c = mesh.connections[i];
			auto it_a = std::find(edge_list.begin(), edge_list.end(), c.a);
			auto it_b = std::find(edge_list.begin(), edge_list.end(), c.b);
			if(it_a != edge_list.end()){
				std::cout << "Flip for connection between " << c.a.face << "/" << c.a.edge << " and " << c.b.face << "/" << c.b.edge << " is " << c.flip << std::endl;
				// c.a->a, b->c.b
				assert(it_b == edge_list.end());
				int idx_a = std::distance(edge_list.begin(), it_a);
				int idx_b = (edge_list.size() - idx_a) + edge_list.size();
				sm::Mesh::FaceEdge a; a.face = fid; a.edge = idx_a;
				sm::Mesh::FaceEdge b; b.face = fid; b.edge = idx_b;

				mesh.connections.emplace_back();
				mesh.connections.back().a = b;
				mesh.connections.back().b = c.b;
				mesh.connections.back().flip = false; // how flip is treated in the ui is a tad confusing, todo revisit
				c.b = a;
				c.flip = false;
			}
			else if(it_b != edge_list.end()){
				//b->a ==> c.b->b,a->c.a
				std::cout << "Flip for connection between " << c.a.face << "/" << c.a.edge << " and " << c.b.face << "/" << c.b.edge << " is " << c.flip << std::endl;
				assert(it_a == edge_list.end());
				int idx_b = std::distance(edge_list.begin(), it_b);
				int idx_a = (edge_list.size() - idx_b) + edge_list.size();
				sm::Mesh::FaceEdge a; a.face = fid; a.edge = idx_a;
				sm::Mesh::FaceEdge b; b.face = fid; b.edge = idx_b;

				mesh.connections.emplace_back();
				mesh.connections.back().a = a;
				mesh.connections.back().b = c.a;
				mesh.connections.back().flip = false;
				c.a = b;
				c.flip = false;
			}
		}

	}


	// add xfer or (mov) instructions based on required_rack using any number of new needles (at added depth).

	// debug print
	{
		std::cout << "creating out_slack_xfer face before moving forward from " << face_id  << std::endl;
		std::cout << "\t required-shift : " << required_shift << std::endl;
		std::cout << "\t loop-out-edge : " << loop_out_edge << " yarn-out-edge : " << yarn_out_edge << std::endl;
		std::cout << "\t active loops : ";
		for( auto fe: active_set){
			std::cout << fe.face << "/" << fe.edge << ", ";
		}
		std::cout << std::endl;
	}


	return true;
}

bool sm::create_in_slack_xfer(uint32_t face_id,  sm::Mesh &mesh, sm::Library &library, sm::Code &code){

	// similar to previous, but created by the "next" stitch to adjust incoming yarn.
	// in general, the created stitch calls for adjustment if its outgoing loop doesn't match incoming loop, so this might not be needed.
	return false;
}
