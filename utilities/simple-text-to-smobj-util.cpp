#include "sm.hpp"

#include <string>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <sstream>
#include <list>
#include <iterator>



// as most crocheters are right handed and they crochet to the left, 
// i will designate left as 0. sorry, left handers
enum Direction{
	left = 0, 
	right = 1,
};

//likewise, the front is the natural state.
//need a better name.... could just make it a bool
enum Flipped{
	front = 0,
	back = 1,
};

//AD added new face with direction struct
struct FaceWithDirection{
	sm::Library::Face face;
	Direction direction;
	Flipped flip; 
};


int main(int argc, char **argv) {
	if (argc != 4) {
		std::cerr <<
			"Usage:\n\t./text-to-smobj <in.txt> <library.sf> <out.smobj>\n"
			" Where in.txt looks something like this:\n"
			"    face-name e1 e2 e3 e4\n"
			"    face-name e2 e3 . .\n"
			"    ...\n"
			" That is, each line has a face name followed by variables indicating\n"
			" where that face's edge connect. '.' is a special variable meaning\n"
			" 'does not connect'"
			<< std::endl;
		return 1;
	}
	std::string in_text = argv[1];
	std::string in_library = argv[2];
	std::string out_smobj = argv[3];
	std::cout << "Will use connections in '" << in_text << "' and faces in library '" << in_library << "' to create smobj '" << out_smobj << "'" << std::endl;

	//-------------------------

	sm::Library library = sm::Library::load(in_library);

	// AD: change the map to work with the new datastruture 
	// map of string to list of FaceWithDirections
	// the string is the input, the list is the potential interpertations.
	std::unordered_map< std::string, std::list<FaceWithDirection> > shorthand_to_face;
	
	for (auto const &f : library.faces) {
		FaceWithDirection face_to_add;
		//add the face, if it wasnt made from a derive, it goes to the right and faces forward
		//if made from a derive, then calculate which side it is.
		face_to_add.face = f;
		if (f.derive.from == "") { //no derivation
			face_to_add.direction = Direction::right;
			face_to_add.flip = Flipped::front;
		}
		else { // there is a derivation
			if (f.derive.by & 001) {
				face_to_add.direction = Direction::left;
				face_to_add.flip = Flipped::front;
			}
			if (f.derive.by & 010) {
				face_to_add.direction = Direction::right;
				face_to_add.flip = Flipped::back;
			}
			if (f.derive.by & 011) {
				face_to_add.direction = Direction::left;
				face_to_add.flip = Flipped::back;
			}
			

		}

		//the first two characters as the alias
		std::string shorthand = face_to_add.face.shorthand;
		
		//

		// if the lookup fails, we will create a list
		
  		std::unordered_map<std::string,std::list<FaceWithDirection>>::iterator got = shorthand_to_face.find(shorthand);
		// std::cout << shorthand << std::endl;
		if (got == shorthand_to_face.end()){
			
			auto ret = shorthand_to_face.emplace(shorthand, std::list<FaceWithDirection>{});
			assert(ret.second && "No duplicate face shorthand names.");

			ret.first->second.emplace_back(face_to_add);

		}
		else { 	// if the lookup is successful, we add it to the list. 
			std::list<FaceWithDirection>  face_list = got->second;
			face_list.emplace_back(face_to_add);
			got->second = face_list;
			std::cout << "face name: " << face_to_add.face.name << "face sh: " << face_to_add.face.shorthand << " amount: " << face_list.size() << std::endl;


		
		}

	}
	std::cout << "library load successful" << std::endl;
  	std::unordered_map<std::string,std::list<FaceWithDirection> >::iterator got = shorthand_to_face.find("dc");
	if (got != shorthand_to_face.end()){
		
		for (auto f : got->second){
			std::cout << f.face.name << std::endl;
			

		}
	}


	//library templates, to be used in layout:
	std::vector< std::vector< glm::vec2 > > templates;

	sm::Mesh mesh;

	

	//----------------------------

	{ //parse text file:
		std::unordered_map< std::string, uint32_t > name_to_L;
		auto get_L = [&mesh,&name_to_L,&templates](sm::Library::Face const &face) -> uint32_t {
			auto f = name_to_L.find(face.name);
			if (f != name_to_L.end()) return f->second;
			name_to_L.insert(std::make_pair(face.name, mesh.library.size()));
			mesh.library.emplace_back(face.key());
			templates.emplace_back();
			for (auto const &e : face.edges) {
				templates.back().emplace_back(e.vertex);
			}
			return mesh.library.size()-1;
		};

		std::unordered_map< std::string, uint32_t > open_edges;
		
		
		// the first line is the last set of stitches,
		// so reverse the line order for building from the bottom up
		std::ifstream text(in_text);
		std::string line;
		std::list<std::list<std::string>> reverse_lines;
		// std::vector< std::string > toks;
		while (std::getline(text,line)){
			{ //trim comments:
				auto idx = line.find('#');
				if (idx != std::string::npos) line = line.substr(0,idx);
			}
			std::vector< std::string > toks;
			std::list<std::string> l;
			{
				std::string tok;
				std::istringstream str(line);
				while (str >> tok) {
					l.emplace_back(tok);
				}
			}
			if (l.size() == 0) continue; //row is empty, so we skip
			// for (std::string tok : toks)
			reverse_lines.emplace_front(l);
		}
		for (std::list<std::string> line: reverse_lines){
			for (std::string sh : line){
				std::cout << sh << std::endl;

			}
		}

		
		unsigned int rowIndex = 1; 
		std::list< sm::Mesh::Connection> previous_row_top_connections;
		std::list< sm::Mesh::Connection> current_row_top_connections;
		for (; rowIndex <= reverse_lines.size(); rowIndex++) {
			auto l = std::next(reverse_lines.begin(),rowIndex-1);
			std::list<std::string> line(*l);
			std::cout << "line size " << line.size() << std::endl; 
			uint32_t previous_face = -1U;
			uint32_t previous_right_edge = -1U;
			uint32_t previous_left_edge = -1U;
			current_row_top_connections.clear();
			



			for(unsigned int i = 0; i < l->size(); i++) {
				//capture the stitch
				std::string stitch_token = *line.begin();
				line.erase(line.begin());

				auto f = shorthand_to_face.find(stitch_token);
				if (f == shorthand_to_face.end()) {
					std::cerr << "ERROR: face '" << stitch_token << "' does not appear in library." << std::endl;
					return 1;
				}
				// sm::Library::Face const &face = *f->second;
				// AD: from the face with direction list, find the right face based on the state
				// if a right handed crocheter: 
				// the 
				std::list<FaceWithDirection> const face_with_direction_list = f->second;

				//we will do only one type of face for now
				sm::Library::Face const face = f->second.front().face;

				mesh.faces.emplace_back();
				sm::Mesh::Face &m_face = mesh.faces.back();
				m_face.type = get_L(face);

				for (uint32_t i = 0; i < face.edges.size(); ++i) {
					m_face.emplace_back(mesh.vertices.size());
					mesh.vertices.emplace_back(std::numeric_limits< float >::quiet_NaN());
				}
				uint32_t current_face = mesh.faces.size()-1;
				//edge 1 is right
				//edge 3 is left
				//maybe theres a better way to get this? 
				uint32_t left_edge = 3;
				uint32_t right_edge = 1; 
				uint32_t top_edge = 2;
				uint32_t bottom_edge = 0;


				if (rowIndex % 2 == 1){ // odd
					if (previous_face != -1U) {
						mesh.connections.emplace_back();
						sm::Mesh::Connection &con = mesh.connections.back();

						con.a.face = previous_face;
						con.a.edge = previous_right_edge;
						con.b.face = current_face;
						con.b.edge = left_edge;
						con.flip = true;
					}
					previous_face = current_face;
					previous_right_edge = right_edge; 

				}
				else { //even
					if (previous_face != -1U) {
						mesh.connections.emplace_back();
						sm::Mesh::Connection &con = mesh.connections.back();

						con.a.face = previous_face;
						con.a.edge = previous_left_edge;
						con.b.face = current_face;
						con.b.edge = right_edge;

						con.flip = true;
					}
					previous_face = current_face;
					previous_left_edge = left_edge; 
				}
				current_row_top_connections.emplace_back();
				sm::Mesh::Connection &con_top = current_row_top_connections.back();
				con_top.a.face = current_face;
				con_top.a.edge = top_edge;

				con_top.flip = true;

				if (rowIndex >= 2 && !previous_row_top_connections.empty()){ //its not the first row so we should connect to the row below. 

	
					sm::Mesh::Connection &con_bot = previous_row_top_connections.front();

					con_bot.b.face = current_face;
					con_bot.b.edge = bottom_edge;
					// con.flip = true;
					std::cout << "A: "<< con_bot.a.face << " B: "<<con_bot.b.face <<std::endl;
					mesh.connections.emplace_back(con_bot);
					previous_row_top_connections.erase(previous_row_top_connections.begin());

				}
			}
			previous_row_top_connections = current_row_top_connections;
	
			
		}

		if (!open_edges.empty()) {
			std::cerr << "WARNING: the following open edges are not connected:\n";
			for (auto const &oe : open_edges) {
				std::cerr << "   " << oe.first << "\n";
			}
			std::cerr.flush();
		}
		uint32_t removed = 0;
		for (uint32_t i = 0; i < mesh.connections.size(); /* later */) {
			if (mesh.connections[i].b.face == -1U) {
				std::swap(mesh.connections[i], mesh.connections.back());
				mesh.connections.pop_back();
				++removed;
			} else {
				++i;
			}
		}
		if (removed != 0) {
			std::cout << "WARNING: trimmed " << removed << " incomplete connections." << std::endl;
		}
		assert(removed == open_edges.size());
	}

	std::cout << "Created mesh with " << mesh.vertices.size() << " vertices, " << mesh.faces.size() << " faces, and " << mesh.connections.size() << " connections." << std::endl;

	//----------------------------

	//center templates:
	for (auto &t : templates) {
		glm::vec2 avg = glm::vec2(0.0f);
		for (auto const &p : t) {
			avg += p;
		}
		avg /= t.size();
		for (auto &p : t) {
			p -= avg;
		}
	}


		{ //basic (2d) position/orientation relaxation for positions:
		assert(templates.size() == mesh.library.size());

		std::vector< glm::vec2 > pts;
		struct Face : std::vector< uint32_t > {
			uint32_t type = -1U;
			float angle = 0.0f;
			glm::vec2 center = glm::vec2(0.0f);
		};
		std::vector< Face > faces;

		{ //set up pts/faces from mesh vertices and connections:
			std::vector< uint32_t > verts_to_ids;
			verts_to_ids.reserve(mesh.vertices.size());
			for (uint32_t i = 0; i < mesh.vertices.size(); ++i) {
				verts_to_ids.emplace_back(i);
			}

			auto unify = [&verts_to_ids](uint32_t i, uint32_t j){
				assert(i < verts_to_ids.size());
				assert(j < verts_to_ids.size());
				uint32_t a = i;
				while (verts_to_ids[a] != a) a = verts_to_ids[a];
				uint32_t b = j;
				while (verts_to_ids[b] != b) b = verts_to_ids[b];
				if (a < b) {
					verts_to_ids[b] = a;
				} else if (b < a) {
					verts_to_ids[a] = b;
				}
				//NOTE: this is not efficient; efficient code would clean up the other pointers during unification.
			};

			for (auto const &c : mesh.connections) {
				assert(c.a.face < mesh.faces.size());
				assert(c.a.edge < mesh.faces[c.a.face].size());
				uint32_t a0 = mesh.faces[c.a.face][c.a.edge];
				uint32_t a1 = mesh.faces[c.a.face][(c.a.edge + 1)%mesh.faces[c.a.face].size()];
				assert(c.b.face < mesh.faces.size());
				assert(c.b.edge < mesh.faces[c.b.face].size());
				uint32_t b0 = mesh.faces[c.b.face][c.b.edge];
				uint32_t b1 = mesh.faces[c.b.face][(c.b.edge + 1)%mesh.faces[c.b.face].size()];
				if (c.flip) {
					unify(a0, b1); unify(a1, b0);
				} else {
					unify(a0, b0); unify(a1, b1);
				}
			}

			std::vector< uint32_t > ids_to_pts(verts_to_ids.size(), -1U);
			std::vector< uint32_t > verts_to_pts;
			verts_to_pts.reserve(mesh.vertices.size());
			for (uint32_t i = 0; i < mesh.vertices.size(); ++i) {
				uint32_t id = i;
				while (verts_to_ids[id] != id) id = verts_to_ids[id];
				assert(id < ids_to_pts.size());
				if (ids_to_pts[id] == -1U) {
					ids_to_pts[id] = pts.size();
					pts.emplace_back(glm::vec2(0.0f));
				}
				verts_to_pts.emplace_back(ids_to_pts[id]);
			}

			faces.reserve(mesh.faces.size());
			for (auto const &f : mesh.faces) {
				faces.emplace_back();
				faces.back().type = f.type;
				for (auto const &v : f) {
					faces.back().emplace_back(verts_to_pts[v]);
				}
			}
		} //end of setup

		//brain-dead iterative solve:
		for (uint32_t pass = 0; pass < 2; ++pass) {
			std::vector< glm::vec2 > old_pts = pts;
			bool converged = false;

			for (uint32_t iter = 0; iter < 100000; ++iter) {

				//update point positions:
				std::vector< glm::vec3 > votes(pts.size(), glm::vec3(0.0f));
				for (auto &f : faces) {
					//regularization:
					if (pass == 0) {
						f.angle = 0.0f; //<-- HACK: ignore rotation during first pass
					} else {
						f.angle *= 0.999f;
					}
	
					glm::vec2 right = glm::vec2(std::cos(f.angle), std::sin(f.angle));
					glm::vec2 up = glm::vec2(-right.y, right.x);
					std::vector< glm::vec2 > const &t = templates[f.type];
					assert(t.size() == f.size());
					for (uint32_t i = 0; i < t.size(); ++i) {
						votes[f[i]] += glm::vec3(right * t[i].x + up * t[i].y + f.center, 1.0f);
					}
				}

				std::swap(old_pts, pts); //store old point positions

				glm::vec2 pts_center = glm::vec2(0.0f);
				for (uint32_t i = 0; i < pts.size(); ++i) {
					pts[i] = glm::vec2(votes[i]) / votes[i].z;
					pts_center += pts[i];
				}
				pts_center /= pts.size();
				for (auto &p : pts) {
					p -= pts_center;
				}
				double absolute_movement = 0.0;
				for (uint32_t i = 0; i < pts.size(); ++i) {
					absolute_movement += std::abs(pts[i].x - old_pts[i].x);
					absolute_movement += std::abs(pts[i].y - old_pts[i].y);
				}
				if (absolute_movement < 1e-4) {
					std::cout << "Pass " << (pass+1) << " converged at iteration " << iter << "." << std::endl;
					converged = true;
					break;
				}
	
				for (auto &f : faces) {
					//update face centers:
					glm::vec2 avg = glm::vec2(0.0f);
					for (auto i : f) {
						avg += pts[i];
					}
					f.center = avg /= f.size();
	
					//determine angles:
					//Port of (javascript) code I used in a lunch talk once:
					float coefCos = 0.0;
					float coefSin = 0.0;
					float coefSinCos = 0.0;
					float coefCos2 = 0.0;
					float coefSin2 = 0.0;
					auto addDeriv = [&](float c, float s, float o) {
						//add derivative of (c * cos(t) + s * sin(t) + o * 1)^2
						//==  2 * (c * cos(t) + s * sin(t) + o * 1) * (c *-sin(t) + s * cos(t))
						coefCos += o * s;
						coefSin += o * -c;
						coefCos2 += c * s;
						coefSin2 += s * -c;
						coefSinCos += (c * -c + s * s);
					};
					std::vector< glm::vec2 > const &t = templates[f.type];
					for (uint32_t i = 0; i < f.size(); ++i) {
						// (offsets[i].x * cos(t) + offsets[i].y *-sin(t) - [real offset].x)^2
						//+(offsets[i].x * sin(t) + offsets[i].y * cos(t) - [real offset].y)^2
						addDeriv(t[i].x, -t[i].y, -pts[f[i]].x + f.center.x);
						addDeriv(t[i].y,  t[i].x, -pts[f[i]].y + f.center.y);
					}
					//end up with *just* coefCos, coefSin.
					if (std::abs(coefSinCos) > 1e-6 || std::abs(coefCos2) > 1e-6 || std::abs(coefSin2) > 1e-6) {
						std::cout << "Unexpected coefs: " << coefCos << " " << coefSin << " " << coefSinCos << " " << coefCos2 << " " << coefSin2 << std::endl;
					}
					if (coefCos == 0.0 && coefSin == 0.0) {
						f.angle = 0.0f;
					} else {
						//Want:
						//coefCos * cos + coefSin * sin = 0.0 [deriv == 0]
						//also:
						//coefSin * cos + -coefCos * sin > 0 [second deriv > 0]
		
						//let len = Math.sqrt(coefCos * coefCos + coefSin * coefSin);
						//f.right.x = coefSin / len;
						//f.right.y =-coefCos / len;
						f.angle = std::atan2(-coefCos, coefSin);
					}
				}
			} //for (iteration)
			if (!converged) {
				std::cout << "Pass " << (pass+1) << " did not converge; continuing anyway." << std::endl;
			}
		} //for (pass)

		//copy results back to mesh:
		assert(faces.size() == mesh.faces.size());
		for (uint32_t f = 0; f < faces.size(); ++f) {
			assert(mesh.faces[f].size() == faces[f].size());
			for (uint32_t i = 0; i < faces[f].size(); ++i) {
				mesh.vertices[mesh.faces[f][i]] = glm::vec3(pts[faces[f][i]], 0.0f);
			}
		}

	}


	std::cout << "Writing '" << out_smobj << "'." << std::endl;
	mesh.save(out_smobj);

	return 0;
}
