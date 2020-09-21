#include "sm.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <set>

// assumes that the given smobj is a single connected component...
int main(int argc, char* argv[]){
	if(argc < 5){
		std::cout << "Usage: ./merge-faces lib.sf face1 face2 out.sf" << std::endl;
		return 1;
	}
	std::string in_sf = argv[1];
	std::string in_smobj  = argv[2];
	std::string out_sf = argv[3];
	std::string out_smobj = argv[4];
	
	//load in_sf library, load in_smobj, merge all faces in it based on connections..
	// along any connection, turn an edge point into a 3D spline position
	// update the end points' edges
	//generate a new library face and a new smobj
	
	sm::Library library = sm::Library::load(in_sf);
	sm::Mesh smobj = sm::Mesh::load(in_smobj);
	std::cout << "Library has " << library.faces.size() << " faces." << std::endl;
	std::cout << "Smobj has " << smobj.faces.size() << " faces." << std::endl;

	std::map< sm::Mesh::FaceEdge, sm::Mesh::FaceEdge> old_fe_to_new_fe;
	std::set<uint32_t> deleted_faces;

	std::map<std::string, uint32_t> library_name_to_index;
	for(auto const &l : library.faces){
		library_name_to_index[l.key()] = &l - &library.faces[0];
	}


	// keep track of connections per face for cheaper handling + dealing with multi connection faces.

	for(auto const &c : smobj.connections){
		// merge along connection c

		sm::Mesh::FaceEdge a = c.a;
		sm::Mesh::FaceEdge b = c.b;
		if(smobj.faces[a.face].type < 0 || smobj.faces[a.face].type >= smobj.library.size()){
			std::cout << "face a's type is out of bounds. " << std::endl;
		}
		if(smobj.faces[b.face].type < 0 || smobj.faces[b.face].type >= smobj.library.size()){
			std::cout << "face b's type is out of bounds. " << std::endl;
		}
		std::string a_name = smobj.library[smobj.faces[a.face].type];
		std::string b_name = smobj.library[smobj.faces[b.face].type];
		if(library_name_to_index.count(a_name) == 0){
			std::cerr << "Face a's lib name not found in library " << a_name << std::endl;
		}
		if(library_name_to_index.count(b_name) == 0){
			std::cerr << "Face b's lib name not found in library " << b_name << std::endl;
		}
		sm::Library::Face la = library.faces[library_name_to_index[a_name]];
		sm::Library::Face lb = library.faces[library_name_to_index[b_name]];
		
		if(a.face == b.face) {
			std::cout << "----- face has multiple connects, TODO resolve this ---- " << std::endl;
			continue;
		}

		std::cout << "Merging face " << b.face << " into " << a.face << std::endl;

		// merge library faces
		{
			sm::Library::Face l;
			l.name = "combo_"+la.name+"_"+lb.name;
			// the face and library share the same pattern, so create edges based on the new edges
			// need to translate based on the previous vertex
			assert(la.edges.size() == smobj.faces[a.face].size());
			assert(lb.edges.size() == smobj.faces[b.face].size());
			// l.edges.emplace_back()
			glm::vec2 translate = 0.5f*(la.edges[a.edge].vertex + la.edges[(a.edge+1)%la.edges.size()].vertex) - 0.5f*(lb.edges[b.edge].vertex + lb.edges[(b.edge+1)%lb.edges.size()].vertex);
			std::map<uint32_t, uint32_t> a_edge_map, b_edge_map;
			{

				for(uint32_t k = 0; k < a.edge; ++k){
					l.edges.emplace_back(la.edges[k]);
					a_edge_map[k] = l.edges.size()-1;
				}
				for(uint32_t k = 0; k < smobj.faces[b.face].size()-1; ++k){
					l.edges.emplace_back(lb.edges[(b.edge + k + 1)%smobj.faces[b.face].size()]);
					l.edges.back().vertex += translate;
					b_edge_map[ (b.edge + k + 1)%smobj.faces[b.face].size()] = l.edges.size()-1;
				}
				for(uint32_t k = a.edge + 1; k < smobj.faces[a.face].size(); ++k){
					l.edges.emplace_back(la.edges[k]);
					a_edge_map[k] = l.edges.size()-1;
				}
			}
			// translate the midpoints for all the spline points in face b:
			for(auto &y : lb.yarns){
				for(auto &p : y.middle){
					p += glm::vec3(translate, 0.f);
				}
			}
			// any yarn that does not end or begin on a.edge or b.edge, directly copy them
			for(auto ya : la.yarns){
				if(ya.begin.edge != a.edge && ya.end.edge != a.edge){ 
					l.yarns.emplace_back(ya);
					//map the edge to the correct new index
					//l.yarns.back().edge
					l.yarns.back().begin.edge = a_edge_map[l.yarns.back().begin.edge];
					l.yarns.back().end.edge = a_edge_map[l.yarns.back().end.edge];
				}
			}
			for(auto yb : lb.yarns){
				if(yb.begin.edge != b.edge && yb.end.edge != b.edge){
					l.yarns.emplace_back(yb);
					//map the edge to the correct new index
					//l.yarns.back().edge
					l.yarns.back().begin.edge = b_edge_map[l.yarns.back().begin.edge];
					l.yarns.back().end.edge = b_edge_map[l.yarns.back().end.edge];
					
				}
			}
			// merge yarns based on location 
			// std::cout << "Trying to merge yarns between 2 faces. YA = " << la.yarns.size() << " YB = " << lb.yarns.size() << std::endl;
			for(auto ya : la.yarns){
				//std::cout << "ya begin " << ya.begin.edge << " end " << ya.end.edge << std::endl;
				for(auto yb : lb.yarns){
					//std::cout << "\tyb begin " << yb.begin.edge << " end " << yb.end.edge << std::endl;
					if(ya.end.edge == a.edge && yb.begin.edge == b.edge){
						l.yarns.emplace_back(yb);
						l.yarns.back().begin.edge = b_edge_map[l.yarns.back().begin.edge];
						l.yarns.back().end = ya.end;
						l.yarns.back().end.edge = a_edge_map[ya.end.edge];
						// add middle point 3D 
						for(auto m : ya.middle) l.yarns.back().middle.emplace_back(m);
						
						std::cout << "Merge yarn point between yb and ya (1)" << std::endl;
					}
					else if(ya.begin.edge == a.edge && yb.end.edge == b.edge){
						l.yarns.emplace_back(ya);
						l.yarns.back().begin.edge = a_edge_map[l.yarns.back().begin.edge];
						l.yarns.back().end = yb.end;
						l.yarns.back().end.edge = b_edge_map[yb.end.edge];
						// add middle point 3D 
						for(auto m : yb.middle) l.yarns.back().middle.emplace_back(m);
						
						std::cout << "Merge yarn point between ya and yb (2)" << std::endl;
					}
					else if(ya.begin.edge == a.edge && yb.begin.edge == b.edge){
						l.yarns.emplace_back();
						l.yarns.back().begin = ya.end;
						l.yarns.back().begin.edge = a_edge_map[l.yarns.back().begin.edge];
						l.yarns.back().end = yb.end;
						l.yarns.back().end.edge = b_edge_map[l.yarns.back().end.edge];
						
						glm::vec3 p1 = glm::vec3(la.edges[ya.begin.edge].vertex + ya.begin.along * la.edges[(ya.begin.edge+1)%la.edges.size()].vertex, 0.f); 
						glm::vec3 p2 = glm::vec3((lb.edges[yb.begin.edge].vertex + yb.begin.along * lb.edges[(yb.begin.edge+1)%lb.edges.size()].vertex) + translate, 0.f); 
						// add middle point 3D 
						
						for(auto a_it = ya.middle.rbegin(); a_it != ya.middle.rend(); ++a_it){
							l.yarns.back().middle.emplace_back(*a_it);
						}

						l.yarns.back().middle.emplace_back(p1);
						//l.yarns.back().middle.emplace_back(p2);

						for(auto b_it = yb.middle.begin(); b_it != yb.middle.end(); ++b_it){
							l.yarns.back().middle.emplace_back(*b_it);
						}
						std::cout << "Merge yarn point between ya and yb (3)" << std::endl;
					}
					else if(ya.end.edge == a.edge && yb.end.edge == b.edge){
						l.yarns.emplace_back();
						l.yarns.back().begin = ya.begin;
						l.yarns.back().begin.edge = a_edge_map[l.yarns.back().begin.edge];
						l.yarns.back().end = yb.begin;
						l.yarns.back().end.edge = b_edge_map[l.yarns.back().end.edge];
						
						glm::vec3 p1 = glm::vec3(la.edges[ya.end.edge].vertex + ya.end.along * la.edges[(ya.end.edge+1)%la.edges.size()].vertex, 0.f); 
						glm::vec3 p2 = glm::vec3((lb.edges[yb.end.edge].vertex + yb.end.along * lb.edges[(yb.end.edge+1)%lb.edges.size()].vertex) + translate, 0.f); 
						// add middle point 3D 
						
						for(auto a_it = ya.middle.begin(); a_it != ya.middle.end(); ++a_it){
							l.yarns.back().middle.emplace_back(*a_it);
						}

						l.yarns.back().middle.emplace_back(p1);
						//l.yarns.back().middle.emplace_back(p2);

						for(auto b_it = yb.middle.rbegin(); b_it != yb.middle.rend(); ++b_it){
							l.yarns.back().middle.emplace_back(*b_it);
						}					
						std::cout << "Merge yarn point between ya and yb (4)" << std::endl;

					} 
					
				}
			}
			library.faces.emplace_back(l);
		}
		// TODO handle the case when multiple connections between two faces exist...

		// merge smobj faces
		{
			// a <-- b
			deleted_faces.insert(b.face);
			// update all the fe in a.face and b.face
			// 5-6
			// | |
			// .=.
			// | |
			// 1-2
			sm::Mesh::Face face;
			face.type = smobj.library.size();
			uint32_t expected_total = smobj.faces[a.face].size() + smobj.faces[b.face].size()-2;
			old_fe_to_new_fe.clear();
			// a first half:
			for(uint32_t k = 0; k < a.edge; ++k){
				face.emplace_back(smobj.faces[a.face][k]);
				// replace any connection with a.face, k with a.face, face.size()-1
				sm::Mesh::FaceEdge old_fe, new_fe;
				old_fe.face = a.face; new_fe.face = a.face;
				old_fe.edge = k; new_fe.edge = face.size()-1;
				old_fe_to_new_fe[old_fe] = new_fe;
			}
			// b:
			for(uint32_t k = 0; k < smobj.faces[b.face].size()-1; ++k){
				face.emplace_back(smobj.faces[b.face][(b.edge + k + 1)%smobj.faces[b.face].size()]);
				// replace any connection with b.face, (b.edge + k + 1)%smobj.faces[b.face].size()
				// with a.face, face.size()-1
				sm::Mesh::FaceEdge old_fe, new_fe;
				old_fe.face = b.face; new_fe.face = a.face;
				old_fe.edge = (b.edge + k + 1)%smobj.faces[b.face].size();
				new_fe.edge = face.size()-1;
				old_fe_to_new_fe[old_fe] = new_fe;
			}
			// a second half:
			for(uint32_t k = a.edge + 1; k < smobj.faces[a.face].size(); ++k){
				face.emplace_back(smobj.faces[a.face][k]);
				// replace any connection with a.face, k with a.face, face.size()-1
				sm::Mesh::FaceEdge old_fe, new_fe;
				old_fe.face = a.face; new_fe.face = a.face;
				old_fe.edge = k; new_fe.edge = face.size()-1;
				old_fe_to_new_fe[old_fe] = new_fe;
			}
			for(auto &cc : smobj.connections){
				if(old_fe_to_new_fe.count(cc.a)) cc.a = old_fe_to_new_fe[cc.a];
				if(old_fe_to_new_fe.count(cc.b)) cc.b = old_fe_to_new_fe[cc.b];
			}
			std::cout << "total " << face.size() << " expected-total " << expected_total << std::endl;
			assert(face.size() == expected_total);
			smobj.faces[a.face] = face; // update "a"
		}

		// any new 
		// create a new library face and assign the face type of a.face as that library face
		library_name_to_index[library.faces.back().key()] = library.faces.size()-1;
		smobj.library.emplace_back(library.faces.back().key());

		// update library_name_to_index[]

		// create a  combined library face and add it to the end of the list
		// add this face name to smobj.library and update face a's type to point to this type
		// std::cout <<"... done merging two faces. " << std::endl;
	}
	std::cout << "fixing the rest of the mesh.. " << std::endl;
	// get rid of all connections
	smobj.connections.clear();
	// faces don't get used anywhere else, so just delete faces to be deleted:
	smobj.faces.erase(std::remove_if(smobj.faces.begin(), 
				smobj.faces.end(),
				[&deleted_faces, &smobj](sm::Mesh::Face const &f){return deleted_faces.count(&f - &smobj.faces[0]);}),
			smobj.faces.end());

	// ensure that there is only one face left
	{
		if(smobj.faces.size() != 1){
			std::cerr << "Faces have not been reduced to one face! Total faces = " << smobj.faces.size() << std::endl;
		}
	}

	smobj.faces[0].type = 0;
	smobj.library.clear();
	std::string label = library.faces.back().key();
	smobj.library.emplace_back(label);

	// get rid of any unused vertices and re-index
	{
		
		auto vertices = smobj.vertices;
		smobj.vertices.clear();
		
		for(uint32_t x = 0; x < smobj.faces[0].size(); ++x){
			smobj.vertices.emplace_back(vertices[smobj.faces[0][x]]);
			smobj.faces[0][x] = x;
		}

	}

	// TODO make sure hints are updated, for now erased..
	smobj.hints.clear();

	// generate a template face based on smobj face, keep track edge symbols and yarn values
	smobj.save(out_smobj);
	library.faces.erase(library.faces.begin(), library.faces.begin() + library.faces.size()-1);
	library.save(out_sf);
	
	std::cout << "Input .sf " << in_sf << " smobj " << in_smobj << " Output .sf " << out_sf << " smobj " << out_smobj << std::endl;
	
	return 0;
}
