#include "sm.hpp"
#include "hinters.hpp"

// generate hints assuming:
// 1. Patch has only short-rows <-- verify face signatures only consume/produce 1 loop
// 2. Is entirely on (same) bed <-- verify code faces are the same (all "f" or all "b")
//
sm::Mesh sm::hint_shortrow_only_patch(sm::Mesh const &mesh, sm::Library const &face_library, sm::Code const &code_library){

	sm::Mesh out = mesh;
	
	std::map<std::string, uint32_t> name_to_code_idx;
	std::map<std::string, uint32_t> name_to_lib_idx; // is this even necessary?
	for(auto const &c : code_library.faces){
		name_to_code_idx[c.key()] = &c - &code_library.faces[0];
	}
	for(auto const &c : face_library.faces){
		name_to_lib_idx[c.key()] = &c - &face_library.faces[0];
	}
	int needle = 0;
	char bed = 'f';
	
	for(auto const &f : mesh.faces){
		uint32_t fid = &f - &mesh.faces[0];

		if(!name_to_code_idx.count(mesh.library[f.type])){
			std::cerr << "Face signature does not exist in code library " << mesh.library[f.type] << std::endl;
			break;
		}
		if(!name_to_lib_idx.count(mesh.library[f.type])){
			std::cerr << "Face signature does not exist in face library " << mesh.library[f.type] << std::endl;
			break; 
		}
		auto const &code = code_library.faces[name_to_code_idx[mesh.library[f.type]]];
		auto const &face = face_library.faces[name_to_lib_idx[mesh.library[f.type]]];
		// get edges from the library
		// get template from code
		// if edge is loop edge
		for(auto &e : code.edges){
			auto eid = &e - &code.edges[0];
			//assign hints
			sm::Mesh::FaceEdge fe;
			fe.face = fid; fe.edge = eid;
			sm::Mesh::Hint hint;
			hint.fe = fe;
			hint.bn.bed = e.bn.bed;
			hint.bn.needle = e.bn.needle + needle;
			out.location_hints.emplace_back(hint);
		}
		// update needle based on direction of the yarn
		// TODO (stopped here..)
	}
	return out;
}
