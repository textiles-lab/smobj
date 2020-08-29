#include "sm.hpp"
#include "hinters.hpp"
#include <iostream>

// generate hints assuming:
// 1. Patch has only short-rows <-- verify face signatures only consume/produce 1 loop
// 2. Is entirely on (same) bed <-- verify code faces are the same (all "f" or all "b")
//
sm::Mesh sm::hint_shortrow_only_patch(sm::Mesh const &mesh, sm::Code const &code_library){

	sm::Mesh out = mesh;
	
	std::map<std::string, uint32_t> name_to_code_idx;
	for(auto const &c : code_library.faces){
		name_to_code_idx[c.key()] = &c - &code_library.faces[0];
	}
	{
		// check that all faces have an associated code face
		for(auto const &f : mesh.faces){
			auto name = mesh.library[f.type];
			if(!name_to_code_idx.count(name)){
				std::cerr << "Face with name " << name << " does not exist in the code library, cannot hint. " << std::endl;
				return out;
			}
		}
	}
	std::map<sm::Mesh::FaceEdge, sm::BedNeedle> hints;
	char bed = 'x';
	{
		// assign hints for face 0, from the template directly
		auto l = code_library.faces[name_to_code_idx[mesh.library[mesh.faces[0].type]]];
		for(auto &e : l.edges){
			sm::Mesh::FaceEdge fe;
			fe.face = 0;
			fe.edge = &e - &l.edges[0];
			hints[fe] = e.bn;
			if(e.bn.bed != 'x') bed = e.bn.bed;
		}
		//std::cout << "Starting bed " << bed << std::endl;
		for(auto &h : hints) h.second.bed = bed;
	}

	// verify that pattern is singe-bed based
	for(auto const &f : mesh.faces){
		
		auto l = code_library.faces[name_to_code_idx[mesh.library[mesh.faces[&f - &mesh.faces[0]].type]]];
		for(auto const &e : l.edges){
			if(e.bn.bed != 'x' && e.bn.bed != bed){
				std::cerr << "Faces don't all lie on the same bed, use a different hinter" << std::endl;
				return out;
			}
		}
	}

	bool did_update = true;
	while( did_update ){
		did_update = false;
		for(auto const &c : mesh.connections){
			if(hints.count(c.a) ){
				if(!hints.count(c.b)){
					hints[c.b] = hints[c.a];
					// now assign the rest of c.b
					auto l = code_library.faces[name_to_code_idx[mesh.library[mesh.faces[c.b.face].type]]];
					int offset = -l.edges[c.b.edge].bn.location() + hints[c.b].location();
					for(uint32_t i = 0; i < l.edges.size(); ++i){
						sm::Mesh::FaceEdge fe;
						fe.face = c.b.face;
						fe.edge = i;
						hints[fe] = l.edges[i].bn;
						hints[fe].needle += offset;
						hints[fe].bed = bed;
					}
					
					did_update = true;
				}
				else{
					//std::cout << "B " << hints[c.b].to_string() << " nudge " << int(hints[c.b].nudge) << " A " << hints[c.a].to_string() <<  " nudge " << int(hints[c.a].nudge) << std::endl;
				//assert(hints[c.b] == hints[c.a] &&  "Connected edges should agree if no pentagons exist");
				}
			}
			else if(hints.count(c.b)){
				if(!hints.count(c.a)){
					did_update = true;
					hints[c.a] = hints[c.b];
					// now assign the rest of c.a
					auto l = code_library.faces[name_to_code_idx[mesh.library[mesh.faces[c.a.face].type]]];
					int offset = -l.edges[c.a.edge].bn.location() + hints[c.a].location();
					for(uint32_t i = 0; i < l.edges.size(); ++i){
						sm::Mesh::FaceEdge fe;
						fe.face = c.a.face;
						fe.edge = i;
						hints[fe] = l.edges[i].bn;
						hints[fe].needle += offset;
						hints[fe].bed = bed;
					}
				}
				else{
					//std::cout << "A " << hints[c.a].to_string() << " nudge " << int(hints[c.a].nudge) << " B " << hints[c.b].to_string() << " nudge " << int(hints[c.b].nudge) << std::endl;
				//assert(hints[c.b] == hints[c.a] &&  "Connected edges should agree if no pentagons exist");
				}
			}
		}
		if(!did_update){
			break;
		}
	}
	// clear old hints??
	out.hints.clear();
	for(auto h : hints){
		sm::BedNeedle bn;
		out.hints.emplace_back();
		out.hints.back().lhs = h.first;
		bn.bed = h.second.bed;
		bn.needle = h.second.needle;
		bn.nudge = h.second.nudge;
		out.hints.back().rhs = bn;
		out.hints.back().type = sm::Mesh::Hint::Resource;
		out.hints.back().src = sm::Mesh::Hint::User;
	}
	
	return out;
}
// generate hints assuming:
// 1. Patch has only short-rows <-- verify face signatures only consume/produce 1 loop
// 2. Is entirely on (same) bed <-- verify code faces are the same (all "f" or all "b")
//
sm::Mesh sm::hint_shortrow_only_tubes(sm::Mesh const &mesh, sm::Code const &code_library){
	assert(false && "TODO IMPLEMENT");
	sm::Mesh m = mesh;
	return m;
}
