#include "sm.hpp"
#include "hinters.hpp"
#include <iostream>

// generate hints assuming:
// 1. Patch has only short-rows <-- verify face signatures only consume/produce 1 loop
// 2. Is entirely on (same) bed <-- verify code faces are the same (all "f" or all "b")
//

sm::Mesh sm::constraint_assign_frontface_variant(sm::Mesh const &mesh, sm::Code const &code_library){
	sm::Mesh out = mesh;

	for(auto &f : out.faces){
		std::string name = out.library[f.type];
		std::string variant = "";
		bool found_variant = false;
		for(auto c : code_library.faces){
				if(c.key_library() == name){
					bool has_back = false;
					for(auto e : c.edges){
						if(e.bn.bed == 'b') has_back = true;
					}
					if(!has_back){
						variant = c.variant;
						found_variant = true;
					}
				}
		}
		if(found_variant){
			sm::Mesh::Hint h;
			h.type = sm::Mesh::Hint::Variant;
			h.lhs.face = &f - &out.faces[0];
			h.rhs = variant;
			h.src = sm::Mesh::Hint::Heuristic;
			// TODO if hint already exists, replace else append...

			out.hints.emplace_back(h);
		}
	}
	return out;
}

sm::Mesh sm::hint_shortrow_only_patch(sm::Mesh const &mesh, sm::Code const &code_library){
	// for every face, assign a variant hint that is the first "all front bed " variant, and use that variant

    sm::Mesh out = mesh;
    out.hints.clear();
    out = constraint_assign_frontface_variant(out, code_library);
    std::map<std::string, uint32_t> name_to_code_idx;

	for(auto const &c : code_library.faces){
        // key -> code library name
        // value -> index
		name_to_code_idx[c.key()] = &c - &code_library.faces[0];
	}

	auto find_face_variant = [&](uint32_t fidx)->std::string{
		// go through hints of type Variants
		// if hint.lhs.face == fidx, return std::get<std::string>(h.rhs)
		for(auto &h : out.hints) {
			if (h.type == sm::Mesh::Hint::Variant && h.lhs.face == fidx)
		    return std::get<std::string>(h.rhs);
		}
		return "";
	};

	std::map<sm::Mesh::FaceEdge, sm::BedNeedle> hints;
	char bed = 'x';
	{
		// assign hints for face 0, from the template directly
		std::string lib_name = mesh.library[out.faces[0].type];
		std::string variant = find_face_variant(0);
		std::string code_name = lib_name + " " + variant;
		// ideally
		std::cout << "Code name " << code_name << std::endl;
		assert(name_to_code_idx.count(code_name) && "code name not in name_to_code_idx");
		auto l = code_library.faces[name_to_code_idx[code_name]];
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
		// mesh.library <-- has a list of "strings" that indicate the sf library face signatures of all the faces used in smobj
		// Smobj::Mesh::Face::type is an index into mesh.library
		// mesh.library[f.type] gives you the geometry sf library's face signature for face 'f'
		// the code signature for face f is mesh.library[f.type] + " " + variant_name
		std::string lib_name = mesh.library[mesh.faces[&f - &mesh.faces[0]].type];
		// check if there is a variant associated with this face
		// if there is no variant, find a front bed variant and then hook into that..
		std::string variant = find_face_variant(&f - &mesh.faces[0]);
		std::string code_name = lib_name + " " + variant;
 		auto l = code_library.faces[name_to_code_idx[code_name]];

        for(auto const &e : l.edges){
            // Why face beds are assigned even if there are no constraints?
            std::cout << "code name " << l.key() << " lib name " << l.key_library() << std::endl;
            std::cout << "bed: " << bed << std::endl;
            std::cout << "(code template) e.bn.bed: " << e.bn.bed << std::endl;
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
					did_update = true;
					hints[c.b] = hints[c.a];
					// now assign the rest of c.b
					std::string lib_name = mesh.library[mesh.faces[c.b.face].type];
					std::string variant = find_face_variant(c.b.face);
					std::string code_name = lib_name + " " + variant;
					auto l = code_library.faces[name_to_code_idx[code_name]];

					int offset = -l.edges[c.b.edge].bn.location() + hints[c.b].location();
					for(uint32_t i = 0; i < l.edges.size(); ++i){
						sm::Mesh::FaceEdge fe;
						fe.face = c.b.face;
						fe.edge = i;
						hints[fe] = l.edges[i].bn;
						hints[fe].needle += offset;
						hints[fe].bed = bed;
					}
				}
			}
			else if(hints.count(c.b)){
				if(!hints.count(c.a)){
					did_update = true;
					hints[c.a] = hints[c.b];
					// now assign the rest of c.a
					std::string lib_name = mesh.library[mesh.faces[c.b.face].type];
					std::string variant = find_face_variant(c.b.face);
					std::string code_name = lib_name + " " + variant;
					auto l = code_library.faces[name_to_code_idx[code_name]];

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
			}
		}
		if(!did_update){
			break;
		}
	}
	// clear old hints??
	//out.hints.clear();
	for(auto h : hints){
		sm::BedNeedle bn;
		out.hints.emplace_back();
		out.hints.back().lhs = h.first;
		bn.bed = h.second.bed;
		bn.needle = h.second.needle;
		bn.nudge = h.second.nudge;
		out.hints.back().rhs = bn;
		out.hints.back().type = sm::Mesh::Hint::Resource;
		out.hints.back().src = sm::Mesh::Hint::Inferred;
	}

	// If order was correct, then for this style of transfer free smobjs, the face order is the instruction ordering
	//	std::vector< std::pair<uint32_t, uint32_t> > total_order;

	{
		out.total_order.clear();
		for(uint32_t i = 0; i < out.faces.size(); ++i){
			std::string lib_name = mesh.library[mesh.faces[i].type];
			std::string variant = find_face_variant(i);
			std::string code_name = lib_name + " " + variant;
			auto l = code_library.faces[name_to_code_idx[code_name]];
			for(uint32_t j = 0; j < l.instrs.size(); ++j){
				if(l.instrs[j].op == sm::Code::Face::Instr::Operation::Xfer){
					std::cerr << "This constraint generator should not be used with code that introduces transfers." << std::endl;
					// todo throw
				}
				out.total_order.emplace_back(std::make_pair(i,j));
			}
		}
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


// (a) if all edge in a face has a resource hints and the only reasonable code in .code has one variant, then assign that
bool sm::constraint_assign_variant_from_resource(sm::Mesh &mesh, sm::Code const &code_library){
    // return this
    bool is_changed = false;

	auto find_edge_resource = [&](uint32_t fidx, uint32_t eidx)->sm::BedNeedle{
		// go through resource hints
		for(auto &h : mesh.hints) {
			if (h.type == sm::Mesh::Hint::Resource && h.lhs.face == fidx && h.lhs.edge == eidx)
                return std::get<sm::BedNeedle>(h.rhs);
		}
        sm::BedNeedle bn;
		return bn;
	};

    std::map<std::string, std::vector<uint32_t> > name_to_code_idx;
	for(auto const &c : code_library.faces){
        // What to do if there is a same library name with different variants?
        std::vector<uint32_t> &vec = name_to_code_idx[c.key_library()];
		vec.push_back(&c - &code_library.faces[0]);
	}

    // If all edge in a face has resource hints
    // Check if resource hints are consistant
	for(auto &f : mesh.faces){
		std::string lib_name = mesh.library[f.type];
        uint32_t face_id = &f - &mesh.faces[0];

        // If this face already has variant hint, return
        bool has_variant = false;
		for(auto &h : mesh.hints) {
			if (h.type == sm::Mesh::Hint::Variant && h.lhs.face == face_id) {
                has_variant = true;
            }
		}
        if (has_variant) continue;

        // If the code of this face has more than two variants, return
        std::vector<uint32_t> code_idx = name_to_code_idx[lib_name];
        if (code_idx.size() > 1) continue;

        // TODO: Also check needle number???
        // Assign first edge bed to `first`
        sm::BedNeedle first = find_edge_resource(face_id, 0);
        if (first.bed == 'x') continue;
        bool is_same = true;

        // Iterate through edges of this face and is_same = false if
        // the bed doesn't match with `first`
		auto l = code_library.faces[code_idx[0]];
		for(auto &e : l.edges){
			uint32_t eidx = &e - &l.edges[0];
            sm::BedNeedle bn = find_edge_resource(face_id, eidx);
            if (first.bed != bn.bed) {
                is_same = false;
                break;
            }
		}

        // If is_same == true, it means that all edges of f had resource hints and its bed was `first`.
        // So assign the only variant in the code as variant hint and mark it has Inferred
        if (is_same) {
            is_changed = true;

            sm::Mesh::Hint h;
			h.type = sm::Mesh::Hint::Variant;
			h.lhs.face = &f - &mesh.faces[0];
            // Get variant from the only face that has a variant
			h.rhs = l.variant;
			h.src = sm::Mesh::Hint::Inferred;
			mesh.hints.emplace_back(h);
        }
    }

    return is_changed;
}

bool sm::constraint_extend_resource_from_resource(sm::Mesh  &mesh, sm::Code const &code_library){
    return false;
	//Inferred

	// Inference rules
	// - Resource
	//   (a) one resource hints exist for one edge and you want to infer everything else for that face
	//   (b) one edge A has a resource constraint and this edge has a connection to some other edge B, then assign B the constraint of A
	//   (c) if there is a resource conflict, then remove every infered constraints along that path, and warn the user that something's wrong
}

bool sm::constraint_face_instruction_order(sm::Mesh &mesh, sm::Code const &code_library){
    return false;
// Eventually, something here
// Maybe from a known database ? From instruction ordering
}

sm::Mesh sm::infer_constraints(sm::Mesh &mesh, sm::Code const &code_library) {
    bool flag = true;
    while (flag) {
        flag = false;
		flag = constraint_extend_resource_from_resource(mesh, code_library);
        flag = constraint_assign_variant_from_resource(mesh, code_library);
		flag =  constraint_face_instruction_order(mesh, code_library);
	}

    return mesh;
}
