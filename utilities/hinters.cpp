#include "sm.hpp"
#include "hinters.hpp"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <z3++.h>

// generate hints assuming:
// 1. Patch has only short-rows <-- verify face signatures only consume/produce 1 loop
// 2. Is entirely on (same) bed <-- verify code faces are the same (all "f" or all "b")
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
        assert(name_to_code_idx.count(code_name) && "code name not in name_to_code_idx");
        auto l = code_library.faces[name_to_code_idx[code_name]];
        for(auto &e : l.edges){
            sm::Mesh::FaceEdge fe;
            fe.face = 0;
            fe.edge = &e - &l.edges[0];
            hints[fe] = e.bn;
            if(e.bn.bed != 'x') bed = e.bn.bed;
        }
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
                if(l.instrs[j].op == sm::Instr::Operation::Xfer){
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
sm::Mesh sm::hint_shortrow_only_tubes(sm::Mesh const &mesh, sm::Code const &code_library){
    assert(false && "TODO IMPLEMENT");
    sm::Mesh m = mesh;
    return m;
}

std::pair<sm::BedNeedle, sm::Mesh::Hint::HintSource> find_edge_resource
                                    (sm::Mesh const &mesh, uint32_t fidx, uint32_t eidx) {
    // go through resource hints
    for(auto &h : mesh.hints) {
        if (h.type == sm::Mesh::Hint::Resource && h.lhs.face == fidx && h.lhs.edge == eidx)
            return std::make_pair(std::get<sm::BedNeedle>(h.rhs), h.src);
    }

    sm::BedNeedle bn;
    return std::make_pair(bn, sm::Mesh::Hint::Heuristic);
};

// In the constraint_assign_variant_from_resource, we made sure that if we can infer
// any variant from any resource constraints, then that is already assigned.
// So at this point, we can assume that if a face doesn't have a variant, that
// means this face has no compatible code library code.
// 1. If a variant constraint doesn't exist for this particular face, then continue;
// 2. If it has a variant,
//    Make sure we pick `l` the right face with the right variants
//    Rest is the same

// Inference rules
// - Resource
//   (a) one resource hints exist for one edge and you want to infer everything else for that face
//   (b) one edge A has a resource constraint and this edge has a connection to some other edge B, then assign B the constraint of A
//   (c) if there is a resource conflict, then remove every infered constraints along that path, and warn the user that something's wrong
bool sm::constraint_extend_resource_from_resource(sm::Mesh &mesh, sm::Code const &code_library){

    std::map<std::string, std::vector<uint32_t> > name_to_code_idx;
    for(auto const &c : code_library.faces){
        // What to do if there is a same library name with different variants?
        std::vector<uint32_t> &vec = name_to_code_idx[c.key_library()];
        vec.push_back(&c - &code_library.faces[0]);
    }

    // Accumulate all changes here and assign afterwards
    std::map<sm::Mesh::FaceEdge, std::pair<sm::BedNeedle, sm::Mesh::Hint::HintSource>> constraints;

    // (a)
    // Iterate through edges and obtain BedNeedle
    for(auto &f : mesh.faces){
        std::string lib_name = mesh.library[f.type];
        uint32_t face_id = &f - &mesh.faces[0];

        // If this face doesn't have a variant, continue
        bool has_variant = false;
        std::string face_variant = "";
        for(auto &h : mesh.hints) {
            if (h.type == sm::Mesh::Hint::Variant && h.lhs.face == face_id) {
                 face_variant = std::get<std::string>(h.rhs);
                 has_variant = true;
                 break;
            }
        }
        if (!has_variant){
            continue;
        }

        // Make sure we pick the code face which has the same variant
        // as the face's variant constraint.
        sm::Code::Face l;
        for (auto f_id : name_to_code_idx[lib_name]) {
            l = code_library.faces[f_id];
            if (l.variant == face_variant) break;
        }

        //NOTE: (Maybe) If code template bedneedle assignment can be different across variants, also look for variant if available.
        // Obtain the first edge that resource hint is not 'x' and store it in edge_constraint
        bool found_resource = false;
        std::pair<sm::BedNeedle, uint32_t> edge_constraint;
        int offset = 0;
        auto f_src = sm::Mesh::Hint::Inferred;
        for(auto &e : l.edges) {
            uint32_t eidx = &e - &l.edges[0];
            auto pair = find_edge_resource(mesh, face_id, eidx);
            sm::BedNeedle bn = pair.first;
            sm::Mesh::Hint::HintSource src = pair.second;

            if (bn.bed != 'x') {
                // If this resource constraint is Heuristic, propagated edge's resource constraint should
                // also be Heuristic.
                if (src == sm::Mesh::Hint::Heuristic)
                    f_src = sm::Mesh::Hint::Heuristic;

                // sanity check that e.bn.bed == bn.bed
                edge_constraint.first = bn;
                edge_constraint.second = eidx;
                // offset between e.bn and bn should be the same for all edges in this face
                if(found_resource && offset != (e.bn.location() - bn.location())){
					std::cerr << "Inconsistent resource allocation for face " << face_id <<
                        " edge " << eidx << ", old offset = " << offset << " current = " <<
                        e.bn.location() - bn.location() << std::endl; 

					return false;
				}

				//assert((!found_resource || offset == e.bn.location() - bn.location())
                //        && "Offset is constant across all edge resources, right?");
                offset = e.bn.location() - bn.location();
                found_resource = true;
            }
        }

        // If no resource constraints are found, we can't propagate
        if (!found_resource) continue;

        // Propagate constraints to other edges in this face.
        // Raise error if there is a resource conflict and return false
        for(auto &e : l.edges) {
            uint32_t eidx = &e - &l.edges[0];
            // skip if this edge is edge_constraint edge
            if (eidx == edge_constraint.second) continue;

            sm::BedNeedle bn = find_edge_resource(mesh, face_id, eidx).first;
            float eps = 1e-3f;
            // If there is already a resource assigned to this edge, make sure that it's not conflicting
            if (!bn.dontcare()  && !e.bn.dontcare()) {
                // needle 1, nudge -1 and needle 0, nudge 1 --> compatible since location is 0.5, location returns a float
                // needle check: std::abs(bn.location()- other.location()) < eps)
                float e_location = e.bn.location() - offset; 
                if (bn.bed != e.bn.bed || 
                        std::abs(e_location - bn.location()) > eps) {
                    // Resource conflict!!
                    std::cerr << "Resource conflict while running verifier." << std::endl;
                    std::cerr << "Check resource constraint of edge " << eidx << " in face " << face_id << std::endl;;

                    return false;
                }
            } else {
                // Assign resource constraint to this edge
                sm::Mesh::FaceEdge fe;
                fe.face = face_id;
                fe.edge = eidx;
                // get the bed-needle value from e.bn and add offset to it
                sm::BedNeedle new_bn;
                new_bn.bed = e.bn.bed;
                float location = e.bn.location() - offset;
                new_bn.needle = location;
                if (std::abs(location - new_bn.needle) > eps)
                    new_bn.nudge = (location - new_bn.needle) > eps ? 1 : -1;

                constraints[fe] = std::make_pair(new_bn, f_src);
            }
        }
    }

	// (b)
	for(auto const &c : mesh.connections){
		auto const &lib_str = mesh.library[mesh.faces[c.a.face].type];
		auto const &l = name_to_code_idx[lib_str].front();
		auto const &edge_type = code_library.faces[l].edges[c.a.edge].type;
		
		sm::Mesh::Hint::HintSource src = sm::Mesh::Hint::Inferred;

		if(edge_type[0] == 'y'){ // todo maybe maintian 'l' vs 'y' as an actual enum type.
			// don't propagate over yarn edges since they can go across beds.
			// so always heuristic
			src = sm::Mesh::Hint::Heuristic;
		}

		// TODO apply any xfers associated with this connections in the "out" direction edge to "in" direction edge
		if(constraints.count(c.a) && !constraints.count(c.b)){
			constraints[c.b] = constraints[c.a];
			if(src == sm::Mesh::Hint::Heuristic) constraints[c.b].second = src;
		} else if (constraints.count(c.b) && !constraints.count(c.a)) {
			constraints[c.a] = constraints[c.b];
			if(src == sm::Mesh::Hint::Heuristic) constraints[c.a].second = src;
		}
	}

    // return this
    bool is_changed = false;
    for (auto h : constraints) {
        // TODO: Need more reasonable fix for infinite loop
        bool flag = true;
        for (auto &hint : mesh.hints) {
            if (hint.type == sm::Mesh::Hint::Resource && hint.lhs == h.first && std::get<sm::BedNeedle>(hint.rhs) == h.second.first)
                flag = false;
        }
        if (!flag) continue;

        is_changed = true;
        mesh.hints.emplace_back();
        mesh.hints.back().lhs = h.first;
        sm::BedNeedle bn = h.second.first;
        mesh.hints.back().src = h.second.second;
        mesh.hints.back().rhs = bn;
        mesh.hints.back().type = sm::Mesh::Hint::Resource;
    }

    return is_changed;
}


// 1. inffered constraints
// When all edges in a face is infered or user constraints and if we have only one compatible face in the
// code library, then assign that variable and call it infered constraints
// 2. heuristic
// When some edges in a face are user constraint, or all edges have constraints but are heuristic, or
// if there are more than 2 compatible faces in code library, then assign the first variant And
// call it heuristic.
bool sm::constraint_assign_variant_from_resource(sm::Mesh &mesh, sm::Code const &code_library){

    // Accumulate all changes here and assign afterwards
    std::vector<sm::Mesh::Hint> constraints;

    std::map<std::string, std::vector<uint32_t> > name_to_code_idx;
    for(auto const &c : code_library.faces){
        std::vector<uint32_t> &vec = name_to_code_idx[c.key_library()];
        vec.push_back(&c - &code_library.faces[0]);
    }

    // If all edge in a face has resource hints
    // Check if resource hints are consistant
    for(auto &f : mesh.faces){
        uint32_t face_id = &f - &mesh.faces[0];

        // If this face already has variant hint, continue
        bool has_variant = false;
        for(auto &h : mesh.hints) {
            if (h.type == sm::Mesh::Hint::Variant && h.lhs.face == face_id) {
                has_variant = true;
            }
        }
        if (has_variant) continue;

        // If resource constraint exists for this face,
        // save that constraint
        std::vector<sm::Mesh::Hint> face_constraints;
        // Obtain constraints of this face
        for (auto &h : mesh.hints) {
            if (h.type == sm::Mesh::Hint::Resource && h.lhs.face == face_id)
                face_constraints.push_back(h);
        }
        // If there is no resource constraint for this face, there is no way we can 
        // infer variant.
        if (face_constraints.empty()) continue;

        // Possible case:
        // 1. All edges have resource constraints and are inferred or user.
        //    Only one compatible face in the code library
        //    -> Assign variable constraint as Inferred
        // 2. Only subset of edges have resource constraints
        //    -> Assign variable constraint as Heuristic
        // 3. All edges have resource constraints but some of them are Heuristic
        //    -> Assign variable constraint as Heuristic
        // 4. There are more than 2 compatible face
        //    -> Assign the first code face variant as Heuristic
        //
        // Go through each variant of this face,
        // and see if edge resource type and variant type are compatible
        // If they are compatible, assign the variant constraint.
        std::string lib_name = mesh.library[f.type];
        sm::Mesh::Hint::HintSource f_src =
            (name_to_code_idx[lib_name].size() == 1) ? sm::Mesh::Hint::Inferred
                                                     : sm::Mesh::Hint::Heuristic;


        /*
        Handle ff and bb mixed in one face case
        face 4 edges, bottom edge: f0 top edge: b0 variants: ff, fb, bb
        Case 1: Only bottom edge has resource hint, that hint is say F123 --> assign "ff" ("fb also compatible but second in list")
        Case 2: Bottom edge has resource hint F123, top edge has resource hint B123 --> assign "fb" ("ff" is not compatible because top edge doesn't match)
        */
        for (auto face_code : name_to_code_idx[lib_name]) {
            // Get code face
            auto l = code_library.faces[face_code];

            // If all edge in this face match with resource constraint
            bool e_resource = true;

            // For each edge in this code face
            for (int i = 0; i < (int)l.edges.size(); i++) {
                sm::Code::Face::Edge e = l.edges[i];
                // If edge code is don't care, continue
                if (e.bn.bed == 'x') continue;
                bool has_constraint = false;

                // For all resource constraints that this face has
                for (auto f_constraint : face_constraints) {
                    // If a resource constraint is for this edge
                    if (i == (int)f_constraint.lhs.edge) {
                        // Get constraint content
                        sm::BedNeedle f_bn = std::get<sm::BedNeedle>(f_constraint.rhs);
                        // Code face and resource constraint matches!
                        if (e.bn.bed == f_bn.bed) {
                            has_constraint = true;
                            // If some edges' resource constraint were Heuristic,
                            // Variant constraint that we're trying to construct should
                            // also be Heuristic.
                            if (f_constraint.src == sm::Mesh::Hint::Heuristic)
                                f_src = sm::Mesh::Hint::Heuristic;

                        } else {
                            // If code edge doesn't match with resource constraint,
                            // Skip this code face.
                            e_resource = false;
                        }
                    }
                }

                if (!has_constraint) f_src = sm::Mesh::Hint::Heuristic;
            }

            if (e_resource) {
                sm::Mesh::Hint face_hint;

                sm::Mesh::FaceEdge fe;
                fe.face = face_id;
                face_hint.lhs = fe;
                face_hint.type = sm::Mesh::Hint::Variant;

                face_hint.rhs = l.variant;
                face_hint.src = f_src;

                constraints.push_back(face_hint);

                // Assign the first compatible code face
                break;
            }
        }
    }

    // return this
    bool is_changed = false;
    for (auto h : constraints) {
        /*
         * We don't need this, but just in case
        bool flag = true;
        for (auto &hint : mesh.hints) {
            if (hint.type == h.type && hint.lhs == h.lhs && hint.rhs == h.rhs)
                flag = false;
        }
        if (!flag) continue;
        */

        is_changed = true;
        mesh.hints.push_back(h);
    }

    return is_changed;
}

// From the face dependency (arising from connections), insert partial orders
// For instructions within the face, insert partial orders
bool sm::constraint_assign_order_from_face_order(sm::Mesh &mesh, sm::Code const &code_library, sm::Library const &face_library){

	
	std::map<uint32_t, std::pair<std::string, sm::Mesh::Hint::HintSource> > face_variants;
	std::vector<sm::Mesh::Hint> candidates;
	std::map<sm::Mesh::FaceEdge, uint32_t> face_instr_idx; 

	for(auto const &h : mesh.hints){
		if(h.type == sm::Mesh::Hint::Variant){
			// TODO deal with multiple possible constraints
			std::string variant = std::get<std::string>(h.rhs);
			face_variants[h.lhs.face] = std::make_pair(variant,
                    (h.src != sm::Mesh::Hint::User ?  h.src : sm::Mesh::Hint::Inferred));
		}
		else if(h.type == sm::Mesh::Hint::Order){
			auto rhs = std::get<sm::Mesh::FaceEdge>(h.rhs);
			if(face_instr_idx.count(h.lhs)){
			}
			else{
				face_instr_idx[h.lhs] = face_instr_idx.size();
			}
			if(face_instr_idx.count(rhs)){
			}
			else{
				face_instr_idx[rhs] = face_instr_idx.size();
			}
		}
	}
	sm::InstrGraph graph;
	compute_instruction_graph(mesh, code_library, &graph);
	//std::cout << "Graph: " << graph.nodes.size() << " edges l(" << graph.edge_loops.size() << ") y(" << graph.edge_yarns.size() << ")" << std::endl;
	for(auto pr : graph.edge_loops){
		sm::Mesh::Hint h;
		auto i1 = graph.nodes[pr.first];
		auto i2 = graph.nodes[pr.second];
		h.lhs.face = i1.face_instr.first; h.lhs.edge = i1.face_instr.second;
		sm::Mesh::FaceEdge rhs;
		rhs.face =i2.face_instr.first; rhs.edge = i2.face_instr.second;
		h.rhs = rhs;
		h.type = sm::Mesh::Hint::Order;
		h.src = sm::Mesh::Hint::Inferred;
		//std::cout << "Loop edge between " << i1.face_instr.first << "/" << i1.face_instr.second << "-> " << i2.face_instr.first << "/" <<i2.face_instr.second << std::endl;
		candidates.emplace_back(h);
	}
	for(auto pr : graph.edge_yarns){
		sm::Mesh::Hint h;
		auto i1 = graph.nodes[pr.first];
		auto i2 = graph.nodes[pr.second];
		h.lhs.face = i1.face_instr.first; h.lhs.edge = i1.face_instr.second;
		sm::Mesh::FaceEdge rhs;
		rhs.face =i2.face_instr.first; rhs.edge = i2.face_instr.second;
		h.rhs = rhs;
		h.type = sm::Mesh::Hint::Order;
		h.src = sm::Mesh::Hint::Inferred;
		//std::cout << "Yarn edge between " << i1.face_instr.first << "/" << i1.face_instr.second << "-> " << i2.face_instr.first << "/" <<i2.face_instr.second << std::endl;
		candidates.emplace_back(h);
	}

	// check if these hints are really new
	bool did_update = false;
	std::set<std::pair<uint32_t, uint32_t>> partials; // note: this should not be necessary here, but in general for checking ordering hints will be useful
	for(auto h : mesh.hints){
		if(h.type == sm::Mesh::Hint::Order){
			assert(face_instr_idx.count(h.lhs));
			auto rhs = std::get<sm::Mesh::FaceEdge>(h.rhs);
			assert(face_instr_idx.count(rhs));
			partials.insert(std::make_pair(face_instr_idx[h.lhs], face_instr_idx[rhs]));
		}
	}
	for(auto h : candidates){
		auto it = std::find_if(mesh.hints.begin(), mesh.hints.end(), [&h](sm::Mesh::Hint const &hh)->bool{
				return (hh.type == sm::Mesh::Hint::Order && hh.lhs == h.lhs && std::get<sm::Mesh::FaceEdge>(hh.rhs) == std::get<sm::Mesh::FaceEdge>(h.rhs));
				});

		if(it == mesh.hints.end()){

			auto rhs = std::get<sm::Mesh::FaceEdge>(h.rhs);
			if(!face_instr_idx.count(h.lhs)){
				face_instr_idx[h.lhs] = face_instr_idx.size();
			}
			if(!face_instr_idx.count(rhs)){
				face_instr_idx[rhs] = face_instr_idx.size();
			}
			auto candidate = std::make_pair(face_instr_idx[h.lhs], face_instr_idx[rhs]);
			partials.insert(candidate);
			std::vector<uint32_t> sequence;
			// does including this hint introduce a cyclic constraint?
			if(partial_order_to_sequence(partials, &sequence)){

				mesh.hints.emplace_back(h);
				did_update = true;
			}
			else{
				partials.erase(candidate);
			}
		}

	}

	return did_update;
}

bool sm::constraint_face_instruction_order(sm::Mesh &mesh, sm::Code const &code_library){

// calling some sort of transfer planner that is "verified"
// Transfer planning problem:
// Going from infinite gauge m/c to discrete m/c by introducing xfer/xfer pairs and moving them around
//
  /*
  verifier:
  1) Code-based
    a) come up with a default schedule for a fully connected smobj using real-valued needles (infinite machine gauge)
    b) given a face/instruction order, verify that by introducing paired xfer instructions the face/instruction order can be implemented
    from the infinite schedule
  2) Geometry-based
    a) construct a new face object with yarns from the face/instruction reordered code
    b) if there is no collision between original yarns and new yarns after some relaxation, then okay
  */

    return false;
    // Eventually, something here
    // Maybe from a known database ? From instruction ordering
}

// To make test cases for infer_constraints
// Add some example smobj, 3x3 and test following
// 1. No valid variant exists for a particular resource assignment
//    -> Code does not try to assign random variant
// 2. Edge resource have combination of front and back
//    code:
//    knit f0 
//    xfer f0 b0
//    bottom loop will have f0, top will have b0
//

// Error reporting
// We want to tell users which constraints doesn't make sense to users
// Get a list of incompatible constraints, and somehow report it to users
// Maybe implement a better API which can be called from EditMesh.cpp??
//
// Have selected_faces like infrastructure in EditMesh.cpp
// Add problematic faces/edges to that whenever we call infer_constraints
// Highlight it with red?
//
// If constraint 0 and 1 are incompatible, then 1 is incompatibe with 0


// Helper function that checks if code-template-face is compatible with edge-resources
//{

//}
//

static void generate_smt_phi(sm::Mesh &mesh, sm::Code const &code_library) {

    z3::context context;
    z3::solver solver(context);

    std::vector<std::vector<std::pair<z3::expr, std::string>>> variants; 
    std::vector<std::pair<sm::Mesh::FaceEdge, std::pair<z3::expr, z3::expr>>> resources; 
    std::vector<z3::expr> face_expr;
    for(auto const &f : mesh.faces) {
        std::vector<sm::Code::Face> c_faces;
        uint32_t face_num = &f - &mesh.faces[0];
        std::string face_name = mesh.library[mesh.faces[face_num].type];

        // Obtain all constraints for this face here

        for(auto c : code_library.faces){
            if(c.key_library() == face_name){
                c_faces.push_back(c);
            }
        }

        // Declare variables
        // Assume that edge names are consistant within variants
        std::vector<z3::expr> e_beds;
        std::vector<z3::expr> e_needles;
        for (auto &e : c_faces[0].edges) {
            std::string sig;
            if (e.direction == sm::Code::Face::Edge::Out) {
                sig = "_out_";
            } else {
                sig = "_in_";
            }

            std::string v = "f" + std::to_string(face_num) + sig + e.type;

            std::string b_str = v + "_b";
            z3::expr b_expr = context.bool_const(b_str.c_str());
            e_beds.push_back(b_expr);

            std::string n_str = v + "_n";
            z3::expr n_expr = context.real_const(n_str.c_str());
            e_needles.push_back(n_expr);

            sm::Mesh::FaceEdge fe;
            // TODO: Is this correct?
            fe.face = face_num;
            fe.edge = &e - &c_faces[0].edges[0];
            resources.push_back(std::make_pair(fe, std::make_pair(b_expr, n_expr)));
        }

        std::vector<std::pair<z3::expr, std::string>> variant_s;
        // Face Or
        z3::expr face = context.bool_val(false);
        // For each variants
        for (auto &c : c_faces) {
            z3::expr variant = context.bool_val(true);
            // bed formula
            for (int i = 0; i < (int)e_beds.size(); i++) {
                sm::BedNeedle bn = c.edges[i].bn;
                z3::expr tf = (bn.bed == 'f') ? context.bool_val(true) : context.bool_val(false);
                variant = variant && (e_beds[i] == tf);
            }

            // needle formula
            z3::expr base = context.real_val(std::to_string(c.edges[0].bn.location()).c_str());
            z3::expr e_base = e_needles[0];
            for (int i = 1; i < (int)e_needles.size(); i++) {
                z3::expr bn = context.real_val(std::to_string(c.edges[i].bn.location()).c_str());
                variant = variant && (e_base == e_needles[i] + base - bn);
            }

            variant_s.push_back(std::make_pair(variant, c.variant));
            face = face || variant;
        }

        face_expr.push_back(face);
        variants.push_back(variant_s);
    }

	for(auto const &c : mesh.connections){
        auto face_a = face_expr[c.a.face];
        auto face_b = face_expr[c.b.face];
        solver.add(face_a);
        solver.add(face_b);
        z3::expr a_b_expr = context.bool_const("a_b_expr");
        z3::expr b_b_expr = context.bool_const("b_b_expr");
        z3::expr a_n_expr = context.real_const("a_n_expr");
        z3::expr b_n_expr = context.real_const("b_n_expr");
        // Compare edge
        for (auto const &p : resources) {
            sm::Mesh::FaceEdge fe = p.first;
            z3::expr b_expr = p.second.first;
            z3::expr n_expr = p.second.second;
            if (fe.face == c.a.face && fe.edge == c.a.edge) {
                a_b_expr = b_expr;
                a_n_expr = n_expr;
            } else if (fe.face == c.b.face && fe.edge == c.b.edge) {
                b_b_expr = b_expr;
                b_n_expr = n_expr;
            }
        }
        solver.add(a_b_expr == b_b_expr);
        solver.add(a_n_expr == b_n_expr);
    }

    // Run solver here!
    if (solver.check() == z3::sat) {
        z3::model model = solver.get_model();
        for (auto const &variant : variants) {
            for (auto const &p : variant) {
                z3::expr e = model.eval(p.first);
                if (e.is_true()) {
                    sm::Mesh::Hint h;
                    sm::Mesh::FaceEdge fe;
                    fe.face = &variant - &variants[0];
                    h.lhs = fe;
                    h.type = sm::Mesh::Hint::Variant;
                    h.rhs = p.second;
                    h.src = sm::Mesh::Hint::Inferred;
                    std::cout << fe.face << " "  << p.second << std::endl;

                    mesh.hints.push_back(h);
                }
            }
        }

        for (auto const &p : resources) {
            sm::Mesh::FaceEdge fe = p.first;
            z3::expr bed = model.eval(p.second.first);
            z3::expr needle = model.eval(p.second.second);
            std::cout << fe.face << std::endl;
            float needle_num = std::stof(needle.get_decimal_string(100));

            sm::Mesh::Hint h;
            h.lhs = fe;
            h.type = sm::Mesh::Hint::Resource;
            h.src = sm::Mesh::Hint::Inferred;

            sm::BedNeedle bn;
            bn.bed = (bed.is_true()) ? 'f' : 'b';
            bn.needle = needle_num;
            bn.nudge = std::floor((needle_num - bn.needle) * 2);
            h.rhs = bn;
            mesh.hints.push_back(h);
        }
    }

    // Connection
    // Check src and target of connections
    /*
    std::map<uint32_t, std::vector<std::pair<sm::BedNeedle, sm::BedNeedle>>> c_xfer;
    for (auto const &mc : mesh.move_connections) {
        sm::BedNeedle src = mesh.move_instructions[mc.i_idx].src;
        sm::BedNeedle tgt = mesh.move_instructions[mc.i_idx].tgt;
        std::pair p = std::make_pair(src, tgt);

        if (c_xfer.count(mc.c_idx)) {
            c_xfer[mc.c_idx].push_back(p);
        } else {
            std::vector<std::pair<sm::BedNeedle, sm::BedNeedle>> v;
            v.push_back(p);
            c_xfer.insert({mc.c_idx, v});
        }
    }
    */

    return;
}

sm::Mesh sm::infer_constraints(sm::Mesh &mesh, sm::Code const &code_library, sm::Library const &face_library) {
    // Delete all constraints except user constraints when infer_constraints is called
    for(int i = 0;; ) {
        if (i == (int) mesh.hints.size()) break;
        sm::Mesh::Hint h = mesh.hints[i];
        if (h.src != sm::Mesh::Hint::User) {
            mesh.hints.erase(mesh.hints.begin() + i);
        } else {
            i++;
        }
    }

    generate_smt_phi(mesh, code_library);
	bool flag = true;
	while (flag) {
		flag = false;
		flag |= constraint_assign_variant_from_resource(mesh, code_library);
		flag |= constraint_extend_resource_from_resource(mesh, code_library);
		//flag |= constraint_face_instruction_order(mesh, code_library);
		flag |= constraint_assign_order_from_face_order(mesh, code_library, face_library); // update based on variants 
	}

	return mesh;
}
