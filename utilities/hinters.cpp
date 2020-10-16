#include "sm.hpp"
#include "hinters.hpp"
#include <iostream>

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
sm::Mesh sm::hint_shortrow_only_tubes(sm::Mesh const &mesh, sm::Code const &code_library){
    assert(false && "TODO IMPLEMENT");
    sm::Mesh m = mesh;
    return m;
}

sm::BedNeedle find_edge_resource(sm::Mesh const &mesh, uint32_t fidx, uint32_t eidx) {
    // go through resource hints
    for(auto &h : mesh.hints) {
        if (h.type == sm::Mesh::Hint::Resource && h.lhs.face == fidx && h.lhs.edge == eidx)
            return std::get<sm::BedNeedle>(h.rhs);
    }

    sm::BedNeedle bn;
    return bn;
};

// TODO for face:
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

    std::map<std::string, uint32_t > name_to_code_idx;
    // Doesn't work if edges in the code library face are different among variants!
    for(auto const &c : code_library.faces){
        // What to do if there is a same library name with different variants?
        name_to_code_idx[c.key_library()] = &c - &code_library.faces[0];
    }

    // Accumulate all changes here and assign afterwards
    std::map<sm::Mesh::FaceEdge, sm::BedNeedle> constraints;

    // (a)
    // Iterate through edges and obtain BedNeedle
    //
    //
    for(auto &f : mesh.faces){
        std::string lib_name = mesh.library[f.type];
        uint32_t face_id = &f - &mesh.faces[0];
        uint32_t code_idx = name_to_code_idx[lib_name];
        //NOTE: (Maybe) If code template bedneedle assignment can be different across variants, also look for variant if available.
        // Obtain the first edge that resource hint is not 'x' and store it in edge_constraint
        bool found_resource = false;
        auto l = code_library.faces[code_idx]; // TODO: has the right variant given the existing resource hints
        std::pair<sm::BedNeedle, uint32_t> edge_constraint;
        int offset = 0;
        for(auto &e : l.edges) {
            uint32_t eidx = &e - &l.edges[0];
            sm::BedNeedle bn = find_edge_resource(mesh, face_id, eidx);
            if (bn.bed != 'x') {
                // sanity check that e.bn.bed == bn.bed
                found_resource = true;
                edge_constraint.first = bn;
                edge_constraint.second = eidx;
                // offset between e.bn and bn should be the same for all edges in this face
                offset = e.bn.location() - bn.location();
                break;
            }
        }

        // Propagate constraints to other edges in this face.
        // Raise error if there is a resource conflict and return false
        if (found_resource) {
            for(auto &e : l.edges) {
                uint32_t eidx = &e - &l.edges[0];
                // skip if this edge is edge_constraint edge
                if (eidx == edge_constraint.second) continue;

                sm::BedNeedle bn = find_edge_resource(mesh, face_id, eidx);
                float eps = 1e-3f;
                // If there is already a resource assigned to this edge, make sure that it's not conflicting
                if (bn.bed != 'x') {
                    // needle 1, nudge -1 and needle 0, nudge 1 --> compatible since location is 0.5, location returns a float
                    // needle check: std::abs(bn.location()- other.location()) < eps)
                    if (bn.bed != edge_constraint.first.bed ||
                            std::abs(bn.location()- edge_constraint.first.location()) < eps) {
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
                    float location = e.bn.location() + offset;
                    new_bn.needle = location;
                    new_bn.nudge = (location - new_bn.needle) > eps ? 1 : -1;

                    constraints[fe] = new_bn;
                }
            }
        }
    }

    // (b)
    for(auto const &c : mesh.connections){
        if(constraints.count(c.a) && !constraints.count(c.b)){
            constraints[c.b] = constraints[c.a];
        } else if (constraints.count(c.b) && !constraints.count(c.a)) {
            constraints[c.a] = constraints[c.b];
        }
    }

    // return this
    bool is_changed = false;
    for(auto h : constraints){
        is_changed = true;
        sm::BedNeedle bn;
        mesh.hints.emplace_back();
        mesh.hints.back().lhs = h.first;
        bn.bed = h.second.bed;
        bn.needle = h.second.needle;
        bn.nudge = h.second.nudge;
        mesh.hints.back().rhs = bn;
        mesh.hints.back().type = sm::Mesh::Hint::Resource;
        mesh.hints.back().src = sm::Mesh::Hint::Inferred; // TODO: This might need fixing
    }

    return is_changed;
}

// TODO: Error reporting
// We want to tell users which constraints doesn't make sense to users
// Get a list of incompatible constraints, and somehow report it to users
// Maybe implement a better API which can be called from EditMesh.cpp??

// TODO: Also Error reporting
// If constraint 0 and 1 are incompatible, then 1 is incompatibe with 0


// Helper function that checks if code-template-face is compatible with edge-resources
//{

//}


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
        // What to do if there is a same library name with different variants?
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
        for (auto face_code : name_to_code_idx[lib_name]) {
            // Get code face
            auto l = code_library.faces[face_code];
            bool one_edge_has_resource = false;
            std::string f_variant_constraint = "";

            // For each edge in this code face
            for (int i = 0; i < (int)l.edges.size(); i++) {
                bool e_resource = false;
                sm::Code::Face::Edge e = l.edges[i];
                // For all resource constraints that this face has
                for (auto f_constraint : face_constraints) {
                    // If a resource constraint is for this edge
                    if (i == (int)f_constraint.lhs.edge) {
                        // Get constraint content
                        sm::BedNeedle f_bn = std::get<sm::BedNeedle>(f_constraint.rhs);
                        // Code face and resource constraint matches!
                        if (e.bn.bed == f_bn.bed) {
                            // If some edges' resource constraint were Heuristic,
                            // Variant constraint that we're trying to construct should
                            // also be Heuristic.
                            if (f_constraint.src == sm::Mesh::Hint::Heuristic)
                                f_src = sm::Mesh::Hint::Heuristic;
                            // We have found one resource constraint for this edge
                            e_resource = true;
                            // Assign variant
                            f_variant_constraint = l.variant;
                        }
                    }
                }

                if (e_resource)
                    one_edge_has_resource = true;
                else
                    f_src = sm::Mesh::Hint::Heuristic;
            }


            if (one_edge_has_resource) {
                sm::Mesh::Hint face_hint;

                sm::Mesh::FaceEdge fe;
                fe.face = face_id;
                face_hint.lhs = fe;
                face_hint.type = sm::Mesh::Hint::Variant;

                face_hint.rhs = f_variant_constraint;
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
        is_changed = true;
        mesh.hints.push_back(h);
    }

    return is_changed;
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

sm::Mesh sm::infer_constraints(sm::Mesh &mesh, sm::Code const &code_library) {
    bool flag = true;
    while (flag) {
        flag = false;
        flag |= constraint_assign_variant_from_resource(mesh, code_library);
        flag |= constraint_extend_resource_from_resource(mesh, code_library);
        //flag |= constraint_face_instruction_order(mesh, code_library);
    }

    return mesh;
}
