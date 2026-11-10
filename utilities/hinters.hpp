// helpers to generate hinted smobjs
#include "sm.hpp"
// add to the sm namespace
namespace sm{
	sm::Mesh constraint_assign_frontface_variant(sm::Mesh const &mesh, sm::Code const &code_library);
	sm::Mesh hint_shortrow_only_patch(sm::Mesh const &in, sm::Code const &code_library);
	sm::Mesh hint_shortrow_only_tubes(sm::Mesh const &in, sm::Code const &code_library);
	sm::Mesh infer_constraints(sm::Mesh &in, sm::Code const &code_library, sm::Library const &face_library);
	bool constraint_assign_order_from_face_order(sm::Mesh &mesh, sm::Code const &code_library, sm::Library const &face_library);
	bool constraint_assign_variant_from_resource(sm::Mesh  &mesh, sm::Code const &code_library);
    bool constraint_extend_resource_from_resource(sm::Mesh  &mesh, sm::Code const &code_library);
    bool constraint_face_instruction_order(sm::Mesh &mesh, sm::Code const &code_library);

};
