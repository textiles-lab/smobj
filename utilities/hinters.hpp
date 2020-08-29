// helpers to generate hinted smobjs 
#include "sm.hpp"
// add to the sm namespace
namespace sm{
	sm::Mesh hint_shortrow_only_patch(sm::Mesh const &in, sm::Code const &code_library);
	sm::Mesh hint_shortrow_only_tubes(sm::Mesh const &in, sm::Code const &code_library);
	sm::Mesh hint_extend(sm::Mesh const &in, sm::Code const &code_library);
};
