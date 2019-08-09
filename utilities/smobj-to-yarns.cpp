
#include "sm.hpp"

#include <string>
#include <iostream>

int main(int argc, char **argv) {

	if (argc != 4) {
		std::cerr << "Usage:\n\t./smobj-to-yarns <in.smobj> <input-library.sf> <out.yarns>" << std::endl;
		return 1;
	}

	std::string in_smobj = argv[1];
	std::string in_library = argv[2];
	std::string out_yarns = argv[3];
	std::cout << "Will elaborate faces of smobj '" << in_smobj << "' using yarn paths from library '" << in_library << "' and write resulting yarns to '" << out_yarns << "'" << std::endl;

	sm::Mesh mesh = sm::Mesh::load(in_smobj);

	sm::Library library = sm::Library::load(in_library);


	sm::Yarns yarns;

	sm::mesh_and_library_to_yarns(mesh, library, &yarns);

	yarns.save(out_yarns);

}

