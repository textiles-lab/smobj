#include<iostream>
#include "sm.hpp"

int main(int argc, char *argv[]){
	if(argc < 5){
		std::cout << "Usage: ./test-code input.code input.sf input.smobj out-text" << std::endl;
		return 0;
	}
	sm::Code code;
	std::cout << "Loading code  from " << argv[1] << std::endl;
	std::cout << "Loading library  from " << argv[2] << std::endl;
	std::cout << "Loading mesh  from " << argv[3] << std::endl;
	code = sm::Code::load(argv[1]);
	

	std::cout << "Saving code to " << argv[4] << std::endl;
	code.save(argv[4]);
	
	sm::Library lib = sm::Library::load(argv[2]);
	sm::Mesh mesh = sm::Mesh::load(argv[3]);

	// Todo load library (.sf) and mesh (.smobj)
	// verify mesh is hinted
	//bool res = sm::verify_hinted_schedule(mesh, lib, code);
	//std::cout << "verifier result " << res << std::endl;

	sm::Mesh out = sm::order_faces(mesh, lib);
	std::string knitout = sm::knitout(out, code);

	std::cout << knitout << std::endl;

	return 0;
} 
