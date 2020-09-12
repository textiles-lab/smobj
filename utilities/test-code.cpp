#include<iostream>
#include <fstream>
#include "sm.hpp"
#include "hinters.hpp"

int main(int argc, char *argv[]){
	if(argc < 4){
		std::cout << "Usage: ./test-code input.code input.sf input.smobj" << std::endl;
		return 0;
	}
	sm::Code code;
	std::cout << "Loading code  from " << argv[1] << std::endl;
	std::cout << "Loading library  from " << argv[2] << std::endl;
	std::cout << "Loading mesh  from " << argv[3] << std::endl;
	code = sm::Code::load(argv[1]);
	code.save("out.code");

	
	sm::Library lib = sm::Library::load(argv[2]);
	sm::Mesh mesh = sm::Mesh::load(argv[3]);

	// Todo load library (.sf) and mesh (.smobj)
	// verify mesh is hinted
	//bool res = sm::verify_hinted_schedule(mesh, lib, code);
	//std::cout << "verifier result " << res << std::endl;

	sm::Mesh out = sm::order_faces(mesh, lib);

	out = sm::hint_shortrow_only_patch(out, code); // <-- add hints 

	
	out.save("hinted.smobj");
	std::cout << "Saved hinted smobj as ./hinted.smobj" << std::endl; 

	std::string knitout = sm::knitout(out, code);
	std::ofstream kw("out.knitout");
	kw << knitout;
	kw.close();
	std::cout << knitout << std::endl;

	return 0;
} 
