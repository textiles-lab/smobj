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

	std::cout << "mesh xfer instructions:" << mesh.move_instructions.size() << std::endl;
	for(auto op : mesh.move_instructions){
		std::cout << "\t" << op.to_string() << std::endl;
	}
	std::cout << "move connections: " << mesh.move_connections.size()  << std::endl;
	for(auto ci : mesh.move_connections){
		std::cout << "\t Instruction " << ci.i_idx << " is associated with connection " << ci.c_idx << std::endl;
	}
	// Todo load library (.sf) and mesh (.smobj)
	// verify mesh is hinted
	//bool res = sm::verify_hinted_schedule(mesh, lib, code);
	//std::cout << "verifier result " << res << std::endl;

	sm::Instr xop; xop.op = sm::Instr::Xfer; xop.src.bed = 'f'; xop.src.needle = 1000; xop.tgt.bed = 'b'; xop.tgt.needle = 1000;
	sm::Mesh out = mesh;
	out.move_instructions.emplace_back(xop);

	out.move_connections.emplace_back();
	out.move_connections.back().i_idx = out.move_instructions.size()-1;
	out.move_connections.back().c_idx = out.connections.size()-1;
	out.move_connections.back().connection = out.connections.back();

	out.save("out.smobj");
	std::cout << "Saved hinted smobj as ./hinted.smobj" << std::endl; 

	std::vector<sm::Mesh::Hint> offenders;
	if(sm::verify(out, lib, code, &offenders)){
		std::string knitout = sm::knitout(out, lib, code);
		std::ofstream kw("out.knitout");
		kw << knitout;
		kw.close();
		std::cout << knitout << std::endl;
	}
	return 0;
} 
