#include<iostream>
#include <fstream>
#include "sm.hpp"
#include "hinters.hpp"

int main(int argc, char *argv[]){
	if(argc < 4){
		std::cout << "Usage: ./test-code input.code input.sf input.smobj" << std::endl;
		return 0;
	}
	std::cout << "Loading code  from " << argv[1] << std::endl;
	std::cout << "Loading library  from " << argv[2] << std::endl;
	std::cout << "Loading mesh  from " << argv[3] << std::endl;
	sm::Code code;
	code = sm::Code::load(argv[1]);
	sm::Library lib = sm::Library::load(argv[2]);
	sm::Mesh mesh = sm::Mesh::load(argv[3]);

	std::cout << "------------- loaded files --------------- " << std::endl;
	if(!sm::compute_code_graph(code, lib)){
		return 0;
	}
	/*
	for(auto const &f : code.faces){
		std::cout << "f: " << f.key_library() << " has e->e " << f.edge_to_edge_connections.size() << " e->i " << f.edge_to_instruction_connections.size() << " i->i " << f.instruction_to_instruction_connections.size() << " i-> e " << f.instruction_to_edge_connections.size() << std::endl;
	}*/
	std::cout << "------------- code graph computed  --------------- " << std::endl;
	//return 0;

	/*
	std::cout << "mesh xfer instructions:" << mesh.move_instructions.size() << std::endl;
	for(auto op : mesh.move_instructions){
		std::cout << "\t" << op.to_string() << std::endl;
	}
	std::cout << "move connections: " << mesh.move_connections.size()  << std::endl;
	for(auto ci : mesh.move_connections){
		std::cout << "\t Instruction " << ci.i_idx << " is associated with connection " << ci.c_idx << std::endl;
	}*/
	
	std::vector<sm::Mesh::Hint> offenders;
	//if(sm::verify(mesh, lib, code, &offenders) && sm::compute_total_order(mesh, code)){
	if( sm::compute_total_order(mesh, code)  && sm::verify(mesh, lib, code, &offenders, true)){
		std::string knitout = sm::knitout(mesh, lib, code);
		std::ofstream kw("out.knitout");
		kw << knitout;
		kw.close();
		
		std::cout << knitout << std::endl;
		std::cout << "Successfully generated knitout instructions. Total order: " << std::endl;
		for(auto const &fi : mesh.total_order){
			std::cout << fi.first << "/" << fi.second << std::endl; 
		}
	}
	return 0;
} 
