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
	if(!sm::compute_code_graph(code)){
		return 0;
	}
	std::cout << "------------- code graph computed  --------------- " << std::endl;
	
	
	std::vector<sm::Mesh::Hint> offenders;
	//if(sm::verify(mesh, lib, code, &offenders) && sm::compute_total_order(mesh, code)){
	if( sm::verify(mesh, lib, code, &offenders, true)){
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
	{
		if(!sm::compute_total_order(mesh, code, lib)){
			std::cerr << "\tTotal order doesn't exist. " << std::endl;
		}
		/*
		if(!sm::verify(mesh, lib, code, &offenders, true)){
			std::cerr << "\tCould not verify mesh. " << std::endl;
			std::cerr << "\tOffenders: " << offenders.size() << std::endl;
		}*/
	}
	return 0;
} 
