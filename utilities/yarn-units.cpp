#include "sm.hpp"

#include <glm/gtx/hash.hpp>
#include <glm/gtx/norm.hpp>

#include <iostream>
#include <unordered_map>
#include <map>

int main(int argc, char **argv) {
	bool usage = false;
	std::string in_yarns = "";
	std::string out_yarns = "";
	std::map< std::string, float > set;
	bool sum = false;
	bool print = false;

	for (int argi = 1; argi < argc; ++argi) {
		std::string arg = argv[argi];
		if (arg == "--help") {
			usage = true;
		} else if (arg == "--set") {
			if (argi + 2 < argc) {
				argi += 1;
				std::string name = argv[argi];
				argi += 1;
				float length = std::stoi(argv[argi]);
				if (set.count(name)) {
					std::cerr << "ERROR: attempting to '--set " << name << "' twice." << std::endl;
					usage = true;
				}
				set[name] = length;
			} else {
				std::cerr << "ERROR: '--set' should be followed by a unit name and a length" << std::endl;
				usage = true;
			}
		} else if (arg == "--sum") {
			sum = true;
		} else if (arg == "--print") {
			print = true;
		} else if (arg == "--write") {
			if (argi + 1 < argc) {
				if (out_yarns != "") {
					std::cerr << "ERROR: multiple output files given." << std::endl;
					usage = true;
				}
				argi += 1;
				out_yarns = argv[argi];
			} else {
				std::cerr << "Expecting output filename after '--write'" << std::endl;
				usage = true;
			}
		} else if (arg == "--") {
			if (argi + 2 != argc) {
				std::cerr << "ERROR: expecting '--' to be followed by exactly the input filename and nothing more." << std::endl;
				usage = true;
			} else {
				if (in_yarns != "") {
					std::cerr << "ERROR: multiple input files ('" << in_yarns << "', '" << argv[argi+1] << "') given." << std::endl;
					usage = true;
				}
				in_yarns = argv[argi+1];
			}
			break; //no argument parsing past '-- input.yarns'
		} else {
			if (in_yarns != "") {
				std::cerr << "ERROR: multiple input files ('" << in_yarns << "', '" << arg << "') given." << std::endl;
				usage = true;
			}
			in_yarns = arg;
		}
	}
	if (in_yarns == "") {
		std::cerr << "ERROR: no input file specified." << std::endl;
		usage = true;
	}
	if (!(print || out_yarns != "" || usage)) {
		std::cerr << "ERROR: utility asked to produce no output." << std::endl;
		usage = true;
	}
	if (usage) {
		std::cerr << "Usage:\n\t./yarn-units [--set unit length] [--set unit length] [--sum] [--print] [--write out.yarns] [--] <in.yarns>\n"
		<< " --set unit length     set the length factor for 'unit' to length\n"
		<< " --sum                 sum all lengths so file only contains checkpoints in terms of the '1' unit\n"
		<< " --print               print out all units, their default values, and the total of that unit in the file\n"
		<< " --write out.yarns     write output file\n"
		<< " --help                print this text\n"
		<< std::endl;
		return 1;
	}

	std::cerr << "Loading '" << in_yarns << "' ..."; std::cerr.flush();
	sm::Yarns yarns = sm::Yarns::load(in_yarns);
	std::cerr << " done." << std::endl;

	//----------------------------------------------
	for (auto const &[name, length] : set) {
		bool found = false;
		for (auto &unit : yarns.units) {
			if (unit.name == name) {
				std::cerr << "set '" << name << "' to " << length << "\n";
				if (unit.name == "1" && length != 1.0) {
					std::cerr << "WARNING: changing the length of the '1' unit to anything other than 1.0 is not allowed or expected by the specification (but doing it anyway because you asked)." << std::endl;
				}
				unit.length = length;
				found = true;
			}
		}
		if (!found) {
				std::cerr << "WARNING: cannot set '" << name << "' because it doesn't appear in the file." << std::endl;
		}
	}

	//----------------------------------------------
	if (sum) {
		std::cerr << "Summing all checkpoints using current unit values..."; std::cerr.flush();
		for (auto &yarn : yarns.yarns) {
			std::map< uint32_t, float > new_lengths;
			for (auto const &cp : yarn.checkpoints) {
				new_lengths.emplace(cp.point, 0.0f).first->second += cp.length * yarns.units.at(cp.unit).length;
			}
			yarn.checkpoints.clear();
			yarn.checkpoints.reserve(new_lengths.size());
			for (auto const &[point, length] : new_lengths) {
				yarn.checkpoints.emplace_back();
				yarn.checkpoints.back().point = point;
				yarn.checkpoints.back().length = length;
				yarn.checkpoints.back().unit = 0;
			}
		}
		yarns.units.clear();
		yarns.units.emplace_back();
		yarns.units.back().name = "1";
		yarns.units.back().length = 1.0f;
		std::cerr << " done." << std::endl;
	}

	//----------------------------------------------
	if (print) {
		std::vector< std::vector< float > > amounts(yarns.units.size());
		for (auto const &yarn : yarns.yarns) {
			for (auto const &cp : yarn.checkpoints) {
				assert(cp.unit < amounts.size());
				amounts[cp.unit].emplace_back(cp.length);
			}
		}
		std::vector< std::vector< std::string > > table;
		std::vector< float > totals;
		for (uint32_t i = 0; i < yarns.units.size(); ++i) {

			//could be a lot of small values, so do a careful sum by adding the two smallest lengths first:
			std::make_heap(amounts[i].begin(), amounts[i].end(), std::greater< float >());
			while (amounts[i].size() > 1) {
				std::pop_heap(amounts[i].begin(), amounts[i].end(), std::greater< float >());
				std::pop_heap(amounts[i].begin(), amounts[i].end()-1, std::greater< float >());
				*(amounts[i].end()-2) += *(amounts[i].end()-1);
				amounts[i].pop_back();
			}
			assert(amounts.size() <= 1);

			float total = (amounts[i].empty() ? 0.0f : amounts[i][0]);

			totals.emplace_back(total * yarns.units[i].length);
			table.emplace_back();
			table.back().emplace_back(yarns.units[i].name);
			table.back().emplace_back(": ");
			table.back().emplace_back(std::to_string(total));
			table.back().emplace_back(" * ");
			table.back().emplace_back(std::to_string(yarns.units[i].length));
			table.back().emplace_back(" = ");
			table.back().emplace_back(std::to_string(total * yarns.units[i].length));
			//std::cout << yarns.units[i].name << ": " <<  total << " (* " << yarns.units[i].length << ")" << "\n";
		}
		{ //overall total:
			std::sort(totals.begin(), totals.end());
			float total = 0.0f;
			for (float a : totals) {
				total += a;
			}
			table.emplace_back();
			table.back().emplace_back("total");
			table.back().emplace_back("");
			table.back().emplace_back("");
			table.back().emplace_back("");
			table.back().emplace_back("");
			table.back().emplace_back("");
			table.back().emplace_back(std::to_string(total));
		}

		for (uint32_t col = 0; col < table[0].size(); ++col) {
			size_t max = 0;
			for (uint32_t row = 0; row < table.size(); ++row) {
				max = std::max(max, table[row].at(col).size());
			}
			for (uint32_t row = 0; row < table.size(); ++row) {
				auto &val = table[row][col];
				if (col == 0) {
					while (val.size() < max) val = ' ' + val;
				} else {
					while (val.size() < max) val = ' ' + val;
				}
			}
		}
		for (auto const &row : table) {
			for (auto const &val : row) {
				std::cout << val;
			}
			std::cout << '\n';
		}
		std::cout.flush();
	}

	if (out_yarns != "") {
		std::cerr << "INFO: writing '" << out_yarns << "'..."; std::cerr.flush();
		yarns.save(out_yarns);
		std::cerr << " done." << std::endl;
	}

	return 0;
}
