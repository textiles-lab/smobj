#include "sm.hpp"

#include <iostream>
#include <fstream>
#include <cstring>

int main(int argc, char **argv) {
	if (argc != 3) {
		std::cerr << "Usage:\n\t./yarns-to-bccx <in.yarns> <out.bccx>" << std::endl;
		return 1;
	}
	std::string in_yarns = argv[1];
	std::string out_bccx = argv[2];

	std::cout << "Will convert yarns in '" << in_yarns << "' to a bccx-formatted file '" << out_bccx << "'." << std::endl;

	sm::Yarns yarns = sm::Yarns::load(in_yarns);

	{
		uint32_t total_points = 0;
		uint32_t total_checkpoints = 0;
		for (auto const &yarn : yarns.yarns) {
			total_points += yarn.points.size();
			total_checkpoints += yarn.checkpoints.size();
		}
		std::cout << "Input has " << yarns.yarns.size() << " yarns with a total of " << total_points << " points and " << total_checkpoints << " checkpoints." << std::endl;
		if (total_checkpoints != 0) {
			std::cerr << "WARNING: checkpoint conversion code not yet written." << std::endl;
		}
	}


	std::ofstream out(out_bccx, std::ios::binary);

	//BCC(x) format is based on 'sim.h' and 'bcc.cpp' from YarnBallSim
	// as well as http://www.cemyuksel.com/cyCodeBase/soln/using_bcc_files.html
	struct BCCHeader {
		char magic[3] = {'B', 'C', 'C'};
		uint8_t bytes = 0x44; //4-byte integers, 4-byte floating point
		char format[2] = {'C', '0'}; //uniform catmull-rom splines
		uint8_t dimensions = 3; //3d points
		uint8_t up_direction = 2; //1 is y-up, 2 is z-up
		uint64_t num_yarns;
		uint64_t num_nodes;
		char info[40];
	} header;
	static_assert(sizeof(BCCHeader) == 64, "BCC header should be 64 bytes.");

	header.num_yarns = yarns.yarns.size();
	header.num_nodes = 0;
	for (auto const &yarn : yarns.yarns) {
		header.num_nodes += yarn.points.size();
	}
	{ //copy source filename into info:
		memset(header.info, '\0', sizeof(header.info));
		std::string info = "from '" + in_yarns + "'";
		if (info.size() > 40) info = info.substr(0, 40);
		memcpy(header.info, info.data(), info.size());
	}
	out.write(reinterpret_cast< const char * >(&header), sizeof(header));

	std::cerr << "WARNING: copying control points directly from yarns to bccx file -- this is technically incorrect in that yarns uses quadratic bezier curves between midpoints of yarns polylines, while bccx files use uniform cubic catmull-rom splines through yarn points." << std::endl;
	for (auto &yarn : yarns.yarns) {
		int32_t count = yarn.points.size();
		//NOTE: if count is negative, yarn is a closed loop; yarns format doesn't store closed loops
		//though, if yarns ever were used to store closed loops, it would be by having the first and last point in the yarn occupy exactly the same position:
		if (!yarn.points.empty() && yarn.points[0] == yarn.points.back()) {
			count = -count;
		}
		out.write(reinterpret_cast< const char * >(&count), sizeof(count));
		static_assert(sizeof(yarn.points[0]) == 12, "yarn points are 3-vectors of floating point values");
		out.write(reinterpret_cast< const char * >(yarn.points.data()), yarn.points.size() * sizeof(yarn.points[0]));
	}

	{ //bccx stuff:
		//yarn radius:
		float radius = 0.1f;
		if (!yarns.yarns.empty()) {
			radius = yarns.yarns[0].radius;
			for (auto const &yarn : yarns.yarns) {
				if (yarn.radius != radius) {
					std::cerr << "WARNING: bccx supports only a single radius, but yarns[0] has radius " << radius << " and yarns[" << (&yarn - &yarns.yarns[0]) << "] has radius " << yarn.radius << ". (Using yarns[0] radius.)" << std::endl;
				}
			}
		}
		out.write(reinterpret_cast< const char * >(&radius), sizeof(radius));
		//yarn point velocities:
		for (auto &yarn : yarns.yarns) {
			for ([[maybe_unused]] auto const &point : yarn.points) {
				glm::vec3 velocity = glm::vec3(0.0f, 0.0f, 0.0f);
				out.write(reinterpret_cast< const char * >(&velocity), sizeof(velocity));
			}
		}
		//yarn point rest lengths:
		// (YarnBallSim's bcc.cpp initializes these to distance-to-next-control-point, so that's what I will use for now; could use checkpoints for this in the future.)
		for (auto &yarn : yarns.yarns) {
			if (yarn.checkpoints.size()) {
				std::cerr << "WARNING: ignoring checkpoints in yarn[" << (&yarn - &yarns.yarns[0]) << "] -- reading lengths from checkpoints is not yet supported by this conversion code." << std::endl;
			}
			for (auto const &point : yarn.points) {
				//negative value is a flag for "not computed":
				float l = glm::length(point - yarn.points[(&point - &yarn.points[0] + 1) % yarn.points.size()]);
				out.write(reinterpret_cast< const char * >(&l), sizeof(l));
			}
		}
	}

	return 0;
}
