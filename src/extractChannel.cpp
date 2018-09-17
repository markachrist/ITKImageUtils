/*
 * extractChannel.cpp
 *
 *  Created on: Jul 28, 2015
 *      Author: mchristopher
 */

#include "ImageUtils.h"
#include "RGBUtils.h"

/**
 * Print usage message and exit.
 */
void usage(int argc, char **argv){

	std::cout << std::endl << "Usage:" << std::endl << std::endl;
	std::cout << "%" << argv[0] << " <input file> <output file> <channel>" << std::endl << std::endl;
	std::cout << "Extract single channel (R = 0, G = 1, B = 2) from <input file> and save to <output file>.";

	exit(-1);
}

/**
 * Tool for extracting a single channel from an RGB image. See usage() for usage information.
 */
int main(int argc, char **argv){

	//Get cl args
	int c = atoi(argv[3]);
	std::string inPath = argv[1];
	std::string outPath = argv[2];

	//Extract channel
	RGBUtils rgb(inPath);
	ByteUtils channel(ByteUtils::floatToByteImage(rgb.getChannel(c)));

	//Write output
	channel.writeToPath(outPath.c_str());

	return 0;
}
