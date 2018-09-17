/*
 * makeTiff.cpp
 *
 *  Created on: Oct 19, 2015
 *      Author: mchristopher
 */

#include <iostream>
#include "ImageUtils.h"

void usage(int argc, char **argv){

	std::cout << std::endl << "Usage:" << std::endl << std::endl;
	std::cout << "%" << argv[0] << " <output type> <output path>" << std::endl << std::endl;
	std::cout << "Converts image at <input path> of <input type> to <output type> and writes it to <output path>.";

	std::cout << std::endl;

	exit(-1);
}


/**
 *
 */
int main(int argc, char **argv){

	std::string inPath = argv[1];
	std::string outPath = argv[2];

	if(argc < 3){
		std::cout << "ERROR - missing required arguments!" << std::endl;
		usage(argc, argv);
	}

	FloatUtils vtk(inPath);

	ByteImageType::Pointer byte = FloatUtils::floatToByteImage(vtk.getImage());
	ByteUtils out(byte);

	out.writeToPath(outPath);

	return 0;
}

