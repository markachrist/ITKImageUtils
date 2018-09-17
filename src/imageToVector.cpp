/*
 * imageToVector.cpp
 *
 *  Created on: Oct 26, 2015
 *      Author: mchristopher
 */

#include <iostream>
#include "ImageUtils.h"
#include "CVMatUtils.h"

void usage(int argc, char **argv){

	std::cout << std::endl << "Usage:" << std::endl << std::endl;
	std::cout << "%" << argv[0] << " <input image> <output txt file>" << std::endl << std::endl;
	std::cout << "Flattens image and writes pixel vlaues to text.";

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

	ByteUtils reader(inPath);
	FloatImageType::Pointer image = FloatUtils::byteToFloatImage(reader.getImage());
	cv::Mat vector;

	CVMatUtils::imageToVector(image, vector);
	vector = vector.t();
	CVMatUtils::writeCSV(vector, outPath.c_str());

	return 0;
}
