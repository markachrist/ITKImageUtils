/*
 * vectorToImage.cpp
 *
 *  Created on: Oct 26, 2015
 *      Author: mchristopher
 */

#include <iostream>
#include "ImageUtils.h"
#include "CVMatUtils.h"

void usage(int argc, char **argv){

	std::cout << std::endl << "Usage:" << std::endl << std::endl;
	std::cout << "%" << argv[0] << " <vector txt file> <width> <height> <output path>" << std::endl << std::endl;
	std::cout << "Converts list of text numbers to image.";

	std::cout << std::endl;

	exit(-1);
}


/**
 *
 */
int main(int argc, char **argv){

	std::string inPath = argv[1];
	int w = atoi(argv[2]);
	int h = atoi(argv[3]);
	std::string outPath = argv[4];

	if(argc < 3){
		std::cout << "ERROR - missing required arguments!" << std::endl;
		usage(argc, argv);
	}

	cv::Mat vector;
	FloatUtils image(w, h, 0);
	FloatImageType::Pointer ptr = image.getImage();

	CVMatUtils::readCSV(vector, inPath.c_str());

	std::cout << "Vector input: " << vector.rows << " x " << vector.cols << std::endl;
	std::cout << "Image output: " << w << " x " << h << std::endl;

	vector.convertTo(vector, CV_32FC1);

	CVMatUtils::vectorToImage(vector, w, h, ptr);

	ByteImageType::Pointer byte = FloatUtils::floatToByteImage(ptr);
	ByteUtils out(byte);

//	FloatUtils out(ptr);

	out.writeToPath(outPath);

	return 0;
}
