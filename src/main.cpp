/*
 * main.cpp
 *
 *  Created on: Dec 1, 2012
 *      Author: mchristopher
 */

#include <iostream>
#include "opencv2/ml/ml.hpp"
#include "ImageUtils.h"
#include "RGBUtils.h"
#include "ImageDataset.h"
#include "CVMatUtils.h"
#include "Performance.h"
#include "ImageAligner.h"

void usage(){

	std::cout << std::endl << "Usage:" << std::endl << std::endl;
	std::cout << "%imageUtilsMain <list file> <output dir>" << std::endl << std::endl;
	std::cout << "Converts set of float images to byte images and saves as TIFFs.";

	std::cout << std::endl;

	exit(-1);
}

int main(int argc, char **argv){

	std::string inPath = argv[1];
	std::string outPath = argv[2];

//	FloatUtils in(inPath);
//	ByteImageType::Pointer bite = FloatUtils::floatToByteImage(in.getImage());
//	ByteUtils out(bite);
//
//	out.writeToPath(outPath.c_str());

	RGBUtils in(inPath);

	std::cout << "Size = " << in.getWidth() << " x " << in.getHeight() << std::endl;

	in.writeToPath(outPath);

	return 0;
}
