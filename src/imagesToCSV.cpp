/*
 * imagesToCSV.cpp
 *
 *  Created on: Nov 6, 2015
 *      Author: mchristopher
 */

#include <iostream>
#include "ImageUtils.h"
#include "ImageDataset.h"
#include "CVMatUtils.h"

void usage(int argc, char **argv){

	std::cout << std::endl << "Usage:" << std::endl << std::endl;
	std::cout << "%" << argv[0] << " <image list> <output csv file>" << std::endl << std::endl;
	std::cout << "Flattens images and writes each to single row in CSV file.";

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


	std::vector<std::string> list;
	cv::Mat mat;

	ImageDataset::readLines(inPath, list);
	ImageDataset data(list);
	CVMatUtils::imagesToMat(data, mat);

	CVMatUtils::writeCSV(mat, outPath.c_str());

	return 0;
}

