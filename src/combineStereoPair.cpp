/*
 * combineStereoPair.cpp
 *
 *  Created on: Oct 12, 2015
 *      Author: mchristopher
 */

#include <iostream>

#include "ImageUtils.h"
#include "RGBUtils.h"

#define PAD_PIXELS 5

void usage(int argc, char **argv){

	std::cout << "%" << argv[0] << " <left image> <right image> <out image>" << std::endl << std::endl;
	std::cout << "Combines two images into a single image and stores output at the given destination.";

	std::cout << std::endl << std::endl;

	exit(-1);
}

int main(int argc, char **argv){

	int w, h;
	RGBImageType::IndexType idx;
	std::string leftPath = argv[1];
	std::string rightPath = argv[2];
	std::string outPath = argv[3];

	if(argc <= 3){
		std::cout << "Error: Missing arguments!" << std::endl;
		usage(argc, argv);
	}

	RGBUtils right(rightPath);
	RGBUtils left(leftPath);

	w = right.getWidth() + PAD_PIXELS + left.getWidth();
	h = right.getHeight();
	h = (h >= left.getHeight()) ? h : left.getHeight();

	RGBUtils out(w, h);

	idx[0] = 0;
	idx[1] = 0;
	out.setImage(out.copyImageAt(left.getImage(), idx));

	idx[0] = PAD_PIXELS + left.getWidth();
	idx[1] = 0;
	out.setImage(out.copyImageAt(right.getImage(), idx));

	out.writeToPath(outPath.c_str());

	return 0;
}
