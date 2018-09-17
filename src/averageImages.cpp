/*
 * averageImages.cpp
 *
 *  Created on: Sep 24, 2015
 *      Author: mchristopher
 */

#include <iostream>
#include "ImageUtils.h"

void usage(){

	std::cout << std::endl << "Usage:" << std::endl << std::endl;
	std::cout << "%averageImages <output path> <image 1> [<image 2> <image 3> ...]" << std::endl << std::endl;
	std::cout << "Averages the given set of float-valued images and stores output in <output path>";

	std::cout << std::endl;

	exit(-1);
}

int main(int argc, char **argv){

	int n = 0;

	if(argc < 3){
		std::cout << "ERROR - need atleast one ipnut image!" << std::endl;
		usage();
	}

	std::string outPath = argv[1];

	FloatUtils out;

	for(int i = 2; i < argc; ++i){
		std::cout << "Adding " << argv[i] << " to average." << std::endl;
		FloatUtils current(argv[i]);

		if(i == 2){
			out.setImage(current.getImage());
		}
		else{
			out.setImage(out.add(current.getImage()));
		}
		++n;
	}

	out.setImage(out.mult(1.0 / n));
	out.writeToPath(outPath);

	return 0;
}
