/*
 * convertImage.cpp
 *
 *  Created on: Sep 24, 2015
 *      Author: mchristopher
 */

//Define supported image types
#define BYTE_TYPE 0
#define INT_TYPE 1
#define FLOAT_TYPE 2
#define RGB_TYPE 3

#include <iostream>
#include "ImageUtils.h"
#include "RGBUtils.h"

void usage(){

	std::cout << std::endl << "Usage:" << std::endl << std::endl;
	std::cout << "%convertImage <input type>  <input path> <output type> <output path>" << std::endl << std::endl;
	std::cout << "Converts image at <input path> of <input type> to <output type> and writes it to <output path>.";

	std::cout << std::endl;

	std::cout << "Types:" << std::endl;
	std::cout << "Byte: " << BYTE_TYPE << std::endl;
	std::cout << "Integer: " << INT_TYPE << std::endl;
	std::cout << "Float: " << FLOAT_TYPE << std::endl;
	std::cout << "RGB: " << RGB_TYPE << std::endl;

	exit(-1);
}

/**
 *
 */
int main(int argc, char **argv){

	//Declare/define default parameters
	int inType;
	int outType;

	std::string inPath, outPath;

	if(argc < 5){
		std::cout << "ERROR - missing required arguments!" << std::endl;
		usage();
	}

	//Get types
	inType = atoi(argv[1]);
	inPath = argv[2];
	outType = atoi(argv[3]);
	outPath = argv[4];

	//Perform conversion
	if(inType == BYTE_TYPE){

		ByteUtils in(inPath);

		if(outType == BYTE_TYPE){
			ByteUtils out(in.getImage());
			out.writeToPath(outPath);
		}
		else if(outType == INT_TYPE){
			IntUtils out(ByteUtils::byteToIntImage(in.getImage()));
			out.writeToPath(outPath);
		}
		else if(outType == FLOAT_TYPE){
			FloatUtils out(ByteUtils::byteToFloatImage(in.getImage()));
			out.writeToPath(outPath);
		}
		else{
			RGBUtils out(in.getWidth(), in.getHeight());

			out.setChannel(0, ByteUtils::byteToFloatImage(in.getImage()));
			out.setChannel(1, ByteUtils::byteToFloatImage(in.getImage()));
			out.setChannel(2, ByteUtils::byteToFloatImage(in.getImage()));

			out.writeToPath(outPath);
		}

	}
	else if(inType == INT_TYPE){

		IntUtils in(inPath);

		if(outType == BYTE_TYPE){
			ByteUtils out(IntUtils::intToByteImage(in.getImage()));
			out.writeToPath(outPath);
		}
		else if(outType == INT_TYPE){
			IntUtils out(in.getImage());
			out.writeToPath(outPath);
		}
		else if(outType == FLOAT_TYPE){
			FloatUtils out(IntUtils::intToFloatImage(in.getImage()));
			out.writeToPath(outPath);
		}
		else{
			RGBUtils out(in.getWidth(), in.getHeight());

			out.setChannel(0, IntUtils::intToFloatImage(in.getImage()));
			out.setChannel(1, IntUtils::intToFloatImage(in.getImage()));
			out.setChannel(2, IntUtils::intToFloatImage(in.getImage()));

			out.writeToPath(outPath);
		}

	}
	else if(inType == FLOAT_TYPE){

		FloatUtils in(inPath);

		if(outType == BYTE_TYPE){
			ByteUtils out(FloatUtils::floatToByteImage(in.getImage()));
			out.writeToPath(outPath);
		}
		else if(outType == INT_TYPE){
			IntUtils out(FloatUtils::floatToIntImage(in.getImage()));
			out.writeToPath(outPath);
		}
		else if(outType == FLOAT_TYPE){
			FloatUtils out(in.getImage());
			out.writeToPath(outPath);
		}
		else{
			RGBUtils out(in.getWidth(), in.getHeight());

			out.setChannel(0, in.getImage());
			out.setChannel(1, in.getImage());
			out.setChannel(2, in.getImage());

			out.writeToPath(outPath);
		}

	}
	else{
		RGBUtils in(inPath);
		FloatUtils gray(in.convertToGrayscale());

		if(outType == BYTE_TYPE){
			ByteUtils out(FloatUtils::floatToByteImage(in.convertToGrayscale()));
			out.writeToPath(outPath);
		}
		else if(outType == INT_TYPE){
			IntUtils out(FloatUtils::floatToIntImage(in.convertToGrayscale()));
			out.writeToPath(outPath);
		}
		else if(outType == FLOAT_TYPE){
			FloatUtils out(in.convertToGrayscale());
			out.writeToPath(outPath);
		}
		else{
			RGBUtils out(in.getImage());
			out.writeToPath(outPath);
		}

	}

	return 0;
}
