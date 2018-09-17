/*
 * transformMain.cpp
 *
 *  Created on: May 18, 2015
 *      Author: mchristopher
 */
#include <iostream>

#include "itkAffineTransform.h"
#include "itkCenteredAffineTransform.h"

#include "ImageUtils.h"
#include "RGBUtils.h"
#include "CVMatUtils.h"

#define NUM_PARAMS 6

void usage(int argc, char **argv){

	std::cout << "%" << argv[0] << " <input image> <output image> <a> <b> <c> <d> <tx> <ty>" << std::endl << std::endl;
	std::cout << "Applies an affine transform to the given image and stores output at the given output destination.";

	std::cout << std::endl << std::endl;

	exit(-1);
}

/**
 *
 */
int main(int argc, char **argv){

	//Check cl arguments
	std::string inPath = argv[1];
	std::string outPath = argv[2];

	if(argc <= 8){
		std::cout << "Error: Missing arguments!" << std::endl;
		usage(argc, argv);
	}

	//Create transform
	itk::AffineTransform<double, 2>::Pointer transform = itk::AffineTransform<double, 2>::New();
	itk::AffineTransform<double, 2>::ParametersType params(6);

//	itk::CenteredAffineTransform<double, 2>::Pointer transform = itk::CenteredAffineTransform<double, 2>::New();
//	itk::CenteredAffineTransform<double, 2>::ParametersType params(8);

	for(int i = 3; i < argc; ++i){
		std::cout << "Reading transform value " << argv[i] << " at index " << i << "( " << (i - 3) << " )" << std::endl;
		params[i - 3] = atof(argv[i]);
	}
	transform->SetParameters(params);

	//Now apply transform
	RGBUtils in(inPath);
	std::cout << "Input size = " << in.getWidth() << " x " << in.getHeight() << std::endl;

	itk::ResampleImageFilter<RGBImageType, RGBImageType>::Pointer mapper =
			itk::ResampleImageFilter<RGBImageType, RGBImageType>::New();
	mapper->SetInput(in.getImage());
	mapper->SetTransform(transform);

	RGBImageType::PixelType def;
	def.Fill(0);

	mapper->SetSize(in.getImage()->GetLargestPossibleRegion().GetSize());
	mapper->SetOutputOrigin(in.getImage()->GetOrigin());
	mapper->SetOutputSpacing(in.getImage()->GetSpacing());
	mapper->SetOutputDirection(in.getImage()->GetDirection());
	mapper->SetDefaultPixelValue(def);
	mapper->Update();

	in.setImage(mapper->GetOutput());
	std::cout << "Output size = " << in.getWidth() << " x " << in.getHeight() << std::endl;
	in.writeToPath(outPath);

	return 0;
}
