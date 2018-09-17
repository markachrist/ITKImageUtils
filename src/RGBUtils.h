/*
 * RGBUtils.h
 *
 *  Created on: Oct 28, 2012
 *      Author: mchristopher
 */

#include "ImageUtils.h"
#include "itkRGBToLuminanceImageFilter.h"
#include "itkVectorResampleImageFilter.h"

#ifndef RGBUTILS_H_
#define RGBUTILS_H_

//Constants used to define channels, can be used in calls to getChannel()
#define RGBUTILS_RED 0
#define RGBUTILS_GREEN 1
#define RGBUTILS_BLUE 2
#define RGBUTILS_RGOPP 3
#define RGBUTILS_BYOPP 4
#define RGBUTILS_GRAY 5

/**
 * Image utility object for RGB images. Provides extensions to handle color channel
 * data access. Image operations such as blurring are performed on each channel
 * separately.
 */
class RGBUtils : public ImageUtils<RGBImageType> {

protected:

	RGBImageType::Pointer createBlankImage(int w, int h);

public:

	//Just call superclass constructors in each case
	RGBUtils() : ImageUtils<RGBImageType>(){}
	RGBUtils(RGBImageType::Pointer i) : ImageUtils<RGBImageType>(i){}
	RGBUtils(int w, int h, RGBPixelType p) : ImageUtils<RGBImageType>(w, h, p){};
	RGBUtils(int w, int h);
//	RGBUtils(std::string &path) : ImageUtils<RGBImageType>(path){}
	RGBUtils(std::string path) : ImageUtils<RGBImageType>(path){}
	RGBUtils(const char *path) : ImageUtils<RGBImageType>(path){}

	FloatImageType::Pointer getChannel(int c);
	void setChannel(int c, FloatImageType::Pointer data);

	FloatImageType::Pointer getRGOppChannel();
	FloatImageType::Pointer getBYOppChannel();
	FloatImageType::Pointer convertToGrayscale();

	RGBImageType::Pointer resize(int w, int h);
	RGBImageType::Pointer crop(int x, int y, int w, int h);
	RGBImageType::Pointer cropWithPadding(int x, int y, int w, int h);
	RGBImageType::Pointer blur(int w, double sigma);
	RGBImageType::Pointer average(int w, int h);
	RGBImageType::Pointer convolve(FloatImageType::Pointer filter);
	RGBImageType::Pointer convolve(RGBImageType::Pointer filter);
	RGBImageType::Pointer correlate(FloatImageType::Pointer filter);
	RGBImageType::Pointer correlate(RGBImageType::Pointer filter);
	RGBImageType::Pointer filter(FloatImageType::Pointer filter);
	RGBImageType::Pointer filter(RGBImageType::Pointer filter);
	RGBImageType::Pointer threshold(RGBImageType::PixelType t);
	RGBImageType::Pointer flipX();
	RGBImageType::Pointer flipY();
	RGBImageType::Pointer simpleXGradient();
	RGBImageType::Pointer simpleYGradient();

	RGBImageType::Pointer mult(double s);
	RGBImageType::Pointer mult(RGBImageType::Pointer s);
	RGBImageType::Pointer add(double s);
	RGBImageType::Pointer add(RGBImageType::Pointer s);
	RGBImageType::Pointer round();
	RGBImageType::Pointer mirrorPadImage(int l, int r, int a, int b);

	RGBImageType::PixelType max();
	RGBImageType::PixelType min();
	RGBImageType::IndexType minLocation();
	RGBImageType::IndexType maxLocation();

	RGBImageType::PixelType otsuThreshold();

	virtual ~RGBUtils();
};

#endif /* RGBUTILS_H_ */
