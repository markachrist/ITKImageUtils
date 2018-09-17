/*
 * RGBUtils.cpp
 *
 *  Created on: Oct 28, 2012
 *      Author: mchristopher
 */

#include "RGBUtils.h"

/**
 * Create a new RGB image of size w x h with black (zero) pixels in each channel.
 *
 * @param w
 */
RGBUtils::RGBUtils(int w, int h){

	RGBImageType::PixelType black;
	black.Fill(0);

	itk::Index<2> idx;
	idx.Fill(0);

	itk::Size<2> s;
	s[0] = w; s[1] = h;

	RGBImageType::RegionType region;
	region.SetIndex(idx);
	region.SetSize(s);

	RGBImageType::Pointer img = RGBImageType::New();
	img->SetRegions(region);
	img->Allocate();
	img->FillBuffer(black);

	this->setImage(img);
}

/**
 * Returns the pixel data from the image referenced by this utility object for the
 * indicated component (channel). Pixel values in the channel are converted to doubles.
 *
 * @param c
 *   The component for which data should be retrieved. Values of 0-2 indicate RGB channels, value of
 *   3 indicates RG opponent channel, value of 4 indicates BY opponent channel.
 * @return
 *   The data from the indicated channel (new copy).
 */
FloatImageType::Pointer RGBUtils::getChannel(int c){

	if(c == RGBUTILS_RGOPP){
		return this->getRGOppChannel();
	}
	else if(c == RGBUTILS_BYOPP){
		return this->getBYOppChannel();
	}
	else if(c == RGBUTILS_GRAY){
		return this->convertToGrayscale();
	}

	//Select channel
	typedef itk::NthElementImageAdaptor<RGBImageType, double> Channel;
	Channel::Pointer channel = Channel::New();
	channel->SelectNthElement(c);
	channel->SetImage(this->image);

	//Use cast filter to copy data to output image
	itk::CastImageFilter<Channel, FloatImageType>::Pointer copy =
			itk::CastImageFilter<Channel, FloatImageType>::New();
	copy->SetInput(channel);
	copy->Update();

	return copy->GetOutput();
}

FloatImageType::Pointer RGBUtils::getRGOppChannel(){

	FloatImageType::Pointer red = this->getChannel(0);
	FloatUtils r(red);
	FloatImageType::Pointer green = this->getChannel(1);
	FloatUtils g(green);

	return r.add(g.mult(-1.0));
}


FloatImageType::Pointer RGBUtils::getBYOppChannel(){

	FloatImageType::Pointer red = this->getChannel(0);
	FloatUtils r(red);
	FloatImageType::Pointer green = this->getChannel(1);
	FloatUtils g(green);
	FloatImageType::Pointer blue = this->getChannel(2);
	FloatUtils b(blue);

	FloatImageType::Pointer yellow = r.add(g.getImage());
	FloatUtils y(yellow);

	return y.add(b.mult(-1.0));
}

/**
 * Replaces the pixel data for the given channel in the image referenced by this utility
 * object with the given image data. Pixel values are cast to RGBComponentType.
 *
 * @param c
 *   Index of the channel to replace.
 * @param data
 *   Image with data that should be copied into the component
 */
void RGBUtils::setChannel(int c, FloatImageType::Pointer data){

	itk::ImageRegionIterator<RGBImageType> it(this->image, this->image->GetLargestPossibleRegion());
	itk::ImageRegionIterator<FloatImageType> fit(data, data->GetLargestPossibleRegion());

	RGBPixelType black;
	black[0] = 0; black[1] = 0; black[2] = 0;

	fit.GoToBegin();
	for(it.GoToBegin(); !it.IsAtEnd(); ++it){
		RGBPixelType p = it.Get();
		p[c] = (RGBComponentType)(fit.Get() + 0.5);

		int r = p[0]; int g = p[1]; int b = p[2];
		float f = fit.Get();

		it.Set(p);
		++fit;
	}
}

/**
 * Resizes the image to the given size. Bicubic interpolation is used to compute pixel values.
 *
 * @param w
 *   The new width of the image
 * @param h
 *   The new height of the image
 *
 */
RGBImageType::Pointer RGBUtils::resize(int w, int h){

	RGBImageType::SizeType inputSize = this->image->GetLargestPossibleRegion().GetSize();
	RGBImageType::SizeType outputSize = inputSize;
	outputSize[0] = w;
	outputSize[1] = h;

	RGBImageType::SpacingType outputSpacing;
	outputSpacing[0] = this->image->GetSpacing()[0] * (static_cast<double>(inputSize[0]) / static_cast<double>(outputSize[0]));
	outputSpacing[1] = this->image->GetSpacing()[1] * (static_cast<double>(inputSize[1]) / static_cast<double>(outputSize[1]));

	itk::VectorResampleImageFilter<RGBImageType, RGBImageType>::Pointer rescale =
			itk::VectorResampleImageFilter<RGBImageType, RGBImageType>::New();
	rescale->SetInput(this->image);
	rescale->SetSize(outputSize);
	rescale->SetOutputSpacing(outputSpacing);
	rescale->SetTransform(itk::IdentityTransform<double, 2>::New());
	rescale->UpdateLargestPossibleRegion();
	rescale->Update();

	RGBImageType::Pointer out = rescale->GetOutput();
	out->SetSpacing(this->image->GetSpacing());

	return out;
}

/**
 * Crop the image to the specified rectangle.
 *
 * @param x
 *   x-coord of upper-left hand corner of the rectangle
 * @param y
 *   y-coord of upper-left hand corner of the rectangle
 * @param w
 *   width of the rectangle
 * @param h
 *   height of the rectangle
 * @return
 *   Image representing the region given by the rectangle
 */
RGBImageType::Pointer RGBUtils::crop(int x, int y, int w, int h){

	itk::Index<2> idx;
	idx[0] = x;
	idx[1] = y;

	itk::Size<2> s;
	s[0] = w;
	s[1] = h;

	itk::ImageRegion<2> r;
	r.SetIndex(idx);
	r.SetSize(s);

	itk::RegionOfInterestImageFilter<RGBImageType, RGBImageType>::Pointer cropper =
			itk::RegionOfInterestImageFilter<RGBImageType, RGBImageType>::New();

	cropper->SetRegionOfInterest(r);
	cropper->SetInput(this->image);

	cropper->Update();

	RGBImageType::Pointer out = cropper->GetOutput();
	itk::Point<double, 2> o;
	o.Fill(0);
	out->SetOrigin(o);

	return out;
}

/**
 * Crop the image to the specified rectangle. If the region extends outside the image,
 * the returned image is padded with 0-valued pixels.
 *
 * @param x
 *   x-coord of upper-left hand corner of the rectangle
 * @param y
 *   y-coord of upper-left hand corner of the rectangle
 * @param w
 *   width of the rectangle
 * @param h
 *   height of the rectangle
 * @return
 *   Image representing the region given by the rectangle
 */
RGBImageType::Pointer RGBUtils::cropWithPadding(int x, int y, int w, int h){

	itk::Index<2> idx;
	idx[0] = x;
	idx[1] = y;

	itk::Size<2> s;
	s[0] = w;
	s[1] = h;

	itk::ImageRegion<2> r;
	r.SetIndex(idx);
	r.SetSize(s);

	int px0 = 0, px1 = 0;
	int py0 = 0, py1 = 0;

	if(x < 0){
		px0 = -x;
	}
	if(this->getWidth() - (x + w) < 0){
		px1 = -(this->getWidth() - (x + w));
	}
	if(y < 0){
		py0 = -y;
	}
	if(this->getHeight() - (y + h) < 0){
		py1 = -(this->getHeight() - (y + h));
	}

	itk::RegionOfInterestImageFilter<RGBImageType, RGBImageType>::Pointer cropper =
			itk::RegionOfInterestImageFilter<RGBImageType, RGBImageType>::New();

	if(px0 > 0 || px1 > 0 || py0 > 0 || py1 > 0){

		RGBImageType::PixelType zero;
		zero.Fill(0);

		RGBUtils padded(this->getWidth() + px0 + px1, this->getHeight() + py0 + py1, zero);

		RGBImageType::IndexType idx;
		idx[0] = px0;
		idx[1] = py0;

		RGBImageType::SizeType size;
		size[0] = this->getWidth();
		size[1] = this->getHeight();

		RGBImageType::RegionType region(idx, size);

		itk::ImageRegionIterator<RGBImageType> it(padded.getImage(), region);
		itk::ImageRegionIterator<RGBImageType> origIt(this->getImage(), this->image->GetLargestPossibleRegion());

		origIt.GoToBegin();
		for(it.GoToBegin(); !it.IsAtEnd(); ++it){
			it.Set(origIt.Value());
			++origIt;
		}

		idx[0] = x + px0;
		idx[1] = y + py0;
		r.SetIndex(idx);
		cropper->SetRegionOfInterest(r);
		cropper->SetInput(padded.getImage());
		cropper->Update();
	}
	else{
		cropper->SetRegionOfInterest(r);
		cropper->SetInput(this->image);
		cropper->Update();
	}

	RGBImageType::Pointer out = cropper->GetOutput();
	itk::Point<double, 2> o;
	o.Fill(0);
	out->SetOrigin(o);

	return out;
}

/**
 * Applies Gaussian blurring using the a filter specified by w and sigma to the image.
 * Performs blur in each color channel independently.
 *
 * @param w
 *   Half width of the Gaussian filter. Filter applied to image is square with size 2*w + 1
 * @param sigma
 *   The standard deviation of the Gaussian blurring function
 * @return
 *   Blurred version of the image
 */
RGBImageType::Pointer RGBUtils::blur(int w, double sigma){

	RGBImageType::SizeType s = this->image->GetLargestPossibleRegion().GetSize();
	RGBImageType::Pointer out = (RGBImageType::Pointer) this->createBlankImage(s[0], s[1]);
	RGBUtils rout(out);

	for(int i = 0; i < 3; ++i){
		FloatUtils fc(this->getChannel(i));
		rout.setChannel(i, fc.blur(w, sigma));
	}

	return rout.getImage();
}

/**
 * Converts the RGB image to grayscale, float valued image. Uses following weighting scheme
 * to compute grayscale values:
 *
 * gray = 0.299*r + 0.587*g + 0.114*b
 *
 * @return
 *   Grayscale version of the RGB image
 */
FloatImageType::Pointer RGBUtils::convertToGrayscale(){

	double coef[3] = {0.2989, 0.5870, 0.1140};
	FloatUtils gray(this->getWidth(), this->getHeight(), 0.0);

	for(int c = 0; c < 3; ++c){
		FloatUtils channel(this->getChannel(c));
		gray.setImage(gray.add(channel.mult(coef[c])));
	}

	return gray.getImage();

}

/**
 * Applies rectangular averaging filter to the image. The filter size used is equal to
 * (2w + 1) by (2h + 1). The filter is applied to each channel separately.
 *
 * @param w
 *   Half width of the filter
 * @param h
 *   Half height of the filter
 * @return
 *   The result of applying the averaging filter to each channel the image
 */
RGBImageType::Pointer RGBUtils::average(int w, int h){
	RGBUtils result(this->copyImage());
	FloatUtils avg;

	for(int c = 0; c < 3; ++c){
		avg.setImage(this->getChannel(c));
		result.setChannel(c, avg.average(w, h));
	}

	return result.getImage();
}

/**
 * Convolve this image with the given filter. The convolution is performed in each RGB channel
 * independently by applying the filter to each channel.
 *
 * Mirroring is performed at the image edge to perform the convolution.
 *
 * @param filter
 *   The filter by which to convolve this image
 * @return
 *   The result of the convolution
 */
RGBImageType::Pointer RGBUtils::convolve(FloatImageType::Pointer filter){
	RGBUtils result(this->copyImage());
	FloatUtils chan;

	for(int c = 0; c < 3; ++c){
		chan.setImage(this->getChannel(c));
		chan.setImage(chan.convolve(filter));
		result.setChannel(c, chan.getImage());
	}

	return result.getImage();
}

/**
 * Convolve this image with the given filter. The convolution is performed in each RGB channel
 * independently using the corresponding channels of the filter.
 *
 * Mirroring is performed at the image edge to perform the convolution.
 *
 * @param filter
 *   The filter by which to convolve this image
 * @return
 *   The result of the convolution
 */
RGBImageType::Pointer RGBUtils::convolve(RGBImageType::Pointer filter){
	RGBUtils result(this->copyImage());
	RGBUtils filt(filter);
	FloatUtils chan;

	for(int c = 0; c < 3; ++c){
		chan.setImage(this->getChannel(c));
		chan.setImage(chan.convolve(filt.getChannel(c)));
		result.setChannel(c, chan.getImage());
	}

	return result.getImage();
}

/**
 * Correlate this image with the given filter. The correlation is performed in each RGB channel
 * independently by applying the filter to each channel.
 *
 * @param filter
 *   The filter with which to correlate this image
 * @return
 *   The result of the normalized correlation
 */
RGBImageType::Pointer RGBUtils::correlate(FloatImageType::Pointer filter){
	RGBUtils result(this->copyImage());
	FloatUtils chan;

	for(int c = 0; c < 3; ++c){
		chan.setImage(this->getChannel(c));
		chan.setImage(chan.correlate(filter));
		result.setChannel(c, chan.getImage());
	}

	return result.getImage();
}

/**
 * Correlate this image with the given filter. The correlation is performed in each RGB channel
 * independently using the corresponding channels of the filter.
 *
 * @param filter
 *   The filter with which to correlate this image
 * @return
 *   The result of the normalized correlation
 */
RGBImageType::Pointer RGBUtils::correlate(RGBImageType::Pointer filter){
	RGBUtils result(this->copyImage());
	RGBUtils filt(filter);
	FloatUtils chan;

	for(int c = 0; c < 3; ++c){
		chan.setImage(this->getChannel(c));
		chan.setImage(chan.correlate(filt.getChannel(c)));
		result.setChannel(c, chan.getImage());
	}

	return result.getImage();
}

/**
 * Filter this image with the given image. The correlation is performed in each RGB channel
 * independently by applying the filter to each channel.
 *
 * @param filter
 *   The filter to apply to this image
 * @return
 *   The result of the filtration
 */
RGBImageType::Pointer RGBUtils::filter(FloatImageType::Pointer filter){
	RGBUtils result(this->copyImage());
	FloatUtils chan;

	for(int c = 0; c < 3; ++c){
		chan.setImage(this->getChannel(c));
		chan.setImage(chan.filter(filter));
		result.setChannel(c, chan.getImage());
	}

	return result.getImage();
}

/**
 * Filter this image with the given image. The correlation is performed in each RGB channel
 * independently using the corresponding channels of the filter.
 *
 * @param filter
 *   The filter to apply to this image
 * @return
 *   The result of the filtration
 */
RGBImageType::Pointer RGBUtils::filter(RGBImageType::Pointer filter){
	RGBUtils result(this->copyImage());
	RGBUtils filt(filter);
	FloatUtils chan;

	for(int c = 0; c < 3; ++c){
		chan.setImage(this->getChannel(c));
		chan.setImage(chan.filter(filt.getChannel(c)));
		result.setChannel(c, chan.getImage());
	}

	return result.getImage();
}

/**
 * Apply the given threshold to the image. Produces a binary output with pixels >= t set equal
 * to 1 and all others set to 0. The thresholding is applied independently to each channel, with
 * the RGB value of t applied to corresponding channels of the image.
 *
 * @param t
 *   The threshold to apply to each channel
 * @return
 *   An image with binary RGB channels resulting from the thresholding
 */
RGBImageType::Pointer RGBUtils::threshold(RGBImageType::PixelType t){
	RGBUtils result(this->copyImage());
	FloatUtils chan;

	for(int c = 0; c < 3; ++c){
		chan.setImage(this->getChannel(c));
		chan.setImage(chan.threshold((t[c])));
		result.setChannel(c, chan.getImage());
	}

	return result.getImage();
}

/**
 * Calculates the gradient of the image in the x direction in each channel independently.
 *
 * The 1D gradient is calculated as a two side difference, as if by passing the
 * filter [0.5 0.0 -0.5] over the image. Image boundaries are handled by taking
 * a one sided difference.
 *
 * @return
 *   The gradient of this image in the x direction.
 */
RGBImageType::Pointer RGBUtils::simpleXGradient(){
	RGBUtils result(this->copyImage());
	FloatUtils avg;

	for(int c = 0; c < 3; ++c){
		avg.setImage(this->getChannel(c));
		result.setChannel(c, avg.simpleXGradient());
	}

	return result.getImage();
}

/**
 * Calculates the gradient of the image in the y direction in each channel independently.
 *
 * The 1D gradient is calculated as a two side difference, as if by passing the
 * transpose of the filter [0.5 0.0 -0.5] over the image. Image boundaries are handled by
 * taking a one sided difference.
 *
 * @return
 *   The gradient of this image in the y direction.
 */
RGBImageType::Pointer RGBUtils::simpleYGradient(){
	RGBUtils result(this->copyImage());
	FloatUtils avg;

	for(int c = 0; c < 3; ++c){
		avg.setImage(this->getChannel(c));
		result.setChannel(c, avg.simpleYGradient());
	}

	return result.getImage();
}

/**
 * Multiply each pixel in the image by the given scale factor.
 *
 * Each channel is multiplied by the factor.
 *
 * @param s
 *   Value by which to multiply each pixel
 * @return
 *   Image after performing multiplication
 */
RGBImageType::Pointer RGBUtils::mult(double s){
	RGBUtils result(this->copyImage());
	FloatUtils avg;

	for(int c = 0; c < 3; ++c){
		avg.setImage(this->getChannel(c));
		result.setChannel(c, avg.mult(s));
	}

	return result.getImage();
}

/**
 * Performs element-by-element multiplication between this image and s.
 *
 * Each channel is handled separately.
 *
 * @param s
 *   Image by which to multiply this image
 * @return
 *   Image after performing multiplication
 */
RGBImageType::Pointer RGBUtils::mult(RGBImageType::Pointer s){
	RGBUtils result(this->copyImage());
	RGBUtils chan(s);
	FloatUtils avg;

	for(int c = 0; c < 3; ++c){
		avg.setImage(this->getChannel(c));
		result.setChannel(c, avg.mult(chan.getChannel(c)));
	}

	return result.getImage();
}

/**
 * Add the given value to each pixel in the image.
 *
 * The value is added to each channel.
 *
 * @param s
 *   Value to add to each pixel
 * @return
 *   Image after performing addition
 */
RGBImageType::Pointer RGBUtils::add(double s){
	RGBUtils result(this->copyImage());
	FloatUtils avg;

	for(int c = 0; c < 3; ++c){
		avg.setImage(this->getChannel(c));
		result.setChannel(c, avg.add(s));
	}

	return result.getImage();
}

/**
 * Performs element-by-element multiplication between this image and s.
 *
 * Each channel is handled separately.
 *
 * @param s
 *   Image by which to multiply this image
 * @return
 *   Image after performing multiplication
 */
RGBImageType::Pointer RGBUtils::add(RGBImageType::Pointer s){
	RGBUtils result(this->copyImage());
	RGBUtils chan(s);
	FloatUtils avg;

	for(int c = 0; c < 3; ++c){
		avg.setImage(this->getChannel(c));
		result.setChannel(c, avg.add(chan.getChannel(c)));
	}

	return result.getImage();
}

/**
 * Rounds the value of each pixel in this image.
 *
 * Effectively does nothing for RGB images.
 *
 * @return
 *   Image after performing the rounding
 */
RGBImageType::Pointer RGBUtils::round(){
	RGBUtils result(this->copyImage());
	FloatUtils avg;

	for(int c = 0; c < 3; ++c){
		avg.setImage(this->getChannel(c));
		result.setChannel(c, avg.round());
	}

	return result.getImage();
}

/**
 * Finds pixel value with maximum luminance.
 *
 * @return
 *   Value of the pixel with maximum luminance.
 */
RGBImageType::PixelType RGBUtils::max(){

	FloatUtils f(this->convertToGrayscale());
	itk::Index<2> idx = f.maxLocation();

	return this->image->GetPixel(idx);
}

/**
 * Finds pixel value with minimum luminance.
 *
 * @return
 *   Location of the pixel with minimum luminance.
 */
RGBImageType::IndexType RGBUtils::maxLocation(){

	FloatUtils f(this->convertToGrayscale());
	return f.maxLocation();
}

/**
 * Finds pixel value with minimum luminance.
 *
 * @return
 *   Value of the pixel with minimum luminance.
 */
RGBImageType::PixelType RGBUtils::min(){

	FloatUtils f(this->convertToGrayscale());
	itk::Index<2> idx = f.minLocation();

	return this->image->GetPixel(idx);
}

/**
 * Finds pixel value with minimum luminance.
 *
 * @return
 *   Location of the pixel with minimum luminance.
 */
RGBImageType::IndexType RGBUtils::minLocation(){

	FloatUtils f(this->convertToGrayscale());
	return f.minLocation();
}

/**
 * Finds Otsu threshold for each of the three color channels.
 *
 * @return
 *   RGB pixel with each element having the Otsu threshold for the corresponding channel in this image
 */
RGBImageType::PixelType RGBUtils::otsuThreshold(){

	RGBPixelType threshold;
	threshold.Fill(0);

	for(int c = 0; c < 3; ++c){
		itk::CastImageFilter<FloatImageType, ByteImageType>::Pointer cast =
				itk::CastImageFilter<FloatImageType, ByteImageType>::New();
		cast->SetInput(this->getChannel(c));
		cast->Update();

		ByteUtils channel(cast->GetOutput());

		threshold[c] = channel.otsuThreshold();
	}

	return threshold;
}

/**
 * Creates a blank image of the same type as the image referenced by this utility object.
 *
 * @param w
 *   Width of the created image
 * @param h
 *   Height of the created image
 * @return
 *   Image of the indicated size and the same type as the image referenced by this object,
 *   with a pixel values set to 0
 */
RGBImageType::Pointer RGBUtils::createBlankImage(int w, int h){

	itk::Index<2> idx;
	idx.Fill(0);

	itk::Size<2> s;
	s[0] = w; s[1] = h;

	RGBImageType::RegionType region;
	region.SetIndex(idx);
	region.SetSize(s);

	RGBImageType::Pointer img = RGBImageType::New();
	img->SetRegions(region);
	img->Allocate();

	RGBPixelType black;
	black[0] = 0; black[1] = 0; black[2] = 0;
	img->FillBuffer(black);

	return img;
}

RGBUtils::~RGBUtils() {
	// TODO Auto-generated destructor stub
}
