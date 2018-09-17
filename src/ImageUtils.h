/*
 * ImageUtils.h
 *
 *  Created on: Oct 24, 2012
 *      Author: mchristopher
 */

#include "itkImage.h"
#include "itkRGBPixel.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVectorResampleImageFilter.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkIdentityTransform.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkVectorImage.h"
#include "itkNthElementImageAdaptor.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkFFTConvolutionImageFilter.h"
#include "itkFFTNormalizedCorrelationImageFilter.h"
#include "itkConvolutionImageFilter.h"
#include "itkConstantBoundaryCondition.h"
#include "itkMeanImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkBSplineDownsampleImageFilter.h"
#include "itkBSplineUpsampleImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkDivideOrZeroOutImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkConstantBoundaryCondition.h"
#include "itkImageDuplicator.h"
#include "itkSqrtImageFilter.h"
#include "itkRoundImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkNormalizeImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkLineIterator.h"
#include "itkTileImageFilter.h"

#ifndef IMAGEUTILS_H_
#define IMAGEUTILS_H_

//Define 2D byte image types
typedef unsigned char BytePixelType;
typedef itk::Image<BytePixelType, 2> ByteImageType;

//Define 2D byte image types
typedef int IntPixelType;
typedef itk::Image<IntPixelType, 2> IntImageType;

//Define 2D RGB image types
typedef unsigned char RGBComponentType;
typedef itk::RGBPixel<RGBComponentType> RGBPixelType;
typedef itk::Image<RGBPixelType, 2> RGBImageType;

//Define 2D float grayscale image types
typedef float FloatPixelType;
typedef itk::Image<FloatPixelType, 2> FloatImageType;

/**
 *  Provides interface for performing simple and common image operations such as basic I/O
 *  operations and simple image manipulations (resizing, blurring, cropping, etc.).
 *
 */
template<class TImage>
class ImageUtils {

protected:

	/** Image for which this utility object serves as a wrapper*/
	typename TImage::Pointer image;

	typename TImage::Pointer createBlankImage(int w, int h);
	static FloatImageType::Pointer createGaussianImage(int w, double sigma);

public:

	ImageUtils(){};
	ImageUtils(typename TImage::Pointer i);
	ImageUtils(int w, int h, typename TImage::PixelType p);

	//Some I/O methods
	ImageUtils(std::string path);
	ImageUtils(const char *path);
	void writeToPath(std::string &path);
	void writeToPath(const char *path);
	void printImage();

	void reallocate(int w, int h, typename TImage::PixelType p);

	//Image operations
	typename TImage::Pointer copyImage();
	typename TImage::Pointer resize(int w, int h);
	typename TImage::Pointer downSample();
	typename TImage::Pointer crop(int x, int y, int w, int h);
	typename TImage::Pointer cropWithPadding(int x, int y, int w, int h);
	typename TImage::Pointer blur(int w, double sigma);
	typename TImage::Pointer average(int w, int h);
	typename TImage::Pointer convolve(typename TImage::Pointer filter);
	typename TImage::Pointer correlate(typename TImage::Pointer filter);
	typename TImage::Pointer filter(typename TImage::Pointer filter);
	typename TImage::Pointer threshold(typename TImage::PixelType t);
	typename TImage::Pointer normalize();
	typename TImage::Pointer flipX();
	typename TImage::Pointer flipY();
	typename TImage::Pointer simpleXGradient();
	typename TImage::Pointer simpleYGradient();
	double sumPixels();

	typename TImage::Pointer mult(double s);
	typename TImage::Pointer mult(typename TImage::Pointer s);
	typename TImage::Pointer divide(typename TImage::Pointer s);
	typename TImage::Pointer add(double s);
	typename TImage::Pointer add(typename TImage::Pointer s);
	typename TImage::Pointer sqrt();
	typename TImage::Pointer round();
	typename TImage::Pointer copyImageAt(typename TImage::Pointer s, typename TImage::IndexType loc);
	typename TImage::Pointer mirrorPadImage(int l, int r, int a, int b);
	typename TImage::Pointer stackImages(typename TImage::Pointer s, int dir);
	typename TImage::Pointer drawLine(typename TImage::IndexType start, typename TImage::IndexType end,
			typename TImage::PixelType value);

	typename TImage::PixelType getBicubicInterpolatedValue(double x, double y);

	typename TImage::PixelType max();
	typename TImage::IndexType maxLocation();
	typename TImage::PixelType min();
	typename TImage::IndexType minLocation();

	typename TImage::PixelType otsuThreshold();

	bool isInBounds(typename TImage::IndexType idx);

	void sampleIndices(int n, std::vector<typename TImage::IndexType> &idxs);
	void sampleIndices(int n, std::vector<typename TImage::IndexType> &idxs, ByteImageType::Pointer mask);

	typename TImage::Pointer markLocations(std::vector<typename TImage::IndexType> &idxs,
			typename TImage::PixelType value, bool singlePoint = true);

	//Conversion methods
	static ByteImageType::Pointer floatToByteImage(FloatImageType::Pointer f);
	static FloatImageType::Pointer byteToFloatImage(ByteImageType::Pointer b);
	static ByteImageType::Pointer intToByteImage(IntImageType::Pointer i);
	static IntImageType::Pointer byteToIntImage(ByteImageType::Pointer b);
	static IntImageType::Pointer floatToIntImage(FloatImageType::Pointer f);
	static FloatImageType::Pointer intToFloatImage(IntImageType::Pointer i);

	//Getters/Setters
	inline int getWidth(){
		return this->image->GetLargestPossibleRegion().GetSize()[0];
	}

	inline int getHeight(){
		return this->image->GetLargestPossibleRegion().GetSize()[1];
	}

	inline void setImage(typename TImage::Pointer i){
		this->image = i;
	}

	inline typename TImage::Pointer getImage(){
		return this->image;
	}

	virtual ~ImageUtils();
};

/**
 * Image utility object for images with byte valued pixels.
 */
typedef ImageUtils<ByteImageType> ByteUtils;

/**
 * Image utility object for images with integer valued pixels.
 */
typedef ImageUtils<IntImageType> IntUtils;

/**
 * Image utility object for images with float valued pixels.
 */
typedef ImageUtils<FloatImageType> FloatUtils;

/**
 * Create utility object for the given image.
 *
 * @param i
 *  Image for which utility object should be created
 */
template<class TImage>
ImageUtils<TImage>::ImageUtils(typename TImage::Pointer i) {
	this->setImage(i);
}

/**
 * Create a new utility object and image with given dimensions, filled with the
 * value p.
 *
 * @param w
 *   Width of the image to create
 * @param h
 *   Height of the image to create
 * @param p
 *   Value to assign to each pixel in the new image
 */
template<class TImage>
ImageUtils<TImage>::ImageUtils(int w, int h, typename TImage::PixelType p){

	itk::Index<2> idx;
	idx.Fill(0);

	itk::Size<2> s;
	s[0] = w; s[1] = h;

	typename TImage::RegionType region;
	region.SetIndex(idx);
	region.SetSize(s);

	typename TImage::Pointer img = TImage::New();
	img->SetRegions(region);
	img->Allocate();
	img->FillBuffer(p);

	this->setImage(img);
}

/**
 * Create a new utility object and image with given dimensions, filled with the
 * value p.
 *
 * @param w
 *   Width of the image to create
 * @param h
 *   Height of the image to create
 * @param p
 *   Value to assign to each pixel in the new image
 */
template<class TImage>
void ImageUtils<TImage>::reallocate(int w, int h, typename TImage::PixelType p){

	itk::Index<2> idx;
	idx.Fill(0);

	itk::Size<2> s;
	s[0] = w; s[1] = h;

	typename TImage::RegionType region;
	region.SetIndex(idx);
	region.SetSize(s);

	typename TImage::Pointer img = TImage::New();
	img->SetRegions(region);
	img->Allocate();
	img->FillBuffer(p);

	this->setImage(img);

}

/**
 * Create utility object for the image at the specified file path.
 *
 * @param path
 *   Path to the image file
 */
template<class TImage>
ImageUtils<TImage>::ImageUtils(std::string path){
	typename itk::ImageFileReader<TImage>::Pointer reader = itk::ImageFileReader<TImage>::New();
	reader->SetFileName(path);
	reader->Update();
	this->setImage(reader->GetOutput());
}

/**
 * Create utility object for the image at the specified file path.
 *
 * @param path
 *   Path to the image file
 */
template<class TImage>
ImageUtils<TImage>::ImageUtils(const char *path){
	std::string p(path);
	typename itk::ImageFileReader<TImage>::Pointer reader = itk::ImageFileReader<TImage>::New();
	reader->SetFileName(p);
	reader->Update();
	this->setImage(reader->GetOutput());
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
template<class TImage>
typename TImage::Pointer ImageUtils<TImage>::createBlankImage(int w, int h){

	itk::Index<2> idx;
	idx.Fill(0);

	itk::Size<2> s;
	s[0] = w; s[1] = h;

	typename TImage::RegionType region;
	region.SetIndex(idx);
	region.SetSize(s);

	typename TImage::Pointer img = TImage::New();
	img->SetRegions(region);
	img->Allocate();
	img->FillBuffer((typename TImage::PixelType)0);

	return img;
}

/**
 * Creates a duplicate of the image referenced by this utility object.
 *
 * @return
 *   A copy of the image
 */
template<class TImage>
typename TImage::Pointer ImageUtils<TImage>::copyImage(){

	typename itk::ImageDuplicator<TImage>::Pointer dupe =
			itk::ImageDuplicator<TImage>::New();
	dupe->SetInputImage(this->image);
	dupe->Update();
	typename TImage::Pointer out = dupe->GetOutput();

	return out;
}

/**
 *
 */
template<class TImage>
typename TImage::PixelType ImageUtils<TImage>::getBicubicInterpolatedValue(double x, double y){


	double a = -0.5;
	double z;
	int w = this->getWidth();
	int h = this->getHeight();
	int u0 = (int) floor(x);  //use floor to handle negative coordinates too
	int v0 = (int) floor(y);

	itk::Index<2> idx;
	double q = 0;
	for (int j = 0; j <= 3; j++) {
		int v = v0 - 1 + j;
		double p = 0;
		for (int i = 0; i <= 3; i++) {
			int u = u0 - 1 + i;
			idx[0] = (u >= 0) ? ((u < w) ? u : w - 1) : 0;
			idx[1] = (v >= 0) ? ((v < h) ? v : h - 1) : 0;

			z = x - u;
			z = (z < 0) ? -z : z;
			if(z <= 1.0){
				z = z*z*(z*(-a + 2.0) + (a - 3.0)) + 1.0;
			}
			else if(z <= 2.0){
				z = -a*z*z*z + 5.0*a*z*z - 8.0*a*z + 4.0*a;
			}
			else{
				z = 0;
			}

			p = p + static_cast<double>(this->image->GetPixel(idx)) * z;
		}

		z = y - v;
		z = (z < 0) ? -z : z;
		if(z < 1.0){
			z = z*z*(z*(-a + 2.0) + (a - 3.0)) + 1.0;
		}
		else if(z < 2.0){
			z = -a*z*z*z + 5.0*a*z*z - 8.0*a*z + 4.0*a;
		}
		else{
			z = 0;
		}

		q = q + p * z;
	}

	return static_cast<typename TImage::PixelType>(q);
}

/**
 * Resamples the image to transform it to the given given size using bicubic
 * interpolation to calculate new pixel values.
 *
 * @param w
 *   New width of the image
 * @param h
 *   New height of the image
 * @return
 *   The image resized to the given dimensions
 */
template<class TImage>
typename TImage::Pointer ImageUtils<TImage>::resize(int w, int h){

	typename TImage::SizeType inputSize = this->image->GetLargestPossibleRegion().GetSize();
	typename TImage::SizeType outputSize = inputSize;
	outputSize[0] = w;
	outputSize[1] = h;

	typename itk::BSplineInterpolateImageFunction<TImage, double, double>::Pointer interp =
			itk::BSplineInterpolateImageFunction<TImage, double, double>::New();
	interp->SetSplineOrder(3);
	interp->SetInputImage(this->image);

	typename TImage::SpacingType outputSpacing;
	outputSpacing[0] = this->image->GetSpacing()[0] * (static_cast<double>(inputSize[0]) / static_cast<double>(outputSize[0]));
	outputSpacing[1] = this->image->GetSpacing()[1] * (static_cast<double>(inputSize[1]) / static_cast<double>(outputSize[1]));

	typename itk::ResampleImageFilter<TImage, TImage>::Pointer rescale =
			itk::ResampleImageFilter<TImage, TImage>::New();
	rescale->SetInterpolator(interp);
	rescale->SetInput(this->image);
	rescale->SetSize(outputSize);
	rescale->SetOutputSpacing(outputSpacing);
	rescale->SetTransform(itk::IdentityTransform<double, 2>::New());
	rescale->UpdateLargestPossibleRegion();
	rescale->Update();

	typename TImage::Pointer out = rescale->GetOutput();
	out->SetSpacing(this->image->GetSpacing());

	return out;
}

/**
 * Down samples the image by a factor of two.
 *
 * @return
 *   Downsampled version of this image created using bicubic interpolation.
 */
template<class TImage>
typename TImage::Pointer ImageUtils<TImage>::downSample(){

	typename itk::BSplineDownsampleImageFilter<TImage, TImage>::Pointer down =
			itk::BSplineDownsampleImageFilter<TImage, TImage>::New();
//	down->SetSplineOrder(3);
	down->SetSplineOrder(1);
	down->SetInput(this->image);
	down->Update();

	return down->GetOutput();
}

/**
 * Crop the image to the specified rectangle. No bounds checking is performed.
 * If the specified rectangle extends beyond the image boundaries, an error is
 * generated.
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
template<class TImage>
typename TImage::Pointer ImageUtils<TImage>::crop(int x, int y, int w, int h){

	itk::Index<2> idx;
	idx[0] = x;
	idx[1] = y;

	itk::Size<2> s;
	s[0] = w;
	s[1] = h;

	itk::ImageRegion<2> r;
	r.SetIndex(idx);
	r.SetSize(s);

	typename itk::RegionOfInterestImageFilter<TImage, TImage>::Pointer cropper =
			itk::RegionOfInterestImageFilter<TImage, TImage>::New();

	cropper->SetRegionOfInterest(r);
	cropper->SetInput(this->image);
	cropper->Update();

	typename TImage::Pointer out = cropper->GetOutput();
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
template<class TImage>
typename TImage::Pointer ImageUtils<TImage>::cropWithPadding(int x, int y, int w, int h){

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


	typename itk::RegionOfInterestImageFilter<TImage, TImage>::Pointer cropper =
			itk::RegionOfInterestImageFilter<TImage, TImage>::New();

	if(px0 > 0 || px1 > 0 || py0 > 0 || py1 > 0){

		ImageUtils<TImage> padded(this->getWidth() + px0 + px1, this->getHeight() + py0 + py1, 0);

		typename TImage::IndexType i;
		i[0] = px0;
		i[1] = py0;

		typename TImage::SizeType size;
		size[0] = this->getWidth();
		size[1] = this->getHeight();

		typename TImage::RegionType region(i, size);

		typename itk::ImageRegionIterator<TImage> it(padded.getImage(), region);
		typename itk::ImageRegionIterator<TImage> origIt(this->getImage(), this->image->GetLargestPossibleRegion());

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

	typename TImage::Pointer out = cropper->GetOutput();
	itk::Point<double, 2> o;
	o.Fill(0);
	out->SetOrigin(o);

	return out;
}

/**
 * Applies a Gaussian blurring to the image using the given filter width and sigma.
 *
 * Applies the blur symmetrically (sigma x = sigma y = sigma).
 *
 * @param w
 *   The filter half width (total window size is 2*w + 1 in each direction)
 * @param sigma
 *   Standard deviation to use for Gaussian blurring functions
 * @return
 *   Blurred version of the image represented by this utility object
 */
template<class TImage>
typename TImage::Pointer ImageUtils<TImage>::blur(int w, double sigma){

	FloatImageType::Pointer kernel = ImageUtils<TImage>::createGaussianImage(w, sigma);

	itk::Index<2> center;
	center.Fill(w);

	typename TImage::Pointer mirrored = this->mirrorPadImage(w + 1, w + 1, w + 1, w + 1);

	//Perform convolution in Fourier domain, MUCH faster for most kernel sizes
	typename itk::FFTConvolutionImageFilter<TImage, FloatImageType, TImage>::Pointer conv =
			itk::FFTConvolutionImageFilter<TImage, FloatImageType, TImage>::New();

	typename itk::ConstantBoundaryCondition<TImage> zeroPad;
	zeroPad.SetConstant(0);
	conv->SetBoundaryCondition(&zeroPad);

	conv->SetInput(mirrored);
//	conv->SetInput(this->image);
	conv->SetKernelImage(kernel);
	conv->Update();

	ImageUtils<TImage> out(conv->GetOutput());

	return out.crop(w + 1, w + 1, out.getWidth() - 2*(w + 1), out.getHeight() - 2*(w + 1));
}

/**
 * Applies rectangular averaging filter to the image. The filter size used is equal to
 * (2w + 1) by (2h + 1).
 *
 * @param w
 *   Half width of the filter
 * @param h
 *   Half height of the filter
 * @return
 *   The result of applying the averging filter to the image
 */
template<class TImage>
typename TImage::Pointer ImageUtils<TImage>::average(int w, int h){

	itk::Size<2> r;
	r[0] = w; r[1] = h;

	typename itk::MeanImageFilter<TImage, TImage>::Pointer avg =
			itk::MeanImageFilter<TImage, TImage>::New();
	avg->SetRadius(r);
	avg->SetInput(this->image);
	avg->Update();

	return avg->GetOutput();
}

///**
// * Computes standard deviation over recutangular image neighborhoods. The nhood size used is
// * equal to (2w + 1) by (2h + 1).
// *
// * @param w
// *   Half width of the nhood
// * @param h
// *   Half height of the nhood
// * @return
// *   The neighborhood standard deviation of the image
// */
//template<class TImage>
//typename TImage::Pointer ImageUtils<TImage>::stdev(int w, int h){
//
//	typename itk::VarianceImageFunction<TImage>
//
//}

/**
 * Convolve this image with the given filter.
 *
 * Mirroring is performed at the image edge to perform the convolution.
 *
 * @param filter
 *   The filter by which to convovle this image
 * @return
 *   The result of the convolution
 */
template<class TImage>
typename TImage::Pointer ImageUtils<TImage>::convolve(typename TImage::Pointer filter){

	typename TImage::SizeType fsize = filter->GetLargestPossibleRegion().GetSize();
	typename TImage::Pointer mirrored = this->mirrorPadImage(fsize[0]/2, fsize[0]/2, fsize[1]/2, fsize[1]/2);

	//Convolution in fourier space much faster for most filter/image sizes
	typename itk::FFTConvolutionImageFilter<TImage, TImage>::Pointer conv =
			itk::FFTConvolutionImageFilter<TImage, TImage>::New();
	conv->SetKernelImage(filter);
	conv->SetInput(mirrored);
	conv->Update();

	ImageUtils<TImage> result(conv->GetOutput());

	return result.crop(fsize[0]/2, fsize[1]/2, this->getWidth(), this->getHeight());
}

/**
 * Compute the normalized correlation between each location in this image and
 * the given filter.
 *
 * @param filter
 *   The filter with which to correlate this image
 * @return
 *   The result of the normalized correlation
 */
template<class TImage>
typename TImage::Pointer ImageUtils<TImage>::correlate(typename TImage::Pointer filter){

	typename TImage::SizeType fsize = filter->GetLargestPossibleRegion().GetSize();
	typename TImage::Pointer mirrored = this->mirrorPadImage(fsize[0]/2, fsize[0]/2, fsize[1]/2, fsize[1]/2);

	//Convolution in fourier space much faster for most filter/image sizes
	typename itk::FFTNormalizedCorrelationImageFilter<TImage, TImage>::Pointer conv =
			itk::FFTNormalizedCorrelationImageFilter<TImage, TImage>::New();
	conv->SetInput(mirrored);
	conv->SetFixedImage(this->getImage());
	conv->SetMovingImage(filter);
	conv->Update();

	ImageUtils<TImage> result(conv->GetOutput());

	return result.crop(fsize[0]/2, fsize[1]/2, this->getWidth(), this->getHeight());
}

/**
 * Apply the given filter to this image. The filtration is performed by centering
 * the filter at each image location and summing the filter and image pixels values.
 * Essentially, non-normalized correlation.
 *
 * @param filter
 *   The filter to apply to this image
 * @return
 *   The result of the filtration
 */
template<class TImage>
typename TImage::Pointer ImageUtils<TImage>::filter(typename TImage::Pointer filter){

	typename TImage::SizeType fsize = filter->GetLargestPossibleRegion().GetSize();
	typename TImage::Pointer mirrored = this->mirrorPadImage(fsize[0]/2, fsize[0]/2, fsize[1]/2, fsize[1]/2);

	typename itk::FlipImageFilter<TImage>::Pointer flip =
			itk::FlipImageFilter<TImage>::New();
	flip->SetInput(filter);
	itk::FixedArray<bool, 2> axes;
	axes[0] = true;
	axes[1] = true;
	flip->SetFlipAxes(axes);
	flip->Update();

	typename TImage::Pointer flipped = flip->GetOutput();

	//Convolution in fourier space much faster for most filter/image sizes
	typename itk::FFTConvolutionImageFilter<TImage, TImage>::Pointer conv =
			itk::FFTConvolutionImageFilter<TImage, TImage>::New();
	conv->SetKernelImage(flipped);
	conv->SetInput(mirrored);
	conv->Update();

	ImageUtils<TImage> result(conv->GetOutput());

	return result.crop(fsize[0]/2, fsize[1]/2, this->getWidth(), this->getHeight());
}

/**
 * Apply the given threshold to the image. Produces a binary output with pixels >= t set equal
 * to 1 and all others set to 0.
 *
 * @param t
 *   The threshold to apply to each image
 * @return
 *   A binary image result of the thresholding
 */
template<class TImage>
typename TImage::Pointer ImageUtils<TImage>::threshold(typename TImage::PixelType t){

	typename itk::BinaryThresholdImageFilter<TImage, TImage>::Pointer thresh =
			itk::BinaryThresholdImageFilter<TImage, TImage>::New();
	thresh->SetInput(this->image);
	thresh->SetLowerThreshold(t);
	thresh->SetInsideValue(1);
	thresh->SetOutsideValue(0);

	thresh->Update();

	return thresh->GetOutput();
}

/**
 * Normalizes the pixel values of this image. After this call, the pixel values will have
 * zero mean and unit variance.
 *
 * @return
 *   This image with standardized pixel values
 */
template<class TImage>
typename TImage::Pointer ImageUtils<TImage>::normalize(){

	typename itk::NormalizeImageFilter<TImage, TImage>::Pointer norm =
			itk::NormalizeImageFilter<TImage, TImage>::New();
	norm->SetInput(this->image);

	norm->Update();

	return norm->GetOutput();
}

/**
 * Flips the image in the X direction.
 *
 * @return
 *   This image after taking the mirror image about the x-axis
 */
template<class TImage>
typename TImage::Pointer ImageUtils<TImage>::flipX(){

	typename itk::FlipImageFilter<TImage>::Pointer flip =
			itk::FlipImageFilter<TImage>::New();
	flip->SetInput(this->image);
	itk::FixedArray<bool, 2> axes;
	axes[0] = false;
	axes[1] = true;
	flip->SetFlipAxes(axes);
	flip->Update();

	return flip->GetOutput();
}

/**
 * Flips the image in the Y direction.
 *
 * @return
 *   This image after taking the mirror image about the y-axis
 */
template<class TImage>
typename TImage::Pointer ImageUtils<TImage>::flipY(){

	typename itk::FlipImageFilter<TImage>::Pointer flip =
			itk::FlipImageFilter<TImage>::New();
	flip->SetInput(this->image);
	itk::FixedArray<bool, 2> axes;
	axes[0] = true;
	axes[1] = false;
	flip->SetFlipAxes(axes);
	flip->Update();

	return flip->GetOutput();
}

/**
 * Calculates the gradient of the image in the x direction.
 *
 * The 1D gradient is calculated as a two side difference, as if by passing the
 * filter [0.5 0.0 -0.5] over the image. Image boundaries are handled by taking
 * a one sided difference.
 *
 * @return
 *   The gradient of this image in the x direction.
 */
template<class TImage>
typename TImage::Pointer ImageUtils<TImage>::simpleXGradient(){

	FloatUtils filter(3, 3, 0.0);
	itk::Index<2> idx;
	idx[0] = 0; idx[1] = 1;
	filter.getImage()->SetPixel(idx, 0.5);
	idx[0] = 2; idx[1] = 1;
	filter.getImage()->SetPixel(idx, -0.5);

	typename itk::FFTConvolutionImageFilter<TImage, FloatImageType, TImage>::Pointer conv =
			itk::FFTConvolutionImageFilter<TImage, FloatImageType, TImage>::New();

	//Set boundary condition to zero pad
	typename itk::ConstantBoundaryCondition<TImage> zeroPad;
	zeroPad.SetConstant(0);
	conv->SetBoundaryCondition(&zeroPad);

	conv->SetInput(this->image);
	conv->SetKernelImage(filter.getImage());
	conv->Update();

	typename TImage::Pointer out = conv->GetOutput();

	//Handle image boundaries
	idx.Fill(0);
	itk::Size<2> size;
	size[0] = 1; size[1] = this->getHeight();
	typename TImage::RegionType reg(idx, size);

	typename itk::ImageRegionIterator<TImage> outItL(out, reg);
	typename itk::ImageRegionIterator<TImage> targetL(this->image, reg);
	idx[0] = 1;
	reg.SetIndex(idx);
	typename itk::ImageRegionIterator<TImage> sideL(this->image, reg);

	outItL.GoToBegin();
	sideL.GoToBegin();
	for(targetL.GoToBegin(); !targetL.IsAtEnd(); ++targetL){
		typename TImage::PixelType r = sideL.Get() - targetL.Get();
		outItL.Set(r);
		++sideL;
		++outItL;
	}

	idx[0] = this->getWidth() - 1; idx[1] = 0;
	reg.SetIndex(idx);

	typename itk::ImageRegionIterator<TImage> outItR(out, reg);
	typename itk::ImageRegionIterator<TImage> targetR(this->image, reg);
	idx[0] = this->getWidth() - 2;
	reg.SetIndex(idx);
	typename itk::ImageRegionIterator<TImage> sideR(this->image, reg);

	outItR.GoToBegin();
	sideR.GoToBegin();
	for(targetR.GoToBegin(); !targetR.IsAtEnd(); ++targetR){
		typename TImage::PixelType r = targetR.Get() - sideR.Get();
		outItR.Set(r);
		++sideR;
		++outItR;
	}

	return out;
}

/**
 * Calculates the gradient of the image in the y direction.
 *
 * The 1D gradient is calculated as a two side difference, as if by passing the
 * transpose of the filter [0.5 0.0 -0.5] over the image. Image boundaries are handled by
 * taking a one sided difference.
 *
 * @return
 *   The gradient of this image in the y direction.
 */
template<class TImage>
typename TImage::Pointer ImageUtils<TImage>::simpleYGradient(){

	FloatUtils filter(3, 3, 0.0);
	itk::Index<2> idx;
	idx[0] = 1; idx[1] = 0;
	filter.getImage()->SetPixel(idx, 0.5);
	idx[0] = 1; idx[1] = 2;
	filter.getImage()->SetPixel(idx, -0.5);

	typename itk::FFTConvolutionImageFilter<TImage, FloatImageType, TImage>::Pointer conv =
			itk::FFTConvolutionImageFilter<TImage, FloatImageType, TImage>::New();
	//Set boundary condition to zero pad
	typename itk::ConstantBoundaryCondition<TImage> zeroPad;
	zeroPad.SetConstant(0);
	conv->SetBoundaryCondition(&zeroPad);

	conv->SetInput(this->image);
	conv->SetKernelImage(filter.getImage());
	conv->Update();

	typename TImage::Pointer out = conv->GetOutput();

	//Handle image boundaries
	idx.Fill(0);
	itk::Size<2> size;
	size[0] = this->getWidth(); size[1] = 1;
	typename TImage::RegionType reg(idx, size);

	typename itk::ImageRegionIterator<TImage> outItL(out, reg);
	typename itk::ImageRegionIterator<TImage> targetL(this->image, reg);
	idx[1] = 1;
	reg.SetIndex(idx);
	typename itk::ImageRegionIterator<TImage> sideL(this->image, reg);

	outItL.GoToBegin();
	sideL.GoToBegin();
	for(targetL.GoToBegin(); !targetL.IsAtEnd(); ++targetL){
		typename TImage::PixelType r = sideL.Get() - targetL.Get();
		outItL.Set(r);
		++sideL;
		++outItL;
	}

	idx[0] = 0; idx[1] = this->getHeight() - 1;
	reg.SetIndex(idx);

	typename itk::ImageRegionIterator<TImage> outItR(out, reg);
	typename itk::ImageRegionIterator<TImage> targetR(this->image, reg);
	idx[1] = this->getHeight() - 2;
	reg.SetIndex(idx);
	typename itk::ImageRegionIterator<TImage> sideR(this->image, reg);

	outItR.GoToBegin();
	sideR.GoToBegin();
	for(targetR.GoToBegin(); !targetR.IsAtEnd(); ++targetR){
		typename TImage::PixelType r = targetR.Get() - sideR.Get();
		outItR.Set(r);
		++sideR;
		++outItR;
	}

	return out;
}

/**
 * Sums all pixels in the image.
 *
 * @return
 *   The sum of all pixels in the image.
 */
template<class TImage>
double ImageUtils<TImage>::sumPixels(){

	double sum = 0;

	typename itk::ImageRegionIterator<TImage> it(this->image, this->image->GetLargestPossibleRegion());

	for(it.GoToBegin(); !it.IsAtEnd(); ++it){
		sum += (double)(it.Get());
	}

	return sum;
}

/**
 * Multiply each pixel in the image by the given scale factor.
 *
 * @param s
 *   Value by which to multiply each pixel
 * @return
 *   Image after performing multiplication
 */
template<class TImage>
typename TImage::Pointer ImageUtils<TImage>::mult(double s){

	typename itk::MultiplyImageFilter<TImage, TImage, TImage>::Pointer mlt =
			itk::MultiplyImageFilter<TImage, TImage, TImage>::New();
	mlt->SetInput(this->image);
	mlt->SetConstant(s);
	mlt->Update();

	return mlt->GetOutput();
}

/**
 * Perform pixel-by-pixel multiplication.
 *
 * @param s
 *   Image by which to multiply this image
 * @return
 *   Image after performing multiplication
 */
template<class TImage>
typename TImage::Pointer ImageUtils<TImage>::mult(typename TImage::Pointer s){

	typename itk::MultiplyImageFilter<TImage, TImage, TImage>::Pointer mlt =
			itk::MultiplyImageFilter<TImage, TImage, TImage>::New();
	mlt->SetInput1(this->image);
	mlt->SetInput2(s);
	mlt->Update();

	return mlt->GetOutput();
}

/**
 * Perform pixel-by-pixel division.
 *
 * When the divisor is zero,  the result is set to the maximum number that can
 * be represented by default to avoid exception.
 *
 * @param s
 *   Image by which to divide this image
 * @return
 *   Image after performing division.
 */
template<class TImage>
typename TImage::Pointer ImageUtils<TImage>::divide(typename TImage::Pointer s){

	typename itk::DivideOrZeroOutImageFilter<TImage, TImage, TImage>::Pointer mlt =
			itk::DivideOrZeroOutImageFilter<TImage, TImage, TImage>::New();
	mlt->SetInput1(this->image);
	mlt->SetInput2(s);
	mlt->Update();

	return mlt->GetOutput();
}

/**
 * Add given value to each pixel in the image.
 *
 * @param s
 *   Value to add to each pixel
 * @return
 *   Image after performing addition
 */
template<class TImage>
typename TImage::Pointer ImageUtils<TImage>::add(double s){

	typename itk::AddImageFilter<TImage, TImage, TImage>::Pointer mlt =
			itk::AddImageFilter<TImage, TImage, TImage>::New();
	mlt->SetInput1(this->image);
	mlt->SetConstant2(s);
	mlt->Update();

	return mlt->GetOutput();
}

/**
 * Performs element-by-element addition between this image and s.
 *
 * @param s
 *   Image to add to this
 * @return
 *   Image after performing addition
 */
template<class TImage>
typename TImage::Pointer ImageUtils<TImage>::add(typename TImage::Pointer s){

	this->image->SetOrigin(s->GetOrigin());
	this->image->SetSpacing(s->GetSpacing());

	typename itk::AddImageFilter<TImage, TImage, TImage>::Pointer mlt =
			itk::AddImageFilter<TImage, TImage, TImage>::New();
	mlt->SetInput1(this->image);
	mlt->SetInput2(s);
	mlt->Update();

	return mlt->GetOutput();
}

/**
 * Computes the square root of each pixel value.
 *
 * @return
 *   Image after applying a square root function to each pixel.
 */
template<class TImage>
typename TImage::Pointer ImageUtils<TImage>::sqrt(){

	typename itk::SqrtImageFilter<TImage, TImage>::Pointer sq =
			itk::SqrtImageFilter<TImage, TImage>::New();
	sq->SetInput(this->image);
	sq->Update();

	return sq->GetOutput();
}

/**
 * Rounds the value of each pixel in this image.
 *
 * @return
 *   Image after performing the rounding
 */
template<class TImage>
typename TImage::Pointer ImageUtils<TImage>::round(){

	typename itk::RoundImageFilter<TImage, TImage>::Pointer round =
			itk::RoundImageFilter<TImage, TImage>::New();
	round->SetInput(this->image);
	round->Update();

	return round->GetOutput();
}

/**
 * Finds maximum value within this image.
 *
 * @return
 *   The maximum pixel value
 */
template<class TImage>
typename TImage::PixelType ImageUtils<TImage>::max(){

	typename itk::MinimumMaximumImageCalculator<TImage>::Pointer mx =
			itk::MinimumMaximumImageCalculator<TImage>::New();
	mx->SetImage(this->image);
	mx->ComputeMaximum();

	return mx->GetMaximum();
}

/**
 * Finds maximum value within this image.
 *
 * @return
 *   The location of the maximum pixel
 */
template<class TImage>
typename TImage::IndexType ImageUtils<TImage>::maxLocation(){

	typename itk::MinimumMaximumImageCalculator<TImage>::Pointer mn =
			itk::MinimumMaximumImageCalculator<TImage>::New();
	mn->SetImage(this->image);
	mn->ComputeMaximum();

	return mn->GetIndexOfMaximum();
}


/**
 * Finds minimum value within this image.
 *
 * @return
 *   The minimum pixel value
 */
template<class TImage>
typename TImage::PixelType ImageUtils<TImage>::min(){

	typename itk::MinimumMaximumImageCalculator<TImage>::Pointer mn =
			itk::MinimumMaximumImageCalculator<TImage>::New();
	mn->SetImage(this->image);
	mn->ComputeMinimum();

	return mn->GetMinimum();
}

/**
 * Finds minimum value within this image.
 *
 * @return
 *   The location of the minimum pixel
 */
template<class TImage>
typename TImage::IndexType ImageUtils<TImage>::minLocation(){

	typename itk::MinimumMaximumImageCalculator<TImage>::Pointer mn =
			itk::MinimumMaximumImageCalculator<TImage>::New();
	mn->SetImage(this->image);
	mn->ComputeMinimum();

	return mn->GetIndexOfMinimum();
}

/**
 * Computes the Otsu threshold value for this image.
 *
 * @return
 *   Pixel value corresponding to the Otsu threshold value
 */
template<class TImage>
typename TImage::PixelType ImageUtils<TImage>::otsuThreshold(){

	typename itk::OtsuThresholdImageFilter<TImage, TImage>::Pointer otsu =
			itk::OtsuThresholdImageFilter<TImage, TImage>::New();

	otsu->SetInput(this->image);
	otsu->Update();

	return otsu->GetThreshold();
}

/**
 * Performs bounds checking on the given index.
 *
 * @return
 *   True if idx is within the bounds of this image, false otherwise.
 */
template<class TImage>
bool ImageUtils<TImage>::isInBounds(typename TImage::IndexType idx){
	return (idx[0] >= 0 && idx[0] < this->getWidth()) && (idx[1] >= 0 && idx[1] < this->getHeight());
}

/**
 * Creates a random sample of n pixel locations within this image.
 *
 * @param n
 *   The number of locations to sample
 * @param idxs
 *   Output destination. Vector is cleared of contents and filled with n random indices on exit.
 */
template<class TImage>
void ImageUtils<TImage>::sampleIndices(int n, std::vector<typename TImage::IndexType> &idxs){

	int w = this->getWidth();
	int h = this->getHeight();
	srand (time(NULL));
	idxs.clear();
	idxs.reserve(n);

	for(int i = 0; i < n; ++i){
		typename TImage::IndexType idx;
//		int xy = rand() % (w*h);
//		idx[0] = xy % w;
//		idx[1] = xy / h;

		idx[0] = rand() % w;
		idx[1] = rand() % h;

		idxs.push_back(idx);
	}
}

/**
 * Creates a random sample of n non-masked pixel locations within this image. The locations
 * are sampled only from pixels with mask values != 0.
 *
 * @param n
 *   The number of locations to sample
 * @param idxs
 *   Output destination. Vector is cleared of contents and filled with n random indices on exit.
 * @param mask
 *   Mask to direct the sampling. Only locations w/ non-zero mask values are returned.
 */
template<class TImage>
void ImageUtils<TImage>::sampleIndices(int n, std::vector<typename TImage::IndexType> &idxs, ByteImageType::Pointer mask){

	int w = this->getWidth();
	int h = this->getHeight();
	srand (time(NULL));
	idxs.clear();
	idxs.reserve(n);

	for(int i = 0; i < n; ){
		typename TImage::IndexType idx;
//		int xy = rand() % (w*h);
//		idx[0] = xy % w;
//		idx[1] = xy / h;

		idx[0] = rand() % w;
		idx[1] = rand() % h;

		if(mask->GetPixel(idx) != 0){
			idxs.push_back(idx);
			++i;
		}
	}
}

/**
 * Sets the value of the indicated locations to the given value.
 *
 * Optionally, creates a distinctive mark around the indicated locations rather than set just the
 * pixel value at the location.
 *
 * @param idxs
 *   Indices of the locations to mark
 * @param value
 *   The value to set in the indicated locations
 * @param singlePoint
 *   If true, sets only the indicated location to the given value. If false, creates a distinctive
 *   marker centered at each of the locations.
 * @return
 *   A new image with the given locations marked
 */
template<class TImage>
typename TImage::Pointer ImageUtils<TImage>::markLocations(std::vector<typename TImage::IndexType> &idxs,
			typename TImage::PixelType value, bool singlePoint){

	typename itk::ImageDuplicator<TImage>::Pointer dupe = itk::ImageDuplicator<TImage>::New();
	dupe->SetInputImage(this->image);
	dupe->Update();
	typename TImage::Pointer result = dupe->GetOutput();

	for(int i = 0; i < idxs.size(); ++i){
		typename TImage::IndexType idx = idxs[i];

		result->SetPixel(idx, value);
		if(!singlePoint){

			typename TImage::OffsetType off;
			off[0] = -2; off[1] = 0;
			if(this->isInBounds(idx + off)){
				result->SetPixel(idx + off, value);
			}

			off[0] = -1; off[1] = 0;
			if(this->isInBounds(idx + off)){
				result->SetPixel(idx + off, value);
			}

			off[0] = 1; off[1] = 0;
			if(this->isInBounds(idx + off)){
				result->SetPixel(idx + off, value);
			}

			off[0] = 2; off[1] = 0;
			if(this->isInBounds(idx + off)){
				result->SetPixel(idx + off, value);
			}

			off[0] = 0; off[1] = -2;
			if(this->isInBounds(idx + off)){
				result->SetPixel(idx + off, value);
			}

			off[0] = 0; off[1] = -1;
			if(this->isInBounds(idx + off)){
				result->SetPixel(idx + off, value);
			}

			off[0] = 0; off[1] = 1;
			if(this->isInBounds(idx + off)){
				result->SetPixel(idx + off, value);
			}

			off[0] = 0; off[1] = 2;
			if(this->isInBounds(idx + off)){
				result->SetPixel(idx + off, value);
			}

		}
	}

	return result;
}


/**
 *
 * Performs no bounds checking and doesn't enforce overlap - if s doesn't overlap this image using
 * the given location, the output is identical to this image.
 *
 * @param l
 *   Image containing
 * @param loc
 *   Position (in this image) where copying region begins
 * @return
 *   Image with data from the specified image replacing the original content
 */
template<class TImage>
typename TImage::Pointer ImageUtils<TImage>::copyImageAt(typename TImage::Pointer s, typename TImage::IndexType loc){

	typename TImage::RegionType region;
	region.SetSize(s->GetLargestPossibleRegion().GetSize());
	region.SetIndex(loc);

	typename TImage::Pointer out = this->copyImage();

	itk::ImageRegionIterator<TImage> from(s, s->GetLargestPossibleRegion());
	itk::ImageRegionIterator<TImage> to(out, region);

	typename TImage::RegionType wholeImage = out->GetLargestPossibleRegion();
	typename TImage::IndexType idx;

	from.GoToBegin();
	for(to.GoToBegin(); !to.IsAtEnd(); ++to){

//		idx = to.GetIndex();

		if(wholeImage.IsInside(to.GetIndex())){
			to.Set(from.Get());
		}
		++from;
	}

	return out;
}

/**
 * Pads this image with extra pixels using mirroring across the input boundary to fill in added
 * pixels.
 *
 * Produces slightly different output than the ITK built-in class MirrorPadImageFilter.
 *
 * @param l
 *   Number of pixels to add to the left side of the image
 * @param r
 *   Number of pixels to add to the right side of the image
 * @param a
 *   Number of pixels to add above the image
 * @param b
 *   Number of pixels to add below the image
 * @return
 *   Image with the specified number of pixels added with values determined by mirroring
 */
template<class TImage>
typename TImage::Pointer ImageUtils<TImage>::mirrorPadImage(int l, int r, int a, int b){

	itk::Size<2> s = this->image->GetLargestPossibleRegion().GetSize();
	int w0 = s[0];
	int h0 = s[1];
	FloatUtils mirror(w0 + r + l, h0 + a + b, 0.0);

	//Copy original image into approriate region
	itk::Index<2> idx;
	idx[0] = l;
	idx[1] = a;
	typename TImage::RegionType reg(idx, s);

	itk::ImageRegionIterator<TImage> *mirrorIt =
			new itk::ImageRegionIterator<TImage>(mirror.getImage(), reg);
	itk::ImageLinearConstIteratorWithIndex<TImage> *it =
			new itk::ImageLinearConstIteratorWithIndex<TImage>(this->image, this->image->GetLargestPossibleRegion());

	it->SetDirection(0);
	it->GoToBegin();
	for(mirrorIt->GoToBegin(); !mirrorIt->IsAtEnd(); ++(*mirrorIt)){
		mirrorIt->Set(it->Get());
		++(*it);

		if(it->IsAtEndOfLine()){
			it->NextLine();
		}
	}

	delete mirrorIt;
	delete it;

	//Above region
	idx[0] = l;
	idx[1] = 0;
	s[0] = w0;
	s[1] = a;

	reg.SetIndex(idx);
	reg.SetSize(s);
	mirrorIt = new itk::ImageRegionIterator<TImage>(mirror.getImage(), reg);

	idx[0] = 0;
	idx[1] = 1;
	reg.SetIndex(idx);
	it = new itk::ImageLinearConstIteratorWithIndex<TImage>(this->image, reg);
	it->SetDirection(0);

	it->GoToReverseBegin();
	it->GoToBeginOfLine();
	for(mirrorIt->GoToBegin(); !mirrorIt->IsAtEnd(); ++(*mirrorIt)){
		mirrorIt->Set(it->Get());
		++(*it);
		if(it->IsAtEndOfLine()){
			it->PreviousLine();
			it->GoToBeginOfLine();
		}
	}

	delete mirrorIt;
	delete it;

	//Below region
	idx[0] = l;
	idx[1] = a + h0;
	s[0] = w0;
	s[1] = b;

	reg.SetIndex(idx);
	reg.SetSize(s);
	mirrorIt = new itk::ImageRegionIterator<TImage>(mirror.getImage(), reg);

	idx[0] = 0;
	idx[1] = h0 - b - 1;
	reg.SetIndex(idx);
	it = new itk::ImageLinearConstIteratorWithIndex<TImage>(this->image, reg);
	it->SetDirection(0);

	it->GoToReverseBegin();
	it->GoToBeginOfLine();
	for(mirrorIt->GoToBegin(); !mirrorIt->IsAtEnd(); ++(*mirrorIt)){
		mirrorIt->Set(it->Get());
		++(*it);
		if(it->IsAtEndOfLine()){
			it->PreviousLine();
			it->GoToBeginOfLine();
		}
	}

	delete mirrorIt;
	delete it;

	//Left region
	idx[0] = 0;
	idx[1] = 0;
	s[0] = l;
	s[1] = a + h0 + b;

	reg.SetIndex(idx);
	reg.SetSize(s);
	mirrorIt = new itk::ImageRegionIterator<TImage>(mirror.getImage(), reg);

	idx[0] = l + 1;
	idx[1] = 0;
	reg.SetIndex(idx);
	it = new itk::ImageLinearConstIteratorWithIndex<TImage>(mirror.getImage(), reg);
	it->SetDirection(0);

	int count = 0;
	it->GoToReverseBegin();
	it->GoToBeginOfLine();
	mirrorIt->GoToEnd();
	for(--(*mirrorIt); count < s[0]*s[1]; --(*mirrorIt)){
		mirrorIt->Set(it->Get());
		++(*it);
		if(it->IsAtEndOfLine()){
			it->PreviousLine();
			it->GoToBeginOfLine();
		}
		++count;
	}

	delete mirrorIt;
	delete it;

	//Right region
	idx[0] = l + w0;
	idx[1] = 0;
	s[0] = r;
	s[1] = a + h0 + b;

	reg.SetIndex(idx);
	reg.SetSize(s);
	mirrorIt = new itk::ImageRegionIterator<TImage>(mirror.getImage(), reg);

	idx[0] = l + w0 - r - 1;
	idx[1] = 0;
	reg.SetIndex(idx);
	it = new itk::ImageLinearConstIteratorWithIndex<TImage>(mirror.getImage(), reg);
	it->SetDirection(0);

	count = 0;
	it->GoToReverseBegin();
	it->GoToBeginOfLine();
	mirrorIt->GoToEnd();
	for(--(*mirrorIt); count < s[0]*s[1]; --(*mirrorIt)){
		mirrorIt->Set(it->Get());
		++(*it);
		if(it->IsAtEndOfLine()){
			it->PreviousLine();
			it->GoToBeginOfLine();
		}
		++count;
	}

	delete mirrorIt;
	delete it;

	return mirror.getImage();
}

/**
 * Create a composite image by combining the given image with this one.
 *
 * @param s
 *   The image to tile next to this one.
 * @param dir
 *   The direction in which the images are tiled. dir <= 0 -> s is placed to the right of this image.
 *   dir >0  -> s is placed below this image.
 * @return
 *   The composite image created by tiling the two images together
 */
template<class TImage>
typename TImage::Pointer ImageUtils<TImage>::stackImages(typename TImage::Pointer s, int dir){

	typename itk::TileImageFilter<TImage, TImage>::Pointer tiler =
			itk::TileImageFilter<TImage, TImage>::New();

	itk::FixedArray< unsigned int, 2> layout;
	layout[0] = (dir <= 0) ? 2 : 0;
	layout[1] = (dir <= 0) ? 0 : 2;

	tiler->SetLayout(layout);

	tiler->SetInput(0, this->image);
	tiler->SetInput(1, s);
	tiler->Update();

	return tiler->GetOutput();
}

/**
 * Draws a line (one pixel thick) from the start to the end index.
 *
 * Performs no bounds checking.
 *
 * @param start
 *   Starting point of the line
 * @param end
 *   Ending point of the line
 * @return
 *   Copy of this image with the specified line drawn
 */
template<class TImage>
typename TImage::Pointer ImageUtils<TImage>::drawLine(typename TImage::IndexType start, typename TImage::IndexType end,
			typename TImage::PixelType value){

	typename itk::ImageDuplicator<TImage>::Pointer dupe = itk::ImageDuplicator<TImage>::New();
	dupe->SetInputImage(this->image);
	dupe->Update();
	typename TImage::Pointer result = dupe->GetOutput();

	typename itk::LineIterator<TImage> it(result, start, end);

	for(it.GoToBegin(); !it.IsAtEnd(); ++it){
		it.Set(value);
	}

	return result;
}

/**
 * Write the image represented by this utility object to the given path.
 *
 * @param path
 *   Path to image destination file
 */
template<class TImage>
void ImageUtils<TImage>::writeToPath(std::string &path){
	typename itk::ImageFileWriter<TImage>::Pointer writer = itk::ImageFileWriter<TImage>::New();

	writer->SetInput(this->image);
	writer->SetFileName(path);
	writer->Update();
}

/**
 * Write the image represented by this utility object to the given path.
 *
 * @param path
 *   Path to image destination file
 */
template<class TImage>
void ImageUtils<TImage>::writeToPath(const char *path){
	std::string p(path);
	this->writeToPath(p);
}

/**
 * Prints the image to stdout.
 */
template<class TImage>
void ImageUtils<TImage>::printImage(){

	int w = this->getWidth();
	int h = this->getHeight();
	int count = 0;

	std::cout << std::endl;

	typename itk::ImageRegionIterator<TImage> it(this->getImage(), this->getImage()->GetLargestPossibleRegion());

	for(it.GoToBegin(); !it.IsAtEnd(); ++it){
		printf("%9.4lf  ", static_cast<double>(it.Get()));
		++count;
		if(count % w == 0){
			std::cout << std::endl;
		}
	}
	std::cout << std::endl;
}

/**
 * Create a square image with pixel values set equal to the values of a symmetric
 * 2D Gaussian function.
 *
 * The output image always has odd dimensions (2*w + 1).
 *
 * @param w
 *   The half width of the output image
 * @param sigma
 *   Standard deviation of the Gaussian used to generate the output
 * @return
 *   Image with float valued pixels set equal to 2D Gaussian function
 */
template<class TImage>
FloatImageType::Pointer ImageUtils<TImage>::createGaussianImage(int w, double sigma){

	double sum = 0;

	itk::Index<2> idx;
	idx.Fill(0);

	itk::Index<2> center;
	center.Fill(w);

	itk::Size<2> s;
	s.Fill(2*w + 1);

	FloatImageType::RegionType region;
	region.SetIndex(idx);
	region.SetSize(s);

	FloatImageType::Pointer kernel = FloatImageType::New();
	kernel->SetRegions(region);
	kernel->Allocate();
	kernel->FillBuffer(0);

	itk::ImageRegionIterator<FloatImageType> it(kernel, region);

	for(it.GoToBegin(); !it.IsAtEnd(); ++it){
		itk::Index<2> cur = it.GetIndex();
		cur[0] = cur[0] - center[0];
		cur[1] = cur[1] - center[1];

		double v = exp(-1.0*(cur[0]*cur[0] + cur[1]*cur[1])/(2*sigma*sigma));
		v /= 2*M_PI*sigma*sigma;
		it.Set(v);
		sum += v;
	}


	FloatUtils mult(kernel);
	kernel = mult.mult(1.0/sum);

	return kernel;
}

/**
 * Convert the image to a byte respresentation. Pixel values are linearly scaled to
 * the range [0, 255].
 *
 * @param f
 *   The float image to convert to byte
 */
template<class TImage>
ByteImageType::Pointer ImageUtils<TImage>::floatToByteImage(FloatImageType::Pointer f){
	itk::RescaleIntensityImageFilter<FloatImageType, FloatImageType>::Pointer rescale =
			itk::RescaleIntensityImageFilter<FloatImageType, FloatImageType>::New();
	rescale->SetInput(f);
	rescale->SetOutputMinimum(0);
	rescale->SetOutputMaximum(255);
	rescale->Update();

	itk::CastImageFilter<FloatImageType, ByteImageType>::Pointer cast =
			itk::CastImageFilter<FloatImageType, ByteImageType>::New();
	cast->SetInput(rescale->GetOutput());
	cast->Update();

	return cast->GetOutput();
}

/**
 * Convert the image to a byte respresentation. Pixel values are linearly scaled to
 * the range [0, 255].
 *
 * @param f
 *   The int image to convert to byte
 */
template<class TImage>
ByteImageType::Pointer ImageUtils<TImage>::intToByteImage(IntImageType::Pointer i){
	itk::RescaleIntensityImageFilter<IntImageType, IntImageType>::Pointer rescale =
			itk::RescaleIntensityImageFilter<IntImageType, IntImageType>::New();
	rescale->SetInput(i);
	rescale->SetOutputMinimum(0);
	rescale->SetOutputMaximum(255);
	rescale->Update();

	itk::CastImageFilter<IntImageType, ByteImageType>::Pointer cast =
			itk::CastImageFilter<IntImageType, ByteImageType>::New();
	cast->SetInput(rescale->GetOutput());
	cast->Update();

	return cast->GetOutput();
}

/**
 * Convert the byte image to a float representation.
 *
 * @param b
 *   The byte image to convert to a float image
 */
template<class TImage>
FloatImageType::Pointer ImageUtils<TImage>::byteToFloatImage(ByteImageType::Pointer b){

	itk::CastImageFilter<ByteImageType, FloatImageType>::Pointer cast =
			itk::CastImageFilter<ByteImageType, FloatImageType>::New();
	cast->SetInput(b);
	cast->Update();

	return cast->GetOutput();
}

/**
 * Convert the byte image to an int representation.
 *
 * @param b
 *   The byte image to convert to an int image
 */
template<class TImage>
IntImageType::Pointer ImageUtils<TImage>::byteToIntImage(ByteImageType::Pointer b){

	itk::CastImageFilter<ByteImageType, IntImageType>::Pointer cast =
			itk::CastImageFilter<ByteImageType, IntImageType>::New();
	cast->SetInput(b);
	cast->Update();

	return cast->GetOutput();
}

/**
 * Convert the int image to a float representation.
 *
 * @param b
 *   The int image to convert to a float image
 */
template<class TImage>
FloatImageType::Pointer ImageUtils<TImage>::intToFloatImage(IntImageType::Pointer i){

	itk::CastImageFilter<IntImageType, FloatImageType>::Pointer cast =
			itk::CastImageFilter<IntImageType, FloatImageType>::New();
	cast->SetInput(i);
	cast->Update();

	return cast->GetOutput();
}

/**
 * Convert the float image to an int representation.
 *
 * @param b
 *   The float image to convert to an int image
 */
template<class TImage>
IntImageType::Pointer ImageUtils<TImage>::floatToIntImage(FloatImageType::Pointer f){

	itk::CastImageFilter<FloatImageType, IntImageType>::Pointer cast =
			itk::CastImageFilter<FloatImageType, IntImageType>::New();
	cast->SetInput(f);
	cast->Update();

	return cast->GetOutput();
}

template<class TImage>
ImageUtils<TImage>::~ImageUtils() {
	// TODO Auto-generated destructor stub
}

#endif /* IMAGEUTILS_H_ */
