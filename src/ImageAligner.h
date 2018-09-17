/*
 * ImageAligner.h
 *
 *  Created on: May 6, 2014
 *      Author: mchristopher
 */

#ifndef IMAGEALIGNER_H_
#define IMAGEALIGNER_H_

#include<vector>
#include<math.h>

#include "itkTranslationTransform.h"
#include "itkSimilarity2DTransform.h"
#include "itkRigid2DTransform.h"
#include "itkEuler2DTransform.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkSingleValuedCostFunction.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkNormalizedCorrelationImageToImageMetric.h"
#include "itkMutualInformationImageToImageMetric.h"
#include "itkOptimizer.h"
#include "itkSingleValuedNonLinearOptimizer.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkConjugateGradientOptimizer.h"
#include "itkLBFGSOptimizer.h"
#include "itkImageRegistrationMethod.h"

#include "ImageUtils.h"

/**
 * Flags to set one of several types of non-deformable transformation types
 * for the image aligner.
 *
 */
enum AlignTransform {
	Translation = 0,	//Shift in x,y
	Rigid,				//Shift and rotation
	Similarity,			//Rigid and homogenous scale
	Affine				//Rigid and non-homogeneous scale (skew)
};

/**
 *
 */
enum AlignMetric {
	MeanSquares = 0,
	NormalizedCorrelation,
	MutualInformation
};

/**
 *
 */
enum AlignOptimizer {
	GradientDescent = 0,
	ConjugateGradient,
	LBFGS
};

/**
 *
 */
template<class TImage>
class ImageAligner {

public:

	typedef itk::Transform<double, 2, 2> TransformType;
	typedef itk::SingleValuedCostFunction MetricType;
	typedef itk::SingleValuedNonLinearOptimizer OptimizerType;

	ImageAligner(typename TImage::Pointer &i);
	ImageAligner(AlignTransform t, AlignMetric m, AlignOptimizer o, typename TImage::Pointer &i);

	TransformType::Pointer alignToTarget(typename TImage::Pointer &query);
	TransformType::Pointer alignToTarget(typename TImage::Pointer &query, itk::Array<double> &params);
	TransformType::Pointer alignToTargetMultiscale(typename TImage::Pointer &query);
	TransformType::Pointer alignToTargetMultiscale(typename TImage::Pointer &query, itk::Array<double> &params);

	void setTransformType(AlignTransform t);
	void setMetricType(AlignMetric t);
	void setOptimizerType(AlignOptimizer t);

	void setTransformParams(itk::Array<double> &params);
	void createIdentityParams(itk::Array<double> &params);

	void setMetricParameters(std::vector<double> &params);
	void createMetricParameters(std::vector<double> &params);

	virtual ~ImageAligner(){};

	//Getters/Setters
	inline int getMinimumSize(){return this->minSize;}
	inline void setMinimumSize(int s){this->minSize = s;}
	inline AlignTransform getTransformType(){return this->ttype;}
	inline AlignMetric getMetricType(){return this->mtype;}
	inline AlignOptimizer getOptimizerType(){return this->otype;}
	inline OptimizerType::Pointer &getOptimizer(){return this->opt;}

private:

	/** */
	AlignTransform ttype;
	TransformType::Pointer trans;

	/** */
	AlignMetric mtype;
	MetricType::Pointer metric;
	std::vector<double> metricParams;

	/** */
	AlignOptimizer otype;
	OptimizerType::Pointer opt;

	/** */
	ImageUtils<TImage> target;

	/** */
	int minSize;

	TransformType::Pointer createTransform();
	typename MetricType::Pointer createMetric();
	OptimizerType::Pointer createOptimizer();

	void init(AlignTransform t, AlignMetric m, AlignOptimizer o, typename TImage::Pointer &i);
	void setToIdentity(TransformType::Pointer &tform);
	void createScales(std::vector<int> &scales);
	void scaleTransform(TransformType::Pointer &tform, double s);

	TransformType::Pointer alignHelper(typename TImage::Pointer &fixed, typename TImage::Pointer &moving,
			TransformType::Pointer &tform);
};

//Typedefs again. Not the cleanest way to do this...
typedef itk::Transform<double, 2, 2> TransformType;
typedef itk::SingleValuedCostFunction MetricType;
typedef itk::SingleValuedNonLinearOptimizer OptimizerType;

/**
 *
 */
template<class TImage>
ImageAligner<TImage>::ImageAligner(typename TImage::Pointer &i){
	this->init(Affine, MutualInformation, ConjugateGradient, i);
}

/**
 *
 */
template<class TImage>
ImageAligner<TImage>::ImageAligner(AlignTransform t, AlignMetric m, AlignOptimizer o, typename TImage::Pointer &i){
	this->init(t, m, o, i);
}

/**
 *
 */
template<class TImage>
void ImageAligner<TImage>::init(AlignTransform t, AlignMetric m, AlignOptimizer o, typename TImage::Pointer &i){

	this->setTransformType(t);
	this->setMetricType(m);
	this->setOptimizerType(o);

	this->target.setImage(i);

	this->setMinimumSize(16);
}


/**
 *
 */
template<class TImage>
TransformType::Pointer ImageAligner<TImage>::alignToTarget(typename TImage::Pointer &query){
	itk::Array<double> params;
	this->createIdentityParams(params);
	return this->alignToTarget(query, params);
}

/**
 *
 */
template<class TImage>
TransformType::Pointer ImageAligner<TImage>::alignToTargetMultiscale(typename TImage::Pointer &query){
	itk::Array<double> params;
	this->createIdentityParams(params);
	return this->alignToTargetMultiscale(query, params);
}

/**
 *
 */
template<class TImage>
TransformType::Pointer ImageAligner<TImage>::alignToTarget(typename TImage::Pointer &query, itk::Array<double> &params){
	typename TImage::Pointer fixed = this->target.getImage();
	this->setTransformParams(params);
	return this->alignHelper(fixed, query, this->trans);
}

/**
 *
 */
template<class TImage>
TransformType::Pointer ImageAligner<TImage>::alignToTargetMultiscale(typename TImage::Pointer &query, itk::Array<double> &params){

	int fw, fh, mw, mh, s;
	std::vector<int> scales;
	typename TImage::Pointer fixed, moving;
	ImageUtils<TImage> queryUtils(query);

	this->createScales(scales);
	this->setTransformParams(params);

	for(int i = 0; i < scales.size(); ++i){

		if(i > 0){
			this->scaleTransform(this->trans, 2);
		}

		s = scales[i];

		fw = this->target.getWidth() / s;
		fh = this->target.getHeight() / s;
		mw = queryUtils.getWidth() / s;
		mh = queryUtils.getHeight() / s;

		fixed = this->target.resize(fw, fh);
		moving = queryUtils.resize(mw, mh);

		this->trans = this->alignHelper(fixed, moving, this->trans);
	}

	return this->trans;
}

/**
 *
 */
template<class TImage>
TransformType::Pointer ImageAligner<TImage>::alignHelper(typename TImage::Pointer &fixed, typename TImage::Pointer &moving,
		TransformType::Pointer &tform){

	//Determine number of samples using metric params
	int numSamples = fixed->GetLargestPossibleRegion().GetSize()[0] * fixed->GetLargestPossibleRegion().GetSize()[1] * this->metricParams[0];

	//set up registraion methd and interpolator
	typename itk::ImageRegistrationMethod<TImage, TImage>::Pointer reg =
			itk::ImageRegistrationMethod<TImage, TImage>::New();
	typename itk::BSplineInterpolateImageFunction<TImage, double>::Pointer interp =
			itk::BSplineInterpolateImageFunction<TImage, double>::New();

	TransformType::ParametersType transParams = tform->GetParameters();

	if(this->mtype == MeanSquares){
		typename itk::MeanSquaresImageToImageMetric<TImage, TImage>::Pointer tempMetric =
				itk::MeanSquaresImageToImageMetric<TImage, TImage>::New();

		tempMetric->SetNumberOfSpatialSamples(numSamples);

		reg->SetMetric(tempMetric);
	}
	else if(this->mtype == NormalizedCorrelation){
		typename itk::NormalizedCorrelationImageToImageMetric<TImage, TImage>::Pointer tempMetric =
				itk::NormalizedCorrelationImageToImageMetric<TImage, TImage>::New();

		tempMetric->SetNumberOfSpatialSamples(numSamples);

		reg->SetMetric(tempMetric);
	}
	else if(this->mtype == MutualInformation){
		typename itk::MutualInformationImageToImageMetric<TImage, TImage>::Pointer tempMetric =
				itk::MutualInformationImageToImageMetric<TImage, TImage>::New();

		tempMetric->SetNumberOfSpatialSamples(numSamples);
		tempMetric->SetFixedImageStandardDeviation(this->metricParams[1]);
		tempMetric->SetMovingImageStandardDeviation(this->metricParams[2]);

		reg->SetMetric(tempMetric);
	}

	reg->SetOptimizer(this->opt);
	reg->SetTransform(tform);
	reg->SetInterpolator(interp);

	reg->SetFixedImage(fixed);
	reg->SetFixedImageRegion(fixed->GetLargestPossibleRegion());
	reg->SetMovingImage(moving);
	reg->SetInitialTransformParameters(transParams);

	reg->Update();

	transParams = reg->GetLastTransformParameters();
	tform->SetParameters(transParams);

	return tform;
}


/**
 *
 */
template<class TImage>
void ImageAligner<TImage>::createScales(std::vector<int> &scales){

	int n, m, w, h, size, nums;

	w = this->target.getWidth();
	h = this->target.getHeight();
	size = (w < h) ? w : h;

	scales.clear();

	n = floor(log2(size));
	m = floor(log2(this->minSize));

	nums = (int)(ceil(log2(size) - log2(minSize))) - 1;

	for(int i = nums; i >= 0; --i){
		scales.push_back(pow(2, i));
	}
}

/**
 *
 */
template<class TImage>
void ImageAligner<TImage>::scaleTransform(TransformType::Pointer &tform, double s){

	TransformType::ParametersType params = this->trans->GetParameters();

	if(this->ttype == Translation){
		//x,y shifting
		params[0] = s*params[0];
		params[1] = s*params[1];
	}
	else if(this->ttype == Rigid){
		//x,y shifting
		params[1] = s*params[1];
		params[2] = s*params[2];
	}
	else if(this->ttype == Similarity){
		//x,y shifting
		params[2] = s*params[2];
		params[3] = s*params[3];
	}
	else if(this->ttype == Affine){
		//x,y shifting
		params[4] = s*params[4];
		params[5] = s*params[5];
	}

	this->trans->SetParameters(params);
}

/**
 *
 */
template<class TImage>
TransformType::Pointer ImageAligner<TImage>::createTransform(){

	TransformType::Pointer trans = NULL;

	if(this->ttype == Translation){
		trans = itk::TranslationTransform<double, 2>::New();
	}
	else if(this->ttype == Rigid){
		trans = itk::Euler2DTransform<double>::New();
	}
	else if(this->ttype == Similarity){
		trans = itk::Similarity2DTransform<double>::New();
	}
	else if(this->ttype == Affine){
		trans = itk::AffineTransform<double, 2>::New();
	}

	itk::Array<double> params;
	this->createIdentityParams(params);
	TransformType::ParametersType ps(trans->GetNumberOfParameters());

	for(int i = 0; i < trans->GetNumberOfParameters(); ++i){
		ps[i] = params[i];
	}
	trans->SetParameters(ps);

	return trans;
}

/**
 * Sets the parameters that are sent to the image metric during registration. The meaning of these parameters depends on the type of
 * metric in use. See the following descriptions for the number/type of parameters expected in each case.
 *
 * MeanSquares
 *   Param 1 = proportion of pixels in the fixed image to use as sample points (see metric->SetNumberOfSpatialSamples()).
 *     A value of 1.0 should be passed to indicate all pixels, 0.5 uses half, etc. Undefined results for values outside (0, 1].
 *
 * NormalizedCorrelation
 *   Param 1 = proportion of pixels in the fixed image to use as sample points (see metric->SetNumberOfSpatialSamples()).
 *     A value of 1.0 should be passed to indicate all pixels, 0.5 uses half, etc. Undefined results for values outside (0, 1].
 *
 * MutualInformation
 *   Param 1 = proportion of pixels in the fixed image to use as sample points (see metric->SetNumberOfSpatialSamples()).
 *     A value of 1.0 should be passed to indicate all pixels, 0.5 uses half, etc. Undefined results for values outside (0, 1].
 *   Param 2 = standard deviation of the fixed image (as passed to metric->SetFixedImageStandardDeviation())
 *   Param 3 = standard deviation of the fixed image (as passed to metric->SetFixedImageStandardDeviation())
 *
 */
template<class TImage>
void ImageAligner<TImage>::setMetricParameters(std::vector<double> &params){
	this->metricParams = params;
}

/**
 *
 */
template<class TImage>
void ImageAligner<TImage>::createMetricParameters(std::vector<double> &params){

	params.clear();

	if(this->mtype == MeanSquares){
		params.push_back(1.0);
	}
	else if(this->mtype == NormalizedCorrelation){
		params.push_back(1.0);
	}
	else if(this->mtype == MutualInformation){
		params.push_back(1.0);
		params.push_back(0.4);
		params.push_back(0.4);
	}
}

/**
 *
 */
template<class TImage>
OptimizerType::Pointer ImageAligner<TImage>::createOptimizer(){

	OptimizerType::Pointer opt = NULL;

	if(this->otype == GradientDescent){
		itk::RegularStepGradientDescentOptimizer::Pointer gopt = itk::RegularStepGradientDescentOptimizer::New();
		gopt->SetMaximumStepLength(0.1);
		gopt->SetMinimumStepLength(0.01);
		gopt->SetNumberOfIterations(100);
		opt = gopt;
	}
	else if(this->otype == ConjugateGradient){
		opt = itk::ConjugateGradientOptimizer::New();
	}
	else if(this->otype == LBFGS){
		opt = itk::LBFGSOptimizer::New();
	}

	return opt;
}

/**
 *
 */
template<class TImage>
void ImageAligner<TImage>::setTransformParams(itk::Array<double> &params){

	TransformType::ParametersType ps(this->trans->GetNumberOfParameters());

	for(int i = 0; i < ps.GetSize(); ++i){
		ps[i] = params[i];
	}
	this->trans->SetParameters(ps);
}

/**
 *
 */
template<class TImage>
void ImageAligner<TImage>::createIdentityParams(itk::Array<double> &ps){

	if(this->ttype == Translation){
		itk::Array<double> params(2);
		//x,y shifting
		params[0] = 0;
		params[1] = 0;

		ps = params;
	}
	else if(this->ttype == Rigid){
		itk::Array<double> params(3);
		//Rotation
		params[0] = 0;
		//x,y shifting
		params[1] = 0;
		params[2] = 0;

		ps = params;
	}
	else if(this->ttype == Similarity){
		itk::Array<double> params(4);
		//Scale
		params[0] = 1.0;
		//Rotation
		params[1] = 0;
		//x,y shifting
		params[2] = 0;
		params[3] = 0;

		ps = params;
	}
	else if(this->ttype == Affine){
		itk::Array<double> params(6);
		//Scale & rotation
		params[0] = 1.0;
		params[1] = 0;
		params[2] = 0;
		params[3] = 1.0;
		//x,y shifting
		params[4] = 0;
		params[5] = 0;

		ps = params;
	}
}

/**
 *
 */
template<class TImage>
void ImageAligner<TImage>::setTransformType(AlignTransform t){
	this->ttype = t;
	this->trans = this->createTransform();
}

/**
 *
 */
template<class TImage>
void ImageAligner<TImage>::setMetricType(AlignMetric t){
	this->mtype = t;
	this->createMetricParameters(this->metricParams);
}

/**
 *
 */
template<class TImage>
void ImageAligner<TImage>::setOptimizerType(AlignOptimizer t){
	this->otype = t;
	this->opt = this->createOptimizer();
}

#endif /* IMAGEALIGNER_H_ */
