/*
 * Performance.h
 *
 *  Created on: Jul 22, 2013
 *      Author: mchristopher
 */

#include "opencv2/core/core.hpp"
#include "opencv2/ml/ml.hpp"

#ifndef PERFORMANCE_H_
#define PERFORMANCE_H_

typedef cv::Mat Mat;

#define DEF_EVAL_POINTS 100

class Performance {


	/** True class values for each point in the evaluated data set */
	Mat truth;

	/** Set of probablities of each point belonging to class 1 - Need not be bound to [0, 1] */
	Mat probs;

	/** Number of points used to construct ROC curve */
	int aucEvalPoints;

public:
	Performance(Mat &t, Mat &p);

	double sensitivity(double t);
	double specificity(double t);

	void roc(Mat &sens, Mat &fpr);
	double auc();

	virtual ~Performance();

	inline void setAUCEvaluationPoints(int n){
		this->aucEvalPoints = (n >= 2) ? n : 2;
	}
	inline int getAUCEvaluationPoints(){
		return this->aucEvalPoints;
	}

};

#endif /* PERFORMANCE_H_ */
