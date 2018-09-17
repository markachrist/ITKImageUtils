/*
 * Performance.cpp
 *
 *  Created on: Jul 22, 2013
 *      Author: mchristopher
 */

#include "Performance.h"

Performance::Performance(Mat &t, Mat &p) {
	this->truth = t;
	this->probs = p;
	this->aucEvalPoints = DEF_EVAL_POINTS;
}

double Performance::sensitivity(double t){

	double tv, pv;
	int tp = 0;
	int fn = 0;

	for(int i = 0; i < this->truth.rows; ++i){

		tv = truth.at<int>(i);
		pv = probs.at<float>(i);

		if(tv > 0){
			if(pv >= t){
				++tp;
			}
			else{
				++fn;
			}
		}
	}
	return 1.0*tp / (tp + fn);
}


double Performance::specificity(double t){
	double tv, pv;
	int tn = 0;
	int fp = 0;

	for(int i = 0; i < this->truth.rows; ++i){

		tv = truth.at<int>(i);
		pv = probs.at<float>(i);

		if(tv <= 0){
			if(pv >= t){
				++fp;
			}
			else{
				++tn;
			}
		}
	}
	return 1.0*tn / (tn + fp);
}

void Performance::roc(Mat &sens, Mat &fpr){

	double min, max;
	double eps = 1e-6;
	int points = this->getAUCEvaluationPoints();
	Mat ts(points, 1, CV_64FC1);

	sens.create(points + 1, 1, CV_64FC1);
	fpr.create(points + 1, 1, CV_64FC1);

	minMaxLoc(probs, &min, &max);
	ts.at<double>(0) = min - eps;
	for(int i = 1; i < points; ++i){
		ts.at<double>(i) = min + ((max - min)/(points - 1))*(i);
	}

	for(int i = 0; i < points; ++i){
		double t = ts.at<double>(i);
		sens.at<double>(i) = this->sensitivity(t);
		fpr.at<double>(i) = 1.0 - this->specificity(t);
	}

	sens.at<double>(points) = 0.0;
	fpr.at<double>(points) = 0.0;
}

double Performance::auc(){

	double ca, area = 0.0;
	Mat sens;
	Mat fpr;
	roc(sens, fpr);

	for(int i = 0; i < sens.rows - 1; ++i){

		double x0 = fpr.at<double>(i);
		double x1 = fpr.at<double>(i + 1);
		double y0 = sens.at<double>(i);
		double y1 = sens.at<double>(i + 1);

		ca = (x0 - x1)*y1 + (x0 - x1)*(y0 - y1)/2;
		area += ca;
	}

	return area;
}


Performance::~Performance() {
	// TODO Auto-generated destructor stub
}
