/*
 * CVMatUtils.h
 *
 *  Created on: Jul 16, 2013
 *      Author: mchristopher
 */

#include <iostream>
#include <fstream>
#include <assert.h>
#include <vector>
#include <opencv2/core/core.hpp>
#include "ImageUtils.h"
#include "ImageDataset.h"

#ifndef CVMATUTILS_H_
#define CVMATUTILS_H_

typedef cv::Mat Mat;

/**
 * Provides utility methods for calculating row- and column-wise statistics from
 * OpenCV Mat matrix representations.
 */
class CVMatUtils {
public:
	CVMatUtils(){}

	static void sumRows(Mat &m, Mat &r);
	static void sumColumns(Mat &m, Mat &r);
	static double rowSum(Mat &m, int r);
	static double columnSum(Mat &m, int c);

	static void avgRows(Mat &m, Mat &r);
	static void avgColumns(Mat &m, Mat &r);
	static double rowAvg(Mat &m, int r);
	static double columnAvg(Mat &m, int c);

	static void varRows(Mat &m, Mat &r);
	static void varColumns(Mat &m, Mat &r);
	static double rowVar(Mat &m, int r);
	static double columnVar(Mat &m, int c);

	static void standardizeCols(cv::Mat data, cv::Mat avgs, cv::Mat stds);

	static bool fitLinearModel(Mat &y, Mat &x, Mat &b, Mat &eps, Mat &stats);

	static void sortRowsWithKeys(Mat &matrix, Mat &keys, Mat &r, Mat &kr);

	static void removeRows(cv::Mat &data, cv::Mat &mask, cv::Mat &r);
	static void removeCols(cv::Mat &data, cv::Mat &mask, cv::Mat &r);

	static void getRow(cv::Mat &data, int row, cv::Mat &r);
	static void getCol(cv::Mat &data, int col, cv::Mat &r);

	static int sampleRows(Mat &m, Mat &rows, Mat &result);

	static void appendRows(cv::Mat &base, cv::Mat &add, cv::Mat &r);
	static void appendCols(cv::Mat &base, cv::Mat &add, cv::Mat &r);

	static void matToImage(Mat &src, FloatImageType::Pointer &dest);
	static void imageToMat(FloatImageType::Pointer &src, Mat &dest);

	static void imageToVector(FloatImageType::Pointer &src, Mat &dest);
	static void vectorToImage(Mat &src, int w, int h, FloatImageType::Pointer &dest);

	static void imageToVector(FloatImageType::Pointer &src, ByteImageType::Pointer &mask, Mat &dest);
	static void vectorToImage(Mat &src, ByteImageType::Pointer &mask, int w, int h, FloatImageType::Pointer &dest);

	static int imagesToMat(ImageDataset images, Mat &mat);
	static int imagesToMat(ImageDataset images, ByteImageType::Pointer &mask, Mat &mat);

	static bool writeCSV(Mat &m, const char *path);
	static bool writeCSV(Mat &m, const char *path, std::vector<std::string> &colnames);
	static void readCSV(Mat &m, const char *path);
	static void readCSV(Mat &m, const char *path, std::vector<std::string> &colnames);

	static void printMatrix(Mat &m, std::string delim = ",");

	virtual ~CVMatUtils(){}
};

#endif /* CVMATUTILS_H_ */
