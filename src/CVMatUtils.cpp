/*
 * CVMatUtils.cpp
 *
 *  Created on: Jul 16, 2013
 *      Author: mchristopher
 */

#include "CVMatUtils.h"

/**
 * Sums along all rows in matrix m. Input is assumed to be a 2D matrix with a
 * scalar data type. On exit, r is a column vector with length equal to the
 * number of rows in m. Output data is of type double.
 *
 * @param m
 *   Input matrix
 * @param r
 *   Output destination. Overwritten on exit.
 */
void CVMatUtils::sumRows(Mat &m, Mat &r){

	int rows = m.rows;

	r.create(rows, 1, CV_64FC1);

	for(int i = 0; i < rows; ++i){
		r.at<double>(i) = rowSum(m, i);
	}
}

/**
 * Sums along all cols in matrix m. Input is assumed to be a 2D matrix with a
 * scalar data type. On exit, r is a row vector with length equal to the number
 * of cols in m. Output data is of type double.
 *
 * @param m
 *   Input matrix
 * @param r
 *   Output destination. Overwritten on exit.
 */
void CVMatUtils::sumColumns(Mat &m, Mat &r){

	int cols = m.cols;

	r.create(1, cols, CV_64FC1);

	for(int j = 0; j < cols; ++j){
		r.at<double>(j) = columnSum(m, j);
	}
}

/**
 * Sums along a single row within a matrix. Matrix is assumed to be 2D with
 * a scalar data type. Regardless of input data type depth, output is a double.
 *
 * @param m
 *   Input matrix
 * @param r
 *   Row of m to sum
 * @return
 *   Sum of the values in row r of m
 */
double CVMatUtils::rowSum(Mat &m, int r){

	int cols = m.cols;
	double sum = 0;

	int d = m.depth();

	switch(d){
		case CV_8U:
			for(int j = 0; j < cols; ++j){
				sum += m.at<unsigned char>(r, j);
			}
			break;
		case CV_8S:
			for(int j = 0; j < cols; ++j){
				sum += m.at<char>(r, j);
			}
			break;
		case CV_16U:
			for(int j = 0; j < cols; ++j){
				sum += m.at<unsigned short>(r, j);
			}
			break;
		case CV_16S:
			for(int j = 0; j < cols; ++j){
				sum += m.at<short>(r, j);
			}
			break;
		case CV_32S:
			for(int j = 0; j < cols; ++j){
				sum += m.at<int>(r, j);
			}
			break;
		case CV_32F:
			for(int j = 0; j < cols; ++j){
				sum += m.at<float>(r, j);
			}
			break;
		case CV_64F:
			for(int j = 0; j < cols; ++j){
				sum += m.at<double>(r, j);
			}
			break;
	}

	return sum;
}

/**
 * Sums along a single column within a matrix. Matrix is assumed to be 2D with
 * a scalar data type. Regardless of input data type depth, output is a double.
 *
 * @param m
 *   Input matrix
 * @param c
 *   Column of m to sum
 * @return
 *   Sum of the values in col c of m
 */
double CVMatUtils::columnSum(Mat &m, int c){

	int rows = m.rows;
	double sum = 0;

	int d = m.depth();

	switch(d){
		case CV_8U:
			for(int i = 0; i < rows; ++i){
				sum += m.at<unsigned char>(i, c);
			}
			break;
		case CV_8S:
			for(int i = 0; i < rows; ++i){
				sum += m.at<char>(i, c);
			}
			break;
		case CV_16U:
			for(int i = 0; i < rows; ++i){
				sum += m.at<unsigned short>(i, c);
			}
			break;
		case CV_16S:
			for(int i = 0; i < rows; ++i){
				sum += m.at<short>(i, c);
			}
			break;
		case CV_32S:
			for(int i = 0; i < rows; ++i){
				sum += m.at<int>(i, c);
			}
			break;
		case CV_32F:
			for(int i = 0; i < rows; ++i){
				sum += m.at<float>(i, c);
			}
			break;
		case CV_64F:
			for(int i = 0; i < rows; ++i){
				sum += m.at<double>(i, c);
			}
			break;
	}

	return sum;
}

/**
 * Averages along all rows in matrix m. Input is assumed to be a 2D matrix with
 * a scalar data type. On exit, r is a column vector with length equal to the
 * number of rows in m. Output data is of type double.
 *
 * @param m
 *   Input matrix
 * @param r
 *   Output destination. Overwritten on exit.
 */
void CVMatUtils::avgRows(Mat &m, Mat &r){

	int rows = m.rows;
	int cols = m.cols;
	sumRows(m, r);

	for(int i = 0; i < rows; ++i){
		r.at<double>(i) = r.at<double>(i) / cols;
	}
}

/**
 * Averages along all cols in matrix m. Input is assumed to be a 2D matrix with
 * a scalar data type. On exit, r is a row vector with length equal to the
 * number of cols in m. Output data is of type double.
 *
 * @param m
 *   Input matrix
 * @param r
 *   Output destination. Overwritten on exit.
 */
void CVMatUtils::avgColumns(Mat &m, Mat &r){

	int rows = m.rows;
	int cols = m.cols;
	sumColumns(m, r);

	for(int j = 0; j < cols; ++j){
		r.at<double>(j) = r.at<double>(j) / rows;
	}
}

/**
 * Computes average along a single row within a matrix. Matrix is assumed
 * to be 2D with a scalar data type. Regardless of input data type depth,
 * output is a double.
 *
 * @param m
 *   Input matrix
 * @param r
 *   Row of m
 * @return
 *   Average of the values in row r of m
 */
double CVMatUtils::rowAvg(Mat &m, int r){

	int cols = m.cols;
	double sum = rowSum(m, r);
	return sum/cols;
}

/**
 * Computes average along a single column within a matrix. Matrix is assumed
 * to be 2D with a scalar data type. Regardless of input data type depth,
 * output is a double.
 *
 * @param m
 *   Input matrix
 * @param c
 *   Column of m
 * @return
 *   Average of the values in col c of m
 */
double CVMatUtils::columnAvg(Mat &m, int c){

	int rows = m.rows;
	double sum = columnSum(m, c);
	return sum/rows;
}

/**
 * Finds variance along all rows in matrix m. Input is assumed to be a 2D matrix
 * with a scalar data type. On exit, r is a col vector with length equal to the
 * number of rows in m. Output data is of type double.
 *
 * @param m
 *   Input matrix
 * @param r
 *   Output destination. Overwritten on exit.
 */
void CVMatUtils::varRows(Mat &m, Mat &r){

	int rows = m.rows;
	r.create(rows, 1, CV_64FC1);

	for(int i = 0; i < rows; ++i){
		r.at<double>(i) = rowVar(m, i);
	}
}

/**
 * Finds variance along all cols in matrix m. Input is assumed to be a 2D matrix
 * with a scalar data type. On exit, r is a row vector with length equal to the
 * number of cols in m. Output data is of type double.
 *
 * @param m
 *   Input matrix
 * @param r
 *   Output destination. Overwritten on exit.
 */
void CVMatUtils::varColumns(Mat &m, Mat &r){

	int cols = m.cols;
	r.create(1, cols, CV_64FC1);

	for(int j = 0; j < cols; ++j){
		r.at<double>(j) = columnVar(m, j);
	}
}

/**
 * Computes variance along a single row within a matrix. Matrix is assumed
 * to be 2D with a scalar data type. Regardless of input data type depth,
 * output is a double.
 *
 * @param m
 *   Input matrix
 * @param r
 *   Row of m
 * @return
 *   Variance of the values in row r of m
 */
double CVMatUtils::rowVar(Mat &m, int r){

	int cols = m.cols;
	double avg = rowAvg(m, r);

	Mat devs = m(cv::Range(r, r + 1), cv::Range(0, cols)) - avg;
	devs = devs.mul(devs);

	avg = rowSum(devs, 0) / cols;

	return avg;
}

/**
 * Computes variance along a single column within a matrix. Matrix is assumed
 * to be 2D with a scalar data type. Regardless of input data type depth,
 * output is a double.
 *
 * @param m
 *   Input matrix
 * @param c
 *   Column of m
 * @return
 *   Variance of the values in col c of m
 */
double CVMatUtils::columnVar(Mat &m, int c){

	int rows = m.rows;
	double avg = columnAvg(m, c);

	Mat devs = m(cv::Range(0, rows), cv::Range(c, c + 1)) - avg;
	devs = devs.mul(devs);

	avg = columnSum(devs, 0) / rows;

	return avg;
}

/**
 * Standardizes the values in each column of the given data matrix. Also,
 * returns the averages and standard deviations used to perform the standardization.
 *
 * @param data
 *   The data to standardize. On exit, each column has a zero mean and unit variance.
 * @param avgs
 *   On exit, a row vector containing the (pre-standardization) averages of each column.
 * @param stds
 *   On exit, a row vector containing the (pre-standardization) std devs of each column.
 */
void CVMatUtils::standardizeCols(cv::Mat data, cv::Mat avgs, cv::Mat stds){

	int n = data.rows;
	int m = data.cols;

	CVMatUtils::avgColumns(data, avgs);
	CVMatUtils::varColumns(data, stds);
	cv::sqrt(stds, stds);

	for(int f = 0; f < m; ++f){
		data(cv::Range(0, n), cv::Range(f, f + 1)) = data(cv::Range(0, n), cv::Range(f, f + 1)) - avgs.at<double>(f);
		data(cv::Range(0, n), cv::Range(f, f + 1)) = data(cv::Range(0, n), cv::Range(f, f + 1)) / stds.at<double>(f);
	}
}

/**
 * Fits a linear model to the given data. Fits model using least squares regression
 * fitting. Returns intercept, slope for p predictor variables, and residuals from the
 * fit model.
 *
 * @param y
 *   Column vector of response variable, length n column vector
 * @param x
 *   Matrix predictor of variables, n x p matrix
 * @param b
 *   Output destination - model intercept (element zero) and variable slope values (elements 1 - p+1)
 * @param eps
 *   Output destination - residuals (model output - y) from the fit model, length n column vector
 * @param stats
 *   Output destination - test statistics (t values) of each fit parameter from the fit model, length p+1 column vector
 * @return
 *   True if a model was successfully fit, false otherwise
 */
bool CVMatUtils::fitLinearModel(Mat &y, Mat &x, Mat &b, Mat &eps, Mat &stats){

	int p, n;
	double sigma = 0.0;
	bool result = true;

	n = x.rows;
	p = x.cols;

	Mat Xt, copy, inv, X(n, p + 1, CV_64FC1, 1.0);

	//Copy x into submatrix of X and get tranpsose
	copy = X.colRange(1, X.cols);
	x.copyTo(copy);
	Xt = X.t();

//	std::cout << std::endl << "X = " << std::endl;
//	std::cout << X << std::endl;
//
//	std::cout << std::endl << "Xt * X = " << std::endl;
//	std::cout << Xt * X << std::endl;

	//Compute coefficients
	inv = (Xt * X).inv();
	b = inv * Xt * y;

	eps = X*b - y;

	//Compute test stats for each model parameter
	for(int i = 0; i < eps.rows; ++i){
		sigma += eps.at<double>(i) * eps.at<double>(i);
	}

	sigma = sigma / (n - p - 1);
	inv = sigma * inv;
	inv.diag(0).copyTo(stats);

	for(int i = 0; i < stats.rows; ++i){
		stats.at<double>(i) = b.at<double>(i) / sqrt(stats.at<double>(i));
	}

	return result;
}

/**
 * Sorts thew rows of matrix based on the corresponding value in keys. Rows are sorted so into
 * ascending order (based on key values) on exit.
 *
 * @param matrix
 *   The matrix to be row-sorted
 * @param keys
 *   Values on which rows are sorted (in ascending order)
 * @param r
 *   Output destination, contains row-sorted version of matrix on exit
 * @param rk
 *   Output destination, contains sorted version of keys on exit
 */
void CVMatUtils::sortRowsWithKeys(Mat &matrix, Mat &keys, Mat &r, Mat &kr){

	cv::Mat idxs;

	cv::sortIdx(keys, idxs, cv::SORT_EVERY_COLUMN + cv::SORT_ASCENDING);
	cv::sort(keys, kr, cv::SORT_EVERY_COLUMN + cv::SORT_ASCENDING);

//	std::cout << idxs << std::endl;
	r.create(matrix.rows, matrix.cols, CV_64FC1);

	for(int i = 0; i < idxs.rows; ++i){
		int idx = idxs.at<int>(i);

//		r.rowRange(cv::Range(i, i + 1));
		matrix(cv::Range(idx, idx + 1), cv::Range::all()).copyTo(r(cv::Range(i, i + 1), cv::Range::all()));
//		r(cv::Range(i, i + 1), cv::Range::all()) = matrix(cv::Range(idx, idx + 1), cv::Range::all());

		std::cout << matrix(cv::Range(idx, idx + 1), cv::Range::all()) << std::endl;

	}
}

/**
 * Applies a mask to the given data to remove a subset of the rows.
 *
 * @param data
 *   The matrix from which to remove rows. Should be double type (CV_64FC1).
 * @param mask
 *   The mask indicating the rows to retain. Any value >0 indicates the row
 *   should be kept, others indicate it should be removed. Should be vector
 *   with length equal to data.rows. Should be int type (CV_32SC1).
 * @param r
 *   Output destination, contains retained rows on exit.
 */
void CVMatUtils::removeRows(cv::Mat &data, cv::Mat &mask, cv::Mat &r){

	int i, j, ii, n, m, count;

	n = data.rows;
	m = data.cols;
	count = 0;

	for(i = 0; i < n; ++i){
		if(mask.at<int>(i) > 0){
			++count;
		}
	}

	r.create(count, m, CV_64FC1);

	ii = 0;
	for(i = 0; i < n; ++i){
		if(mask.at<int>(i) > 0){

			for(j = 0; j < m; ++j){
				r.at<double>(ii, j) = data.at<double>(i, j);
			}
			++ii;
		}
	}
}

/**
 * Applies a mask to the given data to remove a subset of the columns.
 *
 * @param data
 *   The matrix from which to remove cols. Should be double type (CV_64FC1).
 * @param mask
 *   The mask indicating the cols to retain. Any value >0 indicates the col
 *   should be kept, others indicate it should be removed. Should be vector
 *   with length equal to data.cols. Should be int type (CV_32SC1).
 * @param r
 *   Output destination, contains retained cols on exit.
 */
void CVMatUtils::removeCols(cv::Mat &data, cv::Mat &mask, cv::Mat &r){

	int i, j, jj, n, m, count;

	n = data.rows;
	m = data.cols;
	count = 0;

	for(j = 0; j < m; ++j){
		if(mask.at<int>(j) > 0){
			++count;
		}
	}

	r.create(n, count, CV_64FC1);

	jj = 0;
	for(j = 0; j < m; ++j){
		if(mask.at<int>(j) > 0){

			for(i = 0; i < n; ++i){
				r.at<double>(i, jj) = data.at<double>(i, j);
			}
			++jj;
		}
	}
}

/**
 * Creates a new matrix containing only the specified row. The resulting matrix
 * is a row vector.
 *
 * @param data
 *   Matrix from which the row is extracted
 * @param r
 *   The row to extract. No bounds checking is performed.
 * @param r
 *   Output destination
 */
void CVMatUtils::getRow(cv::Mat &data, int row, cv::Mat &r){

	int m = data.cols;

	r.create(1, m, CV_64FC1);

	for(int i = 0; i < m; ++i){
		r.at<double>(0, i) = data.at<double>(row, i);
	}
}

/**
 * Creates a new matrix containing only the specified row. The resulting matrix
 * is a row vector.
 *
 * @param data
 *   Matrix from which the row is extracted
 * @param r
 *   The row to extract. No bounds checking is performed.
 * @param r
 *   Output destination
 */
void CVMatUtils::getCol(cv::Mat &data, int col, cv::Mat &r){

	int n = data.rows;

	r.create(n, 1, CV_64FC1);

	for(int i = 0; i < n; ++i){
		r.at<double>(i, 0) = data.at<double>(i, col);
	}
}

/**
 * Get a sample of the rows in the
 *
 * @param m
 *   Matrix from which to sample
 * @param rows
 *   List of row indices indicating the rows to sample. Should be a vector of ints.
 * @param r
 *   Output destination, contains sampled rows on output
 */
int CVMatUtils::sampleRows(Mat &m, Mat &rows, Mat &r){

	int n, c, nn;

	n = m.rows;
	c = m.cols;
	nn = rows.rows;

	r.create(nn, c, CV_64FC1);

	for(int i = 0; i < nn; ++i){
		int idx = rows.at<int>(i);
		for(int j = 0; j < c; ++j){
			r.at<double>(i, j) = m.at<double>(idx, j);
		}
	}

	return r.rows;
}

/**
 * Create a new matrix by appending all the rows of add to the matrix
 * base. The inputs base and add should have the same number of columns.
 *
 * @param base
 *   Matrix to which the rows of add will be appended
 * @param add
 *   Matrix to append to base
 * @param r
 *   Output destination
 */
void CVMatUtils::appendRows(cv::Mat &base, cv::Mat &add, cv::Mat &r){

	assert(base.cols == add.cols);

	int n, nn, m, i;

	n = base.rows;
	m = base.cols;
	nn = n + add.rows;

	r.create(nn, m, CV_64FC1);

	for(int j = 0; j < m; ++j){

		for(i = 0; i < n; ++i){
			r.at<double>(i, j) = base.at<double>(i, j);
		}
		for(; i < n; ++i){
			r.at<double>(i, j) = add.at<double>(i - n, j);
		}
	}

}

/**
 * Create a new matrix by appending all the columns of add to the matrix
 * base. The inputs base and add should have the same number of rows.
 *
 * @param base
 *   Matrix to which the cols of add will be appended
 * @param add
 *   Matrix to append to base
 * @param r
 *   Output destination
 */
void CVMatUtils::appendCols(cv::Mat &base, cv::Mat &add, cv::Mat &r){

	assert(base.rows == add.rows);

	int n, m, mm, j;

	n = base.rows;
	m = base.cols;
	mm = m + add.cols;

	r.create(n, mm, CV_64FC1);

	for(int i = 0; i < n; ++i){

		for(j = 0; j < m; ++j){
			r.at<double>(i, j) = base.at<double>(i, j);
		}
		for(; j < mm; ++j){
			r.at<double>(i, j) = add.at<double>(i, j - m);
		}
	}
}

/**
 * Writes a matrix to a comma-separated value file.
 *
 * @param m
 *   The matrix to write out
 * @param path
 *   The output file path. It is overwritten, if it already exists.
 * @return
 *   True if matirx is successfully written, false otherwise.
 */
bool CVMatUtils::writeCSV(Mat &m, const char *path){

	int rows = m.rows;
	int cols = m.cols;
	std::ofstream out(path);

	Mat md;
	m.convertTo(md, CV_64FC1);

	for(int i = 0; i < rows; ++i){
		for(int j = 0; j < cols - 1; ++j){
			out << md.at<double>(i,j) << ",";
		}
		out << md.at<double>(i,cols - 1);
		out << std::endl;
	}

	out.close();
	return true;
}

/**
 * Writes a matrix to a comma-separated value file with a header line including the
 * given column names.
 *
 * @param m
 *   The matrix to write out
 * @param path
 *   The output file path. It is overwritten, if it already exists.
 * @param colnames
 *   List of col names to write to header line.
 * @return
 *   True if matirx is successfully written, false otherwise.
 */
bool CVMatUtils::writeCSV(Mat &m, const char *path, std::vector<std::string> &colnames){

	int rows = m.rows;
	int cols = m.cols;
	std::ofstream out(path);

	Mat md;
	m.convertTo(md, CV_64FC1);

	for(int j = 0; j < cols - 1; ++j){
		out << colnames[j] << ",";
	}
	out << colnames[cols - 1] << std::endl;

	for(int i = 0; i < rows; ++i){
		for(int j = 0; j < cols - 1; ++j){
			out << md.at<double>(i,j) << ",";
		}
		out << md.at<double>(i,cols - 1);
		out << std::endl;
	}

	out.close();
	return true;
}

/**
 * Prints the a matrix to stdout using the given delimiter.
 *
 * @param m
 *   The matrix to print
 * @param delim
 *   Delimiter used to separate columns. Rows are separated by new lines.
 */
void CVMatUtils::printMatrix(Mat &m, std::string delim){

	int rows = m.rows;
	int cols = m.cols;

	Mat md;
	m.convertTo(md, CV_64FC1);

	for(int i = 0; i < rows; ++i){
		for(int j = 0; j < cols - 1; ++j){
			std::cout << md.at<double>(i,j) << ",";
		}
		std::cout << md.at<double>(i,cols - 1);
		std::cout << std::endl;
	}
}

/**
 * Reads a comma-separated value file into a matrix.
 *
 * @param m
 *   The output destination. Filled with contents of file oon exit.
 * @param path
 *   Path to the csv file.
 */
void CVMatUtils::readCSV(Mat &m, const char *path){

	int pos, len, rows = 0, cols = 0;

	//Count rows/cols in file
	std::ifstream in(path);
	std::string line;
	while (std::getline(in, line)){
		if(!line.empty()){
			++rows;
		}
		if(rows == 1){
			cols = std::count(line.begin(), line.end(), ',') + 1;
		}
	}

	//Go to begin and read data
	in.clear();
	in.seekg(0, std::ios::beg);

	m.create(rows, cols, CV_64FC1);

	for(int i = 0; i < rows; ++i){

		std::getline(in, line);

		pos = 0;
		len = line.find(",", pos);

		for(int j = 0; j < cols; ++j){
			m.at<double>(i, j) = atof(line.substr(pos, len).c_str());
			pos = pos + len + 1;
			len = (j == cols - 1) ? len = line.length() - pos : len = line.find(",", pos) - pos;
		}
	}

	in.close();
}

/**
 * Reads a comma-separated value file into a matrix. Assumes the csv file has a header line
 * containing column names. These names are read into the given vector.
 *
 * @param m
 *   The output destination. Filled with contents of file on exit.
 * @param colnames
 *   Output destination. Filled with the col names read from the csv.
 * @param path
 *   Path to the csv file.
 */
void CVMatUtils::readCSV(Mat &m, const char *path, std::vector<std::string> &colnames){

	int pos, len, rows = 0, cols = 0;

	//Count rows/cols in file
	std::ifstream in(path);
	std::string line;
	while (std::getline(in, line)){
		if(!line.empty()){
			++rows;
		}
		if(rows == 1){
			cols = std::count(line.begin(), line.end(), ',') + 1;
		}
	}

	//Go to begin and read data
	in.clear();
	in.seekg(0, std::ios::beg);

	colnames.clear();
	std::getline(in, line);
	pos = 0;
	len = line.find(",", pos);

	for(int j = 0; j < cols; ++j){
		colnames.push_back(line.substr(pos, len));
		pos = pos + len + 1;
		len = (j == cols - 1) ? len = line.length() - pos : len = line.find(",", pos) - pos;
	}

	rows--;

	m.create(rows, cols, CV_64FC1);

	for(int i = 0; i < rows; ++i){

		std::getline(in, line);

		pos = 0;
		len = line.find(",", pos);

		for(int j = 0; j < cols; ++j){
			m.at<double>(i, j) = atof(line.substr(pos, len).c_str());
			pos = pos + len + 1;
			len = (j == cols - 1) ? len = line.length() - pos : len = line.find(",", pos) - pos;
		}
	}

	in.close();
}

/**
 * Convert the given OpenCV matrix to an ITK image. Assumes the matrix should
 * have scalar values of type float (CV_32FC1).
 *
 * @param src
 *   The matrix of values to convert to an ITK image
 * @param dest
 *   Output destination. Pixels values are set using the input matrix.
 */
void CVMatUtils::matToImage(Mat &src, FloatImageType::Pointer &dest){

	int rows = src.rows;
	int cols = src.cols;

	FloatImageType::IndexType idx;
	idx.Fill(0);
	FloatImageType::SizeType s;
	s[0] = rows; s[1] = cols;

	FloatImageType::RegionType r;
	r.SetIndex(idx); r.SetSize(s);

//	dest = FloatImageType::New();
//	dest->Delete();
	dest->SetRegions(r);
	dest->Allocate();

	itk::ImageRegionIterator<FloatImageType> it(dest, dest->GetLargestPossibleRegion());

	for(int i = 0; i < rows; ++i){
		for(int j = 0; j < cols; ++j){
			it.Set(src.at<double>(i, j));
			++it;
		}
	}
}

/**
 * Convert the given ITK image to an OpenCV matrix. Creates a matrix with
 * scalar values of type float (CV_32FC1).
 *
 * @param src
 *   The image to convert to a matrix
 * @param dest
 *   Output destination. Matrix elements are set using image pixel values.
 */
void CVMatUtils::imageToMat(FloatImageType::Pointer &src, Mat &dest){

	int rows, cols;

	FloatImageType::RegionType r = src->GetLargestPossibleRegion();
	rows = r.GetSize()[0];
	cols = r.GetSize()[1];

	itk::ImageRegionIterator<FloatImageType> it(src, src->GetLargestPossibleRegion());

	dest.create(rows, cols, CV_32FC1);

	for(int i = 0; i < rows; ++i){
		for(int j = 0; j < cols; ++j){
			dest.at<float>(i, j) = it.Get();
			++it;
		}
	}
}

/**
 * Reads a set of images into a data matrix. The images should all be of the same size.
 * The output matrix is of dimension (number of images) x (width * height).
 *
 * @param images
 *   A set of images to read into a matrix
 * @param mat
 *   Output destination, contains image data (one per row) on exit
 * @return
 *   Returns the width of the images stored in the matrix
 */
int CVMatUtils::imagesToMat(ImageDataset images, Mat &mat){

	int n, w, h, dim;

	Mat curimage, out;

	n = images.getN();

	//Read and vectorized each image
	for(int i = 0; i < n; ++i){

		std::string ipath = images.getPathForImage(i);
		FloatUtils disp(ipath);

		if(i == 0){
			w = disp.getWidth();
			h = disp.getHeight();
			dim = w*h;
			mat.create(n, dim, CV_64FC1);
		}

		FloatImageType::Pointer dispimage = disp.getImage();
		CVMatUtils::imageToVector(dispimage, curimage);

		for(int j = 0; j < dim; ++j){
			mat.at<double>(i, j) = curimage.at<float>(j);
		}
	}

	return w;
}

/**
 * Reads a set of images into a data matrix. The images should all be of the same size.
 * The output matrix is of dimension (number of images) x (number of non-zero values in mask).
 *
 * @param images
 *   A set of images to read into a matrix
 * @param mask
 *   A mask to apply to all images. Only locations with non-zero mask values are extracted from the images.
 * @param mat
 *   Output destination, contains image data (one per row) on exit
 * @return
 *   Returns the width of the images stored in the matrix
 */
int CVMatUtils::imagesToMat(ImageDataset images, ByteImageType::Pointer &mask, Mat &mat){

	int n, w, h, dim;

	Mat curimage, out;

	n = images.getN();

	//Read and vectorized each image
	for(int i = 0; i < n; ++i){

		std::string ipath = images.getPathForImage(i);
		FloatUtils disp(ipath);

		FloatImageType::Pointer dispimage = disp.getImage();
		CVMatUtils::imageToVector(dispimage, mask, curimage);

		if(i == 0){
			w = disp.getWidth();
			h = disp.getHeight();
			dim = curimage.cols;
			mat.create(n, dim, CV_64FC1);
		}

		for(int j = 0; j < dim; ++j){
			mat.at<double>(i, j) = curimage.at<float>(j);
		}
	}

	return w;
}

/**
 * Vectorizes the given image, creating a row vector of length image width times
 * height. Raster scan performed in row major order.
 *
 * @param src
 *   The image to vectorize
 * @param dest
 *    The output destination
 */
void CVMatUtils::imageToVector(FloatImageType::Pointer &src, Mat &dest){

	int i = 0;
	FloatImageType::SizeType s = src->GetLargestPossibleRegion().GetSize();

	dest.create(1, s[0]*s[1], CV_32FC1);

	itk::ImageRegionIterator<FloatImageType> it(src, src->GetLargestPossibleRegion());

	for(it.GoToBegin(); !it.IsAtEnd(); ++it){
		dest.at<float>(i++) = it.Get();
	}
}

/**
 * Vectorizes the given image, creating a row vector with length equal to the number of non-zero
 * values in mask. Raster scan performed in row major order.
 *
 * @param src
 *   The image to vectorize
 * @param mask
 *   Mask to apply to image. Zero-valued mask locations are not copied into the vector.
 * @param dest
 *   The output destination
 */
void CVMatUtils::imageToVector(FloatImageType::Pointer &src, ByteImageType::Pointer &mask, Mat &dest){

	int i = 0, zeros = 0;
	FloatImageType::SizeType s = src->GetLargestPossibleRegion().GetSize();

	itk::ImageRegionIterator<ByteImageType> mIt(mask, mask->GetLargestPossibleRegion());
	for(mIt.GoToBegin(); !mIt.IsAtEnd(); ++mIt){
		if(mIt.Get() == 0){
			++zeros;
		}
	}

	dest.create(1, s[0]*s[1] - zeros, CV_32FC1);

	itk::ImageRegionIterator<FloatImageType> it(src, src->GetLargestPossibleRegion());

	mIt.GoToBegin();
	for(it.GoToBegin(); !it.IsAtEnd(); ++it){
		if(mIt.Get() != 0){
			dest.at<float>(i++) = it.Get();
		}
		++mIt;
	}
}

/**
 * Reshapes/converts vector into a 2D image.
 *
 * @param src
 *   Vector to convert. Should be length w*h.
 * @param w
 *   Width of the output image
 * @param h
 *   Height of the output image
 * @param dest
 *   Output destination, float-valued image of size w x h
 */
void CVMatUtils::vectorToImage(Mat &src, int w, int h, FloatImageType::Pointer &dest){

	int i = 0;

	cv::Mat tempflt;
	src.convertTo(tempflt, CV_32FC1);

	itk::ImageRegionIterator<FloatImageType> it(dest, dest->GetLargestPossibleRegion());

	for(it.GoToBegin(); !it.IsAtEnd(); ++it){
		it.Set(tempflt.at<float>(i++));
	}
}

/**
 * Reshapes/converts vector into a 2D image using the given mask.
 *
 * @param src
 *   Vector to convert. Should have length equal to number of non-zero values in mask.
 * @param mask
 *   Mask indicating locations to ignore. Zero-valued mask locations are set to zero in the output.
 * @param w
 *   Width of the output image
 * @param h
 *   Height of the output image
 * @param dest
 *   Output destination, float-valued image of size w x h
 */
void CVMatUtils::vectorToImage(Mat &src, ByteImageType::Pointer &mask, int w, int h, FloatImageType::Pointer &dest){

	int i = 0;

	cv::Mat tempflt;
	src.convertTo(tempflt, CV_32FC1);

	itk::ImageRegionIterator<FloatImageType> it(dest, dest->GetLargestPossibleRegion());
	itk::ImageRegionIterator<ByteImageType> mIt(mask, mask->GetLargestPossibleRegion());

	mIt.GoToBegin();
	for(it.GoToBegin(); !it.IsAtEnd(); ++it){
		if(mIt.Get() != 0){
			it.Set(tempflt.at<float>(i++));
		} else{
			it.Set(0);
		}
		++mIt;
	}
}
