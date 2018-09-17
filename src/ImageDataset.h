/*
 * ImageDataset.h
 *
 *  Created on: Jun 20, 2013
 *      Author: mchristopher
 */

#ifndef IMAGEDATASET_H_
#define IMAGEDATASET_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <algorithm>

//#include <boost/algorithm/string.hpp>

#define EXTS "tiff,tif,gif,jpg,jpeg,png,vtk"

/**
 * Provides functionality to handle managing the location of the files associated
 * with an image dataset. The dataset can be defined based on a list of images or a
 * data directory storing image files. Provides methods for storing those locations
 * and uses simple conventions to create names for intermediate and final result
 * images.
 *
 * Images are accessed using an integer index.
 *
 * Sets of files to be used as truth or masks for each image can also be set.
 *
 */
class ImageDataset {

private:

	/** List of directories for each image */
	std::vector<std::string> dirs;

	/** List of names for each image */
	std::vector<std::string> names;

	/** List of file extensions for each image */
	std::vector<std::string> exts;

	/** List of files containing truth for each image */
	std::vector<std::string> truths;

	/** List of files containing truth for each image */
	std::vector<std::string> masks;

	/** Directory to store feature images, default = "" */
	std::string featureDir;

	/** Directory to store result images, default = "" */
	std::string resultDir;

	/** Directory to store auxillary images, default = "" */
	std::string auxDir;

	/** The file extension used to create auxillary image file paths */
	std::string auxExt;

	/** List of files to use as auxillary files rather than the standard naming */
	std::vector<std::string> auxList;

	void init();

public:
	ImageDataset(){init();};
	ImageDataset(std::string &file);
	ImageDataset(std::vector<std::string> &files);
	ImageDataset(std::string &dir, const std::string &allowedExts);

	void setFileList(std::vector<std::string> &files);
	void setAuxFileList(std::vector<std::string> &auxFiles);
	void setFile(std::string &file);
	void setAuxFile(std::string &auxFile);

	std::string getPathForImage(int i);
	std::string getBasenameForImage(int i);
	std::string getExtensionForImage(int i);
	std::string getDirForImage(int i);

	std::string getPathForFeatureImage(int i, int f);
	std::string getPathForAuxillaryImage(int i, int a = 0);
	std::string getPathForResultImage(int i);
	std::string getPathForTruth(int i);
	std::string getPathForMask(int i);

	bool imageHasMask(int i);
	bool imageHasTruth(int i);

	void setDataDir(std::string &dirPath, const std::string &allowedExts = EXTS);
	void setTruthDir(std::string &dirPath, const std::string &allowedExts = EXTS);
	void setMaskDir(std::string &dirPath, const std::string &allowedExts = EXTS);

	//Static methods
	static bool directoryExists(std::string &dir);
	static bool makeDirectory(std::string &dir);
	static bool fileExists(std::string &path);
	static bool isRecognizedImage(std::string &path);
	static bool isAllowedExt(std::string ext, std::vector<std::string> list);
	static bool hasExtension(const std::string &file, const std::string &ext);
	static void readLines(std::string &path, std::vector<std::string> &list);
	static void expandList(std::string list, std::vector<std::string> &result);

	//Getters/Setters
	int inline getN(){
		return this->names.size();
	}

	void setFeatureDir(std::string &dir);
	std::string inline getFeatureDir(){
		return this->featureDir;
	}

	void setAuxillaryDir(std::string &dir);
	std::string inline getAuxillaryDir(){
		return this->auxDir;
	}

	void setResultDir(std::string &dir);
	std::string inline getResultDir(){
		return this->resultDir;
	}

	std::string getAuxillaryExtension(){
		return this->auxExt;
	}
	void setAuxillaryExtension(const std::string &ext){
		this->auxExt = ext;
	}

	void setTruthFiles(std::vector<std::string> truthFiles){
		this->truths = truthFiles;
	}

	void setMaskFiles(std::vector<std::string> maskFiles){
		this->masks = maskFiles;
	}

	virtual ~ImageDataset();
};

#endif /* IMAGEDATASET_H_ */
