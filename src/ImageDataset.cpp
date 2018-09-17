/*
 * ImageDataset.cpp
 *
 *  Created on: Jun 20, 2013
 *      Author: mchristopher
 */

#include "ImageDataset.h"

/**
 * Creates a dataset consisting of the given list of images.
 *
 * Attempts to seperate passed file names into directory, base name, and
 * file extension components based on the assumption that they will be in
 * the form: /path/to/dir/basename.ext. The "/path/to/dir/" may be excluded.
 *
 * @param files
 *   A list of files to be represented by this dataset.
 */
ImageDataset::ImageDataset(std::vector<std::string> &files){
	init();
	this->setFileList(files);
}

/**
 * Create a dataset consisting of the single given image.
 *
 * Attempts to separate passed file name into directory, base name, and
 * file extension components based on the assumption that they will be in
 * the form: /path/to/dir/basename.ext. The "/path/to/dir/" may be excluded.
 *
 * @param file
 *   The single file to be represented by this image.
 */
ImageDataset::ImageDataset(std::string &file){
	init();
	this->setFile(file);
}

/**
 * Creates a dataset consisting of the files contained in the given directory.
 *
 * The extList parameter should be a comma separated list of files extensions.The
 * files in the directory are filtered based on the set of given extensions. Only
 * those files with one of the indicated extension are included in the dataset.
 * (Eg. "tif,jpeg" will result in only files with the tif and jpeg extensions within
 * dir being included in the dataset).
 *
 * @param dir
 *   Directory containing the images
 * @param extList
 *   List of file extension to include
 */
ImageDataset::ImageDataset(std::string &dir, const std::string &allowedExts){
	init();
	this->setDataDir(dir, allowedExts);
}

/**
 * Sets default values.
 */
void ImageDataset::init(){
	this->featureDir = "";
	this->resultDir = "";
}

/**
 * Splits a comma separated list into a vector of strings.
 *
 * @param list
 *   Comma separated list of items
 * @param result
 *   Vector in which to store the results
 */
void ImageDataset::expandList(std::string list, std::vector<std::string> &result){

	int p = 0;
	std::string temp(list);
	result.clear();

	while((p = temp.find(",")) != std::string::npos){
		result.push_back(temp.substr(0, p));
		temp.erase(0, p + 1);
	}
	if(!temp.empty()){
		result.push_back(temp);
	}
}

/**
 * Sets the list of files to be included in this dataset to the single image
 * given by file.
 *
 * @param
 *   A single file to be represented by this dataset.
 */
void ImageDataset::setFile(std::string &file){
	std::vector<std::string> files;
	files.push_back(file);
	this->setFileList(files);
}

/**
 * Sets the list of aux files to be included in this dataset to the single given file.
 *
 * @param
 *   A single file to be used as an auxillary file.
 */
void ImageDataset::setAuxFile(std::string &auxFile){
	std::vector<std::string> files;
	files.push_back(auxFile);
	this->setAuxFileList(files);
}

/**
 * Sets list of files to be included in the dataset.
 *
 * @param files
 *   List of files to include in dataset
 */
void ImageDataset::setFileList(std::vector<std::string> &files){
	for(int i = 0; i < files.size(); ++i){
		std::string f = files[i];
		int didx = f.find_last_of("/");
		int eidx = f.find_last_of(".");

		this->names.push_back(f.substr(didx + 1, eidx - didx - 1));
		this->exts.push_back(f.substr(eidx + 1, f.size() - eidx));

		if(didx != std::string::npos){
			this->dirs.push_back(f.substr(0, didx));
		}
		else{
			this->dirs.push_back("");
		}
	}
}

/**
 * Sets list of files to be used as auxillary files. This overrides the standard convention
 * used for naming auxillary files.
 *
 * @param files
 *   List of files to consider the auxillary files for this dataset
 */
void ImageDataset::setAuxFileList(std::vector<std::string> &auxFiles){

	this->auxList.clear();

	for(int i = 0; i < auxFiles.size(); ++i){
		this->auxList.push_back(auxFiles[i]);
	}
}

/**
 * Set the directory containing image data to be represented by the dataset.
 *
 * @param dirPath
 *   Path to data directory
 * @param allowedExts
 *   Comma separated list of file extensions to include in the dataset
 */
void ImageDataset::setDataDir(std::string &dirPath, const std::string &allowedExts){

	DIR *dir = NULL;
	std::vector<std::string> extList;
	ImageDataset::expandList(allowedExts, extList);

	this->names.clear();
	this->exts.clear();
	this->dirs.clear();

	struct dirent *ent;
	if ((dir = opendir (dirPath.c_str())) != NULL) {
		while ((ent = readdir (dir)) != NULL) {

			std::string f(ent->d_name);
			int eidx = f.find_last_of(".");
			std::string ext(f.substr(eidx + 1, f.size() - eidx));

			if(ImageDataset::isAllowedExt(ext, extList)){
				this->dirs.push_back(dirPath);
				this->names.push_back(f.substr(0, eidx));
				this->exts.push_back(ext);
//				std::cout << f.substr(0, eidx) << "." << ext << std::endl; std::cout.flush();
			}
		}
		closedir (dir);
	} else {
		throw std::string("ImageDataset::Couldn't open directory " + dirPath);
	}
}

/**
 * Set the truth files for the images in the dataset based on the contents of the
 * indicated directory.
 *
 * Order of truth files is set by returned order of readdir().
 *
 * @param dirPath
 *   Path to directory containing truth files
 * @param allowedExts
 *   Comma separated list of file extensions to include
 */
void ImageDataset::setTruthDir(std::string &dirPath, const std::string &allowedExts){
	DIR *dir = NULL;
	std::vector<std::string> extList;
	ImageDataset::expandList(allowedExts, extList);

	this->truths.clear();

	struct dirent *ent;
	if ((dir = opendir (dirPath.c_str())) != NULL) {
		while ((ent = readdir (dir)) != NULL) {
			std::string f(ent->d_name);
			int eidx = f.find_last_of(".");
			std::string ext(f.substr(eidx + 1, f.size() - eidx));

			if(ImageDataset::isAllowedExt(ext, extList)){
				this->truths.push_back(dirPath + "/" + f);
			}
		}
		closedir (dir);
	} else {
		throw std::string("ImageDataset::Couldn't open directory " + dirPath);
	}
}

/**
 * Set the mask files for the images in the dataset based on the contents of the
 * indicated directory.
 *
 * Order of mask files is set by returned order of readdir().
 *
 * @param dirPath
 *   Path to directory containing mask files
 * @param allowedExts
 *   Comma separated list of file extensions to include
 */
void ImageDataset::setMaskDir(std::string &dirPath, const std::string &allowedExts){
	DIR *dir = NULL;
	std::vector<std::string> extList;
	ImageDataset::expandList(allowedExts, extList);

	this->masks.clear();

	struct dirent *ent;
	if ((dir = opendir (dirPath.c_str())) != NULL) {
		while ((ent = readdir (dir)) != NULL) {
			std::string f(ent->d_name);
			int eidx = f.find_last_of(".");
			std::string ext(f.substr(eidx + 1, f.size() - eidx));

			if(ImageDataset::isAllowedExt(ext, extList)){
				this->masks.push_back(dirPath + "/" + f);
			}
		}
		closedir (dir);
	} else {
		throw std::string("ImageDataset::Couldn't open directory " + dirPath);
	}
}

/**
 * Performs case-insensitve comparisons ext to determine if ext occurs within list.
 *
 * @param ext
 *   The extension for which to search
 * @param list
 *   List of strings within which to search
 * @return
 *   True if ext is (case-insensitive) equals a string in list, false otherwise
 */
bool ImageDataset::isAllowedExt(std::string ext, std::vector<std::string> list){

	bool allowed = false;
	std::string extLower(ext);
	std::transform(extLower.begin(), extLower.end(), extLower.begin(), ::tolower);

	for(int i = 0; i < list.size(); ++i){

		std::string curLower(list[i].c_str());
		std::transform(curLower.begin(), curLower.end(), curLower.begin(), ::tolower);

		if(curLower == extLower){
			allowed = true;
			break;
		}
	}
	return allowed;
}

/**
 * Performs case insensitive comparison to determine if file has given extension.
 *
 * @param file
 *   File name to check
 * @param ext
 *   File extension
 * @return
 *   True if file ends in ext, false otherwise
 */
bool ImageDataset::hasExtension(const std::string &file, const std::string &ext){
	std::string str(file);
	std::string suffix(ext);

	std::transform(str.begin(), str.end(), str.begin(), ::tolower);
	std::transform(suffix.begin(), suffix.end(), suffix.begin(), ::tolower);

	return str.size() >= suffix.size() && str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

/**
 * Gets the path to the image i file.
 *
 * @param i
 *   Index of the image
 * @result
 *   The path to image i
 */
std::string ImageDataset::getPathForImage(int i){
	return (this->dirs[i].size() != 0) ? this->dirs[i] + "/" + this->names[i] + "." + this->exts[i] :
										 this->names[i] + "." + this->exts[i];
}
std::string ImageDataset::getBasenameForImage(int i){
	return this->names[i];
}
std::string ImageDataset::getExtensionForImage(int i){
	return this->exts[i];
}
std::string ImageDataset::getDirForImage(int i){
	return this->dirs[i];
}

/**
 * Retrieves path to the file containing the truth for image i.
 *
 * List of truth files should be set using setTruthFiles() before calling this
 * function.
 *
 * @param i
 *   Index of the image for which to retrieve the truth file
 * @return
 *   File containing the truth for image i
 */
std::string ImageDataset::getPathForTruth(int i){
	return this->truths[i];
}

/**
 * Determines if image i has a known truth file.
 *
 * @param i
 *   Index of image
 * @return
 *   True if image i has known truth file, false otherwise
 */
bool ImageDataset::imageHasTruth(int i){
	return i < 0 || i >= this->truths.size();
}

/**
 * Retrieves path to the file containing the mask for image i.
 *
 * List of mask files should be set using setMaskFiles() before calling this
 * function.
 *
 * @param i
 *   Index of the image for which to retrieve the mask file
 * @return
 *   File containing the mask for image i
 */
std::string ImageDataset::getPathForMask(int i){
	return this->masks[i];
}

/**
 * Determines if image i has a known mask file.
 *
 * @param i
 *   Index of image
 * @return
 *   True if image i has known mask file, false otherwise
 */
bool ImageDataset::imageHasMask(int i){
	return i < 0 || i >= this->masks.size();
}

/**
 * Creates a path to store a feature f image for image i. Path is in the form:
 *
 * /path/to/featureDir/<image i name>-f<f>.vtk
 *
 * @param i
 *   Index of the image for which a feature file path is created
 * @param f
 *   Inex of the feature for which a feature file path is created
 * @result
 *   Path to feature f image for image i
 */
std::string ImageDataset::getPathForFeatureImage(int i, int f){
	std::stringstream ss;
	ss << this->featureDir << "/" << this->names[i] << "-f" << f << ".vtk";
	return ss.str();
}

/**
 * Creates a path to store an auxillary image for image i. Path is in the form:
 *
 * /path/to/auxDir/<image i name>-a<a>.vtk
 *
 * @param i
 *   Index of the image for which a auxillary file path is created
 * @param a
 *   Inex of the auxillary image for which a file path is created
 * @result
 *   Path to auxillary image a for image i
 */
std::string ImageDataset::getPathForAuxillaryImage(int i, int a){

//	std::cout << "AUX SHIT: " << this->auxList.size() << std::endl; std::cout.flush();

	if(this->auxList.size() > 0 && a <= 0){
		return this->auxList[i];
	}

//	std::cout << "AUX SHIT" << std::endl; std::cout.flush();

	std::stringstream ss;
	ss << this->auxDir << "/" << this->names[i] << "-a" << a << "." << this->auxExt;
	return ss.str();
}

/**
 * Creates a path to store a result image for image i. Path is in the form:
 *
 * /path/to/resutDir/<image i name>-result.vtk
 *
 * @param i
 *   Index of the image for which a result file path is created
 * @result
 *   Path to result image for image i
 */
std::string ImageDataset::getPathForResultImage(int i){
	std::stringstream ss;
	ss << this->resultDir << "/" << this->names[i] << "-result.vtk";
	return ss.str();
}

/**
 * Determines if the given directory exists.
 *
 * @param dir
 *   Path to directory to check
 * @return
 *   True if the directory exists, false otherwise
 */
bool ImageDataset::directoryExists(std::string &dir){
	DIR *d;
	bool exists = false;

	d = opendir(dir.c_str());

	if (d != NULL){
		exists = true;
		closedir(d);
	}

	return exists;
}

/**
 * Create a directory at the given location, if possible.
 *
 * @param dir
 *   Location a which the directory will be created
 * @return
 *   True if dir successfully created, false otherwise
 */
bool ImageDataset::makeDirectory(std::string &dir){
	return mkdir(dir.c_str(), 0755) == 0;
}

/**
 * Checks for file existence.
 *
 * @param filename
 *   Path to the file
 * @return
 *   True if the files exists, false otherwise
 */
bool ImageDataset::fileExists(std::string &path){
	std::ifstream ifile(path.c_str());
	return ifile;
}

/**
 * Checks to see if file is a recognized image file based on extension.
 *
 * @param path
 *   Path to file to be checked
 * @return
 *   True is path points to a recognized image file, false otherwise
 */
bool ImageDataset::isRecognizedImage(std::string &path){

	bool result;
	int eidx;

	std::string allowed(EXTS);
	std::vector<std::string> exts;
	ImageDataset::expandList(allowed, exts);

	eidx = path.find_last_of(".");
	std::string ext(path.substr(eidx + 1, path.size() - eidx));

	return ImageDataset::fileExists(path) && ImageDataset::isAllowedExt(ext, exts);
}

/**
 * Reads each line in a text file and stores the lines the given vector.
 *
 * Useful for reading lists of file names.
 *
 * @param path
 *   Path to the text file to read
 * @param
 *   Output destination, stores one line per element
 */
void ImageDataset::readLines(std::string &path, std::vector<std::string> &list){

	char line[1000];
	std::ifstream in(path.c_str());

	list.clear();

	if(in.is_open()){
		while(!in.eof()){
			in.getline(line, 1000);
			std::string s(line);

//			boost::trim(s);

			if(!s.empty()){
				list.push_back(s);
			}
		}
	}
	in.close();
}

/**
 * Sets the feature directory for this dataset. If dir doesn't exist, it is
 * created during this call.
 *
 * @param dir
 *   Path to directory to store feature images
 */
void ImageDataset::setFeatureDir(std::string &dir){
	this->featureDir = dir;

	if(!ImageDataset::directoryExists(dir)){
		mkdir(dir.c_str(), 0755);
	}
}

/**
 * Sets the auxillary directory for this dataset. If dir doesn't exist, it is
 * created during this call.
 *
 * @param dir
 *   Path to directory to store auxillary images
 */
void ImageDataset::setAuxillaryDir(std::string &dir){
	this->auxDir = dir;

	if(!ImageDataset::directoryExists(dir)){
		mkdir(dir.c_str(), 0755);
	}
}

/**
 * Sets the result directory for this dataset. If dir doesn't exist, it is
 * created during this call.
 *
 * @param dir
 *   Path to directory to store result images
 */
void ImageDataset::setResultDir(std::string &dir){
	this->resultDir = dir;

	if(!ImageDataset::directoryExists(dir)){
		mkdir(dir.c_str(), 0755);
	}
}

ImageDataset::~ImageDataset() {
	// TODO Auto-generated destructor stub
}
