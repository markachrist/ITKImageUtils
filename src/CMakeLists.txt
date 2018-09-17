cmake_minimum_required(VERSION 2.8)

project(ImageUtilities)

# Find ITK.
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

set(OPENCV_INC "/usr/local/Cellar/opencv/2.4.13/include")
include_directories(${OPENCV_INC})

set(OPENCV_LIB "/usr/local/Cellar/opencv/2.4.13/lib")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++ -std=c++03")

add_executable(imageUtilsMain
			   main.cpp
			   RGBUtils.cpp
			   ImageDataset.cpp
			   CVMatUtils.cpp
			   ProcessObserver.cpp
			   Performance.cpp)

TARGET_LINK_LIBRARIES(imageUtilsMain ${ITK_LIBRARIES} ${OPENCV_LIB}/libopencv_core.dylib ${OPENCV_LIB}/libopencv_ml.dylib)
			   
add_executable(transform
			   transformMain.cpp
			   RGBUtils.cpp
			   ImageDataset.cpp
			   CVMatUtils.cpp)

TARGET_LINK_LIBRARIES(transform ${ITK_LIBRARIES} ${OPENCV_LIB}/libopencv_core.dylib)

add_executable(extractChannel
			   extractChannel.cpp
			   RGBUtils.cpp
			   ImageDataset.cpp
			   CVMatUtils.cpp)

TARGET_LINK_LIBRARIES(extractChannel ${ITK_LIBRARIES} ${OPENCV_LIB}/libopencv_core.dylib)

add_executable(averageImages
			   averageImages.cpp)

TARGET_LINK_LIBRARIES(averageImages ${ITK_LIBRARIES})

add_executable(convertImage
			   convertImage.cpp
			   RGBUtils.cpp)

TARGET_LINK_LIBRARIES(convertImage ${ITK_LIBRARIES})

add_executable(combinePair
			   combineStereoPair.cpp
			   RGBUtils.cpp)

TARGET_LINK_LIBRARIES(combinePair ${ITK_LIBRARIES})

add_executable(makeTiff
			   makeTiff.cpp)

TARGET_LINK_LIBRARIES(makeTiff ${ITK_LIBRARIES})

add_executable(vectorToImage
			   vectorToImage.cpp
			   ImageDataset.cpp
			   CVMatUtils.cpp)

TARGET_LINK_LIBRARIES(vectorToImage ${ITK_LIBRARIES} ${OPENCV_LIB}/libopencv_core.dylib)

add_executable(imageToVector
			   imageToVector.cpp
			   ImageDataset.cpp
			   CVMatUtils.cpp)

TARGET_LINK_LIBRARIES(imageToVector ${ITK_LIBRARIES} ${OPENCV_LIB}/libopencv_core.dylib)

add_executable(imagesToCSV
			   imagesToCSV.cpp
			   ImageDataset.cpp
			   CVMatUtils.cpp)

TARGET_LINK_LIBRARIES(imagesToCSV ${ITK_LIBRARIES} ${OPENCV_LIB}/libopencv_core.dylib)

