cmake_minimum_required(VERSION 3.6)
project(KOYO)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_BUILD_TYPE "Debug")


if( CMAKE_BUILD_TYPE STREQUAL "Debug" )
    set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -std=c++14 -g -Wall -Wno-unused-variable -pthread")
else( CMAKE_BUILD_TYPE STREQUAL "Debug" )
    set(CMAKE_BUILD_TYPE "Release")
    set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -std=c++14 -O2 -pthread -fopenmp")
endif( CMAKE_BUILD_TYPE STREQUAL "Debug" )
message("Build Type:" ${CMAKE_BUILD_TYPE}    ${CMAKE_CXX_FLAGS})



set(OpenCV_DIR E://ClionProjects//opencv_install)
find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
include_directories(include)


set(SOURCE_FILES src/KOYOMain.cpp)

add_executable(TestOpenCV ${SOURCE_FILES})
target_link_libraries(TestOpenCV ${OpenCV_LIBS})


add_executable(KOYOMain ${SOURCE_FILES} utils/utils.cpp src/GetUnions.cpp)
target_link_libraries(KOYOMain ${OpenCV_LIBS})

add_executable(UnionTest src/GetUnions.cpp utils/utils.cpp )
target_link_libraries(UnionTest ${OpenCV_LIBS})

add_executable(HuMatch src/Hu.cpp)
target_link_libraries(HuMatch ${OpenCV_LIBS})

add_executable(GeoMatch src/GeoMatchTest.cpp src/GeoMatch.cpp src/CommandParser.cpp src/RotateImage.cpp utils/utils.cpp)
target_link_libraries(GeoMatch ${OpenCV_LIBS})

add_executable(FeaturesMatch src/FeaturesMatch.cpp)
target_link_libraries(FeaturesMatch ${OpenCV_LIBS})

add_executable(RotateImage1 src/RotateImage.cpp src/IMatImpl.cpp src/sobel.cpp)
target_link_libraries(RotateImage1 ${OpenCV_LIBS})



add_executable(TemplateMatch src/RotateImage.cpp src/TemplateMatchTest.cpp src/TemplateMatch.cpp src/TemplateMatchTest.cpp src/IMatImpl.cpp src/sobel.cpp)
target_link_libraries(TemplateMatch ${OpenCV_LIBS})

add_executable(test1 test/test.cpp)
target_link_libraries(test1 ${OpenCV_LIBS})

add_executable(testpyr test/testPyDown.cpp)
target_link_libraries(testpyr ${OpenCV_LIBS})

add_executable(integral Integral/IntegralMain.cpp)
target_link_libraries(integral ${OpenCV_LIBS})

add_executable(fastLocalization src/fastLocalization.cpp)
target_link_libraries(fastLocalization ${OpenCV_LIBS})

add_executable(testTemplateMatch test/testTemplateMatch.cpp)
target_link_libraries(testTemplateMatch ${OpenCV_LIBS})

add_executable(colorDetection src/color_detection/color_detection.cpp)
target_link_libraries(colorDetection ${OpenCV_LIBS})

add_executable(circle_detection src/circle_detection/circle_detection.cpp)
target_link_libraries(circle_detection ${OpenCV_LIBS})

add_executable(getPoint test/getPoint.cpp)
target_link_libraries(getPoint ${OpenCV_LIBS})
