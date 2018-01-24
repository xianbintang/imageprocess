//
// Created by xianb on 2018/1/24.
//
#include <opencv2/core/mat.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv/cv.hpp>
#include <windef.h>
#include <fstream>

#include "circle_detection.h"

int main( int argc, char** argv )
{

    char* imageName = argv[1];

    std::string filename(argv[1]);
    cv::Mat template_image;
    UINT8 *buf = nullptr;
    if(filename.substr(filename.size() - 3, 3) == "yuv") {
        FILE *yuv_file = fopen(filename.c_str(), "rb+");
        buf = new UINT8[WIDTH * HEIGHT];
        fread(buf, WIDTH * HEIGHT, 1, yuv_file);
    } else {
        template_image = cv::imread(argv[1], 0);
        buf = template_image.data;
    }
    Circle circles[20];
    circle_detection_config(buf, circles);
    return 0;

}