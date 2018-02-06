//
// Created by xianb on 2018/1/24.
//
#include <opencv2/core/mat.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv/cv.hpp>
#include <windef.h>
#include <fstream>
#include <c++/iostream>

#include "circle_detection.h"

int main( int argc, char** argv )
{

    std::string filename(argv[1]);
    cv::Mat template_image, rgb_img;
    UINT8 *buf = nullptr;
    if(filename.substr(filename.size() - 3, 3) == "yuv") {
        FILE *yuv_file = fopen(filename.c_str(), "rb+");
        buf = new UINT8[WIDTH * HEIGHT];
        fread(buf, WIDTH * HEIGHT, 1, yuv_file);
    } else {
        template_image = cv::imread(argv[1], 0);
        buf = template_image.data;
        rgb_img = cv::imread(argv[1], 1);
    }
    std::vector<Result_Circle> circles;
    circle_detection_config(buf, circles);

    for (auto iter = circles.begin(); iter != circles.end(); ++iter) {
        int lineThickness = 2;
        int lineType = 10;
        int shift = 0;
        int xCoord = iter->x;
        int yCoord = iter->y;
        int radius = iter->radius;
        cv::Point2f center(xCoord, yCoord);
        if (radius != 0) {
            cv::circle(rgb_img, center, radius - 1, cv::Scalar(0, 0, 255), lineThickness, lineType, shift);
            std::cout << "center: " << center <<  "radius: " << radius << " count: " << std::endl;
        }
    }
    cv::imshow("hahah", rgb_img);
    cv::waitKey(-1);

    return 0;

}