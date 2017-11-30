//
// Created by xianb on 2017/11/30.
//

#include <iostream>
#include <basetsd.h>
#include <opencv/cv.hpp>
#include "../src/create_template_pc/contour_detection.h"

cv::Mat read_y_from_yuv(UINT8 *yuv, const UINT16 width, const UINT16 height)
{
    cv::Mat gray;
    if (!yuv) {
        std::cout << "ERROR: yuv NULL pointer" << std::endl;
    } else {
        gray.create(height, width, CV_8UC1);
        memcpy(gray.data, yuv, width * height * sizeof(unsigned char)); // 只读取y分量上的数据
    }

    return gray;
}

int main(int argc, char **argv)
{
    FILE *yuv_file = fopen(argv[1], "rb+");
    UINT8 *buf = new UINT8[WIDTH * HEIGHT];
    fread(buf, WIDTH * HEIGHT, 1, yuv_file);
    auto gray = get_y_from_yuv(buf, WIDTH, HEIGHT);
    cv::imshow("gray", gray);
    cvWaitKey(0);
    return 0;
}