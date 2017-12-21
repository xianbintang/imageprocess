//
// Created by xianb on 2017/12/20.
//


#include <opencv2/core/mat.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <iostream>
#include <opencv/cv.hpp>
#include <basetsd.h>
#include "../src/create_template_pc/contour_detection.h"

static void saveMat(cv::Mat mat, const char *path) {
    FILE *fp = fopen(path, "w");
    int i,j;
    for (i = 0; i < mat.rows; ++i) {
        for (j = 0; j < mat.cols; ++j) {
//            fprintf(fp, "%d ", (mat.ptr + i * mat.step)[j]);
            fprintf(fp, "%d ", mat.at<uchar>(i, j));
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

static cv::Mat get_y_from_yuv(const UINT8 *yuv, const UINT16 width, const UINT16 height)
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

int main(int argc, char ** argv)
{
    cv::Mat src, dst;
    std::vector<std::vector<float>> k
            {
//                    {1, 4, 7, 4, 1},
//                    {4, 16, 26, 16, 4},
//                    {7, 26, 41, 26, 7},
//                    {4, 16, 26, 16, 4},
//                    {1, 4, 7, 4, 1}

                    {1,2,1},
                    {2,4,2},
                    {1,2,1},

            };

    cv::Mat kernel(3,3,CV_32F);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            kernel.at<float>(i, j) = k[i][j]/16.0;
            std::cout << kernel.at<float>(i, j) << std::endl;
        }
    }

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

    template_image = get_y_from_yuv(buf, WIDTH, HEIGHT);
//    cv::filter2D(template_image, dst, -1, kernel, cv::Point(-1, -1), 0, cv::BORDER_DEFAULT);
    cv::filter2D(template_image, dst, -1, kernel, cv::Point(-1, -1), 0);
//    cv::imshow("hehe", dst);
//    cvWaitKey(0);
    saveMat(dst, "dst.txt");
}


