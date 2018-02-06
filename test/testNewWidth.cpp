//
// Created by xianb on 2018/1/25.
//
#include <iostream>
#include <opencv2/core/mat.hpp>
#include <opencv2/opencv.hpp>

void saveMat(cv::Mat mat, const char *path) {
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

void saveMatf(cv::Mat mat, const char *path) {
    FILE *fp = fopen(path, "w");
    int i,j;
    for (i = 0; i < mat.rows; ++i) {
        for (j = 0; j < mat.cols; ++j) {
//            fprintf(fp, "%d ", (mat.ptr + i * mat.step)[j]);
            fprintf(fp, "%d ", mat.at<short>(i, j));
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}
int main(int argc, char **argv)
{
    cv::Mat img= cv::imread(argv[1], 0);
    cv::Mat img_gray;
    cv::GaussianBlur(img, img_gray, cv::Size(5, 5), 0);
    cv::Mat sdx, sdy;
    cv::Sobel(img_gray, sdx, CV_16S, 1, 0, 3);
    cv::Sobel(img_gray, sdy, CV_16S, 0, 1, 3);
    cv::Mat result;
    result.create(img_gray.rows, img_gray.cols, CV_16SC1);
    int i,j;
    auto color_img = cv::imread(argv[1], -1);
    for (i = 0; i < result.rows; ++i) {
        for (j = 0; j < result.cols; ++j) {
//            fprintf(fp, "%d ", (mat.ptr + i * mat.step)[j]);
//            std::cout << atan2f(sdy.at<short>(i, j), sdx.at<short>(i, j)) << std::endl;;
            result.at<short>(i, j) = (atan2f(sdy.at<short>(i, j), sdx.at<short>(i, j)) * 180.0 / 3.1415926 + 180.0) / 45.0;
            auto val = result.at<short>(i, j);
            std::cout << val << std::endl;
            cv::Scalar clr;
            if(val == 0) clr = cv::Scalar(255,255,255);
            if(val == 1) clr = cv::Scalar(255,0,0);
            if(val == 2) clr = cv::Scalar(0,255,0);
            if(val == 3) clr = cv::Scalar(0,0,255);
            if(val == 4) clr = cv::Scalar(0,255,255);
            if(val == 5) clr = cv::Scalar(255,0,255);
            if(val == 6) clr = cv::Scalar(255,255,0);
            if(val == 7) clr = cv::Scalar(0,0,0);
            cv::circle(color_img, {j, i}, 1, clr);
        }
    }


    cv::imshow("haha", img_gray);
    cv::imshow("hehe", color_img);
    cv::waitKey(-1);

    saveMatf(result, "data//widthang.txt");
    return 0;
}