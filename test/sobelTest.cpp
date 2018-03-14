//
// Created by xianb on 2018/1/19.
//


#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv/cv.hpp>
#include <iostream>

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
    cv::Mat img = cv::imread(argv[1], 0);
    cv::Mat sdx;
    cv::Sobel(img, sdx, CV_16S, 1, 0, 3);
    cv::imshow("haha", img);
    cv::imshow("hehe", sdx);
    std::cout << atan2f(0,0) << std::endl;;
    std::cout << atan2f(1,0) << std::endl;;
    std::cout << atan2f(0,1) << std::endl;;
    cv::waitKey(0);
    saveMatf(sdx, "data//sdx.txt");
    return 0;
}