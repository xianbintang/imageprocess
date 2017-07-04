//
// Created by xianb on 2017/6/30.
//
#include <opencv/cv.h>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include "TimeTracker.h"
#include <opencv/cv.hpp>
#include <iostream>

using namespace cv;
using namespace std;

void dohist(Mat src);
int main(int argc, char **argv)
{
    TimeTracker timeTracker;
    timeTracker.start();
    int j;
    for (int i = 0; i < 20000000; ++i) {
        if (i > 0)
            j = i + j;
    }
    timeTracker.stop();
    std::cout << "time: " << timeTracker.duration() << "i: " << j << std::endl;


    Mat source = imread(argv[1]);
    cvtColor(source, source, CV_RGB2GRAY);
    threshold(source, source, 100, 255, CV_THRESH_BINARY);
//    dohist(source);
    imshow("test", source);
    cvWaitKey(-1);
    return 0;
}

void dohist(Mat src) {
    Mat image = src;
    MatND hist;       // 在cv中用CvHistogram *hist = cvCreateHist
    int dims = 1;
    float hranges[] = {0, 255};
    const float *ranges[] = {hranges};   // 这里需要为const类型
    int size = 256;
    int channels = 0;
    // 计算图像的直方图
    calcHist(&image, 1, &channels, Mat(), hist, dims, &size, ranges);    // cv 中是cvCalcHist
    int scale = 1;
    Mat imageShow(size * scale, size, CV_8U, Scalar(0));
    // 获取最大值和最小值
    double minVal = 0;
    double maxVal = 0;
    minMaxLoc(hist, &minVal, &maxVal, 0, 0);  //  cv中用的是cvGetMinMaxHistValue
    std::cout << "max: " << maxVal << std::endl;
    std::cout << "min: " << minVal << std::endl;
    //显示直方图的图像
    int hpt = saturate_cast<int>(0.9 * size);

    for (int i = 0; i < 256; i++) {
        float value = hist.at<float>(i);           //   注意hist中是float类型    cv中用cvQueryHistValue_1D
        int realValue = saturate_cast<int>(value * hpt / maxVal);
        //line(imageShow, Point(i, size - 1), Point(i, size - realValue), Scalar(0));
        rectangle(imageShow, Point(i * scale, size - 1), Point((i + 1) * scale - 1, size - realValue), Scalar(255));
    }
    namedWindow("showImage");
    imshow("showImage", imageShow);
    waitKey(0);
}

