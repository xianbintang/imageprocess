//
// Created by xianb on 2017/6/30.
//

#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <iostream>
#include <stdio.h>
#include "TimeTracker.h"

using namespace std;
using namespace cv;

/// 全局变量
Mat img; Mat templ; Mat result;
char* image_window = "Source Image";
char* result_window = "Result window";

int match_method;
int max_Trackbar = 5;

/// 函数声明
void MatchingMethod( int, void* );

/** @主函数 */
int main( int argc, char** argv )
{
    /// 载入原图像和模板块
    img = imread( argv[1], 1 );
    templ = imread( argv[2], 1 );

    /// 创建窗口
    namedWindow( image_window, CV_WINDOW_AUTOSIZE );
    namedWindow( result_window, CV_WINDOW_AUTOSIZE );

    /// 创建滑动条
    char* trackbar_label = "M\n0\n1\n2\n3\n4\n5";
    createTrackbar( trackbar_label, image_window, &match_method, max_Trackbar, MatchingMethod );

    MatchingMethod( 0, 0 );

    waitKey(0);
    return 0;
}

/**
 * @函数 MatchingMethod
 * @简单的滑动条回调函数
 */
void MatchingMethod( int, void* )
{
    /// 将被显示的原图像
    Mat img_display;
    img.copyTo( img_display );

    /// 创建输出结果的矩阵
    int result_cols =  img.cols - templ.cols + 1;
    int result_rows = img.rows - templ.rows + 1;

    result.create( result_cols, result_rows, CV_32FC1 );

    /// 进行匹配和标准化

    TimeTracker timeTracker;
    timeTracker.start();
    matchTemplate( img, templ, result, 1);
    timeTracker.stop();
    std::cout << timeTracker.duration() << std::endl;


    timeTracker.start();
    matchTemplate( img, templ, result, 3);
    timeTracker.stop();
    std::cout << timeTracker.duration() << std::endl;


    timeTracker.start();
    matchTemplate( img, templ, result, 4);
    timeTracker.stop();
    std::cout << timeTracker.duration() << std::endl;

    normalize( result, result, 0, 1, NORM_MINMAX, -1, Mat() );

    /// 通过函数 minMaxLoc 定位最匹配的位置
    double minVal; double maxVal; Point minLoc; Point maxLoc;
    Point matchLoc;

    minMaxLoc( result, &minVal, &maxVal, &minLoc, &maxLoc, Mat() );

    /// 对于方法 SQDIFF 和 SQDIFF_NORMED, 越小的数值代表更高的匹配结果. 而对于其他方法, 数值越大匹配越好
    if( match_method  == CV_TM_SQDIFF || match_method == CV_TM_SQDIFF_NORMED )
    { matchLoc = minLoc; }
    else
    { matchLoc = maxLoc; }

    /// 让我看看您的最终结果
    rectangle( img_display, matchLoc, Point( matchLoc.x + templ.cols , matchLoc.y + templ.rows ), Scalar::all(0), 2, 8, 0 );
    rectangle( result, matchLoc, Point( matchLoc.x + templ.cols , matchLoc.y + templ.rows ), Scalar::all(0), 2, 8, 0 );

    imshow( image_window, img_display );
    imshow( result_window, result );

    return;
}

