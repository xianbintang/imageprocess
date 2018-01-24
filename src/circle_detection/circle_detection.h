//
// Created by xianb on 2018/1/24.
//

#ifndef KOYO_CIRCLE_DETECTION_H
#define KOYO_CIRCLE_DETECTION_H


#include <opencv2/core/types.hpp>
#include <opencv2/core/mat.hpp>

using UINT16 = unsigned short;
using INT8 = signed char;
using UINT8 = unsigned char;
using UINT32 = unsigned int;
const int MIN_CONTOUR_PYRA = 80;
const int WIDTH = 640;
const int HEIGHT = 480;

typedef struct Circle {
    cv::Point2f  center;
    int radius;
    double score;
}Circle ;
typedef struct Region{
    cv::Point2f  center;
    int radius;
}Region;

typedef struct _CIRCLE_RUNTIME_PARAM {
    cv::Mat img_data;
    cv::Mat img_gray;
    cv::Mat dist;
    cv::Mat sdx;
    cv::Mat sdy;
    cv::Mat mag;
    cv::Mat edge;
    cv::Mat h_acc;
    cv::Mat coins;
    Region region;
} Circle_runtime_param;

typedef struct KOYO_TOOL_CIRCLE_PARAMETER_
{
    INT8   tool_name[32];         //工具名称

    UINT16 circle_x;
    UINT16 circle_y;
    UINT16 circle_r;
    UINT16 detect_r;
    UINT8  sensitivity;
    UINT8  detect_method;

    UINT16 top_threshold;
    UINT16 bot_threshold;
} Koyo_Tool_Circle_Parameter;


int circle_detection_config(const UINT8 *yuv, Circle circles[]);

#endif //KOYO_CIRCLE_DETECTION_H
