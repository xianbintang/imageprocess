//
// Created by xianb on 2017/11/28.
//


#include <opencv2/core/mat.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv/cv.hpp>
#include "contour_detection.h"

int create_template(const cv::Mat &src)
{
    cv::Mat pyramid_templates[MAX_NUM_PYRAMID];
    pyramid_templates[0] = src;

    // 建立各层金字塔, 并确定最佳金字塔层数
    for (int i = 0; i < MAX_NUM_PYRAMID - 1; ++i) {
        cv::pyrDown(pyramid_templates[i], pyramid_templates[i+1]);
    }

    // 确定角度步长
    double angle_steps[MAX_NUM_PYRAMID] = {1,1,1,1,1};

    // 对每层每个角度建立模板
    for (int i = 0; i < MAX_NUM_PYRAMID; ++i) {
        for (int j = 0; j < MAX_DEGREE; j += angle_steps[i]) {
//            cv::
        }
    }
    cvWaitKey(0);
}
