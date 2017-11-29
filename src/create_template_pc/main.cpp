//
// Created by xianb on 2017/11/28.
//


#include <opencv2/core/mat.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv/cv.hpp>
#include <windef.h>
#include "contour_detection.h"

void init_contour_parameter(Koyo_Tool_Contour_Parameter &koyo_tool_contour_parameter)
{
    koyo_tool_contour_parameter.algo_strategy = 1;

    koyo_tool_contour_parameter.detect_rect_x0 = 196;
    koyo_tool_contour_parameter.detect_rect_y0 = 130;

    koyo_tool_contour_parameter.detect_rect_x1 = 197;
    koyo_tool_contour_parameter.detect_rect_y1 = 427;

    koyo_tool_contour_parameter.detect_rect_x2 = 483;
    koyo_tool_contour_parameter.detect_rect_y2 = 428;

    koyo_tool_contour_parameter.detect_rect_x3 = 483;
    koyo_tool_contour_parameter.detect_rect_y3 = 130;

    koyo_tool_contour_parameter.sensitivity = CONTOUR_ACCURACY_MEDIUM;
    koyo_tool_contour_parameter.angle_range = 180;
}


int main(int argc, char **argv)
{
    cv::Mat template_image = cv::imread(argv[1], 0);


    Koyo_Tool_Contour_Parameter koyo_tool_contour_parameter;
    init_contour_parameter(koyo_tool_contour_parameter);

    cv::Mat template_roi;
    std::vector<cv::Point> rect =  {
            {koyo_tool_contour_parameter.detect_rect_x0, koyo_tool_contour_parameter.detect_rect_y0},
            {koyo_tool_contour_parameter.detect_rect_x1, koyo_tool_contour_parameter.detect_rect_y1},
            {koyo_tool_contour_parameter.detect_rect_x2, koyo_tool_contour_parameter.detect_rect_y2},
            {koyo_tool_contour_parameter.detect_rect_x3, koyo_tool_contour_parameter.detect_rect_y3},
    };
    cutout_template_image(template_image, rect, template_roi);
    create_template(template_roi, koyo_tool_contour_parameter);
//    cv::imshow("eh" ,template_roi);
    cvWaitKey(0);
    return 0;
}


