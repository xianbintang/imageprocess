//
// Created by xianb on 2017/11/28.
//


#include <opencv2/core/mat.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv/cv.hpp>
#include <windef.h>
#include "contour_detection.h"

int main(int argc, char **argv)
{
    cv::Mat template_image = cv::imread(argv[1], 0);
//    cv::imshow("eh" ,template_image);

    double degree = 9;
    double angle = degree  * CV_PI / 180.; // 弧度
    double a = sin(angle), b = cos(angle);
    int width = template_image.cols;
    int height = template_image.rows;
    int width_rotate = int(height * fabs(a) + width * fabs(b));
    int height_rotate = int(width * fabs(a) + height * fabs(b));
//旋转数组map
// [ m0  m1  m2 ] ===>  [ A11  A12   b1 ]
// [ m3  m4  m5 ] ===>  [ A21  A22   b2 ]
    float map[6];
    cv::Mat map_matrix = cv::Mat(2, 3, CV_32F, map);
// 旋转中心
    CvPoint2D32f center = cvPoint2D32f(width / 2, height / 2);
    CvMat map_matrix2 = map_matrix;
    cv2DRotationMatrix(center, degree, 1.0, &map_matrix2);
    map[2] += (width_rotate - width) / 2;
    map[5] += (height_rotate - height) / 2;
    cv::Mat img_rotate;
    cv::warpAffine(template_image, img_rotate, map_matrix, cv::Size(width_rotate, height_rotate), 1, 0, 0);


    create_template(template_image);
//    cvWaitKey(0);
    return 0;
}


