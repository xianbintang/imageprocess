//
// Created by xianb on 2017/11/28.
//


#include <opencv2/core/mat.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv/cv.hpp>
#include <iostream>
#include <TimeTracker.h>
#include "contour_detection.h"

static std::vector<float> rotate_image(const cv::Mat &src, cv::Mat &dst, cv::Point centerP, float degree)
{
    int width = src.cols;
    int height = src.rows;
    double angle = degree  * CV_PI / 180.; // 弧度
    double a = sin(angle), b = cos(angle);
    int width_rotate = int(height * fabs(a) + width * fabs(b));
    int height_rotate = int(width * fabs(a) + height * fabs(b));
    float map[6];
    cv::Mat map_matrix = cv::Mat(2, 3, CV_32F, map);
// 旋转中心
    CvPoint2D32f center = cvPoint2D32f(centerP.x, centerP.y);
    CvMat map_matrix2 = map_matrix;
    cv2DRotationMatrix(center, degree, 1.0, &map_matrix2);
    map[2] += (width_rotate - width) / 2;
    map[5] += (height_rotate - height) / 2;
    cv::warpAffine(src, dst, map_matrix, cv::Size(width_rotate, height_rotate), 1, 0, 0);
    std::vector<float> rotate_matrix;
    rotate_matrix.push_back(map[0]);
    rotate_matrix.push_back(map[1]);
    rotate_matrix.push_back(map[2]);
    rotate_matrix.push_back(map[3]);
    rotate_matrix.push_back(map[4]);
    rotate_matrix.push_back(map[5]);
    return rotate_matrix;
}

int cutout_template_image(const cv::Mat &template_image, std::vector<cv::Point> rect, cv::Mat &interesting_template)
{
    auto degree = 360 - cv::fastAtan2(rect[1].y - rect[2].y, rect[2].x - rect[1].x);
    cv::Mat img_rotate;
    int width = template_image.cols;
    int height = template_image.rows;
    std::vector<float> rotate_matrix;
    rotate_matrix = rotate_image(template_image, img_rotate, cv::Point(width / 2, height / 2), degree);
    std::cout << degree << std::endl;

    int plx = rect[0].x, ply = rect[0].y, prx = rect[3].x, pry = rect[3].y;
    int plxb = rect[1].x, plyb = rect[1].y, prxb = rect[2].x, pryb = rect[2].y;

    rect[0].x = rotate_matrix[0] * plx + rotate_matrix[1] * ply + rotate_matrix[2];
    rect[0].y = rotate_matrix[3] * plx + rotate_matrix[4] * ply + rotate_matrix[5];

    rect[1].x = rotate_matrix[0] * plxb + rotate_matrix[1] * plyb + rotate_matrix[2];
    rect[1].y = rotate_matrix[3] * plxb + rotate_matrix[4] * plyb + rotate_matrix[5];

    rect[2].x = rotate_matrix[0] * prxb + rotate_matrix[1] * pryb + rotate_matrix[2];
    rect[2].y = rotate_matrix[3] * prxb + rotate_matrix[4] * pryb + rotate_matrix[5];

    rect[3].x = rotate_matrix[0] * prx + rotate_matrix[1] * pry + rotate_matrix[2];
    rect[3].y = rotate_matrix[3] * prx + rotate_matrix[4] * pry + rotate_matrix[5];


//    cv::circle(img_rotate, rect[0], 1, cv::Scalar(255,255,255));
//    cv::circle(img_rotate, rect[1], 2, cv::Scalar(255,255,255));
//    cv::circle(img_rotate, rect[2], 3, cv::Scalar(255,255,255));
//    cv::circle(img_rotate, rect[3], 4, cv::Scalar(255,255,255));

    int rect_width = rect[3].x - rect[0].x;
    int rect_height = rect[1].y - rect[0].y;
    interesting_template = img_rotate(cv::Rect(rect[0].x, rect[0].y, rect_width, rect_height));

//    cv::rectangle(img_rotate, rect[0], rect[2], cv::Scalar(255,255,255));
//    cv::imshow("eh" ,interesting_template);
    cvWaitKey(0);
}

/*
 *  输入图像是二值化图像
 * */
static int get_center_numof_contour(const cv::Mat src, cv::Point &center, unsigned int &numofcontour)
{
    double m00, m10, m01;
    auto moments = cv::moments(src, true);
    m10 = moments.m10;
    m01 = moments.m01;
    m00 = moments.m00;
    if (m00 == 0) {
        return -1;
    } else {
        center.x = static_cast<int>(m10/m00);
        center.y = static_cast<int>(m01/m00);
    }
    numofcontour = m00;
}


/*
 * rect是相对640*480图片的坐标
 * */
int do_create_template(TemplateStruct &tpl, const cv::Mat &src, double low_threshold, double high_threshold)
{
    int s32Ret = 0;
    cv::Mat gx;                //Matrix to store X derivative
    cv::Mat gy;                //Matrix to store Y derivative
    // set width and height
    tpl.modelHeight = src.rows;    //Save Template height
    tpl.modelWidth = src.cols;    //Save Template width

    tpl.noOfCordinates = 0;    //initialize
    tpl.cordinates = new cv::Point[tpl.modelWidth * tpl.modelHeight];

    tpl.edgeDerivativeX = new double[tpl.modelWidth * tpl.modelHeight]; //Allocate memory for edge X derivative for selected points
    tpl.edgeDerivativeY = new double[tpl.modelWidth * tpl.modelHeight]; ////Allocate memory for edge Y derivative for selected points

    cv::Sobel(src, gx, CV_16S, 1,0,3);        //gradient in X direction
    cv::Sobel(src, gy, CV_16S, 0,1,3);        //gradient in Y direction

    cv::Mat binaryContour;
    cv::Canny(src, binaryContour, low_threshold, high_threshold);

    const short *_sdx;
    const short *_sdy;
    double fdx, fdy;

    int RSum = 0, CSum = 0;
    int curX, curY;
    int flag = 1;

    for (int i = 0; i < tpl.modelHeight; i++) {
        for (int j = 0; j < tpl.modelWidth; j++) {
            short vx = gx.at<short>(i,j);
            short vy = gy.at<short>(i,j);
            unsigned char U8 = binaryContour.at<uchar>(i, j);
            cv::Point p;

            fdx = vx;
            fdy = vy;

            p.x = j;
            p.y = i;
            if (U8) {
                /* 如果梯度都为零，那么不需要计算，因为分数不会有贡献 */
                if (fdx != 0 || fdy != 0) {
                    /* 坐标变换到外接矩形左上角为(0, 0) */
                    RSum = RSum + j;
                    CSum = CSum + i;    // Row sum and column sum for center of gravity
                    tpl.cordinates[tpl.noOfCordinates].x = j;
                    tpl.cordinates[tpl.noOfCordinates].y = i;

                    /* TODO 可以修改成使用查找表的形式 */
                    double vector_length = sqrt(fdx * fdx + fdy * fdy);
                    if (fabs(vector_length - 0.) < 0.00001) {
//                        printf(".............................................\n");
                    }
                    tpl.edgeDerivativeX[tpl.noOfCordinates] =
                            fdx / vector_length;
                    tpl.edgeDerivativeY[tpl.noOfCordinates] =
                            fdy / vector_length;
                    tpl.noOfCordinates++;
                }
            }
        }
    }

    if (tpl.noOfCordinates == 0) {
//        printf(".........................");
        tpl.centerOfGravity.x = tpl.modelWidth / 2;
        tpl.centerOfGravity.y = tpl.modelHeight / 2;
    } else {
        tpl.centerOfGravity.x = RSum / tpl.noOfCordinates;    // center of gravity
        tpl.centerOfGravity.y = CSum / tpl.noOfCordinates;    // center of gravity
    }

    // change coordinates to reflect center of gravity
    /* 将重心变换到坐标原点 */
    int m;
    for (m = 0; m < tpl.noOfCordinates; m++) {
        /*int temp;

        temp = tpl.cordinates[m].x;
        tpl.cordinates[m].x = temp - tpl.centerOfGravity.x;
        temp = tpl.cordinates[m].y;
        tpl.cordinates[m].y = temp - tpl.centerOfGravity.y;*/
        tpl.cordinates[m].x -= tpl.centerOfGravity.x;
        tpl.cordinates[m].y -= tpl.centerOfGravity.y;
    }

    tpl.modelDefined = true;

    return 1;
}

static void draw_template(cv::Mat src, const TemplateStruct &tpl)
{
    for (int i = 0; i < tpl.noOfCordinates; ++i) {
        cv::circle(src, cv::Point(tpl.cordinates[i].x + tpl.centerOfGravity.x, tpl.cordinates[i].y + tpl.centerOfGravity.y), 1, cv::Scalar(255,255,255));
    }
    cv::imshow("hehe", src);
    cvWaitKey(0);
}


int free_tpls(TemplateStruct tpls[][MAX_DEGREE])
{

}


// 为koyo_tool_contour_parameter对应的工具建立模板
int create_template(const cv::Mat &src, Koyo_Tool_Contour_Parameter koyo_tool_contour_parameter)
{
    cv::Mat pyramid_templates[MAX_NUM_PYRAMID];
    pyramid_templates[0] = src;

    UINT8 sensitity_threshold_low, sensitity_threshold_high;
    if (koyo_tool_contour_parameter.sensitivity == CONTOUR_ACCURACY_LOW) {
        sensitity_threshold_low = 30;
        sensitity_threshold_high = 150;
    } else if (koyo_tool_contour_parameter.sensitivity = CONTOUR_ACCURACY_MEDIUM) {
        sensitity_threshold_low = 30;
        sensitity_threshold_high = 150;
    } else if (koyo_tool_contour_parameter.sensitivity = CONTOUR_ACCURACY_HIGH) {
        sensitity_threshold_low = 30;
        sensitity_threshold_high = 150;
    }

    // 建立各层金字塔, 并确定最佳金字塔层数
    int optimal_pyr_level = 0;
    for (int i = 0; i < MAX_NUM_PYRAMID - 1; ++i) {
        cv::pyrDown(pyramid_templates[i], pyramid_templates[i+1]);

    }

    // 对每一层使用高斯滤波处理
    for (auto iter = std::begin(pyramid_templates); iter != std::end(pyramid_templates); ++iter) {
        cv::Mat after_gaus;
        cv::GaussianBlur(*iter, after_gaus, cv::Size(3,3),3,3);
//        *iter = after_gaus;
    }

    // 图像的质心
    std::vector<cv::Point> centers;
    for (int i = 0; i < MAX_NUM_PYRAMID; ++i) {
        cv::Mat cannyResult;
        cv::Canny(pyramid_templates[i], cannyResult, sensitity_threshold_high, sensitity_threshold_low);
        cv::Point center;
        unsigned int num_of_contour;
        get_center_numof_contour(cannyResult, center, num_of_contour);
        centers.push_back(center);
        std::cout << num_of_contour << std::endl;
        if (num_of_contour <= MIN_CONTOUR_PYRA) {
            break;
        }
        ++optimal_pyr_level;
//        cv::imshow(std::string("pyr") + std::string(1, i - '0') , cannyResult);
//        cvWaitKey(0);
    }
    std::cout << optimal_pyr_level << std::endl;

    // 确定角度步长
    double angle_steps[MAX_NUM_PYRAMID] = {1,1,1,1,1};

    // 对每层每个角度建立模板
    // tpls中的内存是动态分配的, 在建立完模板后需要释放所有的内存
    TemplateStruct tpls[optimal_pyr_level][MAX_DEGREE];
    TimeTracker tt;
    tt.start();
    for (int i = 0; i < optimal_pyr_level; ++i) {
        for (int j = 0; j < MAX_DEGREE; j += angle_steps[i]) {
            cv::Mat rotated_image;
            rotate_image(pyramid_templates[i], rotated_image, centers[i], j);
            do_create_template(tpls[i][j], rotated_image, sensitity_threshold_low, sensitity_threshold_high);
            draw_template(rotated_image, tpls[i][j]);
//            cv::imshow(std::string("pyr") + std::string(1, i - '0'), rotated_image);
//            cvWaitKey(0);
        }
    }
    tt.stop();
    std::cout << tt.duration() << "ms" << std::endl;

    //建立完模板需要将模板发送给客户端，需要发送的就是tpls这个数据结构
}
