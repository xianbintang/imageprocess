//
// Created by xianb on 2017/11/28.
//


#include <opencv2/core/mat.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv/cv.hpp>
#include <iostream>
#include <TimeTracker.h>
#include <queue>
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

static int rotate_rect(std::vector<cv::Point> &rect, const std::vector<float> rotate_matrix)
{
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
    rotate_rect(rect, rotate_matrix);

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

// 计算点到直线的距离
static double dist_to_line(cv::Point pt2, cv::Point pt1, cv::Point target)
{
    // 根据中心点与直线的距离 排除干扰直线
    // 点(x0,y0)到直线Ax+By+C=0的距离为d = (A*x0+B*y0+C)/sqrt(A^2+B^2)
    double A, B, C,dis;
    // 化简两点式为一般式
    // 两点式公式为(y - y1)/(x - x1) = (y2 - y1)/ (x2 - x1)
    // 化简为一般式为(y2 - y1)x + (x1 - x2)y + (x2y1 - x1y2) = 0
    // A = y2 - y1
    // B = x1 - x2
    // C = x2y1 - x1y2
    A = pt2.y - pt1.y;
    B = pt1.x - pt2.x;
    C = pt2.x * pt1.y - pt1.x * pt2.y;

    double x,y;
    x = target.x;
    y = target.y;
    // 距离公式为d = |A*x0 + B*y0 + C|/√(A^2 + B^2)
    dis = abs(A * x + B * y + C) / sqrt(A * A + B * B);
    return dis;
}

// 如果target点到rect构成的四条直线的距离都大于min_dist，则返回true，否则返回false。若返回false，说明这个点距离至少一条直线的距离小于min_dist
static bool dist_to_lines_less_than(const std::vector<cv::Point> &rect, cv::Point target, double min_dist)
{
    return !((dist_to_line(rect[0], rect[1], target) < min_dist || dist_to_line(rect[1], rect[2], target) < min_dist || \
    dist_to_line(rect[2], rect[3], target) < min_dist || dist_to_line(rect[3], rect[0], target) < min_dist));
}

/*
 * rect是相对640*480图片的坐标
 * */
int do_create_template(TemplateStruct &tpl, const cv::Mat &src, double low_threshold,\
 double high_threshold, const std::vector<cv::Point> &rect)
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
            short fdx = gx.at<short>(i,j);
            short fdy = gy.at<short>(i,j);
            unsigned char U8 = binaryContour.at<uchar>(i, j);
            cv::Point p;
            p.x = j;
            p.y = i;
            // todo 最小距离是多少还需要斟酌，因为在最小分辨率情况下看到边框还是没有去除掉，在最小分辨情况下这个dist太小了。
            if (U8 && dist_to_lines_less_than(rect, p, (tpl.modelHeight + tpl.modelWidth) / 100.0)) {
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

// TODO 发送给客户端以后需要释放
int free_tpls(TemplateStruct tpls[][MAX_DEGREE])
{

}

/*
 *  @param src必须是二值化轮廓图
 *  @return 返回当前层上最佳的旋转步长
 * */
static double get_angle_step(const cv::Mat &src, cv::Point center)
{
    // 保留几个K，然后求平均值，用来排除外点的影响
    int K = 10;
    std::priority_queue<double> max_dist(K, -1);
    std::vector<cv::Point> points;
    // 在目前没有旋转的图片中不会出现因为旋转导致的白边，所以直接在全图搜索就行了，不用考虑别的
    for (int i = 0; i < src.rows; ++i) {
        for (int j = 0; j < src.cols; ++j) {
            if (src.at<uchar>(i, j)) {
                double dist = -1 * sqrt(pow((i - center.y),2) + pow((j - center.x),2));
                if(dist < max_dist.top()) {
                    max_dist.pop();
                    max_dist.push(dist);
                    points.push_back({i, j});
                }
            }
        }
    }
    std::cout << "max dist....." << std::endl;
    double average_max_dist = 0;
    int i = 0;
    while (!max_dist.empty()) {
        average_max_dist += -max_dist.top();
        max_dist.pop();
    }
    average_max_dist /= K;
    std::cout << "average_max_dist: " << average_max_dist << std::endl;

    std::cout <<"optimal angle step: " << acos(1 - 1 / (2 * average_max_dist * average_max_dist)) / CV_PI * 360 << " ~ " << acos(1 - 1 / (average_max_dist * average_max_dist)) / CV_PI * 360 << std::endl;

#ifdef _DEBUG_
    cv::Mat tmp = src;
    for (int k = points.size() - K; k < points.size(); ++k) {
        cv::circle(tmp, points[k], 20, cv::Scalar(255,255,255));
    }
    cv::imshow("max dist", tmp);
    cvWaitKey(0);
#endif
}

// 为koyo_tool_contour_parameter对应的工具建立模板
int create_template(const cv::Mat &src, Koyo_Tool_Contour_Parameter koyo_tool_contour_parameter)
{
    TimeTracker tt1;
    tt1.start();
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
    std::cout << "size of this level: cols " << pyramid_templates[0].cols << ", rows" << pyramid_templates[0].rows << std::endl;
    for (int i = 0; i < MAX_NUM_PYRAMID - 1; ++i) {
        cv::pyrDown(pyramid_templates[i], pyramid_templates[i+1]);
        std::cout << "size of this level: cols " << pyramid_templates[i+1].cols << ", rows" << pyramid_templates[i+1].rows << std::endl;
//        cv::imshow(std::string("pyr") + std::string(1, i - '0') , pyramid_templates[i+1]);
//        cvWaitKey(0);
    }

    // 对每一层使用高斯滤波处理
    for (auto iter = std::begin(pyramid_templates); iter != std::end(pyramid_templates); ++iter) {
        cv::Mat after_gaus;
        cv::GaussianBlur(*iter, after_gaus, cv::Size(3,3),3,3);
        *iter = after_gaus;
    }

    // 图像的质心
    std::vector<cv::Point> centers;
    std::vector<double> angle_steps = {1,1,1,1,1};
    for (int i = 0; i < MAX_NUM_PYRAMID; ++i) {
        cv::Mat cannyResult;
        cv::Canny(pyramid_templates[i], cannyResult, sensitity_threshold_high, sensitity_threshold_low);

        cv::Point center;
        unsigned int num_of_contour;
        get_center_numof_contour(cannyResult, center, num_of_contour);
        centers.push_back(center);

        // 确定角度步长, 使用Canny的轮廓图来计算最远点
        get_angle_step(cannyResult, center);

//        std::cout << num_of_contour << std::endl;
        if (num_of_contour <= MIN_CONTOUR_PYRA) {
            break;
        }
        ++optimal_pyr_level;
//        cv::imshow(std::string("pyr") + std::string(1, i - '0') , cannyResult);
//        cvWaitKey(0);
    }
    std::cout << "optimal level: " <<  optimal_pyr_level << std::endl;

    tt1.stop();
    std::cout << "first half: " << tt1.duration() << std::endl;


    // 对每层每个角度建立模板
    // tpls中的内存是动态分配的, 在建立完模板后需要释放所有的内存
    TemplateStruct tpls[optimal_pyr_level][MAX_DEGREE];
    TimeTracker tt;
    tt.start();
    for (int i = 0; i < optimal_pyr_level; ++i) {
        std::vector<cv::Point> cur_rect = {{0,0}, {0, pyramid_templates[i].rows - 1}, {pyramid_templates[i].cols - 1, pyramid_templates[i].rows - 1}, {pyramid_templates[i].cols - 1, 0}};
        for (int j = 0; j < MAX_DEGREE; j += angle_steps[i]) {
            auto rect = cur_rect;
            cv::Mat rotated_image;
            auto rotate_matrix = rotate_image(pyramid_templates[i], rotated_image, centers[i], j);
            rotate_rect(rect, rotate_matrix);
            do_create_template(tpls[i][j], rotated_image, sensitity_threshold_low, sensitity_threshold_high, rect);
            draw_template(rotated_image, tpls[i][j]);
//            cv::imshow(std::string("pyr") + std::string(1, i - '0'), rotated_image);
//            cvWaitKey(0);
        }
    }
    tt.stop();
    std::cout << tt.duration() << "ms" << std::endl;

    for (int i = 0; i < optimal_pyr_level; ++i) {
        std::cout << tpls[i][0].noOfCordinates << std::endl;
    }
    //建立完模板需要将模板发送给客户端，需要发送的就是tpls这个数据结构
}