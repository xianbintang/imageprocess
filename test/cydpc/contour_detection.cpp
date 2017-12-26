//
// Created by yudong on 2017/12/16
//

#include <opencv2/core/mat.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv/cv.hpp>
#include <iostream>
#include <TimeTracker.h>
#include <queue>
#include <cstdio>
#include <fstream>
#include "contour_detection.h"


// 几个之后主要函数需要用到的工具函数

// 1. 根据待测图长宽，对应到待测图的坐标，以及传递过来的模板的搜索框边长, 求出相交的矩形
std::vector<cv::Point> get_search_rect(int H, int W, cv::Point mapCenter, UINT16 top_search_rect_width){
    int a = top_search_rect_width;
    int x1 = std::max(0, mapCenter.x - a);
    int x3 = std::min(W, mapCenter.x + a);
    int y1 = std::max(0, mapCenter.y - a);
    int y3 = std::min(H, mapCenter.y + a);

    std::vector<cv::Point> rectPoint = {
        {x1, y1},
        {x1, y3},
        {x3, y3},
        {x3, y1}
    };
    return rectPoint;
}

// 2. 根据矩形，求出此时待测图中的轮廓点数
int get_num_of_rect(std::vector<cv::Point> &rectPoint, std::vector<std::vector<bool>> &hasPoint){
    int x1 = rectPoint[0].x;
    int x3 = rectPoint[2].x;
    int y1 = rectPoint[0].y;
    int y3 = rectPoint[2].y;

    int noOfPoint = 0;
    for(int i = x1; i < x3; ++i){
        for(int j = y1; j < y3; ++j){
            if(hasPoint[j][i]){
                ++noOfPoint;
            }
        }
    }
    return noOfPoint;
}

// 3. 该点是否包括在该矩形区域
bool contain_point_in_rect(cv::Point &point, const std::vector<cv::Point> &rectPoint){
    return (point.x >= rectPoint[0].x && point.x < rectPoint[2].x
            && point.y >= rectPoint[0].y && point.y < rectPoint[2].y);
}

// 4. 计算待测图中相对应的点的个数
int num_of_mapping_point(const std::vector<cv::Point> &rect_point, std::vector<cv::Point> &tpl_point,
                         const std::vector<std::vector<bool>> &hasPoint, const cv::Point &mapCenter){
    int noOfPoint = 0;
    for(auto point : tpl_point){
        point.x += mapCenter.x;
        point.y += mapCenter.y;
        if(contain_point_in_rect(point, rect_point) && hasPoint[point.y][point.x]){
            ++noOfPoint;
        }
    }
    return noOfPoint;
}

// 5. 绘制轮廓的函数
static void draw_template(cv::Mat &src, const TemplateStruct &tpl, cv::Point centerPoint)
{
    for (UINT32 i = 0; i < tpl.noOfCordinates; ++i) {
        cv::circle(src, cv::Point(tpl.cordinates[i].x + centerPoint.x, tpl.cordinates[i].y + centerPoint.y),
                   1, cv::Scalar(255,255,255));
    }
}



// 创建一些全局变量来管理待测图的图像金字塔
std::vector<cv::Mat> sample_pyramids;
std::vector<SampleStruct> sampleStructs;


// 创建待测图金字塔以及填充必要的数据结构，进行高斯滤波，Sobel算子的计算，以及大津阈值二值化
int do_create_sample_pyramid(const cv::Mat &src)
{
    TimeTracker tt1;
    tt1.start();

    sample_pyramids.push_back(src);

    // 建立各层金字塔, 并确定最佳金字塔层数
    int optimal_pyr_level = 0;

    for (int i = 0; i < MAX_NUM_PYRAMID - 1; ++i) {
        cv::Mat next_level;

        cv::pyrDown(sample_pyramids[i], next_level);
        std::cout << "size of this level: cols " << next_level.cols << ", rows" << next_level.rows << std::endl;
        sample_pyramids.push_back(next_level);
//        cv::imshow("sample_pyr" , sample_pyramids[i+1]);
//        cvWaitKey(0);
    }

    //
    for (auto &pyr : sample_pyramids) {

        cv::Mat after_gaus;
        cv::GaussianBlur(pyr, after_gaus, cv::Size(5,5),0,0);

        //width and height
        SampleStruct sampleStruct;
        sampleStruct.height = pyr.rows;
        sampleStruct.width = pyr.cols;

        //sobel
        cv::Mat gx;                //Matrix to store X derivative
        cv::Mat gy;                //Matrix to store Y derivative

        cv::Sobel(after_gaus, gx, CV_16S, 1,0,3);        //gradient in X direction
        cv::Sobel(after_gaus, gy, CV_16S, 0,1,3);        //gradient in Y direction

        //Canny
        cv::Mat sampleContour;
        cv::Canny(after_gaus, sampleContour, LOW_THRESHOLD, HIGH_THRESHOLD);

        //Ostu
        cv::Mat binaryContour;
        threshold(sampleContour, binaryContour, 0, 255, CV_THRESH_OTSU);


        UINT32 numOfPoint = 0;
        // 计算待测图中所有点的x,y方向向量
        for (int i = 0; i < binaryContour.rows; i++) {
            std::vector<float> vecX;
            std::vector<float> vecY;
            std::vector<bool> vecPoint;

            // 将宽度看成x轴，高度看成y轴
            for (int j = 0; j < binaryContour.cols; j++) {
                short fdx = gx.at<short>(i, j);
                short fdy = gy.at<short>(i, j);

                double vector_length = sqrt(fdx * fdx + fdy * fdy);
                if(fdx != 0 || fdy != 0){
                    vecX.push_back(fdx / vector_length);
                    vecY.push_back(fdy / vector_length);
                }else{
                    vecX.push_back(0.0);
                    vecY.push_back(0.0);
                }
                if(binaryContour.at<uchar>(i, j)){
                    ++numOfPoint;
                    vecPoint.push_back(true);
                }else{
                    vecPoint.push_back(false);
                }
            }
            sampleStruct.edgeDerivativeX.push_back(vecX);
            sampleStruct.edgeDerivativeY.push_back(vecY);
            sampleStruct.hasPoint.push_back(vecPoint);
        }

        sampleStructs.push_back(sampleStruct);

//        cv::imshow("info", binaryContour);
//        cv::waitKey(0);

//        std::cout << "the num of points is " << numOfPoint << std::endl;
        if(numOfPoint <= MIN_CONTOUR_PYRA){
            break;
        }

        optimal_pyr_level++;
    }

    std::cout << "sample pyramid optimal level: " <<  optimal_pyr_level << std::endl;

    tt1.stop();
    std::cout << "create sample pyramid cost : " << tt1.duration() << " ms" << std::endl;
    return 0;
}


// 模板顶层对应的待测图金字塔的相应层数进行匹配
int do_top_layer_match(int level, std::vector<TemplateStruct> &top_tpls, float top_search_angel_nstep,
                 UINT16 top_search_rect_width, float threshold_S){
    // 判断一下一些是否是非正常的图片，如果不满足要求退出
    if(level > sampleStructs.size()-1){
        std::cout << "Error, the sample pyramid's isn't tall enough!" << std::endl;
        exit(-1);
    }
    // 准备必要的参数(0度模板时的轮廓点总数)
    int K = top_tpls[0].noOfCordinates;
    SampleStruct &spl = sampleStructs[level];
//    std::cout << spl.height << " " << spl.width << std::endl;
    // 最大相似度
    double max_s = - 1.0;
    for(int i = 0; i < spl.height; ++i){
        for(int j = 0; j < spl.width; ++j){
            // 对应到待测图的坐标
            cv::Point mapCenter;
            mapCenter.x = j;
            mapCenter.y = i;
            // 角度模板的序号置为0
            int l = 0;
            std::vector<cv::Point> rectPoint = get_search_rect(spl.height, spl.width, mapCenter, top_search_rect_width);
//            // 判断搜索框的运行情况
//            cv::Mat rectPic = sample_pyramids[4].clone();
//            cv::rectangle(rectPic,
//                       rectPoint[0],
//                       rectPoint[2],
//                       cv::Scalar( 0, 0, 255 ),
//                       1,
//                       8 );

            int num = get_num_of_rect(rectPoint, spl.hasPoint);
//            std::cout << "this rect has num is " << num << std::endl;
//            cv::imshow("rectangle", rectPic);
//            cv::waitKey(0);
            if(num > K * threshold_S){
//                std::cout << "OK the num of Point more than K*S and now center is ("
//                          << j << "," << i << ")" << std::endl;

                for(double k = 0.0; k < MAX_DEGREE; k += top_search_angel_nstep){
                    TemplateStruct &tpl = top_tpls[l++];
                    int K1 = tpl.cordinates.size();
                    int num_mapping = num_of_mapping_point(rectPoint, tpl.cordinates, spl.hasPoint, mapCenter);
//                    std::cout << "Top level num of mapping point is " << num_mapping << std::endl;
                    if(num_mapping > K * std::max(0.0, threshold_S - 0.35)){
//                        std::cout << "OK the num of Point more than K*S and now center is ("
//                                  << j << "," << i << "), and the angle is " << k << std::endl;
                        // 相似度计算公式
                        int num_T = 0;
                        double sum_of_s = 0.0;
                        int order = 0;
                        for(auto point : tpl.cordinates){
                            point.x += mapCenter.x;
                            point.y += mapCenter.y;
                            if(contain_point_in_rect(point, rectPoint) && spl.hasPoint[point.y][point.x]){
                                ++num_T;
                                sum_of_s += fabs(tpl.edgeDerivativeX[order] * spl.edgeDerivativeX[point.y][point.x]
                                        + tpl.edgeDerivativeY[order] * spl.edgeDerivativeY[point.y][point.x]);
                            }
                            ++order;
                            if(sum_of_s < (threshold_S-1)*K1 + num_T){
                                break;
                            }
                        }
                        double s = sum_of_s / num_T;
//                        std::cout << "OK the num of Point more than K*S and now center is ("
//                                  << j << "," << i << "), and the angle is " << k
//                                  << "the similarity is " << s << std::endl;

                        // 匹配相似度最高的质心和角度记录下来
                        if(s > max_s){
                            max_s = s;
                            spl.angle = k;
                            spl.center.x = j;
                            spl.center.y = i;
                            spl.similarity = s;
//                            std::cout << "now center is (" << j << "," << i << "), and the angle is " << k
//                                      << "the similarity is " << s << std::endl;
                        }
                    }
                }
            }
        }
    }
    return 0;
}


// 进行其他层的匹配，从最上层往下一层匹配
int do_other_layer_match(int cur_level, std::vector<TemplateStruct> &cur_tpls, float cur_angle_step,
                                   UINT16 cur_search_rect_width, ImportantInfo &upper_info, float threshold_S){
    ImportantInfo cur_info;
    cur_info.center.x = 2 * upper_info.center.x;
    cur_info.center.y = 2 * upper_info.center.y;
    cur_info.angle = upper_info.angle;

    // 备选的质心坐标周围的点以及以当前角度左右几个该层角度步长的角度
    std::vector<cv::Point> point_bak;
    for(int i = cur_info.center.x-1; i <= cur_info.center.x+1; ++i){
        for(int j = cur_info.center.y-1; j <= cur_info.center.y+1; ++j){
            point_bak.push_back(cv::Point(i, j));
        }
    }
    int l = 0;
    for(double k = 0.0; k < cur_info.angle; k += cur_angle_step){
        l++;
    }
    std::vector<int> angle_index;
    for(int i = -2; i <= 2; ++i){
        angle_index.push_back((l+i)%cur_tpls.size());
    }


    TimeTracker tt1;
    tt1.start();

    // 进行相似度的计算
    SampleStruct &spl = sampleStructs[cur_level];

    double max_s = -1.0;
    for(auto centerPoint : point_bak){
        std::vector<cv::Point> rectPoint = get_search_rect(spl.height, spl.width, centerPoint, cur_search_rect_width);
        for(auto l : angle_index){
            TemplateStruct &tpl = cur_tpls[l];
            int K1 = tpl.cordinates.size();
            int num_T = 0;
            double sum_of_s = 0.0;
            int order = 0;
            for(auto point : tpl.cordinates){
                point.x += centerPoint.x;
                point.y += centerPoint.y;
                if(contain_point_in_rect(point, rectPoint) && spl.hasPoint[point.y][point.x]){
                    ++num_T;
                    sum_of_s += fabs(tpl.edgeDerivativeX[order] * spl.edgeDerivativeX[point.y][point.x]
                            + tpl.edgeDerivativeY[order] * spl.edgeDerivativeY[point.y][point.x]);
                    if((int)cur_level == 2 && centerPoint.x == 81 && centerPoint.y == 69) {
//                        std::cout << l << " " << tpl.cordinates[order] << " " << point << " "<< tpl.edgeDerivativeX[order] << " "<< spl.edgeDerivativeX[point.y][point.x]<< " "
//                            << tpl.edgeDerivativeY[order] << " "<< spl.edgeDerivativeY[point.y][point.x] << std::endl;
//                    std::cout << tpl.noOfCordinates << std::endl;
                    }
                }
                ++order;
                if(sum_of_s < (threshold_S-1)*K1 + num_T){
                    break;
                }
            }
            double s = sum_of_s / num_T;

            // 匹配相似度最高的质心和角度记录下来
            if(s > max_s){
                max_s = s;
                spl.angle = l * cur_angle_step;
                spl.center = centerPoint;
                spl.similarity = s;
//                std::cout << "now center is (" << spl.center.x << "," << spl.center.y
//                          << "), and the angle is " << spl.angle
//                          << "the similarity is " << s << std::endl;
            }
        }
    }

    tt1.stop();
    std::cout << cur_level << " layer match cost : " << tt1.duration()  << " ms"<< std::endl;
    return 0;
}


static void bitmap2Mat(cv::Mat &dst, std::vector<UINT8> &bitmap, UINT16 width, UINT16 height){
    for(int i = 0; i < height; ++i){
        for(int j = 0; j < width; ++j){
            dst.at<uchar>(i, j) = bitmap[i * width + j];
        }
    }
}

int do_create_sample_pyramid_1(){

    cv::Mat pyramid_zero;
    cv::Mat pyramid_one;
    cv::Mat pyramid_two;

    // 0
    std::string filename1("gray0.txt");
    std::vector<UINT8> values1;
    std::ifstream infile1(filename1.c_str());
    while(infile1){
        int data;
        infile1 >> data;
        UINT8 data1 = static_cast<UINT8>(data);
        values1.push_back(data1);
    }
    pyramid_zero.create(480, 640, CV_8UC1);
    bitmap2Mat(pyramid_zero, values1, 640, 480);

    // 1
    std::string filename2("gray1.txt");
    std::vector<UINT8> values2;
    std::ifstream infile2(filename2.c_str());
    while(infile2){
        int data;
        infile2 >> data;
        UINT8 data1 = static_cast<UINT8>(data);
        values2.push_back(data1);
    }
    pyramid_one.create(240, 320, CV_8UC1);
    bitmap2Mat(pyramid_one, values2, 320, 240);

    // 2
    std::string filename3("gray2.txt");
    std::vector<UINT8> values3;
    std::ifstream infile3(filename3.c_str());
    while(infile3){
        int data;
        infile3 >> data;
        UINT8 data1 = static_cast<UINT8>(data);
        values3.push_back(data1);
    }
    pyramid_two.create(120, 160, CV_8UC1);
    bitmap2Mat(pyramid_two, values3, 160, 120);

    // put in pyramid
    sample_pyramids.push_back(pyramid_zero);
    sample_pyramids.push_back(pyramid_one);
    sample_pyramids.push_back(pyramid_two);

    int optimal_pyr_level = 0;
    for (auto &pyr : sample_pyramids) {

        cv::Mat after_gaus;
        cv::GaussianBlur(pyr, after_gaus, cv::Size(5,5),0,0);

        //width and height
        SampleStruct sampleStruct;
        sampleStruct.height = pyr.rows;
        sampleStruct.width = pyr.cols;

        //sobel
        cv::Mat gx;                //Matrix to store X derivative
        cv::Mat gy;                //Matrix to store Y derivative

        cv::Sobel(after_gaus, gx, CV_16S, 1,0,3);        //gradient in X direction
        cv::Sobel(after_gaus, gy, CV_16S, 0,1,3);        //gradient in Y direction

        //Canny
        cv::Mat sampleContour;
        cv::Canny(after_gaus, sampleContour, LOW_THRESHOLD, HIGH_THRESHOLD);

        //Ostu
        cv::Mat binaryContour;
        threshold(sampleContour, binaryContour, 0, 255, CV_THRESH_OTSU);


        UINT32 numOfPoint = 0;
        // 计算待测图中所有点的x,y方向向量
        for (int i = 0; i < binaryContour.rows; i++) {
            std::vector<float> vecX;
            std::vector<float> vecY;
            std::vector<bool> vecPoint;

            // 将宽度看成x轴，高度看成y轴
            for (int j = 0; j < binaryContour.cols; j++) {
                short fdx = gx.at<short>(i, j);
                short fdy = gy.at<short>(i, j);

                double vector_length = sqrt(fdx * fdx + fdy * fdy);
                if(fdx != 0 || fdy != 0){
                    vecX.push_back(fdx / vector_length);
                    vecY.push_back(fdy / vector_length);
                }else{
                    vecX.push_back(0.0);
                    vecY.push_back(0.0);
                }
                if(binaryContour.at<uchar>(i, j)){
                    ++numOfPoint;
                    vecPoint.push_back(true);
                }else{
                    vecPoint.push_back(false);
                }
            }
            sampleStruct.edgeDerivativeX.push_back(vecX);
            sampleStruct.edgeDerivativeY.push_back(vecY);
            sampleStruct.hasPoint.push_back(vecPoint);
        }

        sampleStructs.push_back(sampleStruct);

        cv::imshow("info", binaryContour);
        cv::waitKey(0);

        std::cout << "the num of points is " << numOfPoint << std::endl;
        if(numOfPoint <= MIN_CONTOUR_PYRA){
            break;
        }

        optimal_pyr_level++;
    }

    std::cout << "sample pyramid optimal level: " <<  optimal_pyr_level << std::endl;

    return 0;
}


// 最主要的函数
int match_sample_template(const cv::Mat &sample_image, Koyo_Contour_Template_Runtime_Param &kctrp)
{

    // 建立待测图的图像金字塔
    do_create_sample_pyramid(sample_image);
//    do_create_sample_pyramid_1();


    TimeTracker tt1;


    // 进行金字塔顶层匹配
    tt1.start();
    int top_level = kctrp.run_time_npyramid-1;
    std::vector<TemplateStruct> &top_tpls = kctrp.tpls[top_level];
    float top_search_angel_nstep = kctrp.search_angel_nstep[top_level];
    UINT16 top_search_rect_width = kctrp.search_rect_width[top_level];
    do_top_layer_match(top_level, top_tpls, top_search_angel_nstep, top_search_rect_width, 0.7);
    tt1.stop();
    std::cout << "top layer match cost : " << tt1.duration()  << " ms"<< std::endl;


    // 继续向下层拓展，计算每层的匹配情况
    tt1.start();
    int cur_level = top_level-1;
    while(cur_level != -1){
        // 上一层所得的位置信息
        ImportantInfo upper_info;
        upper_info.center = sampleStructs[cur_level+1].center;
        upper_info.angle = sampleStructs[cur_level+1].angle;

        std::vector<TemplateStruct> &cur_tpls = kctrp.tpls[cur_level];
        float cur_angle_step = kctrp.search_angel_nstep[cur_level];
        UINT16 cur_search_rect_width = kctrp.search_rect_width[cur_level];
        do_other_layer_match(cur_level, cur_tpls, cur_angle_step, cur_search_rect_width, upper_info, 0.7);
        cur_level--;
        std::cout << std::endl;
    }
    tt1.stop();
    std::cout << "other layer match cost : " << tt1.duration() << " ms" << std::endl;


    // 测试一下每一层的匹配情况
    int level = top_level;
    while(level != -1){
        cv::Point result_Center = sampleStructs[level].center;
        float result_Angle = sampleStructs[level].angle;
        float result_s = sampleStructs[level].similarity;

        std::cout << std::endl << "************************************" << std::endl;
        std::cout << "At level" << level << std::endl;
        std::cout << "the result of center is (" << result_Center.x << "," << result_Center.y << "), "
                  << "the result of angle is " << result_Angle << ", the similarity is "
                  << result_s << "." <<std::endl;
        std::cout << "************************************" << std::endl;

        float cur_angle_step = kctrp.search_angel_nstep[level];
        cv::Mat src = sample_pyramids[level].clone();
        int l = 0;
        for(float k = 0.0; k < result_Angle; k += cur_angle_step){
            l++;
        }
        l = std::max(l-1, 0);
        TemplateStruct &tpl = kctrp.tpls[level][l];
        draw_template(src, tpl, result_Center);

        cv::imshow("contour", src);
        cv::imshow("raw", sample_pyramids[level]);
        cvWaitKey(0);
        level--;
    }

    std::cout << "over" << std::endl;
    return 0;
}
