//
// Created by xianb on 2017/11/28.
//


#include <opencv2/core/mat.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv/cv.hpp>
#include <iostream>
#include <TimeTracker.h>
#include <queue>
#include <cstdio>
#include "contour_detection.h"

static cv::Mat get_y_from_yuv(const UINT8 *yuv, const UINT16 width, const UINT16 height)
{
    cv::Mat gray;
    if (!yuv) {
        std::cout << "ERROR: yuv NULL pointer" << std::endl;
    } else {
        gray.create(height, width, CV_8UC1);
        memcpy(gray.data, yuv, width * height * sizeof(unsigned char)); // 只读取y分量上的数据
    }

    return gray;
}

static std::unique_ptr<char[]> pack_template(const Koyo_Contour_Template_Runtime_Param &koyo_contour_template_runtime_param)
{
    // 先计算所有要使用的内存大小，然后分配空间，最后一点点将数据拷贝过去
    std::size_t buf_size = 0;

    buf_size += sizeof(koyo_contour_template_runtime_param.run_time_npyramid);
    buf_size += sizeof(float) * koyo_contour_template_runtime_param.run_time_npyramid;
    buf_size += sizeof(UINT16) * koyo_contour_template_runtime_param.run_time_npyramid;
    // 要记录每层金字塔上的模板个数

    for (const auto &tpl_arr: koyo_contour_template_runtime_param.tpls) {
        // 每层金字塔上的模板个数
        buf_size += sizeof(UINT16);
        for (const auto &tpl : tpl_arr) {
            buf_size += sizeof (tpl.modelDefined);
            buf_size += sizeof (tpl.noOfCordinates);
            buf_size += sizeof (tpl.modelHeight);
            buf_size += sizeof (tpl.modelWidth);

            buf_size += sizeof(short) * 2;

            buf_size += sizeof(short) * 2 * tpl.noOfCordinates;
            buf_size += sizeof(float) * tpl.noOfCordinates;
            buf_size += sizeof(float) * tpl.noOfCordinates;
        }
    }

    // 分配空间
    std::unique_ptr<char[]> buf(new char[buf_size]);

//    std::cout << "hehe " << koyo_contour_template_runtime_param.search_angel_nstep.size() << std::endl;
    std::size_t index = 0;
    memcpy(&buf[index], &koyo_contour_template_runtime_param.run_time_npyramid, sizeof(koyo_contour_template_runtime_param.run_time_npyramid));
    index += sizeof(koyo_contour_template_runtime_param.run_time_npyramid);

    for (int i = 0; i < koyo_contour_template_runtime_param.run_time_npyramid; ++i) {
        memcpy(&buf[index], &koyo_contour_template_runtime_param.search_angel_nstep[i], sizeof(float));
        index += sizeof(float);
    }


    for (int i = 0; i < koyo_contour_template_runtime_param.run_time_npyramid; ++i) {
        memcpy(&buf[index], &koyo_contour_template_runtime_param.search_rect_width[i], sizeof(UINT16));
        index += sizeof(UINT16);
    }


    // 拷贝模板数据
    for (const auto &tpl_arr: koyo_contour_template_runtime_param.tpls) {
        UINT16 tpl_size = static_cast<UINT16>(tpl_arr.size()); //肯定不会超的
        memcpy(&buf[index], &tpl_size, sizeof(UINT16));
        index += sizeof(UINT16);

        for (const auto &tpl : tpl_arr) {
            memcpy(&buf[index], &tpl.modelDefined, sizeof(tpl.modelDefined));
            index += sizeof(tpl.modelDefined);

            memcpy(&buf[index], &tpl.noOfCordinates, sizeof(tpl.noOfCordinates));
            index += sizeof(tpl.noOfCordinates);

            memcpy(&buf[index], &tpl.modelHeight, sizeof(tpl.modelHeight));
            index += sizeof(tpl.modelHeight);

            memcpy(&buf[index], &tpl.modelWidth, sizeof(tpl.modelWidth));
            index += sizeof(tpl.modelWidth);

            memcpy(&buf[index], &tpl.centerOfGravity.x, sizeof(short));
            index += sizeof(short);

            memcpy(&buf[index], &tpl.centerOfGravity.y, sizeof(short));
            index += sizeof(short);

            for (auto const &coord : tpl.cordinates) {
                memcpy(&buf[index], &coord.x, sizeof(short));
                index += sizeof(short);

                memcpy(&buf[index], &coord.y, sizeof(short));
                index += sizeof(short);
            }
            for (auto const & edgeX : tpl.edgeDerivativeX) {
                memcpy(&buf[index], &edgeX, sizeof(float));
                index += sizeof(float);
            }

            for (auto const & edgeY : tpl.edgeDerivativeY) {
                memcpy(&buf[index], &edgeY, sizeof(float));
                index += sizeof(float);
            }
        }
    }

    std::cout << buf_size << ", in MB: " << 1.0 * buf_size / 1024 / 1024 << "MB" << std::endl;
    float y1 = buf[buf_size - 4];
    float y2 = buf[buf_size - 8];
    printf("y1:%f  y2:%f\n", *((float*)&buf[buf_size - 4]), *((float*)&buf[buf_size - 8]));

#ifdef _DEBUG_LEVEL_HIGH_
    // 测试写入文件
    TimeTracker timeTracker;
    timeTracker.start();
    FILE *outfile = fopen("template_file", "wb+");
    if(!outfile) {
        std::cout << "Error while open template file for write" << std::endl;
        exit(-1);
    }
    int more_bytes = buf_size, count = 0;
    int write_size = buf_size > 2048 ? 2048 : buf_size;
    char *pb = buf.get();
    while (count < buf_size && more_bytes > 0) {
        fwrite(pb + count, write_size, 1, outfile);
        count += 2048;
        more_bytes -= 2048;
//        std::cout << "count: " << count << std::endl;
//        std::cout << "more: " << more_bytes << std::endl;
        write_size = more_bytes > 2048 ? 2048 : more_bytes;
    }
    fclose(outfile);
    timeTracker.stop();
    std::cout << "write time: " << timeTracker.duration() << std::endl;


    //测试读取文件, 再把template文件恢复出来
    timeTracker.start();
    FILE *infile = fopen("template_file", "rb+");
    if(!infile) {
        std::cout << "Error while open template file for write" << std::endl;
        exit(-1);
    }

    pb = buf.get();

    count = 0;
    int rc;
    while (0 < (rc = fread(pb + count, 2048, 1, infile))) {
        count += 2048;
    }

    fclose(infile);
    timeTracker.stop();
    std::cout << "read time: " << timeTracker.duration() << std::endl;
#endif
    return buf;
}

//int unpack_template(Koyo_Contour_Template_Runtime_Param &koyo_contour_template_runtime_param, std::unique_ptr<char[]> buf)
static int unpack_template(Koyo_Contour_Template_Runtime_Param &koyo_contour_template_runtime_param, char* buf)
{

#ifdef _DEBUG_LEVEL_HIGH_
    TimeTracker timeTracker;
    timeTracker.start();
    FILE *infile = fopen("template_file", "rb+");
    if(!infile) {
        std::cout << "Error while open template file for write" << std::endl;
        exit(-1);
    }

    char *pb = (char *)malloc(50 * 1024 * 1024); //预分配30MB大小

    int count = 0;
    int rc;
    while (0 < (rc = fread(pb + count, 2048, 1, infile))) {
        count += 2048;
    }

    fclose(infile);
    timeTracker.stop();
    std::cout << "read time: " << timeTracker.duration() << std::endl;
    buf = pb;
#endif

    std::size_t index = 0;

    koyo_contour_template_runtime_param.run_time_npyramid = *((UINT8 *)&buf[index]);
    index += sizeof(koyo_contour_template_runtime_param.run_time_npyramid);

    for (int i = 0; i < koyo_contour_template_runtime_param.run_time_npyramid; ++i) {
        float angle = *((float*)&buf[index]);
        index += sizeof(float);
        koyo_contour_template_runtime_param.search_angel_nstep.push_back(angle);
    }

    for (int i = 0; i < koyo_contour_template_runtime_param.run_time_npyramid; ++i) {
        UINT16 width= *((UINT16*)&buf[index]);
        index += sizeof(UINT16);
        koyo_contour_template_runtime_param.search_rect_width.push_back(width);
    }

    // 拷贝模板数据
    for (int i = 0; i < koyo_contour_template_runtime_param.run_time_npyramid; ++i) {
        std::vector<TemplateStruct> tpl_arr;

        UINT16 tpl_size = *((UINT16*)&buf[index]);
        index += sizeof(UINT16);

//        std::cout << "tpl_size: " << tpl_size << std::endl;
        for (int j = 0; j < tpl_size; ++j) {
            TemplateStruct tpl;
            tpl.modelDefined = *((UINT8*)&buf[index]);
            index += sizeof(tpl.modelDefined);

            tpl.noOfCordinates= *((UINT32*)&buf[index]);
            index += sizeof(tpl.noOfCordinates);

            tpl.modelHeight= *((UINT16*)&buf[index]);
            index += sizeof(tpl.modelHeight);

            tpl.modelWidth= *((UINT16*)&buf[index]);
            index += sizeof(tpl.modelWidth);

            tpl.centerOfGravity.x = *((UINT16*)&buf[index]);
            index += sizeof(short);

            tpl.centerOfGravity.y = *((UINT16*)&buf[index]);
            index += sizeof(short);

            // 拷贝特征数据
//            std::cout << "tpl no ofcoordinate: " << tpl.noOfCordinates << std::endl;
            for (std::size_t i = 0; i < tpl.noOfCordinates; ++i) {

                cv::Point coord;
                coord.x = *((short *) &buf[index]);
                index += sizeof(short);

                coord.y = *((short *) &buf[index]);
                index += sizeof(short);
                tpl.cordinates.push_back(coord);
            }

            for (std::size_t i = 0; i < tpl.noOfCordinates; ++i) {
                float edgeX = *((float*)&buf[index]);
                index += sizeof(float);
                tpl.edgeDerivativeX.push_back(edgeX);

                float edgeY = *((float*)&buf[index]);
                index += sizeof(float);
                tpl.edgeDerivativeY.push_back(edgeY);


            }
            tpl_arr.push_back(tpl);
        }
        koyo_contour_template_runtime_param.tpls.push_back(tpl_arr);
    }

    return 0;
}

static std::ostream &print_tpl(std::ostream &os, const TemplateStruct &templateStruct)
{
    std::cout << "template size, width: " << templateStruct.modelWidth <<\
     ", height: " << templateStruct.modelHeight << std::endl;
    std::cout << "Number of coordinate: " << templateStruct.noOfCordinates << std::endl;
    std::cout << "center of gravity: " << templateStruct.centerOfGravity << std::endl;
}
static std::ostream &print_tpls(std::ostream &os, const Koyo_Contour_Template_Runtime_Param & rhl)
{
    os << "runtime npyramid: " << static_cast<int>(rhl.run_time_npyramid) << std::endl;
    os << "angle steps each level: " << std::endl;
    int k = 0;
    for (auto angle : rhl.search_angel_nstep) {
        os << "level " << k++ << ": " << angle << std::endl;
    }

    k = 0;
    os << "templates each level: " << std::endl;
    for (auto tpl: rhl.tpls) {
        os << "level " << k++ << ": " << tpl.size() << ", Num of coordinate this level: " << tpl[0].noOfCordinates << std::endl;
    }
    return os;
}



static std::vector<float> rotate_image(const cv::Mat &src, cv::Mat &dst, cv::Point centerP, float degree)
{
    int width = src.cols;
    int height = src.rows;
    double angle = degree  * CV_PI / 180.; // 弧度
    double a = sin(angle), b = cos(angle);
    // 适当增大一点宽高，防止像素不在图像内
    int width_rotate = int(height * fabs(a) + width * fabs(b));
    int height_rotate = int(width * fabs(a) + height * fabs(b));
    float map[6];
    cv::Mat map_matrix = cv::Mat(2, 3, CV_32F, map);
// 旋转中心
    CvPoint2D32f center = cvPoint2D32f(centerP.x, centerP.y);
    CvMat map_matrix2 = map_matrix;
    cv2DRotationMatrix(center, degree, 1.0, &map_matrix2);
    map[2] += (width_rotate - width) / 2.0;
    map[5] += (height_rotate - height) / 2.0;
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
    return 0;
}

static int cutout_template_image(const cv::Mat &template_image, std::vector<cv::Point> rect, cv::Mat &interesting_template)
{
    auto degree = 360 - cv::fastAtan2(rect[1].y - rect[2].y, rect[2].x - rect[1].x);
    cv::Mat img_rotate;
    int width = template_image.cols;
    int height = template_image.rows;
    std::vector<float> rotate_matrix;
    rotate_matrix = rotate_image(template_image, img_rotate, cv::Point(width / 2, height / 2), degree);
    std::cout << degree << std::endl;
    rotate_rect(rect, rotate_matrix);

#ifdef DEBUG
    cv::circle(img_rotate, rect[0], 1, cv::Scalar(255,255,255));
    cv::circle(img_rotate, rect[1], 2, cv::Scalar(255,255,255));
    cv::circle(img_rotate, rect[2], 3, cv::Scalar(255,255,255));
    cv::circle(img_rotate, rect[3], 4, cv::Scalar(255,255,255));
#endif

    int rect_width = rect[3].x - rect[0].x;
    int rect_height = rect[1].y - rect[0].y;
    interesting_template = img_rotate(cv::Rect(rect[0].x, rect[0].y, rect_width, rect_height));

//    cv::rectangle(img_rotate, rect[0], rect[2], cv::Scalar(255,255,255));
//    cv::imshow("eh" ,interesting_template);
    cvWaitKey(0);
    return 0;
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
    return 0;
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
static int do_create_template(TemplateStruct &tpl, const cv::Mat &src, double low_threshold,\
 double high_threshold, const std::vector<cv::Point> &rect)
{
    int s32Ret = 0;
    cv::Mat gx;                //Matrix to store X derivative
    cv::Mat gy;                //Matrix to store Y derivative
    // set width and height
    tpl.modelHeight = src.rows;    //Save Template height
    tpl.modelWidth = src.cols;    //Save Template width

    tpl.noOfCordinates = 0;    //initialize

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
                    tpl.cordinates.push_back(p);
//                    tpl.cordinates[tpl.noOfCordinates].x = j;
//                    tpl.cordinates[tpl.noOfCordinates].y = i;

                    /* TODO 可以修改成使用查找表的形式 */
                    double vector_length = sqrt(fdx * fdx + fdy * fdy);
                    if (fabs(vector_length - 0.) < 0.00001) {
//                        printf(".............................................\n");
                    }
                    tpl.edgeDerivativeX.push_back(fdx / vector_length);
                    tpl.edgeDerivativeY.push_back(fdy / vector_length);
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
    UINT32 m;
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
    for (UINT32 i = 0; i < tpl.noOfCordinates; ++i) {
        cv::circle(src, cv::Point(tpl.cordinates[i].x + tpl.centerOfGravity.x, tpl.cordinates[i].y + tpl.centerOfGravity.y), 1, cv::Scalar(255,255,255));
    }
    cv::imshow("hehe", src);
    cvWaitKey(0);
}

// TODO 发送给客户端以后需要释放
//int free_tpls(TemplateStruct tpls[][MAX_DEGREE])
//{
//    return 0;
//}

/*
 *  @param src必须是二值化轮廓图
 *  @return 返回当前层上最佳的旋转步长
 * */
static float get_angle_step_and_search_width(const cv::Mat &src, cv::Point center, std::vector<UINT16> &search_rect_width)
{
    // 保留几个K，然后求平均值，用来排除外点的影响
    int K = 50;
    std::priority_queue<float> max_dist(K, -1);
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
//    std::cout << "max dist....." << std::endl;
    float average_max_dist = 0;
    int i = 0;
    while (!max_dist.empty()) {
        average_max_dist += -max_dist.top();
        max_dist.pop();
    }
    average_max_dist /= K;
//    std::cout << "average_max_dist: " << average_max_dist << std::endl;

    search_rect_width.push_back(static_cast<UINT16>(average_max_dist));
    auto range_low = acos(1 - 1 / (2 * average_max_dist * average_max_dist)) / CV_PI * 360;
    auto range_high = acos(1 - 1 / (average_max_dist * average_max_dist)) / CV_PI * 360;
//    std::cout <<"optimal angle step: " << range_low << " ~ " << range_high << std::endl;

#ifdef _DEBUG_
    cv::Mat tmp = src;
    for (unsigned int k = points.size() - K; k < points.size(); ++k) {
        cv::circle(tmp, points[k], 20, cv::Scalar(255,255,255));
    }
//    cv::imshow("max dist", tmp);
//    cvWaitKey(0);
#endif
    return std::max((range_low + range_high) / 2 ,1.0) ;
}

/*
 * 为koyo_tool_contour_parameter对应的工具建立模板, 同时计算出下列参数
 * 1. UINT8 run_time_npyramid;
 * 2. double search_angel_nstep[MAX_NUM_PYRAMID];
 * 3. TemplateStruct tpls[MAX_NUM_PYRAMID][MAX_DEGREE];
 */

//#define _DEBUG_
#ifdef _DEBUG_
std::vector<cv::Point> centers;
std::vector<float> angle_steps;
std::vector<cv::Mat> pyramid_templates;
std::vector<UINT16> search_rect_width;
#else

#endif

static int do_create_template(const cv::Mat &src, Koyo_Tool_Contour_Parameter koyo_tool_contour_parameter, Koyo_Contour_Template_Runtime_Param &koyo_contour_template_runtime_param)
{
    // 运行完成后需要将这个内容发送给嵌入式

    TimeTracker tt1;
    tt1.start();
    // todo 重构成vector版本
//    cv::Mat pyramid_templates[MAX_NUM_PYRAMID];
#ifndef _DEBUG_
    std::vector<cv::Mat> pyramid_templates;
#endif
    pyramid_templates.push_back(src);
//    pyramid_templates[0] = src;

    UINT8 sensitity_threshold_low, sensitity_threshold_high;
    if (koyo_tool_contour_parameter.sensitivity == CONTOUR_ACCURACY_LOW) {
        sensitity_threshold_low = 10;
        sensitity_threshold_high = 80;
    } else if (koyo_tool_contour_parameter.sensitivity == CONTOUR_ACCURACY_MEDIUM) {
        sensitity_threshold_low = 30;
        sensitity_threshold_high = 150;
    } else if (koyo_tool_contour_parameter.sensitivity == CONTOUR_ACCURACY_HIGH) {
        sensitity_threshold_low = 60;
        sensitity_threshold_high = 220;
    }

    // 建立各层金字塔, 并确定最佳金字塔层数
    int optimal_pyr_level = 0;
//    std::cout << "size of this level: cols " << pyramid_templates[0].cols << ", rows" << pyramid_templates[0].rows << std::endl;
    for (int i = 0; i < MAX_NUM_PYRAMID - 1; ++i) {
        cv::Mat next_level;
        cv::pyrDown(pyramid_templates[i], next_level);
//        std::cout << "size of this level: cols " << next_level.cols << ", rows" << next_level.rows << std::endl;
        pyramid_templates.push_back(next_level);
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
#ifndef _DEBUG_
    std::vector<cv::Point> centers;
    std::vector<float> angle_steps;
    std::vector<UINT16> search_rect_width;
#endif
    for (auto &pyr : pyramid_templates) {
//    for (int i = 0; i < MAX_NUM_PYRAMID; ++i) {
        cv::Mat cannyResult;
        cv::Canny(pyr, cannyResult, sensitity_threshold_high, sensitity_threshold_low);

        cv::Point center;
        unsigned int num_of_contour;
        get_center_numof_contour(cannyResult, center, num_of_contour);
        centers.push_back(center);

        // 确定角度步长, 使用Canny的轮廓图来计算最远点
        auto step = get_angle_step_and_search_width(cannyResult, center, search_rect_width);
        angle_steps.push_back(step);
#ifdef _DEBUG_
        std::cout << "num of coordinate this level: " << num_of_contour << std::endl;
#endif
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
    // todo 重构成vector版本的
//    TemplateStruct tpls[optimal_pyr_level][MAX_DEGREE];
    std::vector<std::vector<TemplateStruct>> tpls;
    TimeTracker tt;
    tt.start();
    // optimal_pyr_level肯定小于pyramid_templates的size
    for (int i = 0; i < optimal_pyr_level; ++i) {
        std::vector<TemplateStruct> cur_level_tpl;
        int k = 0;
        std::vector<cv::Point> cur_rect = {{0,0}, {0, pyramid_templates[i].rows - 1}, {pyramid_templates[i].cols - 1, pyramid_templates[i].rows - 1}, {pyramid_templates[i].cols - 1, 0}};
        for (double j = 0.0; j < MAX_DEGREE; j += angle_steps[i]) {
            TemplateStruct tpl;
            auto rect = cur_rect;
            cv::Mat rotated_image;
            // 还是无法保证完全在图片框内
            auto rotate_matrix = rotate_image(pyramid_templates[i], rotated_image, centers[i], j);
            rotate_rect(rect, rotate_matrix);
            do_create_template(tpl, rotated_image, sensitity_threshold_low, sensitity_threshold_high, rect);
            cur_level_tpl.push_back(tpl);
//            draw_template(rotated_image, tpl);
//            cv::imshow(std::string("pyr") + std::string(1, i - '0'), rotated_image);
//            cvWaitKey(0);
        }
        tpls.push_back(cur_level_tpl);
    }
    tt.stop();
    std::cout << tt.duration() << "ms" << std::endl;

#ifdef _DEBUG_
    for (auto iter = tpls.cbegin(); iter != tpls.end(); ++iter) {
//        std::cout << iter->at(0).noOfCordinates << std::endl;
    }
#endif
    //建立完模板需要将模板发送给客户端，需要发送的就是tpls这个数据结构

    koyo_contour_template_runtime_param.run_time_npyramid = optimal_pyr_level;
    koyo_contour_template_runtime_param.search_angel_nstep = angle_steps;
    // todo 换成move操作会好一些吧
    koyo_contour_template_runtime_param.tpls = tpls;
    koyo_contour_template_runtime_param.search_rect_width = search_rect_width;
    return 0;
}


static void print_debug_info(const std::vector<cv::Mat> &pyramid_template, char* template_data)
{
    std::cout << std::endl << "--------test unpack template-----------" << std::endl;
    Koyo_Contour_Template_Runtime_Param kctrp;
    // 要换成使用C语言的解析
    unpack_template(kctrp, template_data);
//    print_tpls(std::cout, koyo_contour_template_runtime_param);
    std::cout << "unpacked information" << std::endl;
    print_tpls(std::cout, kctrp);

    // 测试恢复出来的模板
    for (int i = 0; i < kctrp.run_time_npyramid; ++i) {
        int k = 0;
        std::vector<cv::Point> cur_rect = {{0,0}, {0, pyramid_templates[i].rows - 1}, {pyramid_templates[i].cols - 1, pyramid_templates[i].rows - 1}, {pyramid_templates[i].cols - 1, 0}};
        for (double j = 0.0; j < MAX_DEGREE; j += angle_steps[i]) {
            TemplateStruct tpl = kctrp.tpls[i][k++];
            auto rect = cur_rect;
            cv::Mat rotated_image;
            // 还是无法保证完全在图片框内
            auto rotate_matrix = rotate_image(pyramid_templates[i], rotated_image, centers[i], j);
            rotate_rect(rect, rotate_matrix);
//            draw_template(rotated_image, tpl);
//            cv::imshow(std::string("pyr") + std::string(1, i - '0'), rotated_image);
//            cvWaitKey(0);
        }
    }
    // 打印倒数二层上的金字塔的模板图片以及相应信息
    cv::imshow("origin", pyramid_templates[0]);
    cvWaitKey(0);
    // 就算0度，pyramid_template和原来的图片还是不一样的。
    std::cout << "optimal level: " << (int)kctrp.run_time_npyramid << std::endl;
    for (auto iter = kctrp.tpls.crbegin(); iter != kctrp.tpls.crbegin() + 3; ++iter) {
        auto target = *iter->cbegin();
        auto level = kctrp.tpls.crend() - iter - 1;
        draw_template(pyramid_template[level], target);
        std::cout << "angle step:" << angle_steps[level] << std::endl;
        print_tpl(std::cout,  target) << std::endl;
    }
}


/*
 *  提供给客户端的接口
 * */
char *create_template(const UINT8 *yuv, Koyo_Tool_Contour_Parameter koyo_tool_contour_parameter)
{
    // 获取灰度图
    auto template_image = get_y_from_yuv(yuv, WIDTH, HEIGHT);
    // 设置参数截取模板图片
    cv::Mat template_roi;
    std::vector<cv::Point> rect =  {
            {koyo_tool_contour_parameter.detect_rect_x0, koyo_tool_contour_parameter.detect_rect_y0},
            {koyo_tool_contour_parameter.detect_rect_x1, koyo_tool_contour_parameter.detect_rect_y1},
            {koyo_tool_contour_parameter.detect_rect_x2, koyo_tool_contour_parameter.detect_rect_y2},
            {koyo_tool_contour_parameter.detect_rect_x3, koyo_tool_contour_parameter.detect_rect_y3},
    };
    cutout_template_image(template_image, rect, template_roi);

    // 使用截取出来的图片进行轮廓建立
    Koyo_Contour_Template_Runtime_Param koyo_contour_template_runtime_param;
    do_create_template(template_roi, koyo_tool_contour_parameter, koyo_contour_template_runtime_param);

    // 打包后的template_data是unique_ptr上的指针，调用release来获取原始指针，但是要记得delete []这个内存
    std::cout << "test pack template" << std::endl;
    auto template_data = pack_template(koyo_contour_template_runtime_param);

#ifdef  _DEBUG_
    print_debug_info(pyramid_templates, template_data.get());
#endif

//    cv::imshow("eh" ,template_roi);
//    cvWaitKey(0);
    return template_data.release();
}
