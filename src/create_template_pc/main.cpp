//
// Created by xianb on 2017/11/28.
//


#include <opencv2/core/mat.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv/cv.hpp>
#include <windef.h>
#include <fstream>
#include "contour_detection.h"

void saveMat(cv::Mat mat, const char *path) {
    FILE *fp = fopen(path, "w");
    int i,j;
    for (i = 0; i < mat.rows; ++i) {
        for (j = 0; j < mat.cols; ++j) {
//            fprintf(fp, "%d ", (mat.ptr + i * mat.step)[j]);
            fprintf(fp, "%d ", mat.at<uchar>(i, j));
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}
void saveMatf(cv::Mat mat, const char *path) {
    FILE *fp = fopen(path, "w");
    int i,j;
    for (i = 0; i < mat.rows; ++i) {
        for (j = 0; j < mat.cols; ++j) {
//            fprintf(fp, "%d ", (mat.ptr + i * mat.step)[j]);
            fprintf(fp, "%d ", mat.at<short>(i, j));
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

void init_contour_parameter(Koyo_Tool_Contour_Parameter &koyo_tool_contour_parameter)
{
    koyo_tool_contour_parameter.algo_strategy = 1;


#define _CPU_TEST45_BAD_

#ifdef _VI42_
    koyo_tool_contour_parameter.detect_rect_x0 = 279;
    koyo_tool_contour_parameter.detect_rect_y0 = 291;

    koyo_tool_contour_parameter.detect_rect_x1 = 367;
    koyo_tool_contour_parameter.detect_rect_y1 = 418;

    koyo_tool_contour_parameter.detect_rect_x2 = 515;
    koyo_tool_contour_parameter.detect_rect_y2 = 314;

    koyo_tool_contour_parameter.detect_rect_x3 = 427;
    koyo_tool_contour_parameter.detect_rect_y3 = 189;
#endif

#ifdef _VI7_
    koyo_tool_contour_parameter.detect_rect_x0 = 176;
    koyo_tool_contour_parameter.detect_rect_y0 = 125;

    koyo_tool_contour_parameter.detect_rect_x1 = 176;
    koyo_tool_contour_parameter.detect_rect_y1 = 339;

    koyo_tool_contour_parameter.detect_rect_x2 = 403;
    koyo_tool_contour_parameter.detect_rect_y2 = 339;

    koyo_tool_contour_parameter.detect_rect_x3 = 403;
    koyo_tool_contour_parameter.detect_rect_y3 = 125;
#endif

#ifdef _VI6_
    koyo_tool_contour_parameter.detect_rect_x0 = 366;
    koyo_tool_contour_parameter.detect_rect_y0 = 259;

    koyo_tool_contour_parameter.detect_rect_x1 = 366;
    koyo_tool_contour_parameter.detect_rect_y1 = 459;

    koyo_tool_contour_parameter.detect_rect_x2 = 579;
    koyo_tool_contour_parameter.detect_rect_y2 = 459;

    koyo_tool_contour_parameter.detect_rect_x3 = 579;
    koyo_tool_contour_parameter.detect_rect_y3 = 259;
#endif

#ifdef _VI5_
    koyo_tool_contour_parameter.detect_rect_x0 = 265;
    koyo_tool_contour_parameter.detect_rect_y0 = 211;

    koyo_tool_contour_parameter.detect_rect_x1 = 265;
    koyo_tool_contour_parameter.detect_rect_y1 = 421;

    koyo_tool_contour_parameter.detect_rect_x2 = 485;
    koyo_tool_contour_parameter.detect_rect_y2 = 421;

    koyo_tool_contour_parameter.detect_rect_x3 = 485;
    koyo_tool_contour_parameter.detect_rect_y3 = 211;
#endif

#ifdef _CPU_TEST_
    koyo_tool_contour_parameter.detect_rect_x0 = 165;
    koyo_tool_contour_parameter.detect_rect_y0 = 85;

    koyo_tool_contour_parameter.detect_rect_x1 = 165;
    koyo_tool_contour_parameter.detect_rect_y1 = 445;

    koyo_tool_contour_parameter.detect_rect_x2 = 514;
    koyo_tool_contour_parameter.detect_rect_y2 = 445;

    koyo_tool_contour_parameter.detect_rect_x3 = 514;
    koyo_tool_contour_parameter.detect_rect_y3 = 85;
#endif

#ifdef _CPU_TEST45_
    koyo_tool_contour_parameter.detect_rect_x0 = 134;
    koyo_tool_contour_parameter.detect_rect_y0 = 246;

    koyo_tool_contour_parameter.detect_rect_x1 = 361;
    koyo_tool_contour_parameter.detect_rect_y1 = 471;

    koyo_tool_contour_parameter.detect_rect_x2 = 578;
    koyo_tool_contour_parameter.detect_rect_y2 = 255;

    koyo_tool_contour_parameter.detect_rect_x3 = 350;
    koyo_tool_contour_parameter.detect_rect_y3 = 29;
#endif

#ifdef _CPU_TEST45_BAD_
    koyo_tool_contour_parameter.detect_rect_x0 = 200;
    koyo_tool_contour_parameter.detect_rect_y0 = 200;

    koyo_tool_contour_parameter.detect_rect_x1 = 200;
    koyo_tool_contour_parameter.detect_rect_y1 = 400;

    koyo_tool_contour_parameter.detect_rect_x2 = 400;
    koyo_tool_contour_parameter.detect_rect_y2 = 400;

    koyo_tool_contour_parameter.detect_rect_x3 = 400;
    koyo_tool_contour_parameter.detect_rect_y3 = 200;
#endif

#ifdef _KOYO_TEST_
    koyo_tool_contour_parameter.detect_rect_x0 = 150;
    koyo_tool_contour_parameter.detect_rect_y0 = 182;

    koyo_tool_contour_parameter.detect_rect_x1 = 134;
    koyo_tool_contour_parameter.detect_rect_y1 = 255;

    koyo_tool_contour_parameter.detect_rect_x2 = 393;
    koyo_tool_contour_parameter.detect_rect_y2 = 313;

    koyo_tool_contour_parameter.detect_rect_x3 = 407;
    koyo_tool_contour_parameter.detect_rect_y3 = 238;
#endif

#ifdef _GEAR_TEST_
    koyo_tool_contour_parameter.detect_rect_x0 = 142;
    koyo_tool_contour_parameter.detect_rect_y0 = 51;

    koyo_tool_contour_parameter.detect_rect_x1 = 143;
    koyo_tool_contour_parameter.detect_rect_y1 = 402;

    koyo_tool_contour_parameter.detect_rect_x2 = 491;
    koyo_tool_contour_parameter.detect_rect_y2 = 402;

    koyo_tool_contour_parameter.detect_rect_x3 = 491;
    koyo_tool_contour_parameter.detect_rect_y3 = 53;
#endif


#ifdef _IRON_TEST_
    koyo_tool_contour_parameter.detect_rect_x0 = 176;
    koyo_tool_contour_parameter.detect_rect_y0 = 109;

    koyo_tool_contour_parameter.detect_rect_x1 = 176;
    koyo_tool_contour_parameter.detect_rect_y1 = 346;

    koyo_tool_contour_parameter.detect_rect_x2 = 413;
    koyo_tool_contour_parameter.detect_rect_y2 = 346;

    koyo_tool_contour_parameter.detect_rect_x3 = 413;
    koyo_tool_contour_parameter.detect_rect_y3 = 109;
#endif

#ifdef _TIMG_TEST_
    koyo_tool_contour_parameter.detect_rect_x0 = 148;
    koyo_tool_contour_parameter.detect_rect_y0 = 27;

    koyo_tool_contour_parameter.detect_rect_x1 = 148;
    koyo_tool_contour_parameter.detect_rect_y1 = 381;

    koyo_tool_contour_parameter.detect_rect_x2 = 548;
    koyo_tool_contour_parameter.detect_rect_y2 = 381;

    koyo_tool_contour_parameter.detect_rect_x3 = 548;
    koyo_tool_contour_parameter.detect_rect_y3 = 26;
#endif
#ifdef _GEAR2_TEST_
    koyo_tool_contour_parameter.detect_rect_x0 = 225;
    koyo_tool_contour_parameter.detect_rect_y0 = 100;

    koyo_tool_contour_parameter.detect_rect_x1 = 225;
    koyo_tool_contour_parameter.detect_rect_y1 = 351;

    koyo_tool_contour_parameter.detect_rect_x2 = 519;
    koyo_tool_contour_parameter.detect_rect_y2 = 351;

    koyo_tool_contour_parameter.detect_rect_x3 = 519;
    koyo_tool_contour_parameter.detect_rect_y3 = 100;
#endif
#ifdef _CHUANKOU_TEST_
    koyo_tool_contour_parameter.detect_rect_x0 = 284;
    koyo_tool_contour_parameter.detect_rect_y0 = 136;

    koyo_tool_contour_parameter.detect_rect_x1 = 224;
    koyo_tool_contour_parameter.detect_rect_y1 = 389;

    koyo_tool_contour_parameter.detect_rect_x2 = 478;
    koyo_tool_contour_parameter.detect_rect_y2 = 454;

    koyo_tool_contour_parameter.detect_rect_x3 = 540;
    koyo_tool_contour_parameter.detect_rect_y3 = 199;
#endif
#ifdef _LANCHUANKOU_TEST_
    koyo_tool_contour_parameter.detect_rect_x0 = 211;
    koyo_tool_contour_parameter.detect_rect_y0 = 53;

    koyo_tool_contour_parameter.detect_rect_x1 = 165;
    koyo_tool_contour_parameter.detect_rect_y1 = 374;

    koyo_tool_contour_parameter.detect_rect_x2 = 592;
    koyo_tool_contour_parameter.detect_rect_y2 = 434;

    koyo_tool_contour_parameter.detect_rect_x3 = 635;
    koyo_tool_contour_parameter.detect_rect_y3 = 122;
#endif

    koyo_tool_contour_parameter.sensitivity = CONTOUR_ACCURACY_MEDIUM;
    koyo_tool_contour_parameter.angle_range = 180;


#ifdef _CPU_TEST45_BAD_
    koyo_tool_contour_parameter.ext_rect_x = 200;
    koyo_tool_contour_parameter.ext_rect_y = 200;
    koyo_tool_contour_parameter.ext_rect_width = 200;
    koyo_tool_contour_parameter.ext_rect_height = 200;
#endif

#ifdef _CPU_TEST45_
    koyo_tool_contour_parameter.ext_rect_x = 131;
    koyo_tool_contour_parameter.ext_rect_y = 29;
    koyo_tool_contour_parameter.ext_rect_width = 445;
    koyo_tool_contour_parameter.ext_rect_height = 445;
#endif

    UINT8 *bitmap = new UINT8[(sizeof(uchar) * koyo_tool_contour_parameter.ext_rect_height * koyo_tool_contour_parameter.ext_rect_width)];
    int num;
    int k = 0;
    std::ifstream fin("data//contour_erased.txt");
    if(!fin.is_open()) {
        exit(-1);
    }
    while(fin >> num && k < koyo_tool_contour_parameter.ext_rect_width * koyo_tool_contour_parameter.ext_rect_height) {
        bitmap[k++] = num;
    }
//     在这里设置bitmap, 这里的bitmap和客户端传来的有区别了，获取到客户端的后就不用原来的了。
    koyo_tool_contour_parameter.bitmaps = bitmap;
}


int main(int argc, char **argv)
{

    std::string filename(argv[1]);
    cv::Mat template_image;
    UINT8 *buf = nullptr;
    if(filename.substr(filename.size() - 3, 3) == "yuv") {
        FILE *yuv_file = fopen(filename.c_str(), "rb+");
        buf = new UINT8[WIDTH * HEIGHT];
        fread(buf, WIDTH * HEIGHT, 1, yuv_file);
    } else {
        template_image = cv::imread(argv[1], 0);
        buf = template_image.data;
    }

    cv::Mat cleanedImage = cv::imread("data//bitmap_erased1.jpg", 0);
    cv::Mat template_roi_ext;
    // todo 换成从bitmap中读取
    // todo 需要删除掉如果是直接从bitmap中取出来的就不做canny了
    cv::Canny(cleanedImage, template_roi_ext, 10, 80);
    saveMat(template_roi_ext, "data//contour_erased.txt");


    Koyo_Tool_Contour_Parameter koyo_tool_contour_parameter;
    init_contour_parameter(koyo_tool_contour_parameter);

    // 返回指向需要被发送的内存缓冲区的指针
    create_template(buf, koyo_tool_contour_parameter);

    return 0;
}


