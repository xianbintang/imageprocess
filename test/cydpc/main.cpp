#include <opencv2/core/mat.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv/cv.hpp>
#include <iostream>
#include <fstream>
#include "contour_detection.h"
#include "TimeTracker.h"

using namespace std;

// unpack模板信息到kctrp变量上
int unpack_template(std::string filename, Koyo_Contour_Template_Runtime_Param &koyo_contour_template_runtime_param)
{

    TimeTracker timeTracker;
    timeTracker.start();
    FILE *infile = fopen(filename.c_str(), "rb+");
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
    char *buf = pb;

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
                float edgeX = *((float *) &buf[index]);
                index += sizeof(float);
                tpl.edgeDerivativeX.push_back(edgeX);
            }

            for (std::size_t i = 0; i < tpl.noOfCordinates; ++i) {
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


static void bitmap2Mat(cv::Mat &dst, std::vector<UINT8> &bitmap, UINT16 width, UINT16 height){
    for(int i = 0; i < height; ++i){
        for(int j = 0; j < width; ++j){
            dst.at<uchar>(i, j) = bitmap[i * width + j];
        }
    }
}


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


#define WIDTH_TXT  640
#define HEIGHT_TXT 480

int main(int argc, char *argv[])
{
    // 导入待测图文件转化为待测图灰度图
    std::string filename(argv[1]);

    cv::Mat sample_image_gray;

    if(filename.substr(filename.size() - 3, 3) == "txt") {
        vector<UINT8> values;
        ifstream infile(filename.c_str());
        while(infile){
            int data;
            infile >> data;
            UINT8 data1 = static_cast<UINT8>(data);
            values.push_back(data1);
        }
        sample_image_gray.create(HEIGHT_TXT, WIDTH_TXT, CV_8UC1);
        bitmap2Mat(sample_image_gray, values, WIDTH_TXT, HEIGHT_TXT);
//        cv::imshow("raw picture", sample_image_gray);
//        cv::waitKey(0);
    } else {
        UINT8 *buf = nullptr;
        cv::Mat sample_image;
        sample_image = cv::imread(filename, 0);
        buf = sample_image.data;
        sample_image_gray = get_y_from_yuv(buf, WIDTH, HEIGHT);
    }

    // 将模板信息导入到变量kctrp上
    Koyo_Contour_Template_Runtime_Param kctrp;
    unpack_template("template_file", kctrp);

    // 进行待测图和模板的匹配
    match_sample_template(sample_image_gray, kctrp);

    return 0;
}
