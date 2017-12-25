//
// Created by yudong on 2017/12/20.
//

#ifndef KOYO_CONTOUR_DETECTION_H
#define KOYO_CONTOUR_DETECTION_H


const int MAX_NUM_PYRAMID = 5;
const int MAX_DEGREE = 360;
using UINT16 = unsigned short;
using INT8 = signed char;
using UINT8 = unsigned char;
using UINT32 = unsigned int;
const int MIN_CONTOUR_PYRA = 100;
const int WIDTH = 640;
const int HEIGHT = 480;


#define CONTOUR_ACCURACY_LOW      0
#define CONTOUR_ACCURACY_MEDIUM   1
#define CONTOUR_ACCURACY_HIGH     2
#define LOW_THRESHOLD 10
#define HIGH_THRESHOLD 80


typedef struct TemplateMatch
{
    /* 模板的由下面的几个性质进行描述：
     * 1. 模板中所有边缘点的位置；
     * 2. 模板中所有边缘点的梯度值；
     * 3. 模板中所有边缘点在X和Y方向上的梯度方向；
     * 4. 模板中边缘点的个数；
     * 5. 模板的高宽；
     * 6. 模板的重心。
     * */
    UINT8 modelDefined;
    UINT32 noOfCordinates;
    UINT16 modelHeight;
    UINT16 modelWidth;
    cv::Point centerOfGravity;
    std::vector<cv::Point> cordinates;
    std::vector<float> edgeDerivativeX;
    std::vector<float> edgeDerivativeY;
} TemplateStruct;


// 把这些东西发送给传感器
typedef struct KOYO_CONTOUR_TEMPLATE_RUNTIME_PARAM{
    UINT8 run_time_npyramid;
    std::vector<float> search_angel_nstep;
    std::vector<UINT16> search_rect_width;
    std::vector<std::vector<TemplateStruct>> tpls;
} Koyo_Contour_Template_Runtime_Param;


// 创建某层待测图金字塔的数据结构
typedef struct SampleMatch{
    UINT16 height;
    UINT16 width;
    std::vector<std::vector<bool>> hasPoint;
    std::vector<std::vector<float>> edgeDerivativeX;
    std::vector<std::vector<float>> edgeDerivativeY;
    cv::Point center;
    float angle;
    float similarity;
}SampleStruct;


// 重要的位置信息数据结构
typedef struct ImportantInfo{
    cv::Point center;
    float angle;
}ImportantInfo;



// 引入所需的函数声明
int match_sample_template(const cv::Mat &sample_image, Koyo_Contour_Template_Runtime_Param &kctrp);

#endif //KOYO_CONTOUR_DETECTION_H
