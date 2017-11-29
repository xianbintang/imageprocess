//
// Created by xianb on 2017/11/28.
//

#ifndef KOYO_CONTOUR_DETECTION_H
#define KOYO_CONTOUR_DETECTION_H

#include <types.h>
#include <opencv2/core/mat.hpp>

const int MAX_NUM_PYRAMID = 5;
const int MAX_DEGREE = 360;
using UINT16 = unsigned short;
using INT8 = signed char;
using UINT8 = unsigned char;
using UINT32 = unsigned int;
const int MIN_CONTOUR_PYRA = 100;


#define CONTOUR_ACCURACY_LOW      0
#define CONTOUR_ACCURACY_MEDIUM   1
#define CONTOUR_ACCURACY_HIGH     2

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
    int				noOfCordinates;		//Number of elements in coordinate array 边缘点的个数
    cv::Point			*cordinates;		//Coordinates array to store mo hjel points	model points 也就是所有的边缘点
    int				modelHeight;		//Template height 模板的高度
    int				modelWidth;			//Template width 模板的宽度
//    double			*edgeMagnitude;		//gradient magnitude 所有边缘点的梯度值
    double			*edgeDerivativeX;	//gradient in X direction
    double			*edgeDerivativeY;	//gradient in Y direction
    cv::Point			centerOfGravity;	//Center of gravity of template 重心
    cv::Point          rect[4];            //相对于左上点的坐标
    bool modelDefined;
} TemplateStruct;

// 工具级别的参数, 客户端传下来的
typedef struct KOYO_TOOL_CONTOUR_PARAMETER_
{
    UINT16 detect_region_type;   //检测区域形状
    INT8   tool_name[32];        //工具名称

    /* 矩形检测框参数 */
    UINT16 detect_rect_x0;        //检测矩形框左上点横坐标
    UINT16 detect_rect_y0;        //检测矩形框左上点纵坐标
    UINT16 detect_rect_x1;        //检测矩形框左下点横坐标
    UINT16 detect_rect_y1;        //检测矩形框左下点纵坐标
    UINT16 detect_rect_x2;        //检测矩形框右下点横坐标
    UINT16 detect_rect_y2;        //检测矩形框右下点纵坐标
    UINT16 detect_rect_x3;        //检测矩形框右上点横坐标
    UINT16 detect_rect_y3;        //检测矩形框右上点纵坐标

    /* 检测圆参数 */
    UINT16 detect_circ_x;        //检测圆心横坐标
    UINT16 detect_circ_y;        //检测圆心纵坐标
    UINT16 detect_circ_radius;   //检测圆半径

    /* 外接矩形区域参数 */
    UINT16 ext_rect_x;	         //检测外接矩形起始点横坐标
    UINT16 ext_rect_y;           //检测外接矩形起始点纵坐标
    UINT16 ext_rect_width;       //检测外接矩形宽度
    UINT16 ext_rect_height;      //检测外接矩形高度

    /* 搜索区域参数 */
    UINT16 search_rect_x;        //搜索区域起始点横坐标
    UINT16 search_rect_y;        //搜索区域起始点纵坐标
    UINT16 search_rect_width;    //搜索区域宽度
    UINT16 search_rect_height;   //搜索区域高度

    UINT8  sensitivity;          //搜索灵敏度
    UINT8  algo_strategy;
    UINT16 angle_range;          //搜索角度范围

    UINT16 top_threshold;        //处理结果阈值上限
    UINT16 bot_threshold;        //处理结果阈值下限

    UINT32 bitmap_size;          //位图大小,单位是bit
    UINT16 reserved;
    INT8 bitmap_path[128];
    UINT8  bitmaps[0];           //检测区域外接矩形位图,
} Koyo_Tool_Contour_Parameter;

// 把这些东西发送给传感器

typedef struct KOYO_CONTOUR_TEMPLATE_RUNTIME_PARAM{
    UINT8 run_time_npyramid;
    std::vector<float> search_angel_nstep;
    std::vector<std::vector<TemplateStruct>> tpls;
} Koyo_Contour_Template_Runtime_Param;

int cutout_template_image(const cv::Mat &template_image, std::vector<cv::Point> rect, cv::Mat &interesting_template);
int create_template(const cv::Mat &src, Koyo_Tool_Contour_Parameter koyo_tool_contour_parameter);
#endif //KOYO_CONTOUR_DETECTION_H
