//
// Created by xianb on 2017/11/28.
//

#ifndef KOYO_CONTOUR_DETECTION_H
#define KOYO_CONTOUR_DETECTION_H

#include <types.h>
#include <opencv2/core/mat.hpp>

const int MAX_NUM_PYRAMID = 5;
const int MAX_DEGREE = 5;

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
    IPoint			*cordinates;		//Coordinates array to store mo hjel points	model points 也就是所有的边缘点
    int				modelHeight;		//Template height 模板的高度
    int				modelWidth;			//Template width 模板的宽度
    double			*edgeMagnitude;		//gradient magnitude 所有边缘点的梯度值
    double			*edgeDerivativeX;	//gradient in X direction
    double			*edgeDerivativeY;	//gradient in Y direction
    IPoint			centerOfGravity;	//Center of gravity of template 重心
    IPoint          rect[4];            //相对于左上点的坐标
    bool modelDefined;


//	void DrawContours(IplImage* pImage,Point COG,Scalar,int);
//	void DrawContours(IplImage* pImage,Scalar,int);
} TemplateStruct;

int create_template(const cv::Mat &src);

#endif //KOYO_CONTOUR_DETECTION_H
