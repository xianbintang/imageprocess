//
// Created by xianb on 2017/5/9.
//

#ifndef KOYO_TEMPLATEMATCH_H
#define KOYO_TEMPLATEMATCH_H

#include "types.h"
#include <stdio.h>
#include <stdbool.h>

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
	IPoint			*cordinates;		//Coordinates array to store model points	model points 也就是所有的边缘点
	int				modelHeight;		//Template height 模板的高度
	int				modelWidth;			//Template width 模板的宽度
	double			*edgeMagnitude;		//gradient magnitude 所有边缘点的梯度值
	double			*edgeDerivativeX;	//gradient in X direction
	double			*edgeDerivativeY;	//gradient in Y direction
	IPoint			centerOfGravity;	//Center of gravity of template 重心

	bool			modelDefined;


//	void DrawContours(IplImage* pImage,Point COG,Scalar,int);
//	void DrawContours(IplImage* pImage,Scalar,int);
} TemplateMatch;

void CreateDoubleMatrix(double ***matrix,ISize size);
void ReleaseDoubleMatrix(double ***matrix,int size);

//double FindGeoMatchModel(const void* srcarr,double minScore,double greediness, IPoint *resultPoint);
int CreateGeoMatchModel( TemplateMatch *tpl,  const IMat *src, double maxContrast,double minContrast);

double FindGeoMatchModel( const TemplateMatch * tpl, const IMat* srcarr,double minScore,double greediness,IPoint *resultPoint,const double degree, double *resultDegree);
#endif //KOYO_TEMPLATEMATCH_H
