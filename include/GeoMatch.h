//***********************************************************************
// Project		    : GeoMatch
// Author           : Shiju P K
// Email			: shijupk@gmail.com
// Created          : 10-01-2010
//
// File Name		: GeoMatch.h
// Last Modified By : Shiju P K
// Last Modified On : 13-07-2010
// Description      : class to implement edge based template matching
//
// Copyright        : (c) . All rights reserved.
//***********************************************************************

#pragma once
#include "opencv/cv.h"
#include <cxcore.h>
#include <highgui.h>
#include <math.h>

class GeoMatch
{
private:
	/* 模板的由下面的几个性质进行描述：
	 * 1. 模板中所有边缘点的位置；
	 * 2. 模板中所有边缘点的梯度值；
	 * 3. 模板中所有边缘点在X和Y方向上的梯度方向；
	 * 4. 模板中边缘点的个数；
	 * 5. 模板的高宽；
	 * 6. 模板的重心。
	 * */
	int				noOfCordinates;		//Number of elements in coordinate array 边缘点的个数
	CvPoint			*cordinates;		//Coordinates array to store model points	model points 也就是所有的边缘点
	int				modelHeight;		//Template height 模板的高度
	int				modelWidth;			//Template width 模板的宽度
	double			*edgeMagnitude;		//gradient magnitude 所有边缘点的梯度值
	double			*edgeDerivativeX;	//gradient in X direction
	double			*edgeDerivativeY;	//gradient in Y direction	
	CvPoint			centerOfGravity;	//Center of gravity of template 重心
	
	bool			modelDefined;
	
	void CreateDoubleMatrix(double **&matrix,CvSize size);
	void ReleaseDoubleMatrix(double **&matrix,int size);
public:
	GeoMatch(void);
	GeoMatch(const void* templateArr);
	~GeoMatch(void);

    double FindGeoMatchModel(const void* srcarr,double minScore,double greediness,CvPoint *resultPoint,const double degree, double *resultDegree);
	int CreateGeoMatchModel(const void* templateArr,double,double);
	double FindGeoMatchModel(const void* srcarr,double minScore,double greediness, CvPoint *resultPoint);
	void DrawContours(IplImage* pImage,CvPoint COG,CvScalar,int);
	void DrawContours(IplImage* pImage,CvScalar,int);
    void saveTpl();
    void ResumeTpl();
};