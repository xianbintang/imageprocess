#include "TemplateMatch.h"
#include "types.h"
#include <math.h>
#include <malloc.h>
#include <stdio.h>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <TimeTracker.h>
#include <fstream>

extern double lookupTableX[1001][1001];
extern double lookupTableY[1001][1001];
void saveResult(IMat mat, const char *path);
void saveResultf(IMat mat, const char *path);
void resumeIMatf(const char *path, IMat *mat);
bool IsPointInMatrix(IPoint p, IPoint *rect);

static const float atan2_p1 = 0.9997878412794807f*(float)(180/M_PI);
static const float atan2_p3 = -0.3258083974640975f*(float)(180/M_PI);
static const float atan2_p5 = 0.1555786518463281f*(float)(180/M_PI);
static const float atan2_p7 = -0.04432655554792128f*(float)(180/M_PI);
float IfastAtan2( float y, float x ) {
    float ax = std::abs(x), ay = std::abs(y);
    float a, c, c2;
    if (ax >= ay) {
        c = ay / (ax + (float) DBL_EPSILON);
        c2 = c * c;
        a = (((atan2_p7 * c2 + atan2_p5) * c2 + atan2_p3) * c2 + atan2_p1) * c;
    } else {
        c = ax / (ay + (float) DBL_EPSILON);
        c2 = c * c;
        a = 90.f - (((atan2_p7 * c2 + atan2_p5) * c2 + atan2_p3) * c2 + atan2_p1) * c;
    }
    if (x < 0)
        a = 180.f - a;
    if (y < 0)
        a = 360.f - a;
    return a;
}

int RotateGrayImage(IMat matSrc, IMat **matDst, const double degree);

void HI_Sobel(const IMat *src, IMat *gx, IMat *gy);
void HI_Canny(const IMat *src, IMat *edge);

int CreateGeoMatchModel(TemplateMatch *tpl,  const IMat *src, double maxContrast,double minContrast, IPoint * rect)
{

    IMat *gx = 0;		//Matrix to store X derivative
    IMat *gy = 0;		//Matrix to store Y derivative
    IMat *nmsEdges = 0;		//Matrix to store temp restult
    ISize Ssize;

    // Convert IplImage to Matrix for integer operations
    if(src->type != U8C1)
    {
        return 0;
    }

    // set width and height
    Ssize.width =  src->width;
    Ssize.height= src->height;
    tpl->modelHeight =src->height;		//Save Template height
    tpl->modelWidth =src->width;			//Save Template width

    tpl->noOfCordinates=0;											//initialize
    tpl->cordinates = (IPoint *)malloc(tpl->modelWidth * tpl->modelHeight * sizeof(IPoint));

    tpl->edgeMagnitude = (double *)malloc(tpl->modelWidth * tpl->modelHeight * sizeof(double));		//Allocate memory for edge magnitude for selected points
    tpl->edgeDerivativeX = 	(double *)malloc(tpl->modelWidth * tpl->modelHeight * sizeof(double));		//Allocate memory for edge X derivative for selected points
    tpl->edgeDerivativeY = (double *)malloc(tpl->modelWidth * tpl->modelHeight * sizeof(double));		////Allocate memory for edge Y derivative for selected points


    // Calculate gradient of Template
    gx = fICreateMat( Ssize.height, Ssize.width, S16C1 );		//create Matrix to store X derivative
    gy = fICreateMat( Ssize.height, Ssize.width, S16C1 );		//create Matrix to store Y derivative
    IMat *edge= ICreateMat( Ssize.height, Ssize.width, U8C1 );		//create Matrix to store Y derivative
    HI_Sobel(src, gx, gy);		//gradient in X direction
    HI_Canny(src, edge);

    const short* _sdx;
    const short* _sdy;
    double fdx,fdy;

    double **magMat;
    CreateDoubleMatrix(&magMat ,Ssize);

    int RSum=0,CSum=0;
    int curX,curY;

    // save selected edge information
    for (int i = 0; i < edge->height; ++i) {
        for (int j = 0; j < edge->width; ++j) {
            curX=i;	curY=j;
            IPoint p;
            p.x = i;
            p.y = j;
            _sdx = (short*)(gx->fptr + gx->step*i);
            _sdy = (short*)(gy->fptr + gy->step*i);
            fdx = _sdx[j]; fdy = _sdy[j];
            char e = *((char *)(edge->ptr + edge->step * i + j));

            if(e !=0 && IsPointInMatrix(p, rect))  {
                if(fdx!=0 || fdy!=0)
                {
                    RSum=RSum+curX;	CSum=CSum+curY; // Row sum and column sum for center of gravity

                    tpl->cordinates[tpl->noOfCordinates].x = curX;
                    tpl->cordinates[tpl->noOfCordinates].y = curY;

                    double vector_length = sqrt(fdx * fdx + fdy * fdy);
                    tpl->edgeDerivativeX[tpl->noOfCordinates] = fdx / vector_length;
                    tpl->edgeDerivativeY[tpl->noOfCordinates] = fdy / vector_length;
                    tpl->noOfCordinates++;
                }
            }
        }
    }


    tpl->centerOfGravity.x = RSum / tpl->noOfCordinates; // center of gravity
    tpl->centerOfGravity.y = CSum/tpl->noOfCordinates ;	// center of gravity

// change coordinates to reflect center of gravity
/* 将重心变换到坐标原点 */
    int m;
    for(m=0;m< tpl->noOfCordinates ;m++)
    {
        int temp;

        temp=tpl->cordinates[m].x;
        tpl->cordinates[m].x=temp-tpl->centerOfGravity.x;
        temp= tpl->cordinates[m].y;
        tpl->cordinates[m].y =temp-tpl->centerOfGravity.y;
    }

    ReleaseDoubleMatrix(&magMat ,Ssize.height);

    tpl->modelDefined=true;
    return 1;
}




//allocate memory for doubel matrix
void CreateDoubleMatrix(double ***matrix,ISize size)
{
//    matrix = new double*[size.height];
//    for(int iInd = 0; iInd < size.height; iInd++)
//        matrix[iInd] = new double[size.width];
    (*matrix) = (double **)malloc(size.height * sizeof(double *));
    int i;
    for (i = 0; i < size.height; ++i) {
        (*matrix)[i] = (double *)malloc(size.width * sizeof(double));
    }
}
// release memory
void ReleaseDoubleMatrix(double ***matrix,int size)
{
//    for(int iInd = 0; iInd < size; iInd++)
//        delete[] matrix[iInd];
    int i = 0;
    for (i = 0; i < size; ++i) {
        free((*matrix)[i]);
    }
}





double FindGeoMatchModelRotateTpl(TemplateMatch * tpls, const IMat * srcarr,double minScore,double greediness,IPoint *resultPoint,const double degree, double *resultDegree, Rect region)
{

    TimeTracker t1;
    double resultScore = 0 ;
    /* rotateIt */
//    double dg = 0;
    TemplateMatch * tpl;
    int dg;
    const IMat *src = srcarr;

    IMat *Sdx = 0, *Sdy = 0;

    double partialSum = 0;
    double sumOfCoords = 0;
    double partialScore;
    const short *_Sdx;
    const short *_Sdy;
    int i, j, m;            // count variables
    double iTx, iTy, iSx, iSy;
    double gradMag;
    int curX, curY;

    double **matGradMag;  //Gradient magnitude matrix

//		IMat srcstub, *src = (IMat *) srcarr;
//		src = cvGetMat(src, &srcstub);
    if (src->type != U8C1) {
        return 0;
    }

    // source image size
    ISize Ssize;
    Ssize.width = src->width;
    Ssize.height = src->height;

    CreateDoubleMatrix(&matGradMag, Ssize); // create image to save gradient magnitude  values

    Sdx = fICreateMat(Ssize.height, Ssize.width, S16C1); // X derivatives
    Sdy = fICreateMat(Ssize.height, Ssize.width, S16C1); // y derivatives


    double **edgeX;  //Gradient magnitude matrix
    double **edgeY;  //Gradient magnitude matrix
    CreateDoubleMatrix(&edgeX, Ssize); // create image to save gradient magnitude  values
    CreateDoubleMatrix(&edgeY, Ssize); // create image to save gradient magnitude  values

//        cvSobel(src, Sdx, 1, 0, 3);  // find X derivatives
//        cvSobel(src, Sdy, 0, 1, 3); // find Y derivatives

    HI_Sobel(src, Sdx, Sdy);		//gradient in X direction
//    resumeIMatf("gxgx.txt", Sdx);
//    resumeIMatf("gygy.txt", Sdy);
#if 0
    if (src->height == 120) {
        resumeIMatf("result/normsdx120.txt", Sdx);
        resumeIMatf("result/normsdy120.txt", Sdy);
        std::cout << "hehe" << std::endl;
    } else if (src->height == 240) {
        resumeIMatf("result/normsdx240.txt", Sdx);
        resumeIMatf("result/normsdy240.txt", Sdy);
        std::cout << "hehe" << std::endl;
    }
#endif
    t1.start();
    for (i = 0; i < Ssize.height; i++) {
        _Sdx = (short *) (Sdx->fptr + Sdx->step * (i));
        _Sdy = (short *) (Sdy->fptr + Sdy->step * (i));

        for (j = 0; j < Ssize.width; j++) {
            double ex = lookupTableX[abs(_Sdx[j])][abs(_Sdy[j])];
            double ey = lookupTableY[abs(_Sdx[j])][abs(_Sdy[j])];
            if (_Sdx[j] > 0)
                edgeX[i][j] = ex;
            else
                edgeX[i][j] = -ex;
            if (_Sdy[j] > 0)
                edgeY[i][j] = ey;
            else
                edgeY[i][j] = -ey;
#if 0
            iSx = _Sdx[j];  // X derivative of Source image
            iSy = _Sdy[j];  // Y derivative of Source image
            double vector_length = sqrt(_Sdx[j] * _Sdx[j] + _Sdy[j] * _Sdy[j]);
            edgeX[i][j] = _Sdx[j] / vector_length;
            edgeY[i][j] = _Sdy[j] / vector_length;

            gradMag = sqrt((iSx * iSx) + (iSy * iSy)); //Magnitude = Sqrt(dx^2 +dy^2)

            if (gradMag != 0) // hande divide by zero
                matGradMag[i][j] = 1 / gradMag;   // 1/Sqrt(dx^2 +dy^2)
            else
                matGradMag[i][j] = 0;

#endif
        }
    }

    t1.stop();
//    std::cout << "normalization time: " << t1.duration() << std::endl;

//    std::ofstream fout("counts.txt");
    bool flag = false;

    for (dg = 0; dg < 60; dg += 1) {
        int count = 0;
        TimeTracker tt;
        tt.start();
        tpl = &tpls[dg];
//        std::cout << "num of cor: " << tpl->noOfCordinates << std::endl;
        // stoping criterias to search for model
        double normMinScore = minScore / tpl->noOfCordinates; // precompute minumum score
        double normGreediness =
                ((1 - greediness * minScore) / (1 - greediness)) / tpl->noOfCordinates; // precompute greedniness

        for (i = 0; i < Ssize.height; i += 1) {

            for (j = 0; j < Ssize.width; j += 1) {
                partialSum = 0; // initilize partialSum measure
                for (m = 0; m < tpl->noOfCordinates; m += 1) {
                    curX = i + tpl->cordinates[m].x;    // template X coordinate
                    curY = j + tpl->cordinates[m].y; // template Y coordinate
                    iTx = tpl->edgeDerivativeX[m];    // template X derivative
                    iTy = tpl->edgeDerivativeY[m];    // template Y derivative

                    if (curX < region.x || curY < region.y || curX > region.x + region.width - 1 || curY > region.y + region.height - 1)
                        continue;

                    _Sdx = (short *) (Sdx->fptr + Sdx->step * (curX));
                    _Sdy = (short *) (Sdy->fptr + Sdy->step * (curX));

                    iSx = _Sdx[curY]; // get curresponding  X derivative from source image
                    iSy = _Sdy[curY];// get curresponding  Y derivative from source image

                    if ((iSx != 0 || iSy != 0) && (iTx != 0 || iTy != 0)) {
                        //partial Sum  = Sum of(((Source X derivative* Template X drivative) + Source Y derivative * Template Y derivative)) / Edge magnitude of(Template)* edge magnitude of(Source))
//                        partialSum =
//                                partialSum + ((iSx * iTx) + (iSy * iTy)) * (tpl->edgeMagnitude[m] * matGradMag[curX][curY]);

                        partialSum =
                                partialSum + ((edgeX[curX][curY]* iTx) + (edgeY[curX][curY]* iTy));
                    }

                    sumOfCoords = m + 1;
                    partialScore = partialSum / sumOfCoords;
                    // check termination criteria
                    // if partial score score is less than the score than needed to make the required score at that position
                    // break serching at that coordinate.
                    if (partialScore < (MIN((minScore - 1) + normGreediness * sumOfCoords, normMinScore * sumOfCoords)))
                        break;

                    count++;
                }
//                fout << count << std::endl;
                if (partialScore > resultScore) {
                    resultScore = partialScore; //  Match score
                    resultPoint->x = i;            // result coordinate X
                    resultPoint->y = j;            // result coordinate Y
                    *resultDegree = dg;
                    if (resultScore > 0.26  && sumOfCoords > tpl->noOfCordinates * 0.9 / 4) {
//                        flag = true;
//                        break;
                    }
                }
            }
            if (flag) {
                break;
            }
        }
        if(flag) {
            break;
        }

//        std::cout << "count: " << count << std::endl;
        tt.stop();
//        printf("time used: %dms\n\n", tt.duration());
        // free used resources and return score
//        cvReleaseMat( &Sdx );
//        cvReleaseMat( &Sdy );
    }
    ReleaseDoubleMatrix(&matGradMag ,Ssize.height);

//    printf("i: %d, j:%d dg: %d\n", i, j, dg);
    return resultScore;
}

double FindGeoMatchModelRotateTplInRange(TemplateMatch * tpls, IPoint where, int angle,  const IMat * srcarr,double minScore,double greediness,IPoint *resultPoint,const double degree, double *resultDegree, Rect region)
{

    double resultScore = 0 ;
    /* rotateIt */
//    double dg = 0;
    TimeTracker timeTracker;
    TemplateMatch * tpl;
    const IMat *src = srcarr;

    IMat *Sdx = 0, *Sdy = 0;

    double partialSum = 0;
    double sumOfCoords = 0;
    double partialScore;
    const short *_Sdx;
    const short *_Sdy;
    int i, j, m;            // count variables
    double iTx, iTy, iSx, iSy;
    double gradMag;
    int curX, curY;

    double **matGradMag;  //Gradient magnitude matrix

//		IMat srcstub, *src = (IMat *) srcarr;
//		src = cvGetMat(src, &srcstub);
    if (src->type != U8C1) {
        return 0;
    }

    // source image size
    ISize Ssize;
    Ssize.width = src->width;
    Ssize.height = src->height;

    CreateDoubleMatrix(&matGradMag, Ssize); // create image to save gradient magnitude  values

    Sdx = fICreateMat(Ssize.height, Ssize.width, S16C1); // X derivatives
    Sdy = fICreateMat(Ssize.height, Ssize.width, S16C1); // y derivatives


    double **edgeX;  //Gradient magnitude matrix
    double **edgeY;  //Gradient magnitude matrix
    CreateDoubleMatrix(&edgeX, Ssize); // create image to save gradient magnitude  values
    CreateDoubleMatrix(&edgeY, Ssize); // create image to save gradient magnitude  values

//        cvSobel(src, Sdx, 1, 0, 3);  // find X derivatives
//        cvSobel(src, Sdy, 0, 1, 3); // find Y derivatives

    HI_Sobel(src, Sdx, Sdy);		//gradient in X direction
//    resumeIMatf("gxgx.txt", Sdx);
//    resumeIMatf("gygy.txt", Sdy);

//    saveResultf(*Sdx, "Sdx.txt");
//    saveResultf(*Sdy, "Sdy.txt");
#if 0
    if (src->height == 120) {
        resumeIMatf("result/normsdx120.txt", Sdx);
        resumeIMatf("result/normsdy120.txt", Sdy);
    } else if (src->height == 240) {
        resumeIMatf("result/normsdx240.txt", Sdx);
        resumeIMatf("result/normsdy240.txt", Sdy);
    }
#endif
//    timeTracker.start();
    for (i = 0; i < Ssize.height; i++) {
        _Sdx = (short *) (Sdx->fptr + Sdx->step * (i));
        _Sdy = (short *) (Sdy->fptr + Sdy->step * (i));

        for (j = 0; j < Ssize.width; j++) {
            double ex = lookupTableX[abs(_Sdx[j])][abs(_Sdy[j])];
            double ey = lookupTableY[abs(_Sdx[j])][abs(_Sdy[j])];
            if (_Sdx[j] > 0)
                edgeX[i][j] = ex;
            else
                edgeX[i][j] = -ex;
            if (_Sdy[j] > 0)
                edgeY[i][j] = ey;
            else
                edgeY[i][j] = -ey;
#if 0
            iSx = _Sdx[j];  // X derivative of Source image
            iSy = _Sdy[j];  // Y derivative of Source image
            double vector_length = sqrt(_Sdx[j] * _Sdx[j] + _Sdy[j] * _Sdy[j]);
            edgeX[i][j] = _Sdx[j] / vector_length;
            edgeY[i][j] = _Sdy[j] / vector_length;

            gradMag = sqrt((iSx * iSx) + (iSy * iSy)); //Magnitude = Sqrt(dx^2 +dy^2)

            if (gradMag != 0) // hande divide by zero
                matGradMag[i][j] = 1 / gradMag;   // 1/Sqrt(dx^2 +dy^2)
            else
                matGradMag[i][j] = 0;

#endif
        }
    }
//    timeTracker.stop();
//    std::cout << "half : " << timeTracker.duration() << std::endl;

        TimeTracker tt;
        tt.start();
        tpl = &tpls[angle];
        // stoping criterias to search for model
        double normMinScore = minScore / tpl->noOfCordinates; // precompute minumum score
        double normGreediness =
                ((1 - greediness * minScore) / (1 - greediness)) / tpl->noOfCordinates; // precompute greedniness

        for (i = where.x - 4; i < where.x + 4; i++) {
            for (j = where.y - 4; j < where.y + 4; j++) {
                partialSum = 0; // initilize partialSum measure
                int count = 0;
                for (m = 0; m < tpl->noOfCordinates; m += 1) {
                    count++;
                    curX = i + tpl->cordinates[m].x;    // template X coordinate
                    curY = j + tpl->cordinates[m].y; // template Y coordinate
                    iTx = tpl->edgeDerivativeX[m];    // template X derivative
                    iTy = tpl->edgeDerivativeY[m];    // template Y derivative

                    if (curX < 0 || curY < 0 || curX > Ssize.height - 1 || curY > Ssize.width - 1)
                        continue;

                    _Sdx = (short *) (Sdx->fptr + Sdx->step * (curX));
                    _Sdy = (short *) (Sdy->fptr + Sdy->step * (curX));

                    iSx = _Sdx[curY]; // get curresponding  X derivative from source image
                    iSy = _Sdy[curY];// get curresponding  Y derivative from source image

                    if ((iSx != 0 || iSy != 0) && (iTx != 0 || iTy != 0)) {
                        //partial Sum  = Sum of(((Source X derivative* Template X drivative) + Source Y derivative * Template Y derivative)) / Edge magnitude of(Template)* edge magnitude of(Source))
//                        partialSum =
//                                partialSum + ((iSx * iTx) + (iSy * iTy)) * (tpl->edgeMagnitude[m] * matGradMag[curX][curY]);

                        partialSum =
                                partialSum + ((edgeX[curX][curY]* iTx) + (edgeY[curX][curY]* iTy));
                    }

                    sumOfCoords = m + 1;
                    partialScore = partialSum / sumOfCoords;
                    // check termination criteria
                    // if partial score score is less than the score than needed to make the required score at that position
                    // break serching at that coordinate.
                    if (partialScore < (MIN((minScore - 1) + normGreediness * sumOfCoords, normMinScore * sumOfCoords)))
                        break;

                }
//                std::cout << "no of cordinates computed: " << count << std::endl;
                if (partialScore > resultScore) {
                    resultScore = partialScore; //  Match score
                    resultPoint->x = i;            // result coordinate X
                    resultPoint->y = j;            // result coordinate Y
                    *resultDegree = angle;
                }
            }
        }

        tt.stop();
//        printf("time used: %dms\n\n", tt.duration());
        // free used resources and return score
//        cvReleaseMat( &Sdx );
//        cvReleaseMat( &Sdy );
    ReleaseDoubleMatrix(&matGradMag ,Ssize.height);

    return resultScore;
}

void FindTemplateInPyramid(TemplateMatch tpls[][60], IMat *images[], IPoint *pos, int *dg)
{
    /* 在第三层中找位置 */
    TimeTracker tt;
    TimeTracker total;
    tt.start();
    total.start();
    IPoint point, where;

    Rect region4;
    region4.x = 50;
    region4.y = 90;
    region4.width = 160;
    region4.height = 100;
    Rect region3;
    region3.x = 25;
    region3.y = 45;
    region3.width = 80;
    region3.height = 50;
    double degree;
    double score = FindGeoMatchModelRotateTpl(tpls[3], images[3], 0.95, 0.999995, &point, 30, &degree, region3);

    tt.stop();

    std::cout <<"pyr level 3" << " found at degree: " << degree << " Points: (" << point.x << "," << point.y << ")" << " score: " << score
              << " Time used: " << tt.duration() << std::endl;
//    where.x = point.x * 2;
//    where.y = point.y * 2;
    double angle = degree;
//
//    tt.start();
//    score = FindGeoMatchModelRotateTplInRange(tpls[3],where, 0, images[3], 0.95, 0.999995, &point, 30, &degree);
//
//    tt.stop();
//    std::cout << "pyr level 4 "<< " found at degree: " << degree << " Points: (" << point.x << "," << point.y << ")" << " score: " << score
//              << " Time used: " << tt.duration() << std::endl;
//
    where.x = point.x * 2;
    where.y = point.y * 2;
    angle = degree;

    tt.start();
    score = FindGeoMatchModelRotateTplInRange(tpls[4],where, degree, images[4], 0.95, 0.999995, &point, 30, &degree, region4);
    tt.stop();
    std::cout << "pyr level gray "<< " found at degree: " << degree << " Points: (" << point.x << "," << point.y << ")" << " score: " << score
              << " Time used: " << tt.duration() << std::endl;



    *pos = point;
    *dg = degree;

    total.stop();
    std::cout << "total time: " << total.duration() << std::endl;
}

void downSampleTpl(TemplateMatch * srcTpl, TemplateMatch * dstTpl)
{
    int i, j, k;
    for (i = 0; i < 60; ++i) {
        int noOfCorod = srcTpl[i].noOfCordinates;
        int dstNoOfCord = 0;
        dstTpl[i].cordinates = (IPoint *)malloc(sizeof(IPoint) * srcTpl[i].modelWidth * srcTpl[i].modelHeight);
//        dstTpl[i].edgeDerivativeX= (double *)malloc(sizeof(double) * (noOfCorod / 2 + 1));
//        dstTpl[i].edgeDerivativeY= (double *)malloc(sizeof(double) * (noOfCorod / 2 + 1));

        dstTpl[i].edgeDerivativeX = (double *)malloc(srcTpl[i].modelWidth * srcTpl[i].modelHeight * sizeof(double));		//Allocate memory for edge X derivative for selected points
        dstTpl[i].edgeDerivativeY = (double *)malloc(srcTpl[i].modelWidth * srcTpl[i].modelHeight * sizeof(double));		////Allocate memory for edge Y derivative for selected points

        int RSum=0,CSum=0;
        j = 0;
        for (k = 0; k < noOfCorod; k += 4) {
            dstTpl[i].cordinates[j] = srcTpl[i].cordinates[k];
            dstTpl[i].edgeDerivativeX[j] = srcTpl[i].edgeDerivativeX[k];
            dstTpl[i].edgeDerivativeY[j] = srcTpl[i].edgeDerivativeY[k];
            dstNoOfCord++;
            j++;
        }

        dstTpl[i].centerOfGravity = srcTpl[i].centerOfGravity;
        dstTpl[i].modelDefined = true;
        dstTpl[i].modelHeight = srcTpl[i].modelHeight;
        dstTpl[i].modelWidth= srcTpl[i].modelWidth;
        dstTpl[i].noOfCordinates = dstNoOfCord;

    }
}

float GetCross(IPoint p1, IPoint p2,IPoint p)
{
    return (p2.x - p1.x) * (p.y - p1.y) -(p.x - p1.x) * (p2.y - p1.y);
}
bool IsPointInMatrix(IPoint p, IPoint *rect)
{
    IPoint p1 = rect[0]; //left up
    IPoint p2 = rect[1]; //left bottom
    IPoint p3 = rect[2]; //right bottom
    IPoint p4 = rect[3]; //right up

    return GetCross(p1,p2,p) * GetCross(p3,p4,p) >= 0 && GetCross(p2,p3,p) * GetCross(p4,p1,p) >= 0;
    //return false;
}
