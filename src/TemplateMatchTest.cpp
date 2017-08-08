//
// Created by xianb on 2017/5/11.
//

#include <TemplateMatch.h>
#include "types.h"
#include <stdlib.h>

#include <highgui.h>
#include <cv.h>
#include <iostream>
#include <TimeTracker.h>
#include <opencv/cv.hpp>
//int CreateGeoMatchModel(TemplateMatch *tpl, const IMat *src, double maxContrast, double minContrast)

/*
 * typedef struct hiIVE_IMAGE_S
{
IVE_IMAGE_TYPE_E enType;
HI_U32 u32PhyAddr[3];
HI_U8 *pu8VirAddr[3];
HI_U16 u16Stride[3];
HI_U16 u16Width;
HI_U16 u16Height;
HI_U16 u16Reserved; /*Can be used such as elemSize
}IVE_IMAGE_S;
*/


void ResumeTpl(TemplateMatch *tpl);
void saveResultf(IMat mat, const char *path) ;
void saveResult(IMat mat, const char *path) ;
void saveTpl(TemplateMatch *tpl);

void saveTpls(TemplateMatch *tpls, int numOfTpls, char**) ;
void ResumeTpls(TemplateMatch *tpls, int numOfTpls, char**);
void test();
void resumeIMatf(const char *path, IMat *mat);
void downSampleTpl(TemplateMatch * srcTpl, TemplateMatch * dstTpl);
void FindTemplateInPyramid(TemplateMatch tpls[][60], IMat *images[], IPoint *where, int *dg);
double FindGeoMatchModelRotateTplInRange(TemplateMatch * tpls, IPoint where, int angle,  const IMat * srcarr,double minScore,double greediness,IPoint *resultPoint,const double degree, double *resultDegree);
double FindGeoMatchModelRotateTpl(TemplateMatch * tpls, const IMat * srcarr,double minScore,double greediness,IPoint *resultPoint,const double degree, double *resultDegree);
int RotateGrayImage(IMat matSrc, IMat **matDst, const double degree, IPoint * rect);
void rotateTpl(TemplateMatch *tpl, int iWidth, int iHeight, double degree);
void Ipl2IMat(IplImage *img, IMat * src);
void IMat2Ipl(const IMat * src, IplImage *img);
void resumeIMat(const char *path, IMat *mat);
int CreateGeoMatchModel(TemplateMatch *tpl,  const IMat *src, double maxContrast,double minContrast, IPoint *rect);
void deflatRect(IPoint *rect, int w, int h);

double lookupTableX[1001][1001];
double lookupTableY[1001][1001];

void TemplateMatchTest()
{
//    TemplateMatch tpl;
//    IMat src;
//    IVE_IMAGE_S imgSrc;
//    imgSrc.u16Width = 139;
//    imgSrc.u16Height = 44;
//    CreateGeoMatchModel(&tpl, &src, 170, 30);
}


void DrawContours(IplImage* source, IPoint COG,TemplateMatch tpl, CvScalar color,int lineWidth)
{
    CvPoint point;
	point.y=COG.x;
	point.x=COG.y;
	for(int i=0; i<tpl.noOfCordinates; i++)
	{
		point.y=tpl.cordinates[i].x + COG.x;
		point.x=tpl.cordinates[i].y + COG.y;
		cvLine(source,point,point,color,lineWidth);
	}

}
int main(int argc, char ** argv) {
    for(int i = 0; i < 1001; i++) {
        for (int j = 0; j < 1001; ++j) {
            double length;
            length = sqrt(i * i + j * j);
            lookupTableX[i][j] = 1.0 * i / length;
            lookupTableY[i][j] = 1.0 * j / length;
        }
    }
    std::cout << "sizeof lookup Table: " << sizeof(lookupTableY) << std::endl;

    IPoint rect[4];
    TemplateMatch tpl;
    TemplateMatch tpls[5][60];
    IplImage *image = cvLoadImage(argv[1], -1);
    IplImage *Tgray= cvLoadImage("images/testImages/Tgray.bmp", CV_LOAD_IMAGE_GRAYSCALE);
    IplImage *Tlevel4= cvLoadImage("images/testImages/Tlevel4.bmp", CV_LOAD_IMAGE_GRAYSCALE);
    IplImage *Tlevel3= cvLoadImage("images/testImages/Tlevel3.bmp", -1);
    IplImage *Tlevel2= cvLoadImage("images/testImages/Tlevel2.bmp", -1);
    IplImage *Tlevel1= cvLoadImage("images/testImages/Tlevel1.bmp", -1);


    cvSmooth(Tgray,Tgray, CV_GAUSSIAN,3,3,2);
    cvSmooth(Tlevel4,Tlevel4, CV_GAUSSIAN,3,3,2);
    cvSmooth(Tlevel3,Tlevel3, CV_GAUSSIAN,3,3,2);
    cvSmooth(Tlevel2,Tlevel2, CV_GAUSSIAN,3,3,2);
    cvSmooth(Tlevel1,Tlevel1, CV_GAUSSIAN,3,3,2);

    IMat * TgrayM = ICreateMat(Tgray->height, Tgray->width, U8C1);
    Ipl2IMat(Tgray, TgrayM);
    saveResult(*TgrayM, "result/Tgray.txt");

    IMat * Tlevel4M = ICreateMat(Tlevel4->height, Tlevel4->width, U8C1);
    Ipl2IMat(Tlevel4, Tlevel4M);
    saveResult(*Tlevel4M, "result/Tlevel4.txt");

#if 0
    IMat * Tlevel3M = ICreateMat(Tlevel3->height, Tlevel3->width, U8C1);
    Ipl2IMat(Tlevel3, Tlevel3M);
    saveResult(*Tlevel3M, "Tlevel3.txt");

    IMat * Tlevel2M = ICreateMat(Tlevel2->height, Tlevel2->width, U8C1);
    Ipl2IMat(Tlevel2, Tlevel2M);
    saveResult(*Tlevel2M, "Tlevel2.txt");
#endif

    IplImage *srcImage;
    srcImage = cvCreateImage(cvSize(image->width, image->height), IPL_DEPTH_8U, 1);
    if (image->nChannels == 3)
        cvCvtColor(image, srcImage, CV_RGB2GRAY);

    IMat *src= ICreateMat(srcImage->height, srcImage->width, U8C1);
    Ipl2IMat(srcImage, src);

    IMat *tg= ICreateMat(Tgray->height, Tgray->width, U8C1);
    Ipl2IMat(Tgray, tg);
//    std::cout << "tg: " << Tgray->width << " " << Tgray->height << std::endl;

    IMat *t4= ICreateMat(Tlevel4->height, Tlevel4->width, U8C1);
    Ipl2IMat(Tlevel4, t4);
//    std::cout << "t4: " << Tlevel4->width << " " << Tlevel4->height << std::endl;

    IMat *t3= ICreateMat(Tlevel3->height, Tlevel3->width, U8C1);
    Ipl2IMat(Tlevel3, t3);
//    std::cout << "t3: " << Tlevel3->width << " " << Tlevel3->height << std::endl;

    IMat *t2= ICreateMat(Tlevel2->height, Tlevel2->width, U8C1);
    Ipl2IMat(Tlevel2, t2);
//    std::cout << "t2: " << Tlevel2->width << " " << Tlevel2->height << std::endl;

    IMat *t1= ICreateMat(Tlevel1->height, Tlevel1->width, U8C1);
    Ipl2IMat(Tlevel1, t1);
//    std::cout << "t1: " << Tlevel1->width << " " << Tlevel1->height << std::endl;

//    tg = ICreateMat(76, 144, U8C1);
//    resumeIMat("images/testImages/Tgray.txt", tg);
    for (int i = -29; i < 30; ++i) {
        IMat *matDst;
        rect[0].x = 8;
        rect[0].y = 8;

        rect[1].x = 8;
        rect[1].y = 40;

        rect[2].x = 135;
        rect[2].y = 40;

        rect[3].x = 135;
        rect[3].y = 8;

        RotateGrayImage(*tg, &matDst, i, rect);
//        deflatRect(rect, 8, 8);
        CreateGeoMatchModel(&tpls[4][i + 29], matDst, 30, 150, rect);
    }
    char *filenames[5] = {
    "tpls/tpls4/cordinates.txt",
    "tpls/tpls4/edgeMag.txt",
    "tpls/tpls4/ex.txt",
    "tpls/tpls4/ey.txt",
    "tpls/tpls4/tplArgs.txt"};

//    saveTpls(tpls[4], 60, filenames);
//    ResumeTpls(tpls[4], 60, filenames);
    //god

    for (int i = -29; i < 30; ++i) {
        IMat *matDst;
        rect[0].x = 3;
        rect[0].y = 3;

        rect[1].x = 3;
        rect[1].y = 18;

        rect[2].x = 75;
        rect[2].y = 18;

        rect[3].x = 75;
        rect[3].y = 3;

        RotateGrayImage(*t4, &matDst, i, rect);
//        deflatRect(rect, 8, 8);
        CreateGeoMatchModel(&tpls[3][i + 29], matDst, 30, 150, rect);
//        std::cout << tpls[3][i].noOfCordinates << std::endl;
    }
    std::cout << tpls[3][0].noOfCordinates <<std::endl;
    char *filenames1[5] = {
    "tpls/tpls3/cordinates.txt",
    "tpls/tpls3/edgeMag.txt",
    "tpls/tpls3/ex.txt",
    "tpls/tpls3/ey.txt",
    "tpls/tpls3/tplArgs.txt"};

//    saveTpls(tpls[3], 60, filenames1);
//    ResumeTpls(tpls[3], 60, filenames1);
//    for (int i = 0; i < 60; ++i) {
//        IMat *matDst;
//        RotateGrayImage(*t3, &matDst, i, rect);
//        deflatRect(rect, 8, 8);
//        CreateGeoMatchModel(&tpls[2][i], matDst, 30, 170, rect);
//    }
//    downSampleTpl(tpls[3], tpls[2]);
//    for (int k = 0; k < 60; ++k) {
//        printf("tpls2 %d: %d\n", k, tpls[2][k].noOfCordinates);
//    }

//    TemplateMatch tmpTpl[60];
//    downSampleTpl(tpls[2], tmpTpl);
//    memcpy(tpls[2], tmpTpl, sizeof(tmpTpl));
//    std::cout << tpls[2][0].noOfCordinates <<std::endl;
#if 0
    for (int i = 0; i < 60; ++i) {
        IMat *matDst;
        RotateGrayImage(*t2, &matDst, i, rect);
        CreateGeoMatchModel(&tpls[1][i], matDst, 30, 170, rect);
    }
    std::cout << tpls[1][0].noOfCordinates <<std::endl;
#endif
//    for (int i = 0; i < 60; ++i) {
//        IMat *matDst;
//        RotateGrayImage(*t1, &matDst, i);
//        CreateGeoMatchModel(&tpls[0][i], matDst, 30, 170);
//    }
//    std::cout << tpls[0][0].noOfCordinates <<std::endl;


#if 0
    IMat *matDst;
    RotateGrayImage(*src, &matDst, 0, rect);
    saveResult(*matDst, std::string("mat.txt").c_str());

    CreateGeoMatchModel(&tpl, matDst, 30, 170);
    std::cout << tpl.noOfCordinates <<std::endl;
#endif

    IplImage *dst = cvCreateImage(cvSize(srcImage->width * 2, srcImage->height * 2), IPL_DEPTH_8U, 1);
    memset(dst->imageData, 0, dst->imageSize);

//    ResumeTpl(&tpl);
//    rotateTpl(&tpl, 142, 66, 3);
//    DrawContours(dst, tpls[3][10], CvScalar(255,255,255), 1);
//
//    cvNamedWindow("tpl", 1);
//    cvShowImage("hehe", dst);
//    cvWaitKey(-1);
    IPoint point;
    double degree;
    image= cvLoadImage(argv[2], -1 );
    IplImage *dstImage = cvCreateImage(cvSize(image->width, image->height), IPL_DEPTH_8U, 1);
    if (image->nChannels == 3)
        cvCvtColor(image, dstImage, CV_RGB2GRAY);
    IMat *searchIMat = ICreateMat(dstImage->height, dstImage->width, U8C1);
    Ipl2IMat(dstImage, searchIMat);

    IplImage *Sg = cvLoadImage("images/testImages/Sgray.bmp", -1);
    IplImage *S4 = cvLoadImage("images/testImages/Slevel4.bmp", -1);
    IplImage *S3 = cvLoadImage("images/testImages/Slevel3.bmp", -1);
    IplImage *S2 = cvLoadImage("images/testImages/Slevel2.bmp", -1);
    IplImage *S1 = cvLoadImage("images/testImages/Slevel1.bmp", -1);


    cvSmooth(Sg,Sg, CV_GAUSSIAN,3,3,2);
    cvSmooth(S4,S4, CV_GAUSSIAN,3,3,2);
    cvSmooth(S3,S3, CV_GAUSSIAN,3,3,2);
    cvSmooth(S2,S2, CV_GAUSSIAN,3,3,2);
    cvSmooth(S1,S1, CV_GAUSSIAN,3,3,2);

    IMat * SgM = ICreateMat(Sg->height, Sg->width, U8C1);
    Ipl2IMat(Sg, SgM);
    saveResult(*SgM, "result/Sg.txt");

    IMat * S4M = ICreateMat(S4->height, S4->width, U8C1);
    Ipl2IMat(S4, S4M);
    saveResult(*S4M, "result/S4.txt");

#if 0
    IMat * S3M = ICreateMat(S3->height, S3->width, U8C1);
    Ipl2IMat(S3, S3M);
    saveResult(*S3M, "result/S3.txt");

    IMat * S2M = ICreateMat(S2->height, S2->width, U8C1);
    Ipl2IMat(S2, S2M);
    saveResult(*S2M, "S2.txt");
#endif

    IMat *SMg = ICreateMat(Sg->height, Sg->width, U8C1);
    Ipl2IMat(Sg, SMg);
//    std::cout << "Sg: " << Sg->width << " " << Sg->height << std::endl;

    IMat *SM4 = ICreateMat(S4->height, S4->width, U8C1);
    Ipl2IMat(S4, SM4);
//    std::cout << "S4: " << S4->width << " " << S4->height << std::endl;

    IMat *SM3 = ICreateMat(S3->height, S3->width, U8C1);
    Ipl2IMat(S3, SM3);
//    std::cout << "S3: " << S3->width << " " << S3->height << std::endl;

    IMat *SM2 = ICreateMat(S2->height, S2->width, U8C1);
    Ipl2IMat(S2, SM2);
//    std::cout << "S2: " << S2->width << " " << S2->height << std::endl;

    IMat *SM1 = ICreateMat(S1->height, S1->width, U8C1);
    Ipl2IMat(S1, SM1);
//    std::cout << "S1: " << S1->width << " " << S1->height << std::endl;
    IMat *mats[5] = {SM1, SM2, SM3, SM4, SMg};


    IPoint here;
    int dg;
    FindTemplateInPyramid(tpls, mats, &here, &dg);

    DrawContours(Sg, here,tpls[4][dg], CvScalar(255,255,255), 1);
//    cvNamedWindow("tpl", 1);
    cvShowImage("hehe", Sg);
    cvWaitKey(-1);


//    test();
    return 0;

    TimeTracker tt;
    tt.start();

    IPoint where;
    where.x = 44;
    where.y = 60;
//    resumeIMat(std::string("SearchArray_vert.txt").c_str(), searchIMat);
//    double score = FindGeoMatchModelRotateTplInRange(tpls[3],where, 0, SM4, 0.95, 0.999995, &point, 30, &degree);
    double score = FindGeoMatchModelRotateTpl(tpls[4],SMg, 0.95, 0.999995, &point, 30, &degree);
    tt.stop();
    std::cout << " found at degree: " << degree << " Points: (" << point.x << "," << point.y << ")" << " score: " << score
              << " Time used: " << tt.duration() << std::endl;
#if 0
    TimeTracker tt;
    tt.start();

//    resumeIMat(std::string("SearchArray_vert.txt").c_str(), searchIMat);
    double score = FindGeoMatchModelRotateTpl(tpls,searchIMat, 0.95, 0.999995, &point, 30, &degree);
    tt.stop();
    std::cout << " found at degree: " << degree << " Points: (" << point.x << "," << point.y << ")" << " score: " << score
              << " Time used: " << tt.duration() << std::endl;
#endif
    return 0;
}

typedef struct TplItemType {
    IPoint point;
    double mag;
    double edgeX;
    double edgeY;
} TplItem;

int compare (const void *elem1, const void *elem2 )
{
    int comparey = ((TplItem *)elem1)->point.y - ((TplItem *)elem2)->point.y;
    int comparex = ((TplItem *)elem1)->point.x - ((TplItem *)elem2)->point.x;
    if (comparey < 0) {
        return -1;
    } else if (comparey == 0 && comparex < 0) {
        return -1;
    } else if (comparey == 0 && comparex > 0) {
        return 1;
    } else if (comparex == 0 && comparey == 0) {
        return 0;
    } else if (comparey > 0) {
        return 1;
    }
}


void init_tplItems(TplItem *items, TemplateMatch * tpl) {
    for (int i = 0; i < tpl->noOfCordinates; ++i) {
        items[i].point = tpl->cordinates[i];
        items[i].mag = tpl->edgeMagnitude[i];
        items[i].edgeX= tpl->edgeDerivativeX[i];
        items[i].edgeY= tpl->edgeDerivativeY[i];
    }
}


void resume_tplItems(TplItem *items, TemplateMatch * tpl) {
    for (int i = 0; i < tpl->noOfCordinates; ++i) {
        tpl->cordinates[i] = items[i].point;
        tpl->edgeMagnitude[i] = items[i].mag;
        tpl->edgeDerivativeX[i] = items[i].edgeX ;
        tpl->edgeDerivativeY[i] = items[i].edgeY;
    }
}

void ResumeTpl(TemplateMatch *tpl)
{
    int i;
    int x, y;
    double d;
    FILE *fp1 = fopen("cordinates.txt", "r");
    FILE *fp2 = fopen("edgeMag.txt", "r");
    FILE *fp3 = fopen("ex.txt", "r");
    FILE *fp4 = fopen("ey.txt", "r");
    tpl->noOfCordinates = 433;
    tpl->modelWidth = 140;
    tpl->modelHeight= 42;
    tpl->centerOfGravity.x = 20;
    tpl->centerOfGravity.y = 71;

    for(i = 0; i < tpl->noOfCordinates; i++) {
        fscanf(fp1, "%d", &x);
        fscanf(fp1, "%d", &y);
        tpl->cordinates[i].x = x;
        tpl->cordinates[i].y = y;
        fscanf(fp2, "%lf", &d);
        tpl->edgeMagnitude[i] = d;
        fscanf(fp3, "%lf", &d);
        tpl->edgeDerivativeX[i] = d;
        fscanf(fp4, "%lf", &d);
        tpl->edgeDerivativeY[i] = d;
    }

}
void saveTpl(TemplateMatch *tpl) {
    int i;
    int x, y;
    double d;
    FILE *fp1 = fopen("cordinates.txt", "w");
    FILE *fp2 = fopen("edgeMag.txt", "w");
    FILE *fp3 = fopen("ex.txt", "w");
    FILE *fp4 = fopen("ey.txt", "w");

    for(i = 0; i < tpl->noOfCordinates; i++) {
        x = tpl->cordinates[i].x;
        y = tpl->cordinates[i].y;
        fprintf(fp1, "%d ", x);
        fprintf(fp1, "%d ", y);

        d = tpl->edgeMagnitude[i];
        fprintf(fp2, "%lf ", d);
        d = tpl->edgeDerivativeX[i];
        fprintf(fp3, "%lf ", d);
        d = tpl->edgeDerivativeY[i];
        fprintf(fp4, "%lf ", d);

    }
    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    fclose(fp4);

}

/*
 * 只对模板进行旋转的话，由于不做插值过程，可能导致模板变形，对应不上。
 * */
void rotateTpl(TemplateMatch *tpl, int iWidth, int iHeight, double degree)
{
    saveTpl(tpl);
    TplItem *TplItems = (TplItem *)malloc(sizeof(TplItem) * tpl->noOfCordinates);
    memset(TplItems, 0, sizeof(TplItem) * tpl->noOfCordinates);


    init_tplItems(TplItems, tpl);

    double angle = degree * CV_PI / 180.;
    double alpha = cos(angle);
    double beta = sin(angle);
    double m[6];
    IPoint *points = (IPoint *)malloc(sizeof(IPoint) * tpl->noOfCordinates);
    memcpy(points, tpl->cordinates, sizeof(IPoint) * tpl->noOfCordinates);

    int iNewWidth = cvRound(iWidth * fabs(alpha) + iHeight * fabs(beta));
    int iNewHeight = cvRound(iHeight * fabs(alpha) + iWidth * fabs(beta));
    iNewHeight = (iNewHeight % 2 == 1) ? iNewHeight + 1 : iNewHeight;
    iNewWidth = (iNewWidth % 2 == 1) ? iNewWidth + 1 : iNewWidth;

    m[0] = alpha;
    m[1] = beta;
    m[2] = (1 - alpha) * (iWidth / 2.) - beta * (iHeight / 2.);
    m[3] = -m[1];
    m[4] = m[0];
    m[5] = beta * (iWidth / 2.) + (1 - alpha) * (iHeight / 2.);
    m[2] += (iNewWidth - iWidth) / 2.0f;	m[5] += (iNewHeight - iHeight) / 2.0f;

    for (int i = 0; i < tpl->noOfCordinates; ++i) {
        double fx = m[0] * tpl->cordinates[i].x + m[1] * tpl->cordinates[i].y + m[2];
        double fy = m[3] * tpl->cordinates[i].x + m[4] * tpl->cordinates[i].y + m[5];
        TplItems[i].point.x = (int)fx;
        TplItems[i].point.y = (int)fy;
    }
    /*不需要对模板进行排序，因为排不排序都是一样的，只需要保证坐标点和其对应的梯度信息是对应的即可，
     * 对模板的匹配不需要完全按照顺序进行，因为在匹配过程中是要用到坐标点在搜索图像中进行匹配。*/
//    qsort(TplItems, tpl->noOfCordinates, sizeof(TplItem), compare);
    resume_tplItems(TplItems, tpl);
    tpl->modelHeight = iNewHeight;
    tpl->modelWidth = iNewWidth;

    int fx = m[0] * tpl->centerOfGravity.x + m[1] * tpl->centerOfGravity.y + m[2];
    int fy = m[3] * tpl->centerOfGravity.x + m[4] * tpl->centerOfGravity.y + m[5];
    tpl->centerOfGravity.x = fx;
    tpl->centerOfGravity.y = fy;
    printf("x: %d, y:%d height: %d, width: %d\n", fx, fy, iNewHeight, iNewWidth);
    saveTpl(tpl);
}





void deflatRect(IPoint *rect, int w, int h)
{
    IPoint center;
    center.x = (rect[0].x + rect[3].x) / 2;
    center.y = (rect[0].y + rect[1].y) / 2;

    IPoint vc, hc;
    vc.x = (rect[0].x + rect[1].x) / 2;
    vc.y = (rect[0].y + rect[1].y) / 2;

    hc.x = (rect[0].x + rect[3].x) / 2;
    hc.y = (rect[0].y + rect[3].y) / 2;

    /*求法向量*/
    IPoint normalVectorH, normalVectorV;

    normalVectorV.x = center.x - vc.x;
    normalVectorV.y = center.y - vc.y;
    double length;
    length = sqrt(normalVectorV.x * normalVectorV.x + normalVectorV.y * normalVectorV.y);
    normalVectorV.x = 1.0 * normalVectorV.x / length;
    normalVectorV.y = 1.0 * normalVectorV.y / length;

    normalVectorH.x = center.x - hc.x;
    normalVectorH.y = center.y - hc.y;

    length = sqrt(normalVectorH.x * normalVectorH.x + normalVectorH.y * normalVectorH.y);
    normalVectorH.x = 1.0 * normalVectorH.x / length;
    normalVectorH.y = 1.0 * normalVectorH.y / length;

    rect[0].x = rect[0].x + w * (normalVectorH.x + normalVectorV.x);
    rect[0].y = rect[0].y + w * (normalVectorH.y + normalVectorV.y);

    rect[1].x = rect[1].x + w * (-normalVectorH.x + normalVectorV.x);
    rect[1].y = rect[1].y + w * (-normalVectorH.y + normalVectorV.y);

    rect[2].x = rect[2].x + h * (-normalVectorH.x - normalVectorV.x);
    rect[2].y = rect[2].y + h * (-normalVectorH.y - normalVectorV.y);

    rect[3].x = rect[3].x + h * (normalVectorH.x - normalVectorV.x);
    rect[3].y = rect[3].y + h * (normalVectorH.y - normalVectorV.y);
}

void test()
{
    TemplateMatch tpl;

    tpl.cordinates = (IPoint *)malloc(500 * sizeof(IPoint));
    tpl.edgeMagnitude = (double *)malloc(500 * sizeof(double));		//Allocate memory for edge magnitude for selected points
    tpl.edgeDerivativeX = 	(double *)malloc(500 * sizeof(double));		//Allocate memory for edge X derivative for selected points
    tpl.edgeDerivativeY = (double *)malloc(500 * sizeof(double));		////Allocate memory for edge Y derivative for selected points

    ResumeTpl(&tpl);

    std::cout << tpl.noOfCordinates <<std::endl;
    std::cout << tpl.modelWidth<<std::endl;
    std::cout << tpl.modelHeight<<std::endl;
    std::cout << tpl.centerOfGravity.x<<std::endl;
    std::cout << tpl.centerOfGravity.y<<std::endl;



    double partialSum;
    double partialScore;
    int curX, curY;
    int i,j,m;
    const short *_Sdx;
    const short *_Sdy;
    double iTx, iTy, iSx, iSy;
    IMat *Sdx = 0, *Sdy = 0;
    int sumOfCoords = 0;

    ISize Ssize;
    Ssize.width = 320;
    Ssize.height = 240;

   double **edgeX;  //Gradient magnitude matrix
    double **edgeY;  //Gradient magnitude matrix
    CreateDoubleMatrix(&edgeX, Ssize); // create image to save gradient magnitude  values
    CreateDoubleMatrix(&edgeY, Ssize); // create image to save gradient magnitude  values

    Sdx = fICreateMat(Ssize.height, Ssize.width, S16C1); // X derivatives
    Sdy = fICreateMat(Ssize.height, Ssize.width, S16C1); // y derivatives

    double minScore = 0.95;
    double greediness = 0.999995;
    double resultScore = 0;
    IPoint resultPoint;


    resumeIMatf("Sdx.txt", Sdx);
    resumeIMatf("Sdy.txt", Sdy);

    for (i = 0; i < Ssize.height; i++) {
        _Sdx = (short *) (Sdx->fptr + Sdx->step * (i));
        _Sdy = (short *) (Sdy->fptr + Sdy->step * (i));

        for (j = 0; j < Ssize.width; j++) {

            iSx = _Sdx[j];  // X derivative of Source image
            iSy = _Sdy[j];  // Y derivative of Source image
            double vector_length = sqrt(_Sdx[j] * _Sdx[j] + _Sdy[j] * _Sdy[j]);
            edgeX[i][j] = _Sdx[j] / vector_length;
            edgeY[i][j] = _Sdy[j] / vector_length;
        }
    }



    double normMinScore = minScore / tpl.noOfCordinates; // precompute minumum score
    double normGreediness =
            ((1 - greediness * minScore) / (1 - greediness)) / tpl.noOfCordinates; // precompute greedniness

    for (i = 0; i < 240; i++) {
        for (j = 0; j < 320; j++) {
            partialSum = 0; // initilize partialSum measure
            int count = 0;
            for (m = 0; m < tpl.noOfCordinates; m++) {
                count++;
                curX = i + tpl.cordinates[m].x;    // template X coordinate
                curY = j + tpl.cordinates[m].y; // template Y coordinate
                iTx = tpl.edgeDerivativeX[m];    // template X derivative
                iTy = tpl.edgeDerivativeY[m];    // template Y derivative

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
                resultPoint.x = i;            // result coordinate X
                resultPoint.y = j;            // result coordinate Y
            }
        }
    }
    printf("score: %lf, where: (%d, %d)", resultScore, resultPoint.x, resultPoint.y);


}

void ResumeTpls(TemplateMatch *tpls, int numOfTpls, char **files)
{
    int i,j;
    int x, y;
    double d;
    FILE *fp1 = fopen(files[0], "r");
    FILE *fp2 = fopen(files[1], "r");
    FILE *fp3 = fopen(files[2], "r");
    FILE *fp4 = fopen(files[3], "r");
    FILE *fp5 = fopen(files[4], "r");

    for (i = 0; i < numOfTpls; ++i) {
        fscanf(fp5, "%d", &(tpls[i].noOfCordinates));
        fscanf(fp5, "%d", &(tpls[i].modelWidth));
        fscanf(fp5, "%d", &(tpls[i].modelHeight));
        fscanf(fp5, "%d", &(tpls[i].centerOfGravity.x));
        fscanf(fp5, "%d", &(tpls[i].centerOfGravity.y));

        tpls[i].cordinates      =(IPoint*)malloc(sizeof(IPoint) * tpls[i].noOfCordinates);
        tpls[i].edgeMagnitude   =(double *)malloc(sizeof(double) * tpls[i].noOfCordinates);
        tpls[i].edgeDerivativeX =(double *)malloc(sizeof(double) * tpls[i].noOfCordinates);
        tpls[i].edgeDerivativeY =(double *)malloc(sizeof(double) * tpls[i].noOfCordinates);
        for(j = 0; j < tpls[i].noOfCordinates; j++) {
            fscanf(fp1, "%d", &x);
            fscanf(fp1, "%d", &y);
            tpls[i].cordinates[j].x = x;
            tpls[i].cordinates[j].y = y;
            fscanf(fp2, "%lf", &d);
            tpls[i].edgeMagnitude[j] = d;
            fscanf(fp3, "%lf", &d);
            tpls[i].edgeDerivativeX[j] = d;
            fscanf(fp4, "%lf", &d);
            tpls[i].edgeDerivativeY[j] = d;
        }

    }

    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    fclose(fp4);
    fclose(fp5);
}
void saveTpls(TemplateMatch *tpls, int numOfTpls, char **files) {
    int i, j;
    int x, y;
    double d;

    FILE *fp1 = fopen(files[0], "w");
    FILE *fp2 = fopen(files[1], "w");
    FILE *fp3 = fopen(files[2], "w");
    FILE *fp4 = fopen(files[3], "w");
    FILE *fp5 = fopen(files[4], "w");

    for (i = 0; i < numOfTpls; ++i) {
        fprintf(fp5, "%d ", tpls[i].noOfCordinates);
        fprintf(fp5, "%d ", tpls[i].modelWidth);
        fprintf(fp5, "%d ", tpls[i].modelHeight);
        fprintf(fp5, "%d ", tpls[i].centerOfGravity.x);
        fprintf(fp5, "%d ", tpls[i].centerOfGravity.y);

        for(j = 0; j < tpls[i].noOfCordinates; j++) {
            x = tpls[i].cordinates[j].x;
            y = tpls[i].cordinates[j].y;
            fprintf(fp1, "%d ", x);
            fprintf(fp1, "%d ", y);

            d = tpls[i].edgeMagnitude[j];
            fprintf(fp2, "%lf ", d);
            d = tpls[i].edgeDerivativeX[j];
            fprintf(fp3, "%lf ", d);
            d = tpls[i].edgeDerivativeY[j];
            fprintf(fp4, "%lf ", d);
        }
    }

    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    fclose(fp4);
    fclose(fp5);

}

