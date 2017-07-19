//
// Created by xianb on 2017/5/17.
//
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include "types.h"
#include <fstream>
#include <iostream>

void Ipl2IMat(IplImage *img, IMat * src) {
    for (int i = 0; i < img->height; ++i) {
        unsigned char * srow = (unsigned char *)(img->imageData + i * img->widthStep);
        for (int j = 0; j < img->width; ++j) {
            (src->ptr + i * src->step)[j] = srow[j];
        }
    }
}

void IMat2Ipl(const IMat * src, IplImage *img) {
    for (int i = 0; i < img->height; ++i) {
        unsigned char * srow = (unsigned char *)(img->imageData + i * img->widthStep);
        for (int j = 0; j < img->width; ++j) {
             srow[j] = (src->ptr + i * src->step)[j];
        }
    }
}

void fIpl2IMat(IplImage *img, IMat * src) {
    for (int i = 0; i < img->height; ++i) {
        const short* srow = (const short*)(img->imageData + i * img->widthStep);
        for (int j = 0; j < img->width; ++j) {
            (src->fptr + i * src->step)[j] = srow[j];
        }
    }
}

void fIMat2Ipl(const IMat * src, IplImage *img) {
    for (int i = 0; i < img->height; ++i) {
        short* srow = (short*)(img->imageData + i * img->widthStep);
        for (int j = 0; j < img->width; ++j) {
             srow[j] = (src->fptr + i * src->step)[j];
        }
    }
}

void saveResultf(IMat mat, const char *path) {
    std::ofstream fout(path);
    for (int i = 0; i < mat.height; ++i) {
        for (int j = 0; j < mat.width; ++j) {
            fout << (mat.fptr + i * mat.step)[j] << " ";
        }
        fout << std::endl;
    }
}
void saveResult(IMat mat, const char *path) {
    std::ofstream fout(path);
    for (int i = 0; i < mat.height; ++i) {
        for (int j = 0; j < mat.width; ++j) {
            fout << (mat.ptr + i * mat.step)[j] << " ";
        }
        fout << std::endl;
    }
}

void HI_Sobel(const IMat *src, IMat *gx, IMat *gy)
{
    IplImage *  img = cvCreateImage(cvSize(src->width, src->height), 8, 1);
    IMat2Ipl(src, img);
//    cvShowImage("gx", img);
//    cvWaitKey(-1);
    IplImage *sdx = cvCreateImage(cvSize(src->width, src->height), IPL_DEPTH_16S, 1);
    IplImage *sdy = cvCreateImage(cvSize(src->width, src->height), IPL_DEPTH_16S, 1);
    cvSobel(img, sdx, 1, 0, 3);
//    cvShowImage("gx", sdx);
//    cvWaitKey(-1);
    cvSobel(img, sdy, 0, 1, 3);
//    cvShowImage("gx", sdy);
//    cvWaitKey(-1);
    fIpl2IMat(sdx, gx);
    fIpl2IMat(sdy, gy);
//    saveResult(*gx, "gx.txt");
//    saveResult(*gy, "gy.txt");

//    fIMat2Ipl(gx, img);
//    cvShowImage("gx", img);
//    cvWaitKey(-1);
//    fIMat2Ipl(gy, img);
//    cvShowImage("gy", img);
//    cvWaitKey(-1);
}
void HI_Canny(const IMat *src, IMat *edge)
{
    IplImage *  img = cvCreateImage(cvSize(src->width, src->height), 8, 1);
    IMat2Ipl(src, img);

    IplImage *edgeIpl = cvCreateImage(cvSize(src->width, src->height), IPL_DEPTH_8U, 1);
    cvCanny(img, edgeIpl, 30, 170, 3);
//    std::cout <<"hehe" <<std::endl;

    Ipl2IMat(edgeIpl, edge);
}
