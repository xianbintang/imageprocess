//
// Created by xianb on 2017/5/24.
//
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <opencv/highgui.h>
#include <opencv/cv.hpp>
#include <bitset>
#include <TimeTracker.h>

void displayTextImagebitMat(const char *path, int width, int height)
{
    IplImage *img = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, 1);
    memset(img->imageData, 0, sizeof(char) * width * height);
    std::ifstream fin(path);
    int num = 0;

    int ct = 0;
    unsigned char* prow = (unsigned char*)(img->imageData);
    while(fin >> num) {
        ct++;
        if(ct >= width * height / 8) {
            std::cout << "out of board" << std::endl;
            break;
        }
//        std::cout << num << std::endl;
        for (int i = 0; i < 8; ++i) {
            auto bit = std::bitset<10>(num);
            *prow = bit[7] * 255;
            ++prow;
            num <<= 1;
        }
    }
    std::cout << "ct: " << ct << std::endl;

    cvShowImage("img", img);
    cv::waitKey(0);
//    cvSaveImage("E://tmp//tmp139.bmp", img);
}
void displayTextImage(const char *path, int width, int height)
{
    IplImage *img = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, 1);
    memset(img->imageData, 0, sizeof(char) * width * height);
    std::ifstream fin(path);
    int num;
    for (int i = 0; i < height; ++i) {
//        unsigned char * prow = (unsigned char *)(img->imageData + i * img->widthStep);
        unsigned char* prow = (unsigned char*)(img->imageData + i * img->widthStep);
        for (int j = 0; j < width; ++j) {
            fin >> num;
//            std::cout << num << std::endl;
            prow[j] = num * 255;
        }
    }
    cvShowImage("img", img);
//    cvSaveImage("E://tmp//tmp.bmp", img);
    cvWaitKey();
}


void displayTextImagef(const char *path, int width, int height)
{
    IplImage *img = cvCreateImage(cvSize(width, height), IPL_DEPTH_16S, 1);
    memset(img->imageData, 0, sizeof(char) * width * height);
    std::ifstream fin(path);
    int num;
    for (int i = 0; i < height; ++i) {
//        unsigned char * prow = (unsigned char *)(img->imageData + i * img->widthStep);
        short* prow = (short *)(img->imageData + i * img->widthStep);
        for (int j = 0; j < width; ++j) {
            fin >> num;
//            std::cout << num << std::endl;
            prow[j] = num *50;
            if(i == 36 && j == 98) {
                std::cout << num << std::endl;
            }
        }
    }
    cvShowImage("img", img);
    cvSaveImage("SearchArray_8U.bmp", img);
    cvWaitKey();
}


void displayTextImageS8(const char *path, int width, int height)
{
    IplImage *img = cvCreateImage(cvSize(width, height), IPL_DEPTH_8S, 1);
    memset(img->imageData, 0, sizeof(char) * width * height);
    std::ifstream fin(path);
    int num;
    for (int i = 0; i < height; ++i) {
//        unsigned char * prow = (unsigned char *)(img->imageData + i * img->widthStep);
        signed char* prow = (signed char*)(img->imageData + i * img->widthStep);
        for (int j = 0; j < width; ++j) {
            fin >> num;
//            std::cout << num << std::endl;
            prow[j] = num;
        }
    }
    cvShowImage("img", img);
//    cvSaveImage("SearchArray_8U.bmp", img);
    cvWaitKey();
}



void Text2HIImage(char *path, int width, int height) {

    int i, j;
    int num;
    printf("\n");
    FILE *fp = fopen(path, "r");
    std::ifstream fin(path);
    for (i = 0; i < height; ++i) {
        for (j = 0; j < width; ++j) {
//            fscanf(fp, "%d", &num);
//            printf("%d ", num);
            fin >> num;
            if(num < -300)
                std::cout << num << " ";
        }
    }
}

int main(int argc, char ** argv)
{
//    Text2HIImage("templateArray.txt", 144, 76);
    if (atoi(argv[4]) == 0)
        displayTextImage(argv[1], atoi(argv[2]), atoi(argv[3]));
    else if(atoi(argv[4]) == 1)
        displayTextImagef(argv[1], atoi(argv[2]), atoi(argv[3]));
    else if(atoi(argv[4]) == 2)
        displayTextImageS8(argv[1], atoi(argv[2]), atoi(argv[3]));
    else if(atoi(argv[4]) == 3)
        displayTextImagebitMat(argv[1], atoi(argv[2]), atoi(argv[3]));
    else if(atoi(argv[4]) == 4) {
        Text2HIImage(argv[1], atoi(argv[2]), atoi(argv[3]));
    }
        return 0;
}