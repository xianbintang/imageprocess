//
// Created by xianb on 2017/4/19.
//
#include<opencv2/core/core.hpp>
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <iostream>
#include <fstream>
#include "../include/common.h"

static void saveMat(int **img, std::string filename, int h, int w)
{
	auto fp = fopen(filename.c_str(), "w");
	for(int i = 0; i < h; i++) {
		for(int j = 0; j < w; j++) {
//			if(img[i][j] == 0)
//				fprintf(fp, "%d ", 0);
//			else
				fprintf(fp, "%d ", img[i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

void SaveImageAsText(IplImage *img, std::string filename)
{
	int **mat = new int*[img->height];
	for(int i = 0; i < img->height; i++)
		mat[i] = new int[img->width];
	for(int i = 0; i < img->height; i++) {
		unsigned char * srow = (unsigned char *)(img->imageData + i * img->widthStep);
		for(int j = 0; j < img->width; j++) {
			mat[i][j] = srow[j];
		}
	}
	saveMat(mat,filename, img->height, img->width);
}

void Text2Ipl(IplImage *img, std::string filename)
{
    std::ifstream fin(filename);
	int num;
	for(int i = 0; i < img->height; i++) {
		unsigned char * srow = (unsigned char *)(img->imageData + i * img->widthStep);
		for(int j = 0; j < img->width; j++) {
            fin >> num;
			srow[j] = num;
		}
	}
}


void displayTextImage(std::vector<std::vector<int>> TextImg, int width, int height)
{
	IplImage *img = cvCreateImage(cvSize(width, height), IPL_DEPTH_16S, 1);
    std::ofstream out("outNum.txt");
	for (const auto &v:TextImg) {
        out << "Union: ";
		for (const auto &p: v) {
//			std::cout << p%width << " " << p/width << std::endl;
//			cvCircle(img, cv::Point(p%width, p/width),1 ,cvScalar(255, 255, 255), -1);
            unsigned char * prow = (unsigned char*)img->imageData + (p/width) * img->widthStep;
//            out << "(" << p/width << ", " << p%width << ")" << " ";
            out << p << " ";
			prow[p%width] = 255;
		}
        out << std::endl;
        cvShowImage("img", img);
		cvWaitKey();
	}
}

