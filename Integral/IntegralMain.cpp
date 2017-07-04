//
// Created by xianb on 2017/6/29.
//


#include <opencv/cv.h>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include "TimeTracker.h"
#include <opencv/cv.hpp>
#include <iostream>

using namespace cv;
#define BLOCKW 4
#define BLOCKH 4
#define DIFFERENCEOFMEAN 10
#define LNE 0

inline int getCurrentMean(const Mat &integral, const Point &posA, const Point &posD, int size)
{
    int A = (posA.y - 1 > 0 && posA.x - 1 > 0) ? integral.at<int>(posA.y - 1, posA.x - 1) : 0;
    int B = (posA.y - 1 > 0) ? integral.at<int>(posA.y - 1, posD.x) : 0;
    int C = (posA.x - 1 > 0) ? integral.at<int>(posD.y, posA.x - 1) : 0;
    int D = integral.at<int>(posD.y, posD.x);
    int mean = (A + D - B - C) / size;
    return mean;
}

int CreateIntegralBinaryPattern(const Mat &mat, Mat &binaryPattern, const Point &start, const Point &end)
{
    int mean = getCurrentMean(mat, start, end, (end.x - start.x + 1) * (end.y - start.y + 1));
//    std::cout << "mean: " << mean <<std::endl;
    for(int i=0;i<binaryPattern.rows;i++)
    {
        for(int j=0;j<binaryPattern.cols;j++)
        {
            /* Point 是(x, y)而不是(y, x)*/
            binaryPattern.at<int>(i, j) = getCurrentMean(mat, Point(start.x + j * BLOCKH, start.y + i * BLOCKW),\
             Point(start.x + (j + 1) * BLOCKH - 1, start.y + (i + 1) * BLOCKW - 1), BLOCKH * BLOCKW) > mean ? 1 : -1;
        }
    }
    return mean;
}

int notEqual(const Mat &sourceIntegral, const Mat &tplBP, const Point &start, Point end, int tplMean)
{
    int count = 0;
    int mean = 0;
    if (end.y > sourceIntegral.rows)
        end.y = sourceIntegral.rows;
    if (end.x > sourceIntegral.cols)
        end.x = sourceIntegral.cols;
    mean = getCurrentMean(sourceIntegral, start, end, (end.x - start.x + 1) * (end.y - start.y + 1));
    if (abs(mean - tplMean) >= DIFFERENCEOFMEAN) {
        return -1;
    }
    static int ct = 0;
    ct++;
//    printf("ct: %d\n", ct);
//    std::cout << "mean: " << mean <<std::endl;
    for(int i=0;i<tplBP.rows;i++)
    {
        for(int j=0;j<tplBP.cols;j++)
        {
            /* Point 是(x, y)而不是(y, x)*/
            int curMean = getCurrentMean(sourceIntegral, Point(start.x + j * BLOCKH, start.y + i * BLOCKW),\
             Point(start.x + (j + 1) * BLOCKH - 1, start.y + (i + 1) * BLOCKW - 1), BLOCKH * BLOCKW) > mean ? 1 : -1;
            if (curMean != tplBP.at<int>(i, j)) {
                count++;
            }
            if (count > LNE)
                return -1;
        }
    }
    return count;
}
void matchIt(const Mat &tplBP, const Mat &sourceIntegral, int tplMean, std::vector<Point> &where)
{
    Mat curBP(tplBP.rows, tplBP.cols, CV_32SC(1));
    Mat dst;
    for (int i = 0; i < sourceIntegral.rows - BLOCKH; ++i) {
        for (int j = 0; j < sourceIntegral.cols - BLOCKW; ++j) {
            /* get current binary test pattern */
            int ne = notEqual(sourceIntegral, tplBP, Point(j, i), Point(j + tplBP.cols * BLOCKW - 1, i + tplBP.rows * BLOCKH - 1), tplMean);
            if (ne != -1) {
                where.push_back(Point(j, i));
            }
//            if (ne != -1)
//                std::cout << "x: " << j << " y: " << i << "sum: " << sum << std::endl;
//            std::cout << sum << std::endl;
        }
    }
}

int main(int argc,char *argv[])
{
    Mat tpl = imread(argv[1]);
    Mat source =imread(argv[2]);

    cvtColor(source,source,CV_RGB2GRAY);
    cvtColor(tpl,tpl,CV_RGB2GRAY);

    Mat sourceIntegral;
    Mat tplIntegral;
    integral(source, sourceIntegral,CV_32S); //计算积分图
    integral(tpl, tplIntegral,CV_32S); //计算积分图

    Mat tplBP;
    tplBP.create(tplIntegral.rows / BLOCKH, tplIntegral.cols / BLOCKW, CV_32SC(1));
    int tplMean = CreateIntegralBinaryPattern(tplIntegral, tplBP, Point(0, 0), Point(tplIntegral.rows - 1, tplIntegral.cols - 1));

    TimeTracker tt;
    std::vector<Point> where;
    tt.start();
    matchIt(tplBP, sourceIntegral, tplMean, where);
    tt.stop();
    std::cout << where << "time used: " << tt.duration() << "ms" << std::endl;
//    std::cout << "M = " << std::endl << " " << tplBP << std::endl << std::endl;

    return 0;
}
