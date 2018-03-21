//
// Created by xianb on 2018/3/19.
//
#include <opencv2/core/mat.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv/cv.hpp>
#include <windef.h>
#include <fstream>
#include <c++/iostream>
#include "TimeTracker.h"

const int WIDTH = 640;
const int HEIGHT = 480;
using namespace std;
using namespace cv;

const int BLOCKH = 8;
const int BLOCKW = 8;
const int NBLOCK = 30;
void saveMat(cv::Mat mat, const char *path);
void saveMatf(cv::Mat mat, const char *path);

cv::Mat colorImg;
void saveMat(cv::Mat mat, const char *path) {
    FILE *fp = fopen(path, "w");
    int i,j;
    for (i = 0; i < mat.rows; ++i) {
        for (j = 0; j < mat.cols; ++j) {
//            fprintf(fp, "%d ", (mat.ptr + i * mat.step)[j]);
            fprintf(fp, "%d ", mat.at<uchar>(i, j));
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}
void saveMatf(cv::Mat mat, const char *path) {
    FILE *fp = fopen(path, "w");
    int i,j;
    for (i = 0; i < mat.rows; ++i) {
        for (j = 0; j < mat.cols; ++j) {
//            fprintf(fp, "%d ", (mat.ptr + i * mat.step)[j]);
            fprintf(fp, "%d ", mat.at<short>(i, j));
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}
static cv::Mat get_y_from_yuv(const UINT8 *yuv, const UINT16 width, const UINT16 height)
{
    cv::Mat gray;
    if (!yuv) {
        std::cout << "ERROR: yuv NULL pointer" << std::endl;
    } else {
        gray.create(height, width, CV_8UC1);
        memcpy(gray.data, yuv, width * height * sizeof(unsigned char)); // 只读取y分量上的数据
    }

    return gray;
}

inline int getCurrentSum(const Mat &integral, Point &posA, Point &posD, int size)
{
    int A = 0, B = 0, C = 0, D = 0;
    if(posD.x > integral.cols) {
        posD.x = integral.cols - 1;
    }

    if(posD.y > integral.rows) {
        posD.y = integral.rows - 1 ;
    }
    A = integral.at<int>(posA.y, posA.x);
    B = integral.at<int>(posA.y, posD.x);
    C = integral.at<int>(posD.y, posA.x);
    D = integral.at<int>(posD.y, posD.x);

//    std::cout << posA << " " << posD << " ";
    int sum = (A + D - B - C) / 255;

    return sum;
}

int compareIntegralSum(cv::Mat &integralSum1, cv::Mat &integralSum2)
{
    if(integralSum1.cols != integralSum2.cols || integralSum1.rows != integralSum2.rows) {
        return -1;
    }
    int height = integralSum1.rows;
    int width = integralSum1.cols;
    int ct = 0;
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            if(integralSum1.at<int>(i, j) == integralSum2.at<int>(i,j)) {
                ct++;
            }
        }
    }
    return ct;
}

Mat CreateIntegralSum(const Mat &mat, Point &start, Point &end, int patternSize)
{
//    int mean = getCurrentSum(mat, start, end, (end.x - start.x + 1) * (end.y - start.y + 1));
    Mat sumPattern;
    sumPattern.create(mat.rows / patternSize, mat.cols / patternSize, CV_32SC(1));

    for(int i=0;i<sumPattern.rows;i++)
    {
        for(int j=0;j<sumPattern.cols;j++)
        {
            /* Point 是(x, y)而不是(y, x)*/
            Point A = Point(start.x + j * patternSize, start.y + i * patternSize), D = Point(start.x + (j + 1) * patternSize, start.y + (i + 1) * patternSize);
            auto curSum = getCurrentSum(mat, A, D , 0);
            sumPattern.at<int>(i, j) = curSum;
//            if(curSum != 0)
//                std::cout << curSum << std::endl;
        }
    }
    return sumPattern;
}


cv::Mat computeIntegralSum(cv::Mat &image)
{
    cv::GaussianBlur(image, image, cv::Size(5,5),0);
    cv::Canny(image, image, 30, 150);
    cv::Mat img_integ, integralSum;
    cv::integral(image, img_integ,CV_32S); //计算积分图
    integralSum.create(img_integ.rows / BLOCKH, img_integ.cols / BLOCKW, CV_32SC(1));
    Point A = Point(0, 0), D = Point(img_integ.cols - 1, img_integ.rows - 1);
//    CreateIntegralSum(img_integ, integralSum, A, D, 3);

    return integralSum;
}

int computePointsInt(cv::Mat &image)
{
    int ct = 0;
    for (int i = 0; i < image.rows; ++i) {
        for (int j = 0; j < image.cols; ++j) {
            ct += image.at<int>(i, j);
        }
    }

    return ct;
}

double computePointsChar(cv::Mat &image)
{
    double ct = 0;
    for (int i = 0; i < image.rows; ++i) {
        for (int j = 0; j < image.cols; ++j) {
            ct += image.at<uchar>(i, j);
        }
    }

    return ct;
}

cv::Mat computeIntegral(cv::Mat &image)
{
    cv::GaussianBlur(image, image, cv::Size(5,5),0);
    cv::Canny(image, image, 30, 150);
    std::cout << "after canny: " << computePointsChar(image) / 255 << std::endl;;
    cv::Mat integ;
    cv::integral(image, integ,CV_32S); //计算积分图

    return integ;
}

bool compareSumPattern(cv::Mat &srcPattern, cv::Mat &targetPattern, cv::Point position)
{
    int ct = 0;
    int ctnozero = 0;
    for (int i = 0; i < targetPattern.rows; ++i) {
        for (int j = 0; j < targetPattern.cols; ++j) {
            if(targetPattern.at<int>(i, j)) {
                ctnozero++;
            }
            if(targetPattern.at<int>(i, j) && abs(srcPattern.at<int>(i + position.y, j + position.x) - targetPattern.at<int>(i, j)) < 16) {
                ct++;
            }
        }
    }
//    std::cout << "compareSumPattern: " << ct << " " << ctnozero << std::endl;
    return 1.0 * ct / ctnozero> 0.70;
}
int findTargetArea(cv::Mat &src, cv::Mat &target)
{

    TimeTracker tt;

    // 返回指向需要被发送的内存缓冲区的指针
    auto integ= computeIntegral(src);
    auto tobe= computeIntegral(target);
    Point A = Point(0, 0), D=  Point(target.cols - 1, target.rows - 1);
    auto tobeSum = getCurrentSum(target, A, D, 0);

    int patternSize = 8;
    // 计算模板的pattern
    auto templatePattern = CreateIntegralSum(tobe, A, D, patternSize);
    std::cout << " ct: " << computePointsInt(templatePattern) << std::endl;

    // 计算图片的pattern

    A = Point(0, 0);
    D=  Point(src.cols - 1, src.rows - 1);
    auto srcPattern = CreateIntegralSum(integ, A, D, patternSize);
    std::cout << " ct: " << computePointsInt(srcPattern) << std::endl;

    std::vector<Point> region;
    tt.start();
    int ct = 0;
    for (int i = 0; i < src.rows; i += BLOCKH) {
        for (int j = 0; j < src.cols; j += BLOCKW) {
            Point A = Point(j, i), D = Point(j + NBLOCK * BLOCKW, i + NBLOCK * BLOCKH);
            int sum = getCurrentSum(integ, A, D, 0);
            if(sum > 1800 ) {
                Point pos = Point(j / patternSize, i / patternSize);
                if(compareSumPattern(srcPattern, templatePattern, pos)) {
                    region.push_back(A);
                    ct++;
                }
            }
//            std::cout << sum << " ";
        }
    }
    tt.stop();
    std::cout << "duration: " << tt.duration() << std::endl;
    for (int i = 0; i < region.size(); ++i) {
//        cv::circle(colorImg, region[i], 1, cv::Scalar(0,0,255));
        cv::circle(colorImg, cv::Point(region[i].x + NBLOCK * BLOCKW / 2, region[i].y + NBLOCK * BLOCKH / 2), 1, cv::Scalar(0,0,255));
//        cv::rectangle(colorImg, region[i],  cv::Point(region[i].x + NBLOCK * BLOCKW, region[i].y + NBLOCK * BLOCKH),  cv::Scalar(0,0,255));
    }
    std::cout << ct << std::endl;
    cv::imshow("integimg", integ);
    cv::imshow("tobe", tobe);
    cv::imshow("colorImg", colorImg);

    cv::waitKey(0);
}

int main(int argc, char **argv)
{

    std::string filename(argv[1]);
    std::string tobeMatch(argv[1]);

    cv::Mat template_image;
    cv::Mat tobeMatchImg;
    UINT8 *buf = nullptr;
    if(filename.substr(filename.size() - 3, 3) == "yuv") {
        FILE *yuv_file = fopen(filename.c_str(), "rb+");
        buf = new UINT8[WIDTH * HEIGHT];
        fread(buf, WIDTH * HEIGHT, 1, yuv_file);
        template_image = get_y_from_yuv(buf, WIDTH, HEIGHT);
    } else {
        template_image = cv::imread(argv[1], 0);
        tobeMatchImg = cv::imread(argv[2], 0);
        colorImg = cv::imread(argv[1], -1);
        buf = template_image.data;
    }

    findTargetArea(template_image, tobeMatchImg);

//    saveMat(integralSum, "data//integralSum1.txt");
//    saveMat(tobeSum, "data//integralSum2.txt");


//    cv::imshow("integral", template_integral);
//    cv::imshow("tobeintegral", tobe_integral);
//    cv::imshow("tobeimg", tobeMatchImg);
//    cv::imshow("template_image", template_image);
//    cv::imshow("sum", integralSum);
//    cv::waitKey(0);
    return 0;
}