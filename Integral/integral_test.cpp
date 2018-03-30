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
const int NBLOCK = 4;
void saveMat(cv::Mat mat, const char *path);
void saveMatf(cv::Mat mat, const char *path);

cv::Mat colorImg;
void saveMat(const cv::Mat mat, const char *path) {
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
void saveMatf(const cv::Mat mat, const char *path) {
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
    if(posD.x >= integral.cols) {
        posD.x = integral.cols - 1;
    }

    if(posD.y >= integral.rows) {
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

cv::Mat computeIntegral(const cv::Mat &image)
{
    cv::Mat imgGB, imgCanny;
    cv::GaussianBlur(image, imgGB, cv::Size(5,5),0);
    cv::Canny(imgGB, imgCanny, 30, 150);
//    std::cout << "after canny: " << computePointsChar(imgCanny) / 255 << std::endl;;
    cv::Mat integ;
    cv::integral(imgCanny, integ,CV_32S); //计算积分图

    return integ;
}

bool compareSumPattern(const cv::Mat &srcPattern, const cv::Mat &targetPattern, const cv::Point position, double thresh)
{
    int ct = 0;
    int ctnozero = 0;
    for (int i = 0; i < targetPattern.rows; ++i) {
        for (int j = 0; j < targetPattern.cols; ++j) {
            int targetCount = targetPattern.at<int>(i, j);
            int srcCount = srcPattern.at<int>(i + position.y, j + position.x);
            int threshold = 1;

            if(targetCount) {
                ctnozero++;
                if(targetCount > 20) {
                    threshold = targetCount * 0.2;
                } else if (targetCount > 10){
                    threshold = targetCount / 3;
                } else if (targetCount > 5){
                    threshold = targetCount / 2;
                } else if(targetCount > 3){
                    threshold = targetCount - 2;
                }
                if(targetCount == 3)
                    threshold = 2;

                if(position.x == 6 && position.y == 6) {
//                    std::cout << "[" <<targetCount << ", " << srcCount << ", " <<  threshold <<  "] ";
                    if(abs(srcCount - targetCount) <= threshold) {
//                        std::cout << "ok------ " << std::endl;;
                    } else {
//                        std::cout << "bad***** " << std::endl;;
                    }
                }


                if(abs(srcCount - targetCount) <= threshold) {
                    ct++;
                }
            }
        }
    }
//    std::cout << "compareSumPattern: " << ct << " " << ctnozero << std::endl;
//    std::cout << "score: " << 1.0 * ct / ctnozero << "ct: " << ct << "ctnozero: " << ctnozero << std::endl;
    if(1.0 * ct / ctnozero> thresh) {
        std::cout << "score: " << 1.0 * ct / ctnozero << "ct: " << ct << "ctnozero: " << ctnozero << "threshold: " << thresh << std::endl;
//        std::cout << "position: " << position << "score: " << std::endl;
    }
    return 1.0 * ct / ctnozero> thresh;
}


static std::vector<float> rotate_image(const cv::Mat &src, cv::Mat &dst, cv::Point centerP, float degree)
{
    int width = src.cols;
    int height = src.rows;
    double angle = degree  * CV_PI / 180.; // 弧度
    double a = sin(angle), b = cos(angle);
    // 适当增大一点宽高，防止像素不在图像内
    int width_rotate = int(height * fabs(a) + width * fabs(b));
    int height_rotate = int(width * fabs(a) + height * fabs(b));
    float map[6];
    cv::Mat map_matrix = cv::Mat(2, 3, CV_32F, map);
//    CvPoint2D32f center = cvPoint2D32f(centerP.x, centerP.y);
// 旋转中心, 以原始图片中心作为旋转中心而不是质心
    CvPoint2D32f center = cvPoint2D32f( width / 2.0, height / 2.0);
    CvMat map_matrix2 = map_matrix;
    cv2DRotationMatrix(center, degree, 1.0, &map_matrix2);
    // 这里不能改
    map[2] += (width_rotate - width) / 2.0;
    map[5] += (height_rotate - height) / 2.0;
    cv::warpAffine(src, dst, map_matrix, cv::Size(width_rotate, height_rotate), 1, 0, 0);
    std::vector<float> rotate_matrix;
    rotate_matrix.push_back(map[0]);
    rotate_matrix.push_back(map[1]);
    rotate_matrix.push_back(map[2]);
    rotate_matrix.push_back(map[3]);
    rotate_matrix.push_back(map[4]);
    rotate_matrix.push_back(map[5]);
    return rotate_matrix;
}
int findTargetArea(cv::Mat &src, cv::Mat &target)
{
    int total = 0;
    cv::Mat origin_src = src;
    cv::Mat origin_target = target;
    TimeTracker tt;
    int totalTime = 0;
    for (double d = 0.0; (int)d < 360; d += 1) {
        src = origin_src;
        cv::Mat rotated_image;
        // 还是无法保证完全在图片框内
        auto rotate_matrix = rotate_image(src, rotated_image, cv::Point(src.rows / 2, src.cols / 2), d);
        src = rotated_image;

//        target = origin_target;
//        imshow("tar", target);
//        imshow("sar", src);
//        rotate_image(target, rotated_image, cv::Point(target.rows / 2, target.cols / 2), d);
//        imshow("rotar", rotated_image);

        colorImg = src;
        std::vector<Point> region;
        // 返回指向需要被发送的内存缓冲区的指针
        auto integ= computeIntegral(src);
        auto tobe= computeIntegral(target);
        Point A = Point(0, 0), D=  Point(target.cols - 1, target.rows - 1);
//        auto tobeSum = getCurrentSum(target, A, D, 0);

        int patternSize = 4;
        // 计算模板的pattern
        auto templatePattern8 = CreateIntegralSum(tobe, A, D, 8);
        auto templatePattern4 = CreateIntegralSum(tobe, A, D, 4);
        auto templatePattern2 = CreateIntegralSum(tobe, A, D, 2);
//        std::cout << " ct: " << computePointsInt(templatePattern) << std::endl;

        // 计算图片的pattern

        A = Point(0, 0);
        D=  Point(src.cols - 1, src.rows - 1);
        auto srcPattern8 = CreateIntegralSum(integ, A, D, 8);
        auto srcPattern4 = CreateIntegralSum(integ, A, D, 4);
        auto srcPattern2 = CreateIntegralSum(integ, A, D, 2);
//        std::cout << " ct: " << computePointsInt(srcPattern) << std::endl;

//        cv::imshow("s", srcPattern8);
//        cv::imshow("t", templatePattern8);
//        cv::waitKey(0);
        tt.start();
        int ct = 0;

        for (int i = 0; i < src.rows; i += BLOCKH) {
            for (int j = 0; j < src.cols; j += BLOCKW) {
                Point A = Point(j, i), D = Point(j + NBLOCK * BLOCKW, i + NBLOCK * BLOCKH);
                int sum = getCurrentSum(integ, A, D, 0);
                if(sum > 130 && sum < 180) {
                    // FIXME 这里要修改，不是从右上角进行扩散匹配，而是从中心位置进行扩散匹配。
                    // FIXME 这里都没有进行扩展，而是直接只在目标位置上进行搜索，所以肯定会导致匹配的不准。。。
                    for (int k = -BLOCKW / patternSize; k < BLOCKW / patternSize; ++k) {
                        Point pos = Point(j/patternSize + k, i/patternSize + k);
                        if(compareSumPattern(srcPattern4, templatePattern4, pos, 0.75)) {
                            // pos是目标位置，从目标位置进行扩展
//                            std::cout << "pos: " << pos << std::endl;
                            for (int m = -1; m < 1; ++m) {
                                Point pos1 = Point(j/2 + m, i/2 + m);
//                                std::cout << "pos1: " << pos1 << std::endl;
                                if(compareSumPattern(srcPattern2, templatePattern2, pos1, 0.65)) {
                                    region.push_back(A);
                                    ct++;
                                }
                            }
                        }
                    }

                }
//            std::cout << sum << " ";
            }
        }
    tt.stop();
        totalTime += tt.duration();
//    std::cout << "duration: " << tt.duration() << std::endl;
    for (int i = 0; i < region.size(); ++i) {
//        cv::circle(colorImg, region[i], 1, cv::Scalar(0,0,255));
        cv::circle(colorImg, cv::Point(region[i].x + NBLOCK * BLOCKW / 2, region[i].y + NBLOCK * BLOCKH / 2), 1, cv::Scalar(255,255,255));
        cv::rectangle(colorImg, region[i],  cv::Point(region[i].x + NBLOCK * BLOCKW, region[i].y + NBLOCK * BLOCKH),  cv::Scalar(255,255,255));
        cv::imshow("colorImg", colorImg);
        cv::waitKey(0);
    }
//    std::cout << "targets: " << ct << " degree: " << d << std::endl;
        total += ct;
//    cv::imshow("integimg", integ);
//    cv::imshow("tobe", tobe);
    }
    std::cout << "total: " << total << "time: " << totalTime << std::endl;
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