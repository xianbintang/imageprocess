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
const int NBLOCK = 8;
void saveMat(cv::Mat mat, const char *path);
void saveMatf(cv::Mat mat, const char *path);

cv::Mat colorImg;
cv::Mat colorImg4;
cv::Mat colorImg2;
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
    if(posA.y < integral.rows && posA.x < integral.cols && posD.y < integral.rows && posD.x < integral.cols) {
        A = integral.at<int>(posA.y, posA.x);
        B = integral.at<int>(posA.y, posD.x);
        C = integral.at<int>(posD.y, posA.x);
        D = integral.at<int>(posD.y, posD.x);
    } else {
//        std::cout << "out of border" << std::endl;
        std::cout << "out of boarder: " << __LINE__ << std::endl;
    }

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
            if(i < sumPattern.rows && j < sumPattern.cols) {
                sumPattern.at<int>(i, j) = curSum;
            } else {
//                std::cout << "out of border" << std::endl;
                std::cout << "out of boarder: " << __LINE__ << std::endl;
            }
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
static double dist_to_line(cv::Point pt2, cv::Point pt1, cv::Point target)
{
    // 根据中心点与直线的距离 排除干扰直线
    // 点(x0,y0)到直线Ax+By+C=0的距离为d = (A*x0+B*y0+C)/sqrt(A^2+B^2)
    double A, B, C,dis;
    // 化简两点式为一般式
    // 两点式公式为(y - y1)/(x - x1) = (y2 - y1)/ (x2 - x1)
    // 化简为一般式为(y2 - y1)x + (x1 - x2)y + (x2y1 - x1y2) = 0
    // A = y2 - y1
    // B = x1 - x2
    // C = x2y1 - x1y2
    A = pt2.y - pt1.y;
    B = pt1.x - pt2.x;
    C = pt2.x * pt1.y - pt1.x * pt2.y;

    double x,y;
    x = target.x;
    y = target.y;
    // 距离公式为d = |A*x0 + B*y0 + C|/√(A^2 + B^2)
    dis = abs(A * x + B * y + C) / sqrt(A * A + B * B);
    return dis;
}

static bool dist_to_lines_less_than(const std::vector<cv::Point> &rect, cv::Point target, double min_dist)
{
    return !((dist_to_line(rect[0], rect[1], target) < min_dist || dist_to_line(rect[1], rect[2], target) < min_dist || \
    dist_to_line(rect[2], rect[3], target) < min_dist || dist_to_line(rect[3], rect[0], target) < min_dist));
}
static int rotate_rect(std::vector<cv::Point> &rect, const std::vector<float> rotate_matrix)
{
    int plx = rect[0].x, ply = rect[0].y, prx = rect[3].x, pry = rect[3].y;
    int plxb = rect[1].x, plyb = rect[1].y, prxb = rect[2].x, pryb = rect[2].y;

    rect[0].x = rotate_matrix[0] * plx + rotate_matrix[1] * ply + rotate_matrix[2];
    rect[0].y = rotate_matrix[3] * plx + rotate_matrix[4] * ply + rotate_matrix[5];

    rect[1].x = rotate_matrix[0] * plxb + rotate_matrix[1] * plyb + rotate_matrix[2];
    rect[1].y = rotate_matrix[3] * plxb + rotate_matrix[4] * plyb + rotate_matrix[5];

    rect[2].x = rotate_matrix[0] * prxb + rotate_matrix[1] * pryb + rotate_matrix[2];
    rect[2].y = rotate_matrix[3] * prxb + rotate_matrix[4] * pryb + rotate_matrix[5];

    rect[3].x = rotate_matrix[0] * prx + rotate_matrix[1] * pry + rotate_matrix[2];
    rect[3].y = rotate_matrix[3] * prx + rotate_matrix[4] * pry + rotate_matrix[5];
    return 0;
}
void Dilation(const cv::Mat &src, cv::Mat &dilation_dst, int size )
{
    int dilation_type = cv::MORPH_RECT;

    cv::Mat element = cv::getStructuringElement( dilation_type,
                                         cv::Size( size, size));
    ///膨胀操作
    dilate( src, dilation_dst, element);
}
cv::Mat computeIntegral(const cv::Mat &image, const std::vector<cv::Point> &cur_rect)
{
    cv::Mat imgGB, imgCanny;
    cv::GaussianBlur(image, imgGB, cv::Size(5,5),0);
    cv::Canny(imgGB, imgCanny, 30, 150);

    if(cur_rect.size()) {
//        std::vector<cv::Point> cur_rect = {{0,0}, {0, image.rows - 1}, {image.cols - 1, image.rows - 1}, {image.cols - 1, 0}};
//        rotate_rect(cur_rect, rotate_matrix);
//        cv::rectangle(imgCanny, cur_rect[0],  cur_rect[3],  cv::Scalar(255,255,255));
        for (int i = 0; i < imgCanny.rows; ++i) {
            for (int j = 0; j < imgCanny.cols; ++j) {
                cv::Point p;
                p.x = j;
                p.y = i;
                unsigned char U8 = imgCanny.at<uchar>(i, j);
                double min_dist = (imgCanny.rows + imgCanny.cols) / 50.0;
                //          min_dist = (min_dist > MIN_DIST ? min_dist : MIN_DIST);
                // todo 这里由于使用了位与操作，所以不用再判断距离了
                if (U8 && !dist_to_lines_less_than(cur_rect, p, min_dist)) {
                    imgCanny.at<uchar>(i, j) = 0;
                }
            }
        }
    }


    Dilation(imgCanny, imgCanny, 3);
//    cv::imshow("canny", imgCanny);
//    cv::waitKey(0);
    std::cout << "after canny: " << computePointsChar(imgCanny) / 255 << std::endl;;
    cv::Mat integ;
    cv::integral(imgCanny, integ,CV_32S); //计算积分图

    return integ;
}

bool compareSumPattern(const cv::Mat &srcPattern, const cv::Mat &targetPattern, const cv::Point position, double thresh)
{
    int ct = 0;
    int ctnozero = 0;


    TimeTracker tt;
    tt.start();
    for (int i = 0; i < targetPattern.rows; ++i) {
        for (int j = 0; j < targetPattern.cols; ++j) {
            int targetCount;

            if(i < targetPattern.rows && j < targetPattern.cols) {
                targetCount = targetPattern.at<int>(i, j);
            } else{
//                std::cout << "out of boarder" << std::endl;
                std::cout << "out of boarder: " << __LINE__ << std::endl;
            }

            int srcCount;
            cv::Point A(j + position.x, i + position.y);

            if(A.y >= srcPattern.rows && A.y >= srcPattern.cols) {
                A.y = srcPattern.rows - 1;
                A.x = srcPattern.cols - 1;
            }
            srcCount = srcPattern.at<int>(A.y, A.x);

//            if(thresh > 0.78 && position.x == 47 && position.y == 47) {
//                cv::Mat tm = colorImg2;
//                cv::rectangle(tm, A,  cv::Point(A.x, A.y),  cv::Scalar(255,255,255));
//                cv::imshow("tmp", tm);
//                cv::waitKey(-1);
//            }

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

//                if(position.x == 6 && position.y == 6) {
////                    std::cout << "[" <<targetCount << ", " << srcCount << ", " <<  threshold <<  "] ";
//                    if(abs(srcCount - targetCount) <= threshold) {
////                        std::cout << "ok------ " << std::endl;;
//                    } else {
////                        std::cout << "bad***** " << std::endl;;
//                    }
//                }
//                if(thresh > 0.78)
//                    std::cout << "[" <<targetCount << ", " << srcCount << ", " <<  threshold <<  "] " << std::endl;;
                if(abs(srcCount - targetCount) <= threshold) {
                    ct++;
                }
            }
        }
    }
//    std::cout << "compareSumPattern: " << ct << " " << ctnozero << std::endl;
//    std::cout << "score: " << 1.0 * ct / ctnozero << "ct: " << ct << "ctnozero: " << ctnozero << std::endl;
    if(1.0 * ct / ctnozero> thresh) {
//        std::cout << "score: " << 1.0 * ct / ctnozero << "ct: " << ct << "ctnozero: " << ctnozero << "threshold: " << thresh << std::endl;
//        std::cout << "position: " << position << "score: " << std::endl;
    }
    tt.stop();
//    std::cout << "compare time: " << tt.duration() << std::endl;
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

//        std::cout << " ct: " << computePointsInt(srcPattern) << std::endl;

//        cv::imshow("s", srcPattern8);
//        cv::imshow("t", templatePattern8);
//        cv::waitKey(0);
    Point A, D;
    A = Point(0, 0);
    D=  Point(src.cols - 1, src.rows - 1);
    auto integ= computeIntegral(src, {});
    auto srcPattern8 = CreateIntegralSum(integ, A, D, 8);
    auto srcPattern4 = CreateIntegralSum(integ, A, D, 4);
    auto srcPattern2 = CreateIntegralSum(integ, A, D, 2);


    cv::pyrDown(colorImg, colorImg2);
    cv::pyrDown(colorImg2, colorImg4);

    tt.start();
    int ct = 0;
    Mat template_integs[360][3];
    for (int i = 0; i < 360; ++i) {
        target = origin_target;
        cv::Mat rotated_image;
        // 还是无法保证完全在图片框内
        auto rotate_matrix = rotate_image(target, rotated_image, cv::Point(target.rows / 2, target.cols / 2), i);
        target = rotated_image;

        //        target = origin_target;
        //                    imshow("tar", target);
        //                    waitKey(0);
        //        imshow("sar", src);
        //        rotate_image(target, rotated_image, cv::Point(target.rows / 2, target.cols / 2), d);
        //        imshow("rotar", rotated_image);

        colorImg = src;
        // 返回指向需要被发送的内存缓冲区的指针

        std::vector<cv::Point> cur_rect = {{0,0}, {0, origin_target.rows - 1}, {origin_target.cols - 1, origin_target.rows - 1}, {origin_target.cols - 1, 0}};
        rotate_rect(cur_rect, rotate_matrix);

        auto tobe= computeIntegral(target, cur_rect);
        Point A = Point(0, 0), D=  Point(target.cols - 1, target.rows - 1);
        //        auto tobeSum = getCurrentSum(target, A, D, 0);

        int patternSize = 4;
        // 计算模板的pattern
        auto templatePattern8 = CreateIntegralSum(tobe, A, D, 8);
        auto templatePattern4 = CreateIntegralSum(tobe, A, D, 4);
        auto templatePattern2 = CreateIntegralSum(tobe, A, D, 2);
        template_integs[i][0] = templatePattern8;
        template_integs[i][1] = templatePattern4;
        template_integs[i][2] = templatePattern2;
//        std::cout << " ct: " << computePointsInt(templatePattern) << std::endl;

    }

    std::vector<Point> region;
    for (int i = 0; i < src.rows; i += BLOCKH) {
        for (int j = 0; j < src.cols; j += BLOCKW) {
            Point A = Point(j, i), D = Point(j + NBLOCK * BLOCKW, i + NBLOCK * BLOCKH);
            int sum = getCurrentSum(integ, A, D, 0);
            if(sum > 1300&& sum < 1600) {
                // FIXME 这里要修改，不是从右上角进行扩散匹配，而是从中心位置进行扩散匹配。
                // FIXME 这里都没有进行扩展，而是直接只在目标位置上进行搜索，所以肯定会导致匹配的不准。。。

                // 在这里旋转角度
                for (int d = 0; (int)d < 360; d += 1) {
#if 0
                    src = origin_src;
                    cv::Mat rotated_image;
                    // 还是无法保证完全在图片框内
                    auto rotate_matrix = rotate_image(src, rotated_image, cv::Point(src.rows / 2, src.cols / 2), d);
                    src = rotated_image;
#endif

                    auto &templatePattern8 = template_integs[d][0];
                    auto &templatePattern4 = template_integs[d][1];
                    auto &templatePattern2 = template_integs[d][2];

                    // 计算图片的pattern

//                    for (int k = -BLOCKW / patternSize; k < BLOCKW / patternSize; ++k) {
                    int rg = 4;
                    for (int rangei = -rg; rangei <= rg; rangei += 1) {
                        for (int rangej = -rg; rangej <= rg; rangej += 1) {
                            Point pos = Point(j/4+ rangej, i/4+ rangei);
                            if(pos.x < 0) {
                                pos.x = 1;
                            }
                            if(pos.y < 0) {
                                pos.y = 1;
                            }
//                        cv::rectangle(colorImg4, pos,  cv::Point(pos.x + NBLOCK * BLOCKW / 16, pos.y + NBLOCK * BLOCKW/ 16),  cv::Scalar(255,255,255));
//                        cv::circle(colorImg4, pos, 1, cv::Scalar(255,255,255));
//                        cv::imshow("tmpX", colorImg4);
//                        cv::waitKey(5);
                            if(compareSumPattern(srcPattern4, templatePattern4, pos, 0.70)) {
//                            cv::rectangle(colorImg4, pos,  cv::Point(pos.x + NBLOCK * BLOCKW / 16, pos.y + NBLOCK * BLOCKW/ 16),  cv::Scalar(255,255,255));
//                                cv::circle(colorImg4, pos, 1, cv::Scalar(255,255,255));
//                                cv::imshow("tmpX", colorImg4);
//                                cv::waitKey(0);

                                // pos是目标位置，从目标位置进       行扩展
//                            std::cout << "pos: " << pos << std::endl;
                                int rg = 4;
                                for (int m = -rg; m <= rg; m += 1) {
                                    for (int n = -rg; n <= rg; n += 1) {
                                        Point pos1 = Point(pos.x * 2 + n, pos.y * 2 + m);
                                        if(pos1.x < 0) {
                                            pos1.x = 1;
                                        }
                                        if(pos1.y < 0) {
                                            pos1.y = 1;
                                        }
//                                std::cout << "pos1: " << pos1 << std::endl;
//                                        if(A.x == 96 && A.y == 96) {
//                                            std::cout << "pos1: " << pos1 << std::endl;
//                                        }
                                        if(compareSumPattern(srcPattern2, templatePattern2, pos1, 0.80)) {
                                            region.push_back(cv::Point(pos1.x * 2, pos1.y * 2));
                                            ct++;
                                            std::cout << "degree: " << d << std::endl;
                                        }
                                    }

                                }
                            }
                        }

                    }
                }

//            std::cout << sum << " ";
            }
        }

//    std::cout << "targets: " << ct << " degree: " << d << std::endl;
//    cv::imshow("integimg", integ);
//    cv::imshow("tobe", tobe);
    }

    tt.stop();
    totalTime += tt.duration();
    for (int i = 0; i < region.size(); ++i) {
//        cv::circle(colorImg, region[i], 1, cv::Scalar(0,0,255));
//        if(abs(region[i].x - 32) < 3 && abs(region[i].y - 32) < 3) {
        std::cout << region[i] << std::endl;
        cv::circle(colorImg, cv::Point(region[i].x + NBLOCK * BLOCKW / 2, region[i].y + NBLOCK * BLOCKH / 2), 1, cv::Scalar(255,255,255));
        cv::rectangle(colorImg, region[i],  cv::Point(region[i].x + NBLOCK * BLOCKW, region[i].y + NBLOCK * BLOCKH),  cv::Scalar(255,255,255));
        cv::imshow("colorImg", colorImg);
        cv::waitKey(0);
//        }
    }
    std::cout << "duration: " << tt.duration() << std::endl;
    std::cout << "total: " << ct << "time: " << totalTime << std::endl;
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