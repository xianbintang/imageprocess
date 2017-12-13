//
// Created by xianb on 2017/12/12.
//


#include <vector>
#include <opencv2/core/types.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>
#include <string>
#include <opencv/cv.hpp>

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
// 旋转中心
    CvPoint2D32f center = cvPoint2D32f(centerP.x, centerP.y);
    CvMat map_matrix2 = map_matrix;
    cv2DRotationMatrix(center, degree, 1.0, &map_matrix2);
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

int main(int argc, char **argv)
{
    std::string filename(argv[1]);
    cv::Mat src = cv::imread(filename, -1);
    for (int i = 0; i < 360; ++i) {
        cv::Mat dst;
        rotate_image(src, dst, cv::Point(src.cols / 2, src.rows / 2), i);
        cv::imshow("haha", dst);
        cvWaitKey();
    }
    return 0;
}

