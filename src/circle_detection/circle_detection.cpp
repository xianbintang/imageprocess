//
// Created by xianb on 2017/7/18.
//


#include <opencv2/opencv.hpp>
#include <opencv/cv.h>

int main(int argc, char ** argv)
{
    // 检测图像中的圆形
    cv::Mat image= cv::imread(argv[1], 0);
    if (!image.data)
    {
        std::cout << "error" << std::endl;
        return 0;
    }
    // 在调用cv::HoughCircles函数前对图像进行平滑，减少误差
    cv::GaussianBlur(image,image,cv::Size(7,7),1.5);
    std::vector<cv::Vec3f> circles;
    cv::HoughCircles(image, circles, CV_HOUGH_GRADIENT,
                     2,   // 累加器的分辨率(图像尺寸/2)
                     50,  // 两个圆之间的最小距离
                     200, // Canny中的高阈值
                     200, // 最小投票数
                     400, 500); // 有效半径的最小和最大值

    // 绘制圆圈
    image= cv::imread(argv[1],0);
    if (!image.data)
    {
        return 0;
    }
    // 一旦检测到圆的向量，遍历该向量并绘制圆形
    // 该方法返回cv::Vec3f类型向量
    // 包含圆圈的圆心坐标和半径三个信息
    std::vector<cv::Vec3f>::const_iterator itc= circles.begin();
    while (itc!=circles.end())
    {
        cv::circle(image,
                   cv::Point((*itc)[0], (*itc)[1]), // 圆心
                   (*itc)[2], // 圆的半径
                   cv::Scalar(255), // 绘制的颜色
                   6); // 圆形的厚度
        ++itc;
    }

    cv::namedWindow("Detected Circles");
    cv::imshow("Detected Circles",image);
    cvWaitKey(-1);
    return 0;
}