//
// Created by xianb on 2017/12/20.
//


#include <opencv2/core/mat.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <iostream>
#include <opencv/cv.hpp>

static void saveMat(cv::Mat mat, const char *path) {
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
int main(int argc, char ** argv)
{
    cv::Mat src, dst;
    std::vector<std::vector<float>> k
            {
//                    {1, 4, 7, 4, 1},
//                    {4, 16, 26, 16, 4},
//                    {7, 26, 41, 26, 7},
//                    {4, 16, 26, 16, 4},
//                    {1, 4, 7, 4, 1}

                    {1,2,1},
                    {2,4,2},
                    {1,2,1},

            };

    cv::Mat kernel(3,3,CV_32F);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            kernel.at<float>(i, j) = k[i][j]/16.0;
            std::cout << kernel.at<float>(i, j) << std::endl;
        }
    }

    src = cv::imread(argv[1], -1);
    cv::filter2D(src, dst, -1, kernel, cv::Point(-1, -1), 0, cv::BORDER_DEFAULT);
//    cv::imshow("hehe", dst);
//    cvWaitKey(0);
    saveMat(dst, "dst.txt");
}


