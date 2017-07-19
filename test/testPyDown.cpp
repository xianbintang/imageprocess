//
// Created by xianb on 2017/6/9.
//


//-----------------------------------【头文件包含部分】---------------------------------------
//      描述：包含程序所依赖的头文件
//----------------------------------------------------------------------------------------------
#include <opencv2/opencv.hpp>
#include "TimeTracker.h"
#include <opencv2/imgproc/imgproc.hpp>
#include <fstream>

//-----------------------------------【命名空间声明部分】---------------------------------------
//      描述：包含程序所使用的命名空间
//-----------------------------------------------------------------------------------------------
using namespace cv;
//-----------------------------------【main( )函数】--------------------------------------------
//      描述：控制台应用程序的入口函数，我们的程序从这里开始
//-----------------------------------------------------------------------------------------------
int main(int argc, char ** argv )
{
    /*
#if 0
    TimeTracker tt;
    tt.start();
    double lookupTableX[300][300];
    double lookupTableY[300][300];
    for (int k = 0; k < 10; ++k) {
        for(int i = 0; i < 300; i++) {
            for (int j = 0; j < 300; ++j) {
                double length;
                length = sqrt(i * i + j * j);
                lookupTableX[i][j] = 1.0 * i / length;
                lookupTableY[i][j] = 1.0 * j / length;

                lookupTableX[i][j] = 1.0 * i / length;
                lookupTableY[i][j] = 1.0 * j / length;
                lookupTableX[i][j] = 1.0 * i / length;
                lookupTableX[i][j] = 1.0 * i / length;
                lookupTableY[i][j] = 1.0 * j / length;
                lookupTableY[i][j] = 1.0 * j / length;
                lookupTableY[i][j] = 1.0 * j / length;
                lookupTableY[i][j] = 1.0 * j / length;
            }
        }

    }

    tt.stop();
    std::cout << "time used: "<< tt.duration() << std::endl;

//    FILE *fpx = fopen("lookupTableX.txt", "wb");
//    FILE *fpy = fopen("lookupTableY.txt", "wb");
//    fwrite(lookupTableX, sizeof(lookupTableX), 1, fpx);
//    fwrite(lookupTableY, sizeof(lookupTableY), 1, fpy);
    return 0;
#endif
#if 0
    std::ifstream fin("counts.txt");
    int num;
    int count100= 0;
    int count80 = 0;
    int count60 = 0;
    int count40 = 0;
    int count20 = 0;
    int count10 = 0;
    int count0 = 0;
    while(fin >> num) {
        if(num >= 100)
            count100++;
        else if (num >= 80)
            count80++;
        else if (num >= 60)
            count60++;
        else if (num >= 40)
            count40++;
        else if (num >= 20)
            count20++;
        else if (num >= 10)
            count10++;
        else if (num >= 0)
            count0++;
    }
    std::cout << "100: " << count100 << " " << 1.0 * count100 / 57600* 100 << "%" << std::endl;
    std::cout << "80: " << count80 << " " <<  1.0 * count80 / 57600* 100 << "%" << std::endl;
    std::cout << "60: " << count60 <<  " " <<  1.0 * count60 / 57600* 100 << "%" << std::endl;
    std::cout << "40: " << count40 <<  " " <<  1.0 * count40 / 57600* 100 << "%" << std::endl;
    std::cout << "20: " << count20 <<  " " <<  1.0 * count20 / 57600* 100 << "%" << std::endl;
    std::cout << "10: " << count10 <<  " " <<  1.0 * count10 / 57600* 100 << "%" << std::endl;
    std::cout << "0: " << count0 <<  " " <<  1.0 * count0 / 57600* 100 << "%" << std::endl;

    return 0;
#endif
     */
#if 1
    //载入原始图
    Mat srcImage = imread(argv[1]);  //工程目录下应该有一张名为1.jpg的素材图
    Mat level4, level3, level2, level1;//临时变量和目标图的定义

//    imshow("test", srcImage);
//    cvWaitKey(-1);
    Mat grayImage; //
    cvtColor(srcImage, grayImage, CV_BGR2GRAY); // 转为灰度图像
    //显示原始图

    //进行向下取样操作
    pyrDown( grayImage, level4, Size( grayImage.cols/2, grayImage.rows/2 ) );
    pyrDown( level4, level3, Size( level4.cols/2, level4.rows/2 ) );
    pyrDown( level3, level2, Size( level3.cols/2, level3.rows/2 ) );
    pyrDown( level2, level1, Size( level2.cols/2, level2.rows/2 ) );
    //显示效果图
//    Mat dstImage1;//临时变量和目标图的定义

//    pyrDown( level4, dstImage1, Size( level4.cols/2, level4.rows/2 ) );

//    imshow("level 4", level4);
//    imshow("level 3", level3);
//    imshow("level 2", level2);
//    imshow("level 1", level1);
    imwrite("images/Tgray.bmp", grayImage);
    imwrite("images/Tlevel4.bmp", level4);
    imwrite("images/Tlevel3.bmp", level3);
    imwrite("images/Tlevel2.bmp", level2);
    imwrite("images/Tlevel1.bmp", level1);
//    Mat dstImage2;//临时变量和目标图的定义
//
//    pyrDown( dstImage1, dstImage2, Size( dstImage1.cols/2, dstImage1.rows/2 ) );
//    imshow("pyrdown", dstImage2);

//    waitKey(0);

    return 0;
#endif
}