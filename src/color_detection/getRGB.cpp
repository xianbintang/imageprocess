#include "cv.h"
#include "highgui.h"
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <opencv2/imgproc.hpp>
#include <opencv/cv.hpp>
#include <c++/fstream>
#include "iostream"
using namespace std;

cv::Mat pImg,imgRGB,imgHSV;
int flags = 0;
CvPoint pt;
CvScalar s = {0.0},ss={0.0};

void on_mouse( int event, int x, int y, int flags, void* param )
{

//    std::cout << "haha" << std::endl;
	switch(event)
	{
	case CV_EVENT_LBUTTONDOWN: 
		{
			auto rgbv =imgRGB.at<cv::Vec3b>(y,x);
            auto hsvv =imgHSV.at<cv::Vec3b>(y,x);
			printf("(%d,%d)B = %d,G = %d, R = %d\n",x,y,rgbv[0],rgbv[1],rgbv[2]);
            printf("(%d,%d)H= %d,S = %d, V = %d\n",x,y,hsvv[0],hsvv[1],hsvv[2]);
		}
		break;
	}
}

void printHSV(cv::Mat hsvimage) {
    std::ofstream of;
    of.open("data//csc.txt");
    for (int i = 0; i < hsvimage.rows; ++i) {
        for (int j = 0; j < hsvimage.cols; ++j) {
            auto hsvv =imgHSV.at<cv::Vec3b>(i,j);
            //printf("(%d,%d)H= %d,S = %d, V = %d\n",i,j,hsvv[0],hsvv[1],hsvv[2]);
            of<< "H: " << (int)hsvv[0] << " " << "S: " <<(int)hsvv[1] << " " << "V: " <<(int)hsvv[2] << " i: " << i << " " << "j: " << j << " " << std::endl;;
        }
    }
    of.close();
}
int main( int argc, char** argv )
{
#if 0
	imgRGB = cv::imread( argv[1], -1);
//	imgHSV = cvCreateImage(cvGetSize(imgRGB),8,3);

	cvNamedWindow( "imgRGB", 2);
	cvSetMouseCallback( "imgRGB", on_mouse, 0 );
	cv::imshow( "imgRGB", imgRGB);
	cvtColor(imgRGB,imgHSV,CV_RGB2HSV);

	cvWaitKey();
//	cvDestroyWindow( "imgRGB" );
//	cvReleaseImage( &imgRGB );
#endif
#if 1
    std::string filename(argv[1]);
    unsigned char *buf = nullptr;
    if(filename.substr(filename.size() - 3, 3) == "yuv") {
        FILE *yuv_file = fopen(filename.c_str(), "rb+");

        buf = new unsigned char[640 * 480 + 320 * 240 * 2];
        fread(buf, 640* 480 + 320 * 240 * 2, 1, yuv_file);
    }
    cv::Mat imgYUV(480+ 480/2, 640, CV_8UC1, (void*) buf);
//    imgRGB.create(480, 640, CV_8UC3);
    cv::cvtColor(imgYUV, imgRGB, CV_YUV420sp2BGR, 3);
    cv::cvtColor(imgRGB, imgHSV, CV_BGR2HSV_FULL, 3);

    printHSV(imgHSV);
    cvNamedWindow("imgRGB1", 2);
    cvSetMouseCallback("imgRGB1", on_mouse, 0 );
    cv::imshow( "imgRGB1", imgRGB );

    cvWaitKey();
    return 0;
#endif


	return 0;
}