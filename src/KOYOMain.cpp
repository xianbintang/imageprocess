#include <opencv2/opencv.hpp>
using namespace cv;
using namespace std;
#include "../include/common.h"

IplImage*  getEdges(IplImage *test_img)
{
	IplImage * img_gray = cvCreateImage(cvGetSize(test_img), IPL_DEPTH_8U, 1);
	//IplImage * img_after_filter = cvCreateImage(cvGetSize(test_img), IPL_DEPTH_8U, 1);
	IplImage * img_after_stre = cvCreateImage(cvGetSize(test_img), IPL_DEPTH_8U, 1);
	IplImage * img_after_canny = cvCreateImage(cvGetSize(test_img), IPL_DEPTH_8U, 1);
	IplImage * img_after_filter = cvCreateImage(cvGetSize(test_img), IPL_DEPTH_8U, 1);
	IplImage * img_final = cvCreateImage(cvGetSize(test_img), IPL_DEPTH_8U, 1);
	IplImage * img_dilate = cvCreateImage(cvGetSize(test_img), IPL_DEPTH_8U, 1);
    IplImage * img_erode = cvCreateImage(cvGetSize(test_img), IPL_DEPTH_8U, 1);
	/*一:灰度化 cvCvtColor(IplImage* src, IplImage* dst, int color)*/
	cvCvtColor(test_img, img_gray, CV_RGB2GRAY);
	/*二:灰度拉伸*/
	cvNormalize(img_gray, img_after_stre, 0, 255, CV_MINMAX);
	/*gray_strecth(img_gray, img_after_stre, 255, 0);  有自带的函数比这个好用*/
	/*三:滤波*/
	img_after_stre = img_gray;
	cvSmooth(img_after_stre, img_after_filter, CV_GAUSSIAN);

    IplConvKernel * t =  cvCreateStructuringElementEx( 3,3,1,1, CV_SHAPE_ELLIPSE );
	cvDilate(img_after_filter, img_dilate, t);
	cvErode(img_dilate, img_erode, t);

//    cvSaveImage("gray.bmp", img_erode);

//    cvSaveImage("binPicture.bmp", img_erode);

//	cvMorphologyEx(img_final1,open, temp, t, CV_MOP_OPEN);
	/*四:边缘检测*/ 
	cvCanny(img_erode, img_after_canny, 30, 150, 3);	/*比用sobel函数好,因为这个不会分方向,把所有的边缘都画了出来了*/
	//cvSobel(img_after_stre, img_after_sobel, 1, 0, 3);
	/*五:二值化*/
	cvThreshold(img_after_canny, img_final, 0, 255, CV_THRESH_BINARY| CV_THRESH_OTSU);
//	cvThreshold(img_erode, img_final, 0, 255, CV_THRESH_BINARY| CV_THRESH_OTSU);

	return img_final;
}

int main() 
{
	IplImage * test_img1 = cvLoadImage("E:\\ClionProjects\\KOYO\\images\\KOYO.bmp", -1);
	IplImage * img_final1 = getEdges(test_img1);
	cvShowImage("img_final1", img_final1);
    SaveImageAsText(img_final1, "ImageAsText.txt");

    std::vector<std::vector<int>> unions = GetUnions("ImageAsText.txt", 139, 44);
    std::cout << unions.size() << " components." << std::endl;
    displayTextImage(unions, 139, 44);

	cvWaitKey();
    return 0;
}

