#if 1
#include <opencv2\opencv.hpp>
#include "opencv/cv.h"
#include "opencv/highgui.h"
#include <iostream>
#include <fstream>
#include <string>
#include "plate.h"
using namespace cv;
using namespace std;
//#include "Point.h"
IplImage*  getContourImage(IplImage * test_img) 
{
	IplImage * img_gray = cvCreateImage(cvGetSize(test_img), IPL_DEPTH_8U, 1);
	//IplImage * img_after_filter = cvCreateImage(cvGetSize(test_img), IPL_DEPTH_8U, 1);
	IplImage * img_after_stre = cvCreateImage(cvGetSize(test_img), IPL_DEPTH_8U, 1);
	IplImage * img_after_canny = cvCreateImage(cvGetSize(test_img), IPL_DEPTH_8U, 1);
	IplImage * img_after_filter = cvCreateImage(cvGetSize(test_img), IPL_DEPTH_8U, 1);
	IplImage * img_final = cvCreateImage(cvGetSize(test_img), IPL_DEPTH_8U, 1);
	IplImage * img_mor = cvCreateImage(cvGetSize(test_img), IPL_DEPTH_8U, 1);
	IplImage * tmp = cvCreateImage(cvGetSize(test_img), IPL_DEPTH_8U, 1);
	/*一:灰度化 cvCvtColor(IplImage* src, IplImage* dst, int color)*/
	cvCvtColor(test_img, img_gray, CV_RGB2GRAY);
	/*二:灰度拉伸*/
	cvNormalize(img_gray, img_after_stre, 0, 255, CV_MINMAX);
	/*gray_strecth(img_gray, img_after_stre, 255, 0);  有自带的函数比这个好用*/
	/*三:滤波*/
	img_after_stre = img_gray;
	cvSmooth(img_after_stre, img_after_filter, CV_GAUSSIAN);
	/*四:边缘检测*/ 
	cvCanny(img_after_filter, img_after_canny, 70, 150, 3);	/*比用sobel函数好,因为这个不会分方向,把所有的边缘都画了出来了*/
	//cvSobel(img_after_stre, img_after_sobel, 1, 0, 3);
	/*五:二值化*/
	cvThreshold(img_after_canny, img_final, 0, 255, CV_THRESH_BINARY| CV_THRESH_OTSU);
	return img_final;
}
void getContours(IplImage * img, CvSeq ** contours  )
{
	CvMemStorage * storage = cvCreateMemStorage(0);
	cvFindContours(img, storage, contours);
}
void saveMat(IplImage *img)
{
	auto fp = fopen("image.txt", "w");
	for(int i = 0; i < img->height; i++) {
		unsigned char * srow = (unsigned char *)(img->imageData + i * img->widthStep);
		for(int j = 0; j < img->width; j++) {
			//if(srow[j] == 0)
			//	fprintf(fp, "  ", 0);
			//else
				fprintf(fp, "%d ", srow[j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

void saveMat(int **img, int h, int w)
{
	auto fp = fopen("image.txt", "w");
	for(int i = 0; i < h; i++) {
		for(int j = 0; j < w; j++) {
			if(img[i][j] == 0)
				fprintf(fp, "%d ", 0);
			else
				fprintf(fp, "%d ", 1);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}
int main() 
{
	IplImage * test_img1 = NULL;
	test_img1 = cvLoadImage("test1.bmp", -1);
	IplImage * result = cvCreateImage (cvGetSize(test_img1), IPL_DEPTH_8U, 1);
	IplImage * temp = cvCreateImage (cvGetSize(test_img1), IPL_DEPTH_8U, 1);
	IplImage * open = cvCreateImage (cvGetSize(test_img1), IPL_DEPTH_8U, 1);
	IplImage * img_final1 = getContourImage(test_img1);
	cvShowImage("img_final1", img_final1);
	cvSaveImage("contour.bmp", img_final1);
	saveMat(img_final1);

	IplConvKernel * t =  cvCreateStructuringElementEx(
		2,2,0,0, CV_SHAPE_RECT);
	cvDilate(img_final1, open, NULL, 20);
//	cvErode(open, open, t);
	//cvMorphologyEx(img_final1,open, temp, NULL, CV_MOP_OPEN);

	int **mat = new int*[img_final1->height];
		for(int i = 0; i < img_final1->height; i++)
			mat[i] = new int[img_final1->width];
		for(int i = 0; i < img_final1->height; i++) {
			unsigned char * srow = (unsigned char *)(img_final1->imageData + i * img_final1->widthStep);
			for(int j = 0; j < img_final1->width; j++) {
				mat[i][j] = srow[j];
			}
		}
		saveMat(mat, img_final1->height, img_final1->width);
	try{

//		DFS DFSTest = DFS(img_final1);
//		vector<vector<Point>> vps = DFSTest.getPoints();

		ifstream fin("allpoints.txt");
		vector<Point> v;
		int n;
		while(fin >> n) {
			v.push_back(Point(n % 320, n / 320));
		}
		for(int i = 0; i < v.size(); i++) {
			cvCircle(result, Point(v[i]),1 ,cvScalar(0, 0, 255));
		}
#if 0
		for(int i = 0; i < vps.size(); i++){
			for(int j = 0; j < vps[i].size(); j++) {
				//	cout << vps[i][j].x << " " <<  vps[i][j].y << endl;
				cvCircle(result, Point(vps[i][j].y, vps[i][j].x),1 ,cvScalar(0, 0, 255));
			}
			char  filename[20];
			sprintf(filename, "result%d.bmp", i);
		//	cvSaveImage(filename, result);
		}
#endif
		cvShowImage("result", result);

	} catch(std::exception &e) {
		cout << e.what() << endl;
	}
#if 0
		getContours(img_final1,&contours1);
		//cvDrawContours(test_img1, contours1,CV_RGB(255, 255, 255), CV_RGB(255, 255, 255), 2, CV_FILLED, 8, cvPoint(0, 0));
		//cvShowImage("test", test_img1);
		getContours(img_final2,&contours2);
		for	(CvSeq * c = contours2; c != NULL; c = c->h_next) {
			double result = cvMatchShapes(contours1, c, 1);   // 根据输入的图像或轮廓来计算它们的hu矩的相似度  
			cout << result << endl;  
		}
#endif
#if 0
		CvMemStorage * storage = cvCreateMemStorage(0);
		CvSeq * contours = NULL;
		List rects = create_list();
		List rects_final = create_list();
		get_contour_rect(img_final, rects, storage, contours);/*没什么可以改进的地方*/
		draw_contour_rect(img_final, rects);			/*显示未筛选前的矩形位置,看看到底有没有把车牌的位置找到*/
		//cvNamedWindow("img_final", 1);
		cvShowImage("img_gray", img_gray);
		cvShowImage("img_after_stre", img_after_stre);
		cvShowImage("img_after_canny", img_after_canny);
		cvShowImage("img_after_filter", img_after_filter);
		//cvShowImage("img_final", img_final);
		cvShowImage("img_final", img_final);
		const int params[2]={CV_IMWRITE_JPEG_QUALITY,100};
		cvSaveImage("img_after_filter.bmp",img_after_filter,params);
		IplConvKernel* kernal = cvCreateStructuringElementEx(2,1, 1, 0, CV_SHAPE_RECT);
		cvMorphologyEx(
			img_final,
			img_mor,
			tmp,
			kernal, //default 3*3
			CV_MOP_OPEN,
			3);
		cvShowImage("img_mor", img_mor);
#endif
		waitKey();
    return 0;
}
CvSeq *getImageContours(CvArr *src)    
{    
    cvThreshold(src, src, 100, 255, CV_THRESH_BINARY);    
    CvMemStorage * storage = cvCreateMemStorage(0);    
    CvSeq * contours;    
    cvFindContours(src, storage, &contours);    
    return contours;    
}   
#endif