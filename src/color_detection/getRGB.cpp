#include "cv.h"
#include "highgui.h"
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include "iostream"
using namespace std;

IplImage* pImg,*imgRGB,*imgHSV;
int flags = 0;
CvPoint pt;
CvScalar s = {0.0},ss={0.0};

void on_mouse( int event, int x, int y, int flags, void* param )
{

	if( !imgRGB )
		return;

	switch(event)
	{
	case CV_EVENT_LBUTTONDOWN: 
		{
			s=cvGet2D(imgRGB,y,x);
			printf("(%d,%d)RGB = %f,G = %f, R = %f\n",x,y,s.val[0],s.val[1],s.val[2]);
			ss = cvGet2D(imgHSV,y,x);
			printf("(%d,%d)H = %f,S = %f, V = %f\n",x,y,ss.val[0],ss.val[1],ss.val[2]);
		}
		break;
	}
}
int main( int argc, char** argv )
{
	imgRGB = cvLoadImage( "images/555.jpg", -1);
	imgHSV = cvCreateImage(cvGetSize(imgRGB),8,3);
	cvNamedWindow( "imgRGB", 2);
	cvSetMouseCallback( "imgRGB", on_mouse, 0 );
	cvShowImage( "imgRGB", imgRGB );
	cvCvtColor(imgRGB,imgHSV,CV_RGB2HSV);


	cvWaitKey();
	cvDestroyWindow( "imgRGB" );
	cvReleaseImage( &imgRGB );
	return 0;
}