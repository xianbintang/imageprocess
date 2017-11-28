#include <iostream>
#include <iomanip>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include "opencv/cv.h"
#include "opencv/highgui.h"
#include "TimeTracker.h"
#include "HuMatch.h"

using namespace cv;
using namespace std;

//void invariantMomentHu(const Mat src, vector<double>& result);

CvHuMoments cvHu(char * filename)
{
	IplImage*img = cvLoadImage(filename,0);
	cout << "hu of " << filename << endl;
	 CvMoments moments;
     CvHuMoments hu;
	 cvMoments(img,&moments,0);
     cvGetHuMoments(&moments, &hu);
     cout<<hu.hu1<<endl<<hu.hu2<<endl<<hu.hu3<<endl<<hu.hu4<<endl<<hu.hu5<<endl<<hu.hu6<<endl<<hu.hu7<<endl<<endl<<endl;
	 return hu;
}


int main(int argc, char **argv)
{
    vector<double> TemplateHuMoments;
    IplImage *img = cvLoadImage(argv[1], 0);
    IplImage * templateImg = cvCreateImage(cvGetSize(img), IPL_DEPTH_8U, 1);

//    cvCanny(img, templateImg, 150, 220);
    cvThreshold(img, templateImg, 100, 255, 0);
    HuMatch hm;
    hm.CreateTemplate(templateImg);

    IplImage * SearchImg = cvLoadImage(argv[2], 0);
    IplImage * si = cvCreateImage(cvGetSize(SearchImg), IPL_DEPTH_8U, 1);
//    cvCanny(SearchImg, si, 150, 220);
    cvThreshold(SearchImg, si, 100, 255, 0);
//    cvShowImage("hehe", si);

    hm.printMoments();

    TimeTracker tt;
    tt.start();
    hm.matchThis(si, std::stod(argv[3]));
    tt.stop();

    cvShowImage("result", si);
    std::cout << "Time used: " << tt.duration() << std::endl;
    waitKey(0);
	return 0;
}



