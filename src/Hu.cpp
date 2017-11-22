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




int main(int argc, char **argv)
{
    vector<double> TemplateHuMoments;
    IplImage *img = cvLoadImage(argv[1], 0);
    IplImage * templateImg = cvCreateImage(cvGetSize(img), IPL_DEPTH_8U, 1);

    cvThreshold(img, templateImg, 150, 255, CV_THRESH_BINARY);
    HuMatch hm;
    hm.CreateTemplate(templateImg);
    IplImage * SearchImg = cvLoadImage(argv[2], 0);
    cvThreshold(SearchImg, SearchImg, 150, 255, CV_THRESH_BINARY);

#if 1
    Mat image = imread(argv[1]);
    cvtColor(image, image, CV_BGR2GRAY);
    Moments mts = moments(image);
    double hu[7];
    HuMoments(mts, hu);
    for (int i=0; i<7; i++)
    {
        cout << hu[i] <<endl;
    }

    cout << "my print" << endl;
   hm.printMoments();
   return 0;
#endif




    TimeTracker tt;
    tt.start();
    hm.matchThis(SearchImg, std::stod(argv[3]));
    tt.stop();

    std::cout << "Time used: " << tt.duration() << std::endl;
    waitKey(0);
	return 0;
}



