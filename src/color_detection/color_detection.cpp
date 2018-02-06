#include <iostream>
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

using namespace cv;
using namespace std;

int main( int argc, char** argv )
{

    namedWindow("Control", CV_WINDOW_AUTOSIZE); //create a window called "Control"

    int iLowH = 100;
    int iHighH = 140;

    int iLowS = 90;
    int iHighS = 255;

    int iLowV = 90;
    int iHighV = 255;

    //Create trackbars in "Control" window
    cvCreateTrackbar("LowH", "Control", &iLowH, 179); //Hue (0 - 179)
    cvCreateTrackbar("HighH", "Control", &iHighH, 179);

    cvCreateTrackbar("LowS", "Control", &iLowS, 255); //Saturation (0 - 255)
    cvCreateTrackbar("HighS", "Control", &iHighS, 255);

    cvCreateTrackbar("LowV", "Control", &iLowV, 255); //Value (0 - 255)
    cvCreateTrackbar("HighV", "Control", &iHighV, 255);

#if 0
    Mat imgOriginal = imread(argv[1], 1);
    Mat imgHSV, imageYUV;
    vector<Mat> hsvSplit;

    cvtColor(imgOriginal, imageYUV, CV_BGR2YUV);
    cvtColor(imageYUV, imgHSV, CV_YUV); //Convert the captured frame from BGR to HSV
    cv::Mat mv[3];
    split(imgHSV, mv);
#endif

    std::string filename(argv[1]);
    unsigned char *buf = nullptr;
    if(filename.substr(filename.size() - 3, 3) == "yuv") {
        FILE *yuv_file = fopen(filename.c_str(), "rb+");

        buf = new unsigned char[640 * 480 + 320 * 240 * 2];
        fread(buf, 640* 480 + 320 * 240 * 2, 1, yuv_file);
    }

    Mat mYUV(480+ 480/2, 640, CV_8UC1, (void*) buf);
    Mat mRGB(480, 640, CV_8UC3);
    cvtColor(mYUV, mRGB, CV_YUV420p2BGR, 3);

    //因为我们读取的是彩色图，直方图均衡化需要在HSV空间做

    imshow("rgb", mRGB); //show the thresholded image
    cv::waitKey(-1);
#if 0

    while (true)
    {
        Mat imgOriginal = imread(argv[1], 1);


        Mat imgHSV;
        vector<Mat> hsvSplit;
        cvtColor(imgOriginal, imgHSV, COLOR_RGB2HSV); //Convert the captured frame from BGR to HSV

        //因为我们读取的是彩色图，直方图均衡化需要在HSV空间做

        Mat imgThresholded;

        inRange(imgHSV, Scalar(iLowH, iLowS, iLowV), Scalar(iHighH, iHighS, iHighV), imgThresholded); //Threshold the image
        std::cout << "iLowH: " << iLowH << " iLowS: " << iLowS << " iLowV: " << iLowV << std::endl;
        std::cout << "iHighH: " << iHighH<< " iHighS: " << iHighS<< " iHighV: " << iHighV << std::endl;

        //开操作 (去除一些噪点)
//        Mat element = getStructuringElement(MORPH_RECT, Size(5, 5));
//        morphologyEx(imgThresholded, imgThresholded, MORPH_OPEN, element);

        //闭操作 (连接一些连通域)
//        morphologyEx(imgThresholded, imgThresholded, MORPH_CLOSE, element);

        imshow("Thresholded Image", imgThresholded); //show the thresholded image
        imshow("Original", imgOriginal); //show the original image

        char key = (char) waitKey(300);
        if(key == 27)
            break;
    }
#endif

    return 0;

}
