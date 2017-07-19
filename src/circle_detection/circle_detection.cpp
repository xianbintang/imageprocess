//
// Created by xianb on 2017/7/18.
//


# if 0
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
                     1,   // 累加器的分辨率(图像尺寸/2)
                     10,  // 两个圆之间的最小距离
                     100, // Canny中的高阈值
                     60, // 最小投票数
                     20, 60); // 有效半径的最小和最大值

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
                   2); // 圆形的厚度
        ++itc;
    }

    cv::namedWindow("Detected Circles");
    cv::imshow("Detected Circles",image);
    cvWaitKey(-1);
    return 0;
}
#endif

#if 0
/*
Author: Chun-Wei Chiang
Object: canny detection of the image and find the circle shape via hough transform
*/

#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sys/time.h>
#include <algorithm>
using namespace std;

using namespace cv;

Mat src;
Mat detected_edges;
Mat dst;

int edgeThresh = 1;
/*
If an edge pixel’s gradient value is higher than the high threshold value, it is marked as a strong edge pixel.
If an edge pixel’s gradient value is smaller than the high threshold value and larger than the low threshold value, it is marked as a weak edge pixel.
If an edge pixel's value is smaller than the low threshold value, it will be suppressed.
Normailly the high Threshold is 2-3 times to low threshold.
*/
int lowThreshold = 70;
int ratio = 3;
/*
size of sobel filter size
*/
int kernel_size = 3;
const int HoughTreshold = 75;

const char* window_name = "original image";
const char* window_name2 = "canny edge detect";
const char* window_name3 = "hough detect";
const char* source_address = "/home/chunwei/computerVision/ShapeDetection/3circle.jpg";
static unsigned int*** Accumulator; //(a,b,r)

void edgefinder(int, void*){

    // Reduce noise with a Gaussian kernel 5*5 with sigma 1.3
    GaussianBlur( src, detected_edges, Size(5,5), 1.3, 1.3);
    // Canny detector
    Canny( detected_edges, detected_edges, lowThreshold, lowThreshold*ratio, kernel_size );
    // Show the edge image
    imshow( window_name2, detected_edges );
    imwrite ("edge.jpg ", detected_edges);
    cout << detected_edges.cols<< "   " << detected_edges.rows << endl;

}

void houghTrasform(unsigned char* hough_src){

    cv::cvtColor(src,dst,CV_GRAY2RGB);

    int w = detected_edges.cols;
    int h = detected_edges.rows;
    int radius = min(w,h)/2;

    // //initialize the matrix
    Accumulator = new unsigned int **[h];
    for(int i = 0 ; i < h; i++){
        Accumulator[i] = new unsigned int *[w];
        for(int j = 0 ; j < w ; j++){
            Accumulator[i][j] = new unsigned int [radius];
            for(int k = 0 ; k < radius; k++){
                Accumulator[i][j][k] = 0;
            }
        }
    }

    int r;

    cout << "start hough transform" <<endl;
    clock_t start_time = clock();

    //traditional Hough transform to detect circle.
    for(int i = 0 ; i < h; i++){
        for(int j = 0; j < w ; j++){

            // cout << hough_src [ (i*w) + j] << " " ;

            if(hough_src [ (i*w) + j]>250){ //edge
                for(int a = 0; a < h ; a++){
                    for (int b = 0 ; b < w ; b++){
                        r = (int) sqrt( (i-a)*(i-a) + (j-b)*(j-b) );
                        if( r >0 && r < radius){
                            Accumulator[a][b][r]++;
                        }
                    }
                }
            }
        }
    }

    for(int i = 0; i < h; i++){
        for(int j = 0; j < w; j++){
            for(int k = 0; k < radius; k++){
                if(Accumulator[i][j][k] > HoughTreshold){
                    //To reduce the duplicate cirlce in (9*9)
                    int max = Accumulator[i][j][k];
                    for(int li = -4 ; li <= 4; li++ ){
                        for(int lj = -4 ; lj <= 4; lj++){
                            for(int lk = -4 ; lk <= 4; lk++){

                                if(i+li >= 0 && i+li < h && j+lj >=0 && j+lj < w && k+lk > 0 && k + lk < radius){

                                    if(Accumulator[i+li][j+lj][k+lk] > max ){
                                        max = Accumulator[i+li][j+lj][k+lk];
                                        li = lj = lk = 5;
                                    }
                                }
                            }
                        }
                    }

                    if(max > Accumulator[i][j][k]) continue;

                    //To avoid two points have the same accumulator, and they would draw two circles at the near place.
                    Accumulator[i][j][k]++;
                    cout << "(" << i << "," << j << "," << k << ")" << endl;
                    cv::circle(dst, Point(j,i), k, Scalar(0,0,255), 1);
                }
            }
        }
    }

    clock_t end_time = clock();
    cout << "It takes " << (end_time - start_time) * 1000 / CLOCKS_PER_SEC << " ms to do the hough transform." << endl;
    cout << "end draw circle" << endl;

    imshow( window_name3, dst);
    imwrite ("output1.jpg ", dst);
}

int main(int argc, char **argv){


    // Load an image in gray scale
    src = imread(argv[1],0 );

    if( !src.data ){
        cout << "cannot find the image"<<endl;
        return -1;
    }

    // Create a matrix of the same type and size as src (for dst)

    // Create a window
    namedWindow( window_name, CV_WINDOW_AUTOSIZE );
    namedWindow( window_name2, CV_WINDOW_AUTOSIZE );
    namedWindow( window_name3, CV_WINDOW_AUTOSIZE );

    // Show the original image
    imshow( window_name, src );
    //find the edge
    edgefinder(0, 0);
    houghTrasform(detected_edges.data);
    // Wait until user exit program by pressing a key
    waitKey(0);

    return 0;
}

#endif

#if 1
#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <math.h>

using namespace cv;

void sobel(Mat img, Mat &dx, Mat &dy, Mat &mag, Mat &dist)
{
    float acc_dx = 0, acc_dy = 0;         //accumulators
    float k1 [] = {-1,-2,-1,0,0,0,1,2,1}; //{-2,-4,-2,0,0,0,2,4,2};//{-1,-2,-1,0,0,0,1,2,1};    //sobel kernal dx
    float k2 [] = {-1,0,1,-2,0,2,-1,0,1};//{-2,0,2,-4,0,4,-2,0,2};//{-1,0,1,-2,0,2,-1,0,1};    //sobel kernal dy

    for(int i=0; i<img.rows; i++) {
        for(int j=0; j<img.cols; j++) {
            acc_dx = acc_dy = 0;

            //apply kernel/mask
            for (int nn=-1; nn<2; nn++) {
                for (int mm = -1; mm < 2; mm++) {
                    if (i + nn > 0 && i + nn < img.rows && j + mm > 0 && j + mm < img.cols) {
                        acc_dx += (float)img.at<uchar>(i+nn,j+mm) * k1[((mm+1)*3) + nn + 1];
                        acc_dy += (float)img.at<uchar>(i+nn,j+mm) * k2[((mm+1)*3) + nn + 1];
                    }
                }
            }
            //write final values
            dx.at<float>(i,j) = acc_dx;
            dy.at<float>(i,j) = acc_dy;
            mag.at<float>(i,j) = (sqrtf(acc_dy*acc_dy + acc_dx*acc_dx)) > 100 ? 255 : 0;
            dist.at<float>(i,j) = atan2f(acc_dy,acc_dx);
            // printf("dist : %f \n", dist.at<float>(i,j) / 3.14159265f * 180 );
        }
    }
}



void inc_if_inside(double *** H, int x, int y, int height, int width, int r )
{
    if (x>0 && x<width && y> 0 && y<height)
        H[y][x][r]++;
}


void hough(Mat &img_data, Mat &dist, double threshold, int minRadius, int maxRadius, double distance, Mat &h_acc, Mat &coins){
    int radiusRange = maxRadius - minRadius;
    int HEIGHT = img_data.rows;
    int WIDTH = img_data.cols;
    int DEPTH = radiusRange;

    double ***H;

    // Allocate memory
    H = new double**[HEIGHT];
    for (int i = 0; i < HEIGHT; ++i) {
        H[i] = new double*[WIDTH];

        for (int j = 0; j < WIDTH; ++j)
            H[i][j] = new double[DEPTH];
    }

    for(int y=0;y<img_data.rows;y++)
    {
        for(int x=0;x<img_data.cols;x++)
        {
            // printf("data point : %f\n", img_data.at<float>(y,x));
            if( (float) img_data.at<float>(y,x) > 250.0 )  //threshold image
            {
                for (int r=minRadius; r<radiusRange; r++)
                {

                    int x0 = round(x + r * cos(dist.at<float>(y,x)) );
                    int x1 = round(x - r * cos(dist.at<float>(y,x)) );
                    int y0 = round(y + r * sin(dist.at<float>(y,x)) );
                    int y1 = round(y - r * sin(dist.at<float>(y,x)) );


                    inc_if_inside(H,x0,y0,HEIGHT, WIDTH, r);
                    // inc_if_inside(H,x0,y1,HEIGHT, WIDTH, r);
                    // inc_if_inside(H,x1,y0,HEIGHT, WIDTH, r);
                    inc_if_inside(H,x1,y1,HEIGHT, WIDTH, r);
                }
            }
        }
    }

    //create 2D image by summing values of the radius dimension
    for(int y0 = 0; y0 < HEIGHT; y0++) {
        for(int x0 = 0; x0 < WIDTH; x0++) {
            for(int r = minRadius; r < radiusRange; r++) {
                h_acc.at<float>(y0,x0) +=  H[y0][x0][r];// > 1 ? 255 : 0;
                // printf("h : %d", H[y0][x0][r]);
            }
        }
    }

    std::vector<Point3f> bestCircles;

    //compute optimal circles
    for(int y0 = 0; y0 < HEIGHT; y0++) {
        for(int x0 = 0; x0 < WIDTH; x0++) {
            for(int r = minRadius; r < radiusRange; r++) {
                if(H[y0][x0][r] > threshold){
                    Point3f circle(x0, y0, r);
                    int i;
                    for(i = 0; i < bestCircles.size(); i++) {
                        int xCoord = bestCircles[i].x;
                        int yCoord = bestCircles[i].y;
                        int radius = bestCircles[i].z;
                        if(abs(xCoord - x0) < distance && abs(yCoord - y0) < distance) {
                            if(H[y0][x0][r] > H[yCoord][xCoord][radius]) {
                                bestCircles.erase(bestCircles.begin()+i);
                                bestCircles.insert(bestCircles.begin(), circle);
                            }
                            break;
                        }
                    }
                    if(i == bestCircles.size()){
                        bestCircles.insert(bestCircles.begin(), circle);
                    }
                }
            }
        }
    }

//    for(int i = 0; i < bestCircles.size(); i++) {
    for(int i = 0; i < 4; i++) {
        int lineThickness = 4;
        int lineType = 10;
        int shift = 0;
        int xCoord = bestCircles[i].x;
        int yCoord = bestCircles[i].y;
        int radius = bestCircles[i].z;
        Point2f center(xCoord, yCoord);
        circle(coins, center, radius-1, Scalar(255,0,0), lineThickness, lineType, shift);
    }
}

int main( int argc, char** argv )
{

    char* imageName = argv[1];

    Mat image, img_grey, img_grey1;     //input mat
    Mat dx,dy,mag,dist;
    Mat dx_out, dy_out, dis_out;  //final output mat
    Mat h_acc, h_out;       //hough space matricies

    image = imread( imageName, 1);
    cvtColor( image, img_grey, COLOR_BGR2GRAY );

    dx.create(img_grey.rows, img_grey.cols, CV_32FC1);
    dy.create(img_grey.rows, img_grey.cols, CV_32FC1);
    mag.create(img_grey.rows, img_grey.cols, CV_32FC1);
    dist.create(img_grey.rows, img_grey.cols, CV_32FC1);

    sobel(img_grey, dx, dy, mag, dist);

    //normalize arrays with max and min values of 255 and 0
    normalize(dx, dx_out, 0, 255, NORM_MINMAX, -1, Mat());
    normalize(dy, dy_out, 0, 255, NORM_MINMAX, -1, Mat());
    normalize(dist, dis_out, 0, 255, NORM_MINMAX, -1, Mat());

    // double H_h = sqrt(2.0) * (double) (mag.rows>mag.cols ? mag.rows : mag.cols); //-r -> +r
    // double H_w = 180;

    h_acc.create(mag.rows, mag.cols, CV_32FC1);
    // threshold(h_acc,h_out,0,255,THRESH_TOZERO);
    // threshold(h_out,h_acc,0,255,THRESH_TOZERO);

    // medianBlur(mag,mag,1);
    //mag,dist,thresh,minRad,maxRad,Dist-circles,output_Hspace, final_result
    hough(mag, dist, 10, 15, 100, 45, h_acc, image);

    // normalize(h_acc, h_out, 0, 255, NORM_MINMAX, -1, Mat());
    // threshold(h_acc,h_out, 200,255,THRESH_TOZERO);

    //save images
    imshow( "dx.jpg", dx_out );
    imshow( "dy.jpg", dy_out );
    imshow( "mag.jpg", mag );
    imshow( "dist.jpg", dis_out );
    imshow( "h_space.jpg", h_acc);
//
    imshow("result.png",image);
    cvWaitKey(-1);
    return 0;
}
#endif

//to install the missing opencl module
//sudo apt-get install ocl-icd-opencl-dev

