//
// Created by xianb on 2017/7/18.
//

#include "TimeTracker.h"
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

        for (int j = 0; j < WIDTH; ++j) {
            H[i][j] = new double[DEPTH];
            for (int k = 0; k < DEPTH; ++k) {
                H[i][j][k] = 0;
            }
        }
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
    for(int y0 = 0; y0 < HEIGHT; y0++) {
        for(int x0 = 0; x0 < WIDTH; x0++) {
            if (h_acc.at<float>(y0,x0) > 130)
                printf("x: %d y: %d h : %lf\n", x0, y0,  h_acc.at<float>(y0,x0));
        }
    }

    std::vector<Point3f> bestCircles;

    //compute optimal circles
    int count = 0;
    for(int y0 = 0; y0 < HEIGHT; y0++) {
        for(int x0 = 0; x0 < WIDTH; x0++) {
            for(int r = minRadius; r < radiusRange; r++) {
                if(H[y0][x0][r] > threshold){
//                    std::cout << H[y0][x0][r]<< std::endl;
                    count++;
                    Point3f circle(x0, y0, r);
                    int i;
                    for(i = 0; i < bestCircles.size(); i++) {
                        int xCoord = bestCircles[i].x;
                        int yCoord = bestCircles[i].y;
                        int radius = bestCircles[i].z;
                        /* 找到领域内得到分数最高的圆 */
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
    std::cout << "count: " << count << std::endl;

    for(int i = 0; i < bestCircles.size(); i++) {
//    for(int i = 0; i < 4; i++) {
        int lineThickness = 1;
        int lineType = 10;
        int shift = 0;
        int xCoord = bestCircles[i].x;
        int yCoord = bestCircles[i].y;
        int radius = bestCircles[i].z;
        Point2f center(xCoord, yCoord);
//        std::cout << H[yCoord][xCoord][radius] << std::endl;
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
    TimeTracker tt;
    tt.start();
    hough(mag, dist, 10, 10, 100, 0, h_acc, image);
    tt.stop();
    std::cout << "time: " << tt.duration() << std::endl;

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


