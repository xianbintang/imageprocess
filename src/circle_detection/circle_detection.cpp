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

bool is_circle(Point p, int radius, Mat mag, Mat dist, Mat dx, Mat dy);
void sobel(Mat img, Mat &sdx, Mat &sdy, Mat &mag, Mat &dist)
{
    short acc_dx = 0, acc_dy = 0;         //accumulators
//    float k1 [] = {-1,-2,-1,0,0,0,1,2,1}; //{-2,-4,-2,0,0,0,2,4,2};//{-1,-2,-1,0,0,0,1,2,1};    //sobel kernal dx
//    float k2 [] = {-1,0,1,-2,0,2,-1,0,1};//{-2,0,2,-4,0,4,-2,0,2};//{-1,0,1,-2,0,2,-1,0,1};    //sobel kernal dy

    cv::Sobel(img, sdx, CV_16S, 1, 0, 3);
    cv::Sobel(img, sdy, CV_16S, 0, 1, 3);

    for(int i=0; i<img.rows; i++) {
        for(int j=0; j<img.cols; j++) {
            acc_dx = (short)sdx.at<short>(i, j);
            acc_dy = (short)sdy.at<short>(i, j);
            mag.at<float>(i,j) = (sqrtf(acc_dy*acc_dy + acc_dx*acc_dx)) > 100 ? 255 : 0;
            dist.at<float>(i,j) = atan2f(acc_dy, acc_dx);
            // printf("dist : %f \n", dist.at<float>(i,j) / 3.14159265f * 180 );
        }
    }
}

void inc_if_inside(double *** H, int x, int y, int height, int width, int r )
{
    if (x>0 && x<width && y> 0 && y<height)
        H[y][x][r]++;
}


void hough(Mat &img_data, Mat &dist, Mat &sdx, Mat &sdy, double threshold, int minRadius, int maxRadius, double distance, Mat &h_acc, Mat &coins)
{
    int radiusRange = maxRadius - minRadius;
    int HEIGHT = img_data.rows;
    int WIDTH = img_data.cols;
    int DEPTH = maxRadius;

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
                for (int r=minRadius; r< maxRadius; r++)
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
            for(int r = minRadius; r < maxRadius; r++) {
                h_acc.at<uchar>(y0,x0) +=  (uchar)H[y0][x0][r];// > 1 ? 255 : 0;
                // printf("h : %d", H[y0][x0][r]);
            }
        }
    }


    Point3f bestCircles[200];
    int number_of_best_cirles = 0;

    //compute optimal circles
    int count = 0;
    for(int y0 = 0; y0 < HEIGHT; y0++) {
        for(int x0 = 0; x0 < WIDTH; x0++) {
            for(int r = minRadius; r < maxRadius; r++) {
                if(H[y0][x0][r] > threshold){
//                    std::cout << H[y0][x0][r]<< std::endl;
                    count++;
                    Point3f circle(x0, y0, r);
                    int i;
                    /*
                     * 检查到一个满足阈值的圆后先与原来的圆判断，是不是在原来的圆附近
                     * 如果在原来的bestcircles附近，那么就更新领域内最大值。
                     * 如果不在附近则退出循环
                     * */
                    for(i = 0; i < number_of_best_cirles; i++) {
                        int xCoord = bestCircles[i].x;
                        int yCoord = bestCircles[i].y;
                        int radius = bestCircles[i].z;
                        /* 找到领域内得到分数最高的圆 */
                        if(abs(xCoord - x0) < distance && abs(yCoord - y0) < distance) {
                            /* 如果在领域附近而且分更高，则更新之 */
                            if(H[y0][x0][r] > H[yCoord][xCoord][radius]) {
                                bestCircles[i] = circle;
                            }
                            /* 如果在领域附近但是分低，则不管，由于i此时不等于bestCircles，所以不会将其加入为新圆 */
                            break;
                        }
                    }
                    /* 如果新的这个圆不是在原来的圆附近，那么将这个新的圆加入vector中 */
                    if(i == number_of_best_cirles){
                        bestCircles[i] = circle;
                        number_of_best_cirles++;
                    }
                }
            }
        }
    }
//    std::cout << "count: " << count << std::endl;

    /* 此时能够保证bestCircles中的是所有在各自领域内的最高分的圆。如何得到 */
    for(int i = 0; i < number_of_best_cirles; i++) {
//    for(int i = 0; i < 4; i++) {
        int lineThickness = 2;
        int lineType = 10;
        int shift = 0;
        int xCoord = bestCircles[i].x;
        int yCoord = bestCircles[i].y;
        int radius = bestCircles[i].z;
        Point2f center(xCoord, yCoord);
//        std::cout << H[yCoord][xCoord][radius] << " radius: " << radius << std::endl;
        if (is_circle(Point(xCoord, yCoord), radius, img_data, dist, sdx, sdy)) {
            circle(coins, center, radius - 1, Scalar(0, 0, 255), lineThickness, lineType, shift);
        }
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

    dx.create(img_grey.rows, img_grey.cols, CV_16SC1);
    dy.create(img_grey.rows, img_grey.cols, CV_16SC1);
    mag.create(img_grey.rows, img_grey.cols, CV_32FC1);
    dist.create(img_grey.rows, img_grey.cols, CV_32FC1);

    sobel(img_grey, dx, dy, mag, dist);

    //normalize arrays with max and min values of 255 and 0
    normalize(dx, dx_out, 0, 255, NORM_MINMAX, -1, Mat());
    normalize(dy, dy_out, 0, 255, NORM_MINMAX, -1, Mat());
    normalize(dist, dis_out, 0, 255, NORM_MINMAX, -1, Mat());

    // double H_h = sqrt(2.0) * (double) (mag.rows>mag.cols ? mag.rows : mag.cols); //-r -> +r
    // double H_w = 180;

    h_acc.create(mag.rows, mag.cols, CV_8UC1);
    // threshold(h_acc,h_out,0,255,THRESH_TOZERO);
    // threshold(h_out,h_acc,0,255,THRESH_TOZERO);

    // medianBlur(mag,mag,1);
    //mag,dist,thresh,minRad,maxRad,Dist-circles,output_Hspace, final_result
    TimeTracker tt;
    tt.start();
    hough(mag, dist, dx, dy, 5, 10, 150, 20, h_acc, image);
    hough(mag, dist, dx, dy, 5, 3, 35, 20, h_acc, image);
//    hough(mag, dist, 10, 20, 28, 20, h_acc, image);
    tt.stop();
    std::cout << "time: " << tt.duration() << std::endl;

    // normalize(h_acc, h_out, 0, 255, NORM_MINMAX, -1, Mat());
    // threshold(h_acc,h_out, 200,255,THRESH_TOZERO);

    //save images
    Mat dis_show;
    Mat dx_show, dy_show;
    convertScaleAbs(dis_out, dis_show);
    convertScaleAbs(dx_out, dx_show);
    convertScaleAbs(dy_out, dy_show);
    imshow( "dx.jpg", dx_show);
    imshow( "dy.jpg", dy_show);
    imshow( "mag.jpg", mag );
    imshow( "dist.jpg", dis_show);
    imshow( "h_space.jpg", h_acc);
//
    imshow("result.png",image);
    cvWaitKey(-1);
    return 0;
}

/* 八个方向上的圆周是不是都是边缘点，他们的梯度是不是都是半径的径向 */
bool is_circle(Point p, int radius, Mat mag,Mat dist,  Mat dx, Mat dy)
{
//    double angles[8] = {3.14159265f * 0 / 180, 3.14159265f * 45 / 180, 3.14159265f * 90 / 180,
//                        3.14159265f * 135 / 180, 3.14159265f * 180 / 180, 3.14159265f * -45 / 180,
//                        3.14159265f * 270 / -90, 3.14159265f * -135/ 180};
    int count = 0;
//    std::cout << "\nis it a circle? " << std::endl;
    for (int i = -180; i < 180; i += 1) {
        double angle = 3.14159265f * i / 180;
        int x0 = (int)round(p.x + radius * cos(angle));
        int y0 = (int)round(p.y + radius * sin(angle));
//        short sdxv, sdyv;
//        sdxv = dx.at<short>(y0, x0);
//        sdyv = dy.at<short>(y0, x0);
//        double radial_direction = atan2f(sdyv, sdxv);
        double radial_direction = dist.at<float>(y0, x0);

//        printf("%lf %lf \n", angle, radial_direction);
        /* 若方向相差角度为30度以内：30 * PI / 180 == 0.52 */
        if (fabs(angle - radial_direction) < 0.52) {
            count++;
        }
    }
//    std::cout << "count: " << count << std::endl;
    /* 至少拟合一半的点数 */
    return count > 180;
}

