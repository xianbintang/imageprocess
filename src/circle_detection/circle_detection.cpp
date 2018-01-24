#include "TimeTracker.h"
#include "circle_detection.h"
#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <math.h>
#include <basetsd.h>

using namespace cv;

/* 八个方向上的圆周是不是都是边缘点，他们的梯度是不是都是半径的径向 */
static bool is_circle(Point p, int radius, Circle_runtime_param &circle_runtime_param, double *score)
{
//    double angles[8] = {3.14159265f * 0 / 180, 3.14159265f * 45 / 180, 3.14159265f * 90 / 180,
//                        3.14159265f * 135 / 180, 3.14159265f * 180 / 180, 3.14159265f * -45 / 180,
//                        3.14159265f * 270 / -90, 3.14159265f * -135/ 180};
    int count = 0;
    int width = circle_runtime_param.mag.cols;
    int height = circle_runtime_param.mag.rows;
//    std::cout << "\nis it a circle? " << std::endl;
    for (int i = -180; i < 180; i += 1) {
        double angle = 3.14159265f * i / 180;
        int x0 = (int)round(p.x + radius * cos(angle));
        int y0 = (int)round(p.y + radius * sin(angle));
        short sdxv, sdyv;
        if (y0 < 0 || x0 < 0 || y0 >= height || x0 >= width)
            continue;
        sdxv = circle_runtime_param.sdx.at<short>(y0, x0);
        sdyv = circle_runtime_param.sdy.at<short>(y0, x0);
        double radial_direction = atan2f(sdyv, sdxv);
        double radial_direction1 = atan2f(-sdyv, -sdxv);

//        printf("%lf %lf \n", angle, radial_direction);
        /* 若方向相差角度为30度以内：30 * PI / 180 == 0.52, 角度如果小的话不能检测出不规整的圆 */
        if (fabs(angle - radial_direction) < 0.26 || fabs(angle - radial_direction1) < 0.26) {
            count++;
        }
    }
    *score = 1.0 * count / 360;
    /* 拟合60%的点 */
    /* 至少拟合一半的点数 */
    std::cout << "count " << count << std::endl;
    return count > 216;
}

/* 判断不同半径情况下检测出来的圆是不是已经检测过了，半径是否类似，圆心位置是否类似 */
static void remove_duplicates(std::vector<Circle> &circles)
{
    for (auto iter = circles.begin(); iter != circles.end(); ++iter) {
        if(iter->radius == 0)
            continue;
        int x0 = iter->center.x, y0 = iter->center.y, r0 = iter->radius;
        double s0 = iter->score;
        for (auto next = iter + 1; next != circles.end(); ++next) {
            int x1 = next->center.x, y1 = next->center.y, r1 = next->radius;
            double s1 = next->score;
            if ((abs(x0 - x1) < 8 && abs(y0 - y1) < 8) && abs(r0 - r1) < 8) {
                if (s0 > s1) {
                    next->center.x = 0;
                    next->center.y = 0;
                    next->radius = 0;
                    next->score = 0;
                } else {
                    iter->center.x = 0;
                    iter->center.y = 0;
                    iter->radius = 0;
                    iter->score = 0;
                }
            }
        }
    }
}

static void prepare_data(Circle_runtime_param &circle_runtime_param)
{
    short acc_dx = 0, acc_dy = 0;         //accumulators
//    float k1 [] = {-1,-2,-1,0,0,0,1,2,1}; //{-2,-4,-2,0,0,0,2,4,2};//{-1,-2,-1,0,0,0,1,2,1};    //sobel kernal dx
//    float k2 [] = {-1,0,1,-2,0,2,-1,0,1};//{-2,0,2,-4,0,4,-2,0,2};//{-1,0,1,-2,0,2,-1,0,1};    //sobel kernal dy

    cv::Sobel(circle_runtime_param.img_gray, circle_runtime_param.sdx, CV_16S, 1, 0, 3);
    cv::Sobel(circle_runtime_param.img_gray, circle_runtime_param.sdy, CV_16S, 0, 1, 3);

    int ct = 0;
    for(int i=0; i<circle_runtime_param.img_gray.rows; i++) {
        for(int j=0; j<circle_runtime_param.img_gray.cols; j++) {
            acc_dx = (short)circle_runtime_param.sdx.at<short>(i, j);
            acc_dy = (short)circle_runtime_param.sdy.at<short>(i, j);
            circle_runtime_param.mag.at<uchar>(i,j) = (sqrt(acc_dy*acc_dy + acc_dx*acc_dx)) > 220? 255 : 0;
            if(sqrt(acc_dy*acc_dy + acc_dx*acc_dx)> 220) {
                ++ct;
            }
            /*
             * TODO dist可以利用查表实现，如果acc_dy和acc_dx能够归一化到8位即可
             * */
            circle_runtime_param.ang.at<float>(i,j) = atan2f(acc_dy, acc_dx);
            // printf("dist : %f \n", dist.at<float>(i,j) / 3.14159265f * 180 );
        }
    }
    std::cout << "ct..:" << ct << std::endl;
}

inline void inc_if_inside(int *** H, int x, int y, int height, int width, int r, int candidates[][3], int *numof_candidates, int threshold, int minRadius)
{
    if (x>0 && x<width && y> 0 && y<height) {
        H[y][x][r]++;
        if (H[y][x][r] > threshold) {
            candidates[*numof_candidates][0] = y;
            candidates[*numof_candidates][1] = x;
            candidates[*numof_candidates][2] = r + minRadius;
            *numof_candidates = *numof_candidates + 1;
        }
    }
}


//void hough(Mat &img_data, Mat &dist, Mat &sdx, Mat &sdy, double threshold, int minRadius, int maxRadius, double distance, Mat &h_acc, Mat &coins, std::vector<Circle> &circles, int *num_circles, Region region)
static void hough(Circle_runtime_param &circle_runtime_param, double threshold, int minRadius, int maxRadius, double distance, std::vector<Circle> &circles, int *num_circles)
{
    Mat img_data = circle_runtime_param.mag, ang= circle_runtime_param.ang, sdx = circle_runtime_param.sdx, sdy = circle_runtime_param.sdy;

    Region region = circle_runtime_param.region;
    int radiusRange = maxRadius - minRadius;
    int regionx = region.center.x - region.radius;
    int regiony = region.center.y - region.radius;
    int region_height = region.radius * 2;
    int region_width = region.radius * 2;

    int HEIGHT = region_height;
    int WIDTH = region_width;
    int DEPTH = radiusRange;

    int ***H;

    // Allocate memory
    /* 这里没必要设置depth为最大半径，只需要设置radiusRange范围这么大的depth即可 */
    /* H的高宽设置成外接矩形的高宽 */
    H = new int**[HEIGHT];
    for (int i = 0; i < HEIGHT; ++i) {
        H[i] = new int*[WIDTH];

        for (int j = 0; j < WIDTH; ++j) {
            H[i][j] = new int[DEPTH];
            for (int k = 0; k < DEPTH; ++k) {
                H[i][j][k] = 0;
            }
        }
    }

    int candidates[1000][3];
    int numof_candidates = 0;
    memset(candidates, 0, sizeof(int) * 1000 * 3);
    TimeTracker t1;
    t1.start();
    TimeTracker t2;
    t2.start();
    int ct = 0, ct2 = 0;
    for(int y=regiony; y > 0 && y < img_data.rows && y<region_height + regiony; y++)
    {
        for(int x=regionx; x > 0 && x<img_data.cols && x < region_width + regionx; x++)
        {
            // printf("data point : %f\n", img_data.at<float>(y,x));
            if( img_data.at<uchar>(y,x) > 250 )  //threshold image
            {
                ct++;
                double theta = ang.at<float>(y,x);
                double cos_theta = cos(theta);
                double sin_theta = sin(theta);
                for (int r=minRadius; r< maxRadius; r++)
                {
                    ++ct2;
                    double r_cos_theta = r * cos_theta;
                    double r_sin_theta = r * sin_theta;
                    int x0 = (x + r_cos_theta);
                    int x1 = (x - r_cos_theta);
                    int y0 = (y + r_sin_theta);
                    int y1 = (y - r_sin_theta);

//                    printf("x0: %d, x1: %d, y0: %d, y1: %d\n", x0, x1, y0, y1);
                    inc_if_inside(H, x0 - regionx, y0 - regiony, HEIGHT, WIDTH, r - minRadius, candidates, &numof_candidates, threshold, minRadius);
                    // inc_if_inside(H,x0,y1,HEIGHT, WIDTH, r);
                    // inc_if_inside(H,x1,y0,HEIGHT, WIDTH, r);
                    inc_if_inside(H, x1 - regionx, y1 - regiony, HEIGHT, WIDTH, r - minRadius, candidates, &numof_candidates, threshold, minRadius);
                }
            }
        }
    }
    std::cout << "ct: " << ct << "ct2: " << ct2 << std::endl;
//    std::out << " numof candidates: " <<  numof_candidates << std::endl;
    t2.stop();
    std::cout << "duration t2: " << t2.duration() << std::endl;
#if 0
    //create 2D image by summing values of the radius dimension
    for(int y0 = 0; y0 < HEIGHT; y0++) {
        for(int x0 = 0; x0 < WIDTH; x0++) {
            for(int r = minRadius; r < maxRadius; r++) {
                int x = x0 + regionx;
                int y = y0 + regiony;
                if (x > 0 && x < img_data.cols && y > 0 && y < img_data.rows)
                    h_acc.at<uchar>(y, x) +=  (uchar)H[y0][x0][r - minRadius];// > 1 ? 255 : 0;
                // printf("h : %d", H[y0][x0][r]);
            }
        }
    }
#endif


    Point3f bestCircles[200];
    int number_of_best_cirles = 0;

    //compute optimal circles
    int count = 0;
    TimeTracker t3;
    t3.start();
    std::cout << "candidates: " << numof_candidates << std::endl;
//    for (int i = 0; i < numof_candidates; ++i) {
//        std::cout << "candidiates " << i << ": " << candidates[i][0] << " " << candidates[i][1] << " " << candidates[i][2] << std::endl;
//    }
    for (int k = 0; k < numof_candidates; ++k) {
        int y0 = candidates[k][0];
        int x0 = candidates[k][1];
        int r = candidates[k][2];

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
//                printf("%d %d\n", H[y0][x0][r - minRadius], H[yCoord][xCoord][radius - minRadius]);
                if(H[y0][x0][r - minRadius] > H[yCoord][xCoord][radius - minRadius]) {
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

    std::cout << "number of best circles" << number_of_best_cirles << std::endl;

    t3.stop();
    std::cout << "duration t3: " << t3.duration() << std::endl;
//    std::cout << "count: " << count << std::endl;
//    std::cout << "num of best: " << number_of_best_cirles << "circle: " << bestCircles[0] << std::endl;
    /* 此时能够保证bestCircles中的是所有在各自领域内的最高分的圆。如何得到 */
    int j = 0;
    for(int i = 0; i < number_of_best_cirles; i++) {
//        int lineThickness = 2;
//        int lineType = 10;
//        int shift = 0;
        int xCoord = bestCircles[i].x + regionx;
        int yCoord = bestCircles[i].y + regiony;
        int radius = bestCircles[i].z;
        Point2f center(xCoord, yCoord);
        double score;
//        std::cout << H[yCoord][xCoord][radius] << " radius: " << radius << std::endl;
//        if (is_circle(Point(xCoord, yCoord), radius, img_data, ang, sdx, sdy, &score)) {
        if (is_circle(Point(xCoord, yCoord), radius, circle_runtime_param, &score)) {
            Circle circle;
            circle.center = center;
            circle.radius = radius;
            circle.score = score;
//            circles[j++] = circle;
            circles.push_back(circle);
            ++j;
//            circle(coins, center, radius - 1, Scalar(0, 0, 255), lineThickness, lineType, shift);
        }
    }

    t1.stop();
    std::cout << "duration: " << t1.duration() << std::endl;
    *num_circles = j;
    /* free H */
    for (int i = 0; i < HEIGHT; ++i) {
        for (int j = 0; j < WIDTH; ++j) {
            delete H[i][j];
        }
        delete H[i];
    }
}
static cv::Mat get_y_from_yuv(const UINT8 *yuv, const UINT16 width, const UINT16 height)
{
    cv::Mat gray;
    if (!yuv) {
        std::cout << "ERROR: yuv NULL pointer" << std::endl;
    } else {
        gray.create(height, width, CV_8UC1);
        memcpy(gray.data, yuv, width * height * sizeof(unsigned char)); // 只读取y分量上的数据
    }

    return gray;
}

int circle_detection_config(const UINT8 *yuv, Circle circles1[])
{
    // 获取灰度图

    Circle_runtime_param circle_runtime_param;

    Mat dx,dy,mag,ang;
//    Mat dx_out, dy_out, dis_out;  //final output mat

    auto image = get_y_from_yuv(yuv, WIDTH, HEIGHT);

#if 1
    circle_runtime_param.mag.create(image.rows, image.cols, CV_8UC1);
    circle_runtime_param.ang.create(image.rows, image.cols, CV_32FC1);
    circle_runtime_param.img_gray = image;
#endif
    TimeTracker tt;
    tt.start();
//    GaussianBlur(img_grey, img_grey, Size(3, 3), 2, 2);
//    sobel(img_grey, dx, dy, mag, ang);
    prepare_data(circle_runtime_param);
//    GaussianBlur(img_grey, img_grey, Size(5, 5), 2, 2);
//    Canny(img_grey, edge, 30, 170, 3);
//    imshow("edge", edge);

    int width = image.cols;
    int height = image.rows;
    tt.stop();
    std::cout << "half: " << tt.duration() << std::endl;
    tt.start();
    /* 需要遍历所有可能的半径组合来检测出所有的圆，这样才能检测出同心圆 */
    /* 半径最大为1/2 * width, 最小为10像素， step为10像素 */
    int step = width / 30;
//    int step = 10;

    std::vector<Circle> total_circles;

    circle_runtime_param.region.center.x = width / 2;
    circle_runtime_param.region.center.y = height / 2;
    circle_runtime_param.region.radius = height / 2.5;


    int num_circles = 0, num_total_circles = 0;
    for (int r = 10; r < (int)(1.0 / 2 * width - 10); r += 20) {
        std::vector<Circle> circles;
        hough(circle_runtime_param, 10, r - 10, r + 10, 20, circles, &num_circles);
        for (auto iter = circles.cbegin(); iter != circles.cend(); ++iter) {
            total_circles.push_back(*iter);
            num_total_circles++;
        }
    }

    remove_duplicates(total_circles);
    tt.stop();

    // 显示效果用的，实际上处理过程在remove_duplicates就已经完成了
    for(int i = 0; i < num_total_circles; i++) {
        int lineThickness = 2;
        int lineType = 10;
        int shift = 0;
        int xCoord = total_circles[i].center.x;
        int yCoord = total_circles[i].center.y;
        int radius = total_circles[i].radius;
        Point2f center(xCoord, yCoord);
        double score = total_circles[i].score;
        if (radius != 0) {
            circle(image, center, radius - 1, Scalar(255, 255, 255), lineThickness, lineType, shift);
            std::cout << "center: " << center <<  "radius: " << radius << " count: " << "score: " << score << std::endl;
        }
    }
    imshow( "gray", image);

    cvWaitKey(-1);
    return 0;
}


