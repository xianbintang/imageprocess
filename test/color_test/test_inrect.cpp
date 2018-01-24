//
// Created by xianb on 2018/1/11.
//

#include <opencv2/core/types.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv/cv.hpp>
#include <cmath>
#include <iostream>
#include <basetsd.h>

typedef struct _LINE_STRUCT_ {
    float k;
    float b;
} line_struct;

typedef struct _CIRCLE_STRUCT_ {
    cv::Point center;
    UINT16 radius;
} circle_struct;


line_struct get_line_equation(cv::Point p1, cv::Point p2) {
    line_struct line;
    line.k = 1.0f * (p1.y - p2.y) / (p1.x - p2.x);
    line.b = p1.y - p1.x * line.k;
    return line;
}

int cmpint( const void *a, const void *b)
{
    return *(int*)a-*(int*)b;
}
// 高度的范围是从点中y最小到y最大
void get_intersections(line_struct lines[], unsigned short intersection[][2], cv::Point minPoint, cv::Point maxPoint)
{
    // 首先判断四条直线中是否有垂直的x轴的，如果有垂直x轴的就直接使用minPoint和maxPoint来计算即可
    int isvertical = false;
    for (int i = 0; i < 4; ++i) {
        if(std::isinf(lines[i].k)) {
            isvertical = true;
        }
    }
    if(isvertical) {
        for (int i = minPoint.y; i < maxPoint.y; ++i) {
            intersection[i][0] = minPoint.x;
            intersection[i][1] = maxPoint.x;
        }
    } else {
        for (int i = minPoint.y; i < maxPoint.y; ++i) {
            // 获取在该处四条直线的交点，之保留在外界矩形内的交点
            int x[4];
            x[0] = 1.0f * (i - lines[0].b) / lines[0].k;
            x[1] = 1.0f * (i - lines[1].b) / lines[1].k;
            x[2] = 1.0f * (i - lines[2].b) / lines[2].k;
            x[3] = 1.0f * (i - lines[3].b) / lines[3].k;
            qsort(x, 4, sizeof(x[0]), cmpint);
            intersection[i][0] = x[1];
            intersection[i][1] = x[2];
        }
    }

}

void get_intersections_circle(circle_struct circle, unsigned short intersection[][2], cv::Point minPoint, cv::Point maxPoint)
{
    UINT16 radius = circle.radius;
    UINT16 cenx = circle.center.x;
    UINT16 ceny = circle.center.y;
    for (int i = minPoint.y; i < maxPoint.y; ++i) {
        int x = sqrt(radius * radius - (ceny - i) * (ceny - i));
        intersection[i][0] = cenx - x;
        intersection[i][1] = cenx + x;
    }
}

int main(int argc, char **argv)
{
    unsigned short intersection[480][2];
    memset(intersection, 0, sizeof(unsigned short) * 480 * 2);
#if 0
    cv::Point p0(150,182);
    cv::Point p1(134,255);
    cv::Point p2(393,313);
    cv::Point p3(407,238);

    cv::Point minPoint(134,182), maxPoint(407,313);
#endif

#if 1
    cv::Point p0(171,106);
    cv::Point p1(150,354);
    cv::Point p2(468,373);
    cv::Point p3(489,134);

    cv::Point maxPoint(150+339,106+266), minPoint(150,106);
#endif

#if 0
    cv::Point p0(10,10);
    cv::Point p1(10,100);
    cv::Point p2(100,100);
    cv::Point p3(100,10);

    cv::Point minPoint(10,10), maxPoint(100,100);
#endif

//    line_struct lines[4];
//    lines[0] = get_line_equation(p0, p1);
//    lines[1] = get_line_equation(p1, p2);
//    lines[2] = get_line_equation(p2, p3);
//    lines[3] = get_line_equation(p3, p0);
//
//    get_intersections(lines, intersection, minPoint, maxPoint);
    circle_struct circleStruct;
    circleStruct.radius = 100;
    circleStruct.center.x = 200;
    circleStruct.center.y = 200;
    maxPoint.x = 300;
    maxPoint.y = 300;

    minPoint.x = 100;
    minPoint.y = 100;
    get_intersections_circle(circleStruct, intersection, minPoint, maxPoint);
    auto img = cv::imread(argv[1], 0);
    for (int i = 0; i < 480; ++i) {
        int startx = intersection[i][0];
        int endx = intersection[i][1];
        int y = i;
        if(startx == 0 && endx == 0) {
            continue;
        }
        std::cout << startx << " " << endx << std::endl;
        for (int j = startx; j < endx; ++j) {
            cv::circle(img, cv::Point(j, i), 1, cv::Scalar(255,255,255), 1);
        }
    }
    cv::imshow("haha", img);
    cv::waitKey(0);
    return 0;
}
