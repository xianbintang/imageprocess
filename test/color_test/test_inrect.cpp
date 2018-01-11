//
// Created by xianb on 2018/1/11.
//

#include <opencv2/core/types.hpp>

typedef struct _LINE_STRUCT_ {
    float k;
    float b;
} line_struct;

line_struct get_line_equation(cv::Point p1, cv::Point p2) {
    line_struct line;
    line.k = 1.0f * (p1.y - p2.y) / (p1.x - p2.x);
    line.b = p1.y - p1.x * line.k;
    return line;
}

// 高度的范围是从点中y最小到y最大
void get_intersections(line_struct lines[], int intersection[][2], cv::Point maxPoint, cv::Point minPoint)
{
    for (int i = minPoint.y; i < maxPoint.y; ++i) {
        // 获取在该处四条直线的交点，之保留在外界矩形内的交点
        int x0 = 1.0f * (i - lines[0].b) / lines[0].k;
        int x1 = 1.0f * (i - lines[1].b) / lines[1].k;
        int x2 = 1.0f * (i - lines[2].b) / lines[2].k;
        int x3 = 1.0f * (i - lines[3].b) / lines[3].k;
    }
}


int main(void)
{
    unsigned short intersection[480][2];
    memset(intersection, 0, sizeof(unsigned short) * 480 * 2);
    cv::Point p0(10,10);
    cv::Point p1(10,100);
    cv::Point p2(100,100);
    cv::Point p3(100,10);

    cv::Point maxPoint(100,100), minPoint(10,10);


    line_struct line1 = get_line_equation(p0, p1);
    line_struct line2 = get_line_equation(p1, p2);
    line_struct line3 = get_line_equation(p2, p3);
    line_struct line4 = get_line_equation(p3, p0);

}
