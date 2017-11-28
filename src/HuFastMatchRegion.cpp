//
// Created by xianb on 2017/11/22.
//

#include <iostream>
#include <iomanip>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include "opencv/cv.h"
#include "opencv/highgui.h"
#include "TimeTracker.h"

using namespace cv;
using namespace std;

void invariantMomentHu(const Mat src, vector<double>& result)
{
    //momen dua dimensi
    double m00=0, m10=0, m01=0;
    int x,y;

    for( int i = 0; i < src.rows; ++i){
        for( int j = 0; j < src.cols; ++j ){
            x = j + 1;
            y = i + 1;
            m00 += src.at<uchar>(i,j);
            m10 += x * src.at<uchar>(i,j);
            m01 += y * src.at<uchar>(i,j);
        }
    }

    //central moment lainnya & kemudian normalized central moments
    double xbar = m10 / m00;
    double ybar = m01 / m00;

    double m11 = 0, m20 = 0, m02 = 0, m30 = 0,
           m03 = 0, m12 = 0, m21 = 0;

    // 这里是针对图像的所有像素，应该改成只用轮廓点的吧
    for( int i = 0; i < src.rows; ++i){
        for( int j = 0; j < src.cols; ++j ){
            if (src.at<uchar>(i,j) >= 100) {

                x = j + 1;
                y = i + 1;
                m11 += ( x - xbar ) * ( y - ybar ) * src.at<uchar>(i,j) / pow( m00, 2 );
                m20 += pow( ( x - xbar ), 2 ) * src.at<uchar>(i,j) / pow( m00, 2 );
                m02 += pow( ( y - ybar ), 2 ) * src.at<uchar>(i,j) / pow( m00, 2 );
                m30 += pow( ( x - xbar ), 3 )  * src.at<uchar>(i,j) / pow( m00, 2.5 );
                m03 += pow( ( y - ybar ), 3 ) * src.at<uchar>(i,j) / pow( m00, 2.5 );
                m12 += ( x - xbar ) * pow( ( y - ybar ), 2 ) * src.at<uchar>(i,j) / pow( m00, 2.5 );
                m21 += pow( ( x - xbar ), 2 ) * ( y - ybar ) * src.at<uchar>(i,j) / pow( m00, 2.5 );
            }
        }
    }

    //ketujuh nilai invariant moments
    double moment;
    //1
    moment = m20 + m02;
    result.push_back(moment);
    //2
    moment = pow( m20 - m02, 2 ) + 4 * pow( m11, 2);
    result.push_back(moment);
    //3
    moment = pow( m30 - 3 * m12, 2 ) + pow( 3 * m21 - m03, 2);
    result.push_back(moment);
    //4
    moment = pow( m30 + m12, 2 ) + pow( m21 + m03, 2);
    result.push_back(moment);
    //5
    moment = ( m30 - 3 * m12 ) * ( m30 + m12 ) * ( pow( m30 + m12, 2 ) - 3* pow( m21 + m03, 2) )
             + ( 3 * m21 - m03 ) * ( m21 + m03 ) * ( 3 * pow( m30 + m12, 2 ) - pow( m21 + m03, 2) );
    result.push_back(moment);
    //6
    moment = ( m20 - m02 )  * ( pow( m30 + m12, 2 ) - pow( m21 + m03, 2) )
             + 4 * m11 * ( m30 + m12 ) * ( m21 + m03 );
    result.push_back(moment);
	/*opencv hu implementing*/
	{
	double t0 = m30 + m12;
    double t1 = m21 + m03;
    double q0 = t0 * t0, q1 = t1 * t1;
    //double n4 = 4 * m.nu11;
    //double s = m20 + m02;
    //double d = m20 - m02;
    t0 *= q0 - 3 * q1;
    t1 *= 3 * q0 - q1;

    q0 = m30 - 3 * m12;
    q1 = 3 * m21 - m03;

    moment = q1 * t0 - q0 * t1;
	}
    //7
    //moment = ( 3 * m21 - m03 )  * ( m21 + m30 )  * ( pow( m30 + m12, 2 ) - 3* pow( m21 + m03, 2) )
     //        - ( m30 + 3 * m12 )  * ( m21 + m03 )  * ( 3 * pow( m30 + m12, 2 ) - pow( m21 + m03, 2) );
    result.push_back(moment);
}

void printMoments(char label[] ,const vector<double> huMoments)
{
    cout<<label<<endl<<setprecision(50);;
    for(unsigned int i = 0; i < huMoments.size(); ++i){
        cout<<i+1<<" = "<<huMoments[i]<<endl;
    }
}

double getScore(const vector<double> &Sa, const vector<double> &Ta) {
    //  计算相似度2
    double dbR2 =0; //相似度
    double temp2 =0;
    double temp3 =0;
    {for(int i=0;i<7;i++)
        {
            temp2 += fabs(Sa[i]-Ta[i]);
            temp3 += fabs(Sa[i]+Ta[i]);
        }}
    return dbR2 =1- (temp2*1.0)/(temp3);
}

int main(int argc, char **argv)
{

    vector<double> huMoments;
    TimeTracker tt;
    Mat g1 = imread(argv[1], 0);
    Mat result;
    Canny(g1, result, 80, 220);
    tt.start();
    invariantMomentHu(result, huMoments);
    tt.stop();
    cout << "time: " << tt.duration() << endl;
//    printMoments("gambar 1 : " ,huMoments);

    namedWindow("gambar1");
    imshow("gambar1", result);

    vector<double> huMoments2;
    Mat g2 = imread(argv[2], 0);
    Canny(g2, result, 80, 220);
    tt.start();
    invariantMomentHu(result, huMoments2);
    tt.stop();
    cout << "time: " << tt.duration() << endl;
//    printMoments("gambar 2 : " ,huMoments2);
	cout << "Score: " << getScore(huMoments, huMoments2) << endl;

    namedWindow("gambar2");
    imshow("gambar2", result);

    waitKey(0);
	return 0;
}
