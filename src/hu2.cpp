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


//#################################################################################//
double M[7] = {0};        //HU不变矩
bool HuMoment(IplImage* img)
{
    int bmpWidth = img->width;
    int bmpHeight = img->height;
    int bmpStep = img->widthStep;
    int bmpChannels = img->nChannels;
    uchar*pBmpBuf = (uchar*)img->imageData;

    double m00=0,m11=0,m20=0,m02=0,m30=0,m03=0,m12=0,m21=0;  //中心矩
    double x0=0,y0=0;    //计算中心距时所使用的临时变量（x-x'）
    double u20=0,u02=0,u11=0,u30=0,u03=0,u12=0,u21=0;//规范化后的中心矩
    //double M[7];    //HU不变矩
    double t1=0,t2=0,t3=0,t4=0,t5=0;//临时变量，
    //double Center_x=0,Center_y=0;//重心
    int Center_x=0,Center_y=0;//重心
    int i,j;            //循环变量

    //  获得图像的区域重心(普通矩)
    double s10=0,s01=0,s00=0;  //0阶矩和1阶矩
    for(j=0;j<bmpHeight;j++)//y
    {
        for(i=0;i<bmpWidth;i++)//x
        {
            s10+=i*pBmpBuf[j*bmpStep+i];
            s01+=j*pBmpBuf[j*bmpStep+i];
            s00+=pBmpBuf[j*bmpStep+i];
        }
    }
    Center_x=(int)(s10/s00+0.5);
    Center_y=(int)(s01/s00+0.5);

    //  计算二阶、三阶矩(中心矩)
    m00=s00;
    for(j=0;j<bmpHeight;j++)
    {
        for(i=0;i<bmpWidth;i++)//x
        {
            x0=(i-Center_x);
            y0=(j-Center_y);
            m11+=x0*y0*pBmpBuf[j*bmpStep+i];
            m20+=x0*x0*pBmpBuf[j*bmpStep+i];
            m02+=y0*y0*pBmpBuf[j*bmpStep+i];
            m03+=y0*y0*y0*pBmpBuf[j*bmpStep+i];
            m30+=x0*x0*x0*pBmpBuf[j*bmpStep+i];
            m12+=x0*y0*y0*pBmpBuf[j*bmpStep+i];
            m21+=x0*x0*y0*pBmpBuf[j*bmpStep+i];
        }
    }

    //  计算规范化后的中心矩: mij/pow(m00,((i+j+2)/2)
    u20=m20/pow(m00,2);
    u02=m02/pow(m00,2);
    u11=m11/pow(m00,2);
    u30=m30/pow(m00,2.5);
    u03=m03/pow(m00,2.5);
    u12=m12/pow(m00,2.5);
    u21=m21/pow(m00,2.5);

    //  计算中间变量
    t1=(u20-u02);
    t2=(u30-3*u12);
    t3=(3*u21-u03);
    t4=(u30+u12);
    t5=(u21+u03);

    //  计算不变矩
    M[0]=u20+u02;
    M[1]=t1*t1+4*u11*u11;
    M[2]=t2*t2+t3*t3;
    M[3]=t4*t4+t5*t5;
    M[4]=t2*t4*(t4*t4-3*t5*t5)+t3*t5*(3*t4*t4-t5*t5);
    M[5]=t1*(t4*t4-t5*t5)+4*u11*t4*t5;
    M[6]=t3*t4*(t4*t4-3*t5*t5)-t2*t5*(3*t4*t4-t5*t5);

    return true;
}

#if 0
void rgb2gray(const Mat src, Mat &result)
{
    CV_Assert(src.depth() != sizeof(uchar)); //harus 8 bit

    result = Mat::zeros(src.rows, src.cols, CV_8UC1); //buat matrik 1 chanel
    uchar data;

    if(src.channels() == 3){

        for( int i = 0; i < src.rows; ++i)
            for( int j = 0; j < src.cols; ++j )
            {
                data = (uchar)(((Mat_<Vec3b>) src)(i,j)[0] * 0.0722 + ((Mat_<Vec3b>) src)(i,j)[1] * 0.7152 + ((Mat_<Vec3b>) src)(i,j)[2] * 0.2126);

                result.at<uchar>(i,j) = data;
            }


    }else{

        result = src;
    }

    threshold(result, result, 120, 255, 0);
}
#endif

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

double getScore(const vector<double> &Sa, const vector<double> &Ta) {
#if 0
	int i, sma, smb;
//	double eps = 1.e-5;
    double eps = 1.e-10;
    double mmm;
    double result = 0;
	for( i = 0; i < 7; i++ ) {
                double ama = fabs( ma[i] );
                double amb = fabs( mb[i] );

                if( ma[i] > 0 )
                    sma = 1;
                else if( ma[i] < 0 )
                    sma = -1;
                else
                    sma = 0;
                if( mb[i] > 0 )
                    smb = 1;
                else if( mb[i] < 0 )
                    smb = -1;
                else
                    smb = 0;

                if( ama > eps && amb > eps )
                {
                    ama = 1. / (sma * log10( ama ));
                    amb = 1. / (smb * log10( amb ));
                    result += fabs( -ama + amb );
                }

	}
	return result;
#endif

#if 0
	//  �������ƶ�1
double dbR =0; //���ƶ�
double dSigmaST =0;
    double dSigmaS =0;
    double dSigmaT =0;
    double temp =0;
    {for(int i=0;i<7;i++)
    {
        temp = fabs(Sa[i]*Ta[i]);
        dSigmaST+=temp;
        dSigmaS+=pow(Sa[i],2);
        dSigmaT+=pow(Ta[i],2);
    }}
    dbR = dSigmaST/(sqrt(dSigmaS)*sqrt(dSigmaT));
	return dbR;


#endif

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

int main(int /*argc*/, char /*argv*/)
{

    vector<double> huMoments;

    TimeTracker tt;
    Mat g1 = imread("images//KOYO.bmp", 0);
    cvHu("images//KOYO.bmp");
//    rgb2gray(g1, g1);
    Mat result;
    Canny(g1, result, 150, 220);
    tt.start();
    invariantMomentHu(result, huMoments);
    tt.stop();
    cout << "time: " << tt.duration() << endl;
    printMoments("gambar 1 : " ,huMoments);

    namedWindow("gambar1");
    imshow("gambar1", result);

    vector<double> huMoments2;
    Mat g2 = imread("images/KOYO2.jpg", 0);
    cvHu("images/KOYO2.jpg");
    Canny(g2, result, 150, 220);
//    rgb2gray(g2, g2);
    tt.start();
    invariantMomentHu(result, huMoments2);
    tt.stop();
    cout << "time: " << tt.duration() << endl;
    printMoments("gambar 2 : " ,huMoments2);

	cout << "Score: " << getScore(huMoments, huMoments2) << endl;

    namedWindow("gambar2");
    imshow("gambar2", result);

    vector<double> huMoments3;
    Mat g3 = imread("images/black640480.bmp", 0);
    Canny(g3, result, 150, 220);
//    rgb2gray(g3, g3);
    invariantMomentHu(result, huMoments3);
//    printMoments("gambar 3 : ", huMoments3);

	cout << "Score: " << getScore(huMoments, huMoments3) << endl;
    namedWindow("gambar3");
    imshow("gambar3", result);

    vector<double> huMoments4;
    Mat g4 = imread("images/region.bmp", 0);
    Canny(g4, result, 150, 220);
//    rgb2gray(g4, g4);
    invariantMomentHu(result, huMoments4);
//    printMoments("gambar 4 : ", huMoments4);

	cout << "Score: " << getScore(huMoments, huMoments4) << endl;
    namedWindow("gambar4");
    imshow("gambar4", result);


    waitKey(0);
	return 0;
}
