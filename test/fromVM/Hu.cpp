#if 0
/*===============================================//
���ܣ������ƥ��
ʱ�䣺3/28/2011 SkySeraph HQU
�ο���
//===============================================*/
#include "iostream"
using namespace std;
#include <opencv2\opencv.hpp>
#include "opencv/cv.h"
#include "opencv/highgui.h"

#include "math.h"

/*
#pragma comment(lib,"highgui.lib")
#pragma comment(lib,"cv.lib")
#pragma comment(lib,"cvaux.lib")
#pragma comment(lib,"cxcore.lib")
*/
//const char* filename ="D:\\My Documents\\My Pictures\\Images\\1.bmp";
//const char* filename2 ="D:\\My Documents\\My Pictures\\Images\\2.bmp";

/*=============================================*/
double M[7] = {0}; //HU�����

bool HuMoment(IplImage* img)
{

	int bmpWidth = img->width;
	int bmpHeight = img->height;
	int bmpStep = img->widthStep;
	int bmpChannels = img->nChannels;
	uchar*pBmpBuf = (uchar*)img->imageData;

	double m00=0,m11=0,m20=0,m02=0,m30=0,m03=0,m12=0,m21=0; //���ľ�
	double x0=0,y0=0; //�������ľ�ʱ��ʹ�õ���ʱ������x-x'��
	double u20=0,u02=0,u11=0,u30=0,u03=0,u12=0,u21=0;//�淶��������ľ�
	//double M[7]; //HU�����
	double t1=0,t2=0,t3=0,t4=0,t5=0;//��ʱ������
	//double Center_x=0,Center_y=0;//����
	int Center_x=0,Center_y=0;//����
	int i,j; //ѭ������

	// ���ͼ�����������
	double s10=0,s01=0,s00=0; //0�׾غ�1�׾� //ע����ֵͼ���0�׾ر�ʾ���
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

	// ������ס����׾�
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

	// ����淶��������ľ�
	u20=m20/pow(m00,2);
	u02=m02/pow(m00,2);
	u11=m11/pow(m00,2);
	u30=m30/pow(m00,2.5);
	u03=m03/pow(m00,2.5);
	u12=m12/pow(m00,2.5);
	u21=m21/pow(m00,2.5);

	// �����м������
	t1=(u20-u02);
	t2=(u30-3*u12);
	t3=(3*u21-u03);
	t4=(u30+u12);
	t5=(u21+u03);

	// ���㲻���
	M[0]=u20+u02;
	M[1]=t1*t1+4*u11*u11;
	M[2]=t2*t2+t3*t3;
	M[3]=t4*t4+t5*t5;
	M[4]=t2*t4*(t4*t4-3*t5*t5)+t3*t5*(3*t4*t4-t5*t5);
	M[5]=t1*(t4*t4-t5*t5)+4*u11*t4*t5;
	M[6]=t3*t4*(t4*t4-3*t5*t5)-t2*t5*(3*t4*t4-t5*t5);


	/*cout<<M[0]<<endl;//<<��"<<M[0]<<"��"<<M[0]<<"��"<<M[0]<<"��"<<M[0]<<"��"<<M[0]<<"��"<<M[0]<<endl;
	cout<<M[1]<<endl;
	cout<<M[2]<<endl;
	cout<<M[3]<<endl;
	cout<<M[4]<<endl;
	cout<<M[5]<<endl;
	cout<<M[6]<<endl;
	cout<<endl;*/
	return true;
}


int main(char argc,char** argv)
{
	int i;
	double Sa[7] = {0},Ta[7] ={0};

	///*Դͼ��
	IplImage*img = cvLoadImage("gp.bmp",0);//�Ҷ�
	HuMoment(img);
	for(i=0;i<7;i++)
	{
		Sa[i] = M[i];
		M[i] =0;
	}
	cout<<Sa[0]<<endl;
	cout<<Sa[1]<<endl;
	cout<<Sa[2]<<endl;
	cout<<Sa[3]<<endl;
	cout<<Sa[4]<<endl;
	cout<<Sa[5]<<endl;
	cout<<Sa[6]<<endl;
	cout<<endl;
	//*/


	///*ģ��ͼ
	IplImage*tpl = cvLoadImage("gpt.bmp",0);//�Ҷ�
	HuMoment(tpl);
	for(i=0;i<7;i++)
	{
		Ta[i] = M[i];
		M[i] =0;
	}
	cout<<Ta[0]<<endl;
	cout<<Ta[1]<<endl;
	cout<<Ta[2]<<endl;
	cout<<Ta[3]<<endl;
	cout<<Ta[4]<<endl;
	cout<<Ta[5]<<endl;
	cout<<Ta[6]<<endl;
	cout<<endl;


	// �������ƶ�
	double dbR =0; //���ƶ�
	double dSigmaST =0;
	double dSigmaS =0;
	double dSigmaT =0;
	double temp =0; 

	for(i=0;i<7;i++)
	{
		temp = Sa[i]*Ta[i];
		dSigmaST+=temp;
		dSigmaS+=pow(Sa[i],2);
		dSigmaT+=pow(Ta[i],2);
	}
	dbR = dSigmaST/(sqrt(dSigmaS)*sqrt(dSigmaT));
	printf("%lf\n",dbR);
	//cout<<dbR<<endl;

	 CvMoments moments;
     CvHuMoments hu;     
	 cvMoments(img,&moments,0); 
     cvGetHuMoments(&moments, &hu);
     cout<<hu.hu1<<"/"<<hu.hu2<<"/"<<hu.hu3<<"/"<<hu.hu4<<"/"<<hu.hu5<<"/"<<hu.hu6<<"/"<<hu.hu7<<"/"<<"/"<<endl;
     cvMoments(tpl,&moments,0); 
     cvGetHuMoments(&moments, &hu);
     cout<<hu.hu1<<"/"<<hu.hu2<<"/"<<hu.hu3<<"/"<<hu.hu4<<"/"<<hu.hu5<<"/"<<hu.hu6<<"/"<<hu.hu7<<"/"<<"/"<<endl;



	cvReleaseImage(&img);
	cvReleaseImage(&tpl);
	system("pause");

	return 0;
}
#endif