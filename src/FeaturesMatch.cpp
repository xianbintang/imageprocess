//
// Created by xianb on 2017/5/2.
//

#include "opencv2/core/core.hpp"
#include "opencv2/calib3d/calib3d.hpp"
#include "opencv2/features2d/features2d.hpp"
//#include "opencv2/nonfree/features2d.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
//#include <opencv2/legacy/legacy.hpp>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>

using namespace cv;
using namespace std;
void getMatchLR(string imgpath);
void getMatchDT(string imgpath);
//void getSURFFeature(string imgpath);
void getORBFeature(string imgpath);
void MatchTemplate(string imgname);
void MatchTemplateORB(string imgname);
void MatchTemplateSIFT(string imgname);
void testMatch();

static string imgdir = "E:\\input_resize\\";
static string Matchsavedir = "E:\\macth_result\\";
static string TemplateDir = "E:\\template_resize\\";
//static string ORBsavedir = "ORB_result\\";
IplImage mat_src_lpl;
IplImage OutImage_lpl;
IplImage TemplateIpl;

IplImage* dst_left;
IplImage* dst_right;
IplImage *src;
IplImage *TemplateIplPtr;

Mat mat_src;
Mat OutImage;
Mat TemplateMat;



int main()
{

    string imgname = "07_in.jpg";
    //testMatch();
    MatchTemplateSIFT(imgname);
    //getMatchLR(imgname);

    return 0;
}

void MatchTemplateSIFT(string imgname)
{
    vector<KeyPoint> LeftKey;
    vector<KeyPoint> RightKey;
    Mat LeftDescriptor;
    Mat RightDescriptor;
    vector<DMatch> Matches;
    //vector<double> meanDisances;
    IplImage* dst;

    ofstream fout("E:\\picin\\rBRIEF_Test\\templateCompareResult.txt");
    int xstep = 2;
    int ystep = 2;
    string TemplateName = "template1.jpg";
    string TemplateImgPath = TemplateDir + TemplateName;

    string imgpath = imgdir + imgname;
    string resultPath = Matchsavedir+imgname;

    TemplateMat = imread(TemplateImgPath,CV_LOAD_IMAGE_GRAYSCALE);
    if(!TemplateMat.data){
        cout<<"no template exist";
        return ;
    }
    int TemplateWidth = TemplateMat.cols;
    int TemplateHeight = TemplateMat.rows;
    std::cout<<"TemplateWidth "<<TemplateWidth<<endl;
    std::cout<<"TemplateHeight "<<TemplateHeight<<endl;

    FeatureDetector *pDetector = new SiftFeatureDetector; 	pDetector->detect(TemplateMat, LeftKey);
    DescriptorExtractor *pExtractor = new SiftDescriptorExtractor;
    pExtractor->compute(TemplateMat, LeftKey, LeftDescriptor);

    DescriptorMatcher *pMatcher = new BruteForceMatcher<L2<float>>;



    mat_src = imread(imgpath, CV_LOAD_IMAGE_GRAYSCALE );
    if(!mat_src.data) {
        cout<<"no src img";
        return ;
    }	mat_src_lpl = IplImage(mat_src);
    src  = &mat_src_lpl;
    long ImgWidth = src->width;
    long ImgHeight = src->height;
    std::cout<<"ImgWidth "<<ImgWidth<<endl;
    std::cout<<"ImgHeight "<<ImgHeight<<endl;

    int x;
    int y;


    //Mat R = Mat(ImgHeight - TemplateHeight, ImgWidth - TemplateWidth, CV_8UC1,255);
    //namedWindow("result", CV_WINDOW_NORMAL );
    //imshow("result",mat_src);
    //uchar *p;
    //while(start_x < ImgWidth - TemplateWidth){
    for(long start_y = 0;start_y <ImgHeight - TemplateHeight;start_y = start_y+ystep){
        for(long start_x = 0;start_x < ImgWidth - TemplateWidth;start_x = start_x+xstep){

            x = start_x;
            y = start_y;


            cvSetImageROI(src,cvRect(x,y,TemplateWidth,TemplateHeight));
            dst = cvCreateImage(cvSize(TemplateWidth,TemplateHeight),
                                IPL_DEPTH_8U,
                                src->nChannels);
            cvCopy(src,dst,0);
            cvResetImageROI(src);

            Mat DstImage = Mat(dst, false); // Do not copy
            pDetector->detect(DstImage, RightKey);

            pExtractor->compute(DstImage, RightKey, RightDescriptor);

            pMatcher->match(LeftDescriptor, RightDescriptor, Matches);
            //double sum = 0;
            double sum = 0;
            //int i = 0;
            for(vector<DMatch>::iterator dite = Matches.begin();dite <Matches.end(); dite++ )
            {
                sum += dite->distance;
            }
            int matchSize = Matches.size()*10;

            if(matchSize>1){

                fout<<exp(-sum/matchSize)<<" ";

            }else{
                fout<<exp(-100.0)<<" ";

            }

        }
    }
    //
    std::cout<<"finish";
    fout.close();

    //destroyWindow("result" );
    //imwrite(resultPath,R);

    delete pDetector;
    delete pExtractor;
    delete pMatcher;

}


void MatchTemplateORB(string imgname)
{
    vector<KeyPoint> LeftKey;
    vector<KeyPoint> RightKey;
    Mat LeftDescriptor;
    Mat RightDescriptor;
    vector<DMatch> Matches;
    //vector<double> meanDisances;
    IplImage* dst;

    ofstream fout("E:\\picin\\rBRIEF_Test\\templateCompareResult.txt");
    int xstep = 2;
    int ystep = 2;
    string TemplateName = "template2.jpg";
    string TemplateImgPath = TemplateDir + TemplateName;

    string imgpath = imgdir + imgname;
    string resultPath = Matchsavedir+imgname;

    TemplateMat = imread(TemplateImgPath,CV_LOAD_IMAGE_GRAYSCALE);
    if(!TemplateMat.data){
        cout<<"no template exist";
        return ;
    }
    int TemplateWidth = TemplateMat.cols;
    int TemplateHeight = TemplateMat.rows;
    std::cout<<"TemplateWidth "<<TemplateWidth<<endl;
    std::cout<<"TemplateHeight "<<TemplateHeight<<endl;
    ORB orb;
    orb(TemplateMat,Mat(),LeftKey,LeftDescriptor);

    DescriptorMatcher *pMatcher = new BruteForceMatcher<HammingLUT>;



    mat_src = imread(imgpath, CV_LOAD_IMAGE_GRAYSCALE );
    if(!mat_src.data) {
        cout<<"no src img";
        return ;
    }
    mat_src_lpl = IplImage(mat_src);
    src  = &mat_src_lpl;
    long ImgWidth = src->width;
    long ImgHeight = src->height;
    std::cout<<"ImgWidth "<<ImgWidth<<endl;
    std::cout<<"ImgHeight "<<ImgHeight<<endl;



    int x;
    int y;


    //Mat R = Mat(ImgHeight - TemplateHeight, ImgWidth - TemplateWidth, CV_8UC1,255);
    //namedWindow("result", CV_WINDOW_NORMAL );
    //imshow("result",mat_src);
    //uchar *p;
    //while(start_x < ImgWidth - TemplateWidth){
    for(long start_y = 0;start_y <ImgHeight - TemplateHeight;start_y = start_y+ystep){
        for(long start_x = 0;start_x < ImgWidth - TemplateWidth;start_x = start_x+xstep){

            x = start_x;
            y = start_y;
            //std::cout<<"<"<<x<<","<<y<<">"<<" ";

            cvSetImageROI(src,cvRect(x,y,TemplateWidth,TemplateHeight));
            dst = cvCreateImage(cvSize(TemplateWidth,TemplateHeight),
                                IPL_DEPTH_8U,
                                src->nChannels);
            cvCopy(src,dst,0);
            cvResetImageROI(src);

            Mat DstImage = Mat(dst, false); // Do not copy
            orb(DstImage,Mat(),RightKey,RightDescriptor);
            //std::cout<<RightDescriptor.size();
            pMatcher->match(LeftDescriptor, RightDescriptor, Matches);
            //double sum = 0;
            double sum = 0;
            //int i = 0;
            for(vector<DMatch>::iterator dite = Matches.begin();dite <Matches.end(); dite++ )
            {
                sum += dite->distance;
            }
            int matchSize = Matches.size()*10;
            //std::cout<<matchSize<<" ";
            if(matchSize>1){
                //int meanDis = sum/matchSize;
                fout<<exp(-sum/matchSize)<<" ";
                //std::cout<<"meanDis"<<meanDis<<" ";
                //R.at<uchar>(x,y) = meanDis;
            }else{
                fout<<exp(-100.0)<<" ";
                //fout<<255<<" ";
                //std::cout<<"meanDis"<<255<<" ";
            }

        }
        //std::cout<<endl;
        fout<<"\n";
        //start_x += step;
    }
    //
    std::cout<<"finish";
    fout.close();

    //destroyWindow("result" );
    //imwrite(resultPath,R);

    //delete pDetector;
    //delete pExtractor;
    delete pMatcher;

}


void getMatchLR(string imgname)
{

    string imgpath = imgdir + imgname;
    mat_src = imread(imgpath, CV_LOAD_IMAGE_GRAYSCALE );
    if(!mat_src.data) {
        cout<<"no img";
        return ;
    }
    mat_src_lpl = IplImage(mat_src);
    src  = &mat_src_lpl;


    //
    cvSetImageROI(src,cvRect(0,0,0.5*src->width,src->height));
    dst_left = cvCreateImage(cvSize(0.5*src->width,src->height),
                             IPL_DEPTH_8U,
                             src->nChannels);
    cvCopy(src,dst_left,0);
    cvResetImageROI(src);

    cvSetImageROI(src,cvRect(0.5*src->width,0,0.5*src->width,src->height));
    dst_right = cvCreateImage(cvSize(0.5*src->width,src->height),
                              IPL_DEPTH_8U,
                              src->nChannels);
    cvCopy(src,dst_right,0);
    cvResetImageROI(src);


    // Convert IplImage to cv::Mat
    Mat matLeftImage = Mat(dst_left, false); // Do not copy
    Mat matRightImage = Mat(dst_right, false);

    // Key point and its descriptor
    vector<KeyPoint> LeftKey;
    vector<KeyPoint> RightKey;
    Mat LeftDescriptor;
    Mat RightDescriptor;
    vector<DMatch> Matches;

    /*
    // Detect key points from image
    FeatureDetector *pDetector = new FastFeatureDetector; //
    pDetector->detect(matLeftImage, LeftKey);
    pDetector->detect(matRightImage, RightKey);
    delete pDetector;

    // Extract descriptors
    DescriptorExtractor *pExtractor = new BriefDescriptorExtractor; //
    pExtractor->compute(matLeftImage, LeftKey, LeftDescriptor);
    pExtractor->compute(matRightImage, RightKey, RightDescriptor);
    delete pExtractor;
    */


    ORB orb;
    orb(matLeftImage,Mat(),LeftKey,LeftDescriptor);
    orb(matRightImage,Mat(),RightKey,RightDescriptor);


    // Matching features
    //DescriptorMatcher *pMatcher = new FlannBasedMatcher; //
    DescriptorMatcher *pMatcher = new BruteForceMatcher<HammingLUT>; //
    pMatcher->match(LeftDescriptor, RightDescriptor, Matches);
    delete pMatcher;


    double max_dist = 0; double min_dist = 200;

    //-- Quick calculation of max and min distances between keypoints
    for( int i = 0; i < LeftDescriptor.rows; i++ )
    { double dist = Matches[i].distance;
        if( dist < min_dist ) min_dist = dist;
        if( dist > max_dist ) max_dist = dist;
    }

    //printf("-- Max dist : %f \n", max_dist );
    //printf("-- Min dist : %f \n", min_dist );

    //-- Draw only "good" matches (i.e. whose distance is less than 2*min_dist )
    //-- PS.- radiusMatch can also be used here.
    std::vector< DMatch > good_matches;

    for( int i = 0; i < LeftDescriptor.rows; i++ )
    { if( Matches[i].distance < 0.5*max_dist )
        { good_matches.push_back( Matches[i]); }
    }

    // Show result
    //drawMatches(matLeftImage, LeftKey, matRightImage, RightKey, Matches, OutImage);
    drawMatches( matLeftImage, LeftKey, matRightImage, RightKey,
                 good_matches, OutImage, Scalar::all(-1), Scalar::all(-1),
                 vector<char>(), DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS );
    OutImage_lpl = IplImage(OutImage);
    cvNamedWindow( "Match features", 1);
    cvShowImage("Match features", &(OutImage_lpl));
    cvWaitKey( 0 );
    cvDestroyWindow( "Match features" );

    string savepath = Matchsavedir + imgname;
    //-- Show detected (drawn) keypoints
    imwrite(savepath, OutImage );//
}

