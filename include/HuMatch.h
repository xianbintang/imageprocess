//
// Created by xianb on 2017/5/3.
//

#ifndef KOYO_HUMATCH_H
#define KOYO_HUMATCH_H
const int HUMOMENTSNUMBIG= 1500;
const int HUMOMENTSNUMSMALL= 500;
const int POINTBIG = 1500;
const int POINTSMALL = 700;


class HuMatch {

private:
    CvSize tSize;
    std::vector<double> templateHuMoments;
    std::vector<cv::Point> searchImgContour;
    int count = 0;
    void DrawContours(IplImage * source, cv::Point p, CvSize size) {
        CvPoint point;
        point.x = p.x + size.width;
        point.y = p.y + size.height;
        cvRectangle(source,p, point, cvScalar(255, 255, 255));
    }


    void invariantMomentHu(const std::vector<cv::Point> &src, std::vector<double>& result)
    {
        //momen dua dimensi
        double m00=0, m10=0, m01=0;
        int x,y;

        for(auto const &p: src) {
            x = p.x + 1;
            y = p.y + 1;

            m00 += 1;
            m10 += x * 1;
            m01 += y * 1;
            count++;
        }
        //central moment lainnya & kemudian normalized central moments
        double xbar = m10 / m00;
        double ybar = m01 / m00;

        double m11 = 0, m20 = 0, m02 = 0, m30 = 0,
                m03 = 0, m12 = 0, m21 = 0;

        for(auto const &p:src) {
            x = p.x + 1;
            y = p.y + 1;
            m11 += ( x - xbar ) * ( y - ybar ) * 1 / pow( m00, 2.0 );
            m20 += pow( ( x - xbar ), 2.0 ) * 1 / pow( m00, 2.0 );
            m02 += pow( ( y - ybar ), 2.0) * 1 / pow( m00, 2.0 );
            m30 += pow( ( x - xbar ), 3.0 )  * 1 / pow( m00, 2.5 );
            m03 += pow( ( y - ybar ), 3.0 ) * 1 / pow( m00, 2.5 );
            m12 += ( x - xbar ) * pow( ( y - ybar ), 2.0 ) * 1 / pow( m00, 2.5 );
            m21 += pow( ( x - xbar ), 2.0 ) * ( y - ybar ) * 1 / pow( m00, 2.5 );
            count++;
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
    /* Time monster */
    void getContourPoints(const IplImage *src, std::vector<cv::Point> &contour) {
        for(int i = 0; i < src->height; i++) {
            unsigned char * srow = (unsigned char *)(src->imageData + i * src->widthStep);
            for(int j = 0; j < src->width; j++) {
                if(srow[j] != 0) {
                    contour.push_back(cv::Point(j, i));
                }
                count++;
            }
        }
    }
public:
    /* the source img should be a binary image*/
    void CreateTemplate(IplImage * templateImg) {
        std::vector<cv::Point> contour;
        getContourPoints(templateImg, contour);

        invariantMomentHu(contour, templateHuMoments);
        tSize = cvGetSize(templateImg);
    }

    double getScore(const std::vector<double> &Sa) {
#if 0
        auto ma = Sa;
        auto mb = templateHuMoments;
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
        auto Ta = templateHuMoments;
double dbR =0;
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
#if 1

        //  计算相似度1
        double dbR =0; //相似度
        double dSigmaST =0;
        double dSigmaS =0;
        double dSigmaT =0;
        double temp =0;
        {
            for(int i=0;i<7;i++)
            {
                temp = fabs(Sa[i]*templateHuMoments[i]);
                dSigmaST+=temp;
                dSigmaS+=pow(Sa[i],2);
                dSigmaT+=pow(templateHuMoments[i],2);
            }
        }
        dbR = dSigmaST/(sqrt(dSigmaS)*sqrt(dSigmaT));
        return dbR;
#endif
    }

    void getROIContour(std::vector<cv::Point> &contour, const cv::Rect &rect) {
        for (auto const &p : searchImgContour) {
            int x = p.x;
            int y = p.y;
            if ((x > rect.x) && x < (rect.x + rect.width) && ((y > rect.y) && y < (rect.y + rect.height))) {
                contour.push_back(p);
            }
            count++;
        }
    }
    /* SearchImg shoud be a binary image */
    void matchThis(IplImage * SearchImg, double threshold) {
        getContourPoints(SearchImg, searchImgContour);
        std::cout << "SearchImgContour: " << searchImgContour.size() << std::endl;
        double maxScore = -99999.9;
        cv::Point p;

        tSize.width += 3;
        tSize.height = tSize.height + 3;
        TimeTracker tt;
        for(int i = 0; i < SearchImg->height - tSize.height; i += 5) {
            for (int j = 0; j < SearchImg->width - tSize.width; j += 5) {
//                tt.start();
                std::vector<double> SearchHuMoments2;
//                SearchHuMoments2.reserve(POINTSMALL);
                /* get ROI */

                count++;
                std::vector<cv::Point> contour;
//                contour.reserve(POINTSMALL);
                getROIContour(contour, cv::Rect(j,i,tSize.width,tSize.height));
                if(contour.size() < searchImgContour.size() / 5) {
                    continue;
                }

//                std::cout << "ROI Contour: " << contour.size() << std::endl;
//            std::cout << "size of Contour: " << contour.size() << std::endl;
                /* getinvariantMomentHu */
                invariantMomentHu(contour, SearchHuMoments2);
//                tt.stop();
//                std::cout << "Debug Time: " << tt.duration() << std::endl;
                /* getScore*/
                double score =  getScore(SearchHuMoments2);
//                std::cout << "Current Score: " << score << std::endl;

                if(score > threshold) {
                    cvCircle(SearchImg,cv::Point(j, i), 1, cvScalar(255, 255, 255));
                    std::cout << "score: " << score << std::endl;
                    p.x = j;
                    p.y = i;
                    DrawContours(SearchImg,p, tSize);
//                    cv::waitKey(0);
                }

                if (score > maxScore) {

                    maxScore = score;
                    p.x = j;
                    p.y = i;
//                    DrawContours(SearchImg,p, tSize);
//                    cv::waitKey(0);
                }
            }

        }
        std::cout << "Point is: " << "x: " << p.x << " y: "<< p.y << std::endl;
        std::cout << "Max Score: " << maxScore << std::endl;
        std::cout << "count: " << count << std::endl;
        DrawContours(SearchImg,p, tSize);
    }

    void printMoments()
    {
//        std::cout<<label<<std::endl<<setprecision(50);;
        for(unsigned int i = 0; i < templateHuMoments.size(); ++i){
            std::cout<<i+1<<" = "<<templateHuMoments[i]<<std::endl;
        }
    }

};


#endif //KOYO_HUMATCH_H
