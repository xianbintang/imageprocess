#include <iostream>  
#include <opencv2/core/core.hpp>  
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/opencv.hpp>
#include <vector>

 typedef uchar UINT8;
const int WIDTH = 640;
const int HEIGHT = 480;

void printHSV(cv::Mat hsvimage) {
	std::ofstream of;
	of.open("data//csc1.txt");
	for (int i = 0; i < hsvimage.rows; ++i) {
		for (int j = 0; j < hsvimage.cols; ++j) {
			auto hsvv =hsvimage.at<cv::Vec3b>(i,j);
			//printf("(%d,%d)H= %d,S = %d, V = %d\n",i,j,hsvv[0],hsvv[1],hsvv[2]);
			of<< "H: " << (int)hsvv[0] << " " << "S: " <<(int)hsvv[1] << " " << "V: " <<(int)hsvv[2] << " i: " << i << " " << "j: " << j << " " << std::endl;;
		}
	}
	of.close();
}
int get_hsv(const UINT8 *yuv, UINT8 *hsv[3] ) {
	cv::Mat imgRGB, imgHSV, img_h, img_s, img_v;

	cv::Mat imgYUV(480 + 480 / 2, 640, CV_8UC1, (void*)yuv);
	// imgRGB.create(480, 640, CV_8UC3);
	cv::cvtColor(imgYUV, imgRGB, CV_YUV420sp2BGR, 3);
	cv::cvtColor(imgRGB, imgHSV, CV_BGR2HSV_FULL, 3);

    printHSV(imgHSV);
	cv::Mat hsv_vec[3];
	cv::split(imgHSV, hsv_vec);
	img_h = hsv_vec[0];
	img_s = hsv_vec[1];
	img_v = hsv_vec[2];

	if (!yuv) {
		std::cout << "No allocate memory!\n" << std::endl;
		return -1;
	}

	cv::imshow("imgRGB", imgRGB);
	cvWaitKey(0);
	cv::imshow("imgHSV_H", img_h);
	cvWaitKey(0);
	cv::imshow("imgHSV_S", img_s);
	cvWaitKey(0);
	cv::imshow("imgHSV_V", img_v);
	cvWaitKey(0);

	for (int i = 0; i < HEIGHT; ++i) {
		for (int j = 0; j < WIDTH; ++j) {
			hsv[0][i * WIDTH + j] = img_h.at<uchar>(i, j);
			hsv[1][i * WIDTH + j] = img_s.at<uchar>(i, j);
			hsv[2][i * WIDTH + j] = img_v.at<uchar>(i, j);
		}
	}

	return 0;
}


int main(int argc, char ** argv) {
	
	std::string filename(argv[1]);
	unsigned char *buf = nullptr;
	uchar *hsv[3];
    hsv[0] = new uchar[640*480];
	hsv[1] = new uchar[640*480];
	hsv[2] = new uchar[640*480];

	if (filename.substr(filename.size() - 3, 3) == "yuv") {
		FILE *yuv_file = fopen(filename.c_str(), "rb+");

		buf = new unsigned char[640 * 480 + 320 * 240 * 2];
		fread(buf, 640 * 480 + 320 * 240 * 2, 1, yuv_file);
		fclose(yuv_file);

		
		get_hsv(buf, hsv);
		std::cout << std::endl;
	}

	return 0;
}