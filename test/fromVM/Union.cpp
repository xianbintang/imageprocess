#if 0
#include<iostream>
#include <string>
#include <list>
#include <vector>
#include <map>
#include <stack>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
using namespace std;

void Seed_Filling(const cv::Mat& binImg, cv::Mat& lableImg)   //种子填充法
{
	// 4邻接方法


	if (binImg.empty() ||
		binImg.type() != CV_8UC1)
	{
		return;
	}

	lableImg.release();
	binImg.convertTo(lableImg, CV_32SC1);

	int label = 1;  

	int rows = binImg.rows - 1;  
	int cols = binImg.cols - 1;
	for (int i = 1; i < rows-1; i++)
	{
		int* data= lableImg.ptr<int>(i);
		for (int j = 1; j < cols-1; j++)
		{
			if (data[j] == 1)
			{
				std::stack<std::pair<int,int>> neighborPixels;   
				neighborPixels.push(std::pair<int,int>(i,j));     // 像素位置: <i,j>
				++label;  // 没有重复的团，开始新的标签
				while (!neighborPixels.empty())
				{
					std::pair<int,int> curPixel = neighborPixels.top(); //如果与上一行中一个团有重合区域，则将上一行的那个团的标号赋给它
					int curX = curPixel.first;
					int curY = curPixel.second;
					lableImg.at<int>(curX, curY) = label;

					neighborPixels.pop();

					if (lableImg.at<int>(curX, curY-1) == 1)
					{//左边
						neighborPixels.push(std::pair<int,int>(curX, curY-1));
					}
					if (lableImg.at<int>(curX, curY+1) == 1)
					{// 右边
						neighborPixels.push(std::pair<int,int>(curX, curY+1));
					}
					if (lableImg.at<int>(curX-1, curY) == 1)
					{// 上边
						neighborPixels.push(std::pair<int,int>(curX-1, curY));
					}
					if (lableImg.at<int>(curX+1, curY) == 1)
					{// 下边
						neighborPixels.push(std::pair<int,int>(curX+1, curY));
					}
				}		
			}
		}
	}
	
}

void Two_Pass(const cv::Mat& binImg, cv::Mat& lableImg)    //两遍扫描法
{
	if (binImg.empty() ||
		binImg.type() != CV_8UC1)
	{
		return;
	}

	// 第一个通路

	lableImg.release();
	binImg.convertTo(lableImg, CV_32SC1);

	int label = 1; 
	std::vector<int> labelSet;
	labelSet.push_back(0);  
	labelSet.push_back(1);  

	int rows = binImg.rows - 1;
	int cols = binImg.cols - 1;
	for (int i = 1; i < rows; i++)
	{
		int* data_preRow = lableImg.ptr<int>(i-1);
		int* data_curRow = lableImg.ptr<int>(i);
		for (int j = 1; j < cols; j++)
		{
			if (data_curRow[j] == 1)
			{
				std::vector<int> neighborLabels;
				neighborLabels.reserve(2);
				int leftPixel = data_curRow[j-1];
				int upPixel = data_preRow[j];
				if ( leftPixel > 1)
				{
					neighborLabels.push_back(leftPixel);
				}
				if (upPixel > 1)
				{
					neighborLabels.push_back(upPixel);
				}

				if (neighborLabels.empty())
				{
					labelSet.push_back(++label);  // 不连通，标签+1
					data_curRow[j] = label;
					labelSet[label] = label;
				}
				else
				{
					std::sort(neighborLabels.begin(), neighborLabels.end());
					int smallestLabel = neighborLabels[0];  
					data_curRow[j] = smallestLabel;

					// 保存最小等价表
					for (size_t k = 1; k < neighborLabels.size(); k++)
					{
						int tempLabel = neighborLabels[k];
						int& oldSmallestLabel = labelSet[tempLabel];
						if (oldSmallestLabel > smallestLabel)
						{							
							labelSet[oldSmallestLabel] = smallestLabel;
							oldSmallestLabel = smallestLabel;
						}						
						else if (oldSmallestLabel < smallestLabel)
						{
							labelSet[smallestLabel] = oldSmallestLabel;
						}
					}
				}				
			}
		}
	}

	// 更新等价对列表
	// 将最小标号给重复区域
	for (size_t i = 2; i < labelSet.size(); i++)
	{
		int curLabel = labelSet[i];
		int preLabel = labelSet[curLabel];
		while (preLabel != curLabel)
		{
			curLabel = preLabel;
			preLabel = labelSet[preLabel];
		}
		labelSet[i] = curLabel;
	}  ;

	for (int i = 0; i < rows; i++)
	{
		int* data = lableImg.ptr<int>(i);
		for (int j = 0; j < cols; j++)
		{
			int& pixelLabel = data[j];
			pixelLabel = labelSet[pixelLabel];	
		}
	}
}
//彩色显示
cv::Scalar GetRandomColor()
{
	uchar r = 255 * (rand()/(1.0 + RAND_MAX));
	uchar g = 255 * (rand()/(1.0 + RAND_MAX));
	uchar b = 255 * (rand()/(1.0 + RAND_MAX));
	return cv::Scalar(b,g,r);
}


void LabelColor(const cv::Mat& labelImg, cv::Mat& colorLabelImg) 
{
	if (labelImg.empty() ||
		labelImg.type() != CV_32SC1)
	{
		return;
	}

	std::map<int, cv::Scalar> colors;

	int rows = labelImg.rows;
	int cols = labelImg.cols;

	colorLabelImg.release();
	colorLabelImg.create(rows, cols, CV_8UC3);
	colorLabelImg = cv::Scalar::all(0);

	for (int i = 0; i < rows; i++)
	{
		const int* data_src = (int*)labelImg.ptr<int>(i);
		uchar* data_dst = colorLabelImg.ptr<uchar>(i);
		for (int j = 0; j < cols; j++)
		{
			int pixelValue = data_src[j];
			if (pixelValue > 1)
			{
				if (colors.count(pixelValue) <= 0)
				{
					colors[pixelValue] = GetRandomColor();
				}

				cv::Scalar color = colors[pixelValue];
				*data_dst++   = color[0];
				*data_dst++ = color[1];
				*data_dst++ = color[2];
			}
			else
			{
				data_dst++;
				data_dst++;
				data_dst++;
			}
		}
	}
}


int main()
{
	int num = 2;
	cv::Mat binImage = cv::imread("contour.bmp", 0);
	cv::threshold(binImage, binImage, 50, 1, CV_THRESH_BINARY_INV);
	cv::Mat labelImg;
	Two_Pass(binImage, labelImg);
	//Seed_Filling(binImage, labelImg);
	//彩色显示
	cv::Mat colorLabelImg;
	LabelColor(labelImg, colorLabelImg);
	cv::imshow("colorImg", colorLabelImg);
	//灰度显示
	cv::Mat grayImg;
	labelImg *= 10;
	labelImg.convertTo(grayImg, CV_8UC1);
	cv::imshow("labelImg", grayImg);


	cv::waitKey(0);
	return 0;
}
#endif