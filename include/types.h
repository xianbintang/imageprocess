//
// Created by xianb on 2017/5/9.
//

#ifndef KOYO_TYPES_H_H
#define KOYO_TYPES_H_H



typedef unsigned char uchar;

typedef struct  IPoint_{

	int x;
	int y;
} IPoint;


typedef struct ISize
{
	int width;
	int height;
} ISize;

#define U8C1 0
#define S16C1 1

typedef struct IMat{
    int width;
	int height;
	int *ptr;
	short *fptr;
    int step;
	int type;
} IMat;


struct Image{
    int width;
    int height;
    char * imageData;
    int imageSize;
    int widthStep;
    int nChannels;
    int depth;
    int roi;
};
typedef struct Image Image;

int rotateImage(const Image imageSrc, Image imageDst, const double degree);

IMat *ICreateMat( const int height, const int width, int type);
IMat *fICreateMat( const int height, const int width, int type);
#define MIN(x, y) (x) < (y) ? (x) : (y)

#endif //KOYO_TYPES_H_H
