//
// Created by xianb on 2017/5/15.
//
#include <malloc.h>
#include <cstdio>
#include "types.h"


void saveResultf(IMat mat, char *path) {
    FILE *fp = fopen(path, "w");
    int i,j;
    for (i = 0; i < mat.height; ++i) {
        for (j = 0; j < mat.width; ++j) {
            fprintf(fp, "%d ", (mat.fptr + i * mat.step)[j]);
        }
        fprintf(fp, "\n");
    }
}
void saveResult(IMat mat, char *path) {
    FILE *fp = fopen(path, "w");
    int i,j;
    for (i = 0; i < mat.height; ++i) {
        for (j = 0; j < mat.width; ++j) {
            fprintf(fp, "%d ", (mat.ptr + i * mat.step)[j]);
        }
        fprintf(fp, "\n");
    }
}


IMat *ICreateMat( const int height, const int width, int type)
{
    IMat * rtv = (IMat *)malloc(sizeof(IMat));
    rtv->height= height;
    rtv->width = width;
    rtv->step = width;
    rtv->type = type;
    rtv->ptr = (int *)malloc(sizeof(int) * height * width);
		int i, j;
    for (i = 0; i < height; ++i) {
        for (j = 0; j < width; ++j) {
            (rtv->ptr + rtv->step*i)[j] = 0;
        }
    }

    return rtv;
}

IMat *fICreateMat( const int height, const int width, int type)
{
    IMat * rtv = (IMat *)malloc(sizeof(IMat));
    rtv->height= height;
    rtv->width = width;
    rtv->step = width;
    rtv->type = type;
    rtv->fptr = (short*)malloc(sizeof(short) * height * width);
		int i, j;
    for (i = 0; i < height; ++i) {
        for (j = 0; j < width; ++j) {
            (rtv->fptr + rtv->step*i)[j] = 0;
        }
    }

    return rtv;
}

void resumeIMat(const char *path, IMat *mat) {
    FILE *fp = fopen(path, "r");
    int i,j;
    int num;
    for (i = 0; i < mat->height; ++i) {
        for (j = 0; j < mat->width; ++j) {
            fscanf(fp, "%d ", &num);
            (mat->ptr + i * mat->step)[j] = num;
        }
    }
}


void resumeIMatf(const char *path, IMat *mat) {
    FILE *fp = fopen(path, "r");
    int i,j;
    int num;
    for (i = 0; i < mat->height; ++i) {
        for (j = 0; j < mat->width; ++j) {
            fscanf(fp, "%d ", &num);
            (mat->fptr + i * mat->step)[j] = num;
        }
    }
}
