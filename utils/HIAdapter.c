//
// Created by xianb on 2017/5/17.
//

#include <types.h>
#include <stdio.h>
#include <malloc.h>
#include "hi_common.h"
#include "hi_comm_video.h"
#include "hi_comm_sys.h"
#include "hi_comm_ive.h"

#include "mpi_vb.h"
#include "mpi_sys.h"
#include "mpi_ive.h"

#include "sample_comm_ive.h"

void HIIMage2IMat(IVE_IMAGE_S *image, IMat *mat) {
	HI_U16 height = image->u16Height;
	HI_U16 width = image->u16Width;
	int i, j;

	HI_U8 *pU8 = NULL;
	HI_S16 *pS16 = NULL;

    if (height != mat->height) {
		printf("ERROR: size of mat and image are not matched!\n");
		return ;
	}
	if (image->enType == IVE_IMAGE_TYPE_U8C1 && mat->type == U8C1) {
		pU8 = image->pu8VirAddr[0];
	} else if (image->enType == IVE_IMAGE_TYPE_S16C1 && mat->type == S16C1) {
		pS16 = image->pu8VirAddr[0];
	} else {
		printf("ERROR: type not matched!\n");
	}

	for (i = 0; i < height; ++i) {
		for (j = 0; j < width; ++j) {
            if (image->enType == IVE_IMAGE_TYPE_U8C1) {
				(mat->ptr + mat->step*i)[j] = *pU8;
				pU8++;
			} else if(image->enType == IVE_IMAGE_TYPE_S16C1){
				(mat->fptr + mat->step*i)[j] = *pS16;
				pS16++;
			}
		}
	}

}

void HIIMat2HIImage(IMat *mat, IVE_IMAGE_S *image) {

	HI_U16 height = image->u16Height;
	HI_U16 width = image->u16Width;
	int i, j;

	HI_U8 *pU8 = NULL;
	HI_S16 *pS16 = NULL;

    if (height != mat->height || width != mat->width) {
		printf("ERROR: size of mat and image are not matched!\n");
		return ;
	}
	if (image->enType == IVE_IMAGE_TYPE_U8C1 && mat->type == U8C1) {
		pU8 = image->pu8VirAddr[0];
	} else if (image->enType == IVE_IMAGE_TYPE_S16C1 && mat->type == S16C1) {
		pS16 = image->pu8VirAddr[0];
	} else {
		printf("ERROR: type not matched!\n");
	}

	for (i = 0; i < height; ++i) {
		for (j = 0; j < width; ++j) {
            if (image->enType == IVE_IMAGE_TYPE_U8C1) {
				*pU8 = (mat->ptr + mat->step*i)[j];
				pU8++;
			} else if(image->enType == IVE_IMAGE_TYPE_S16C1){
				*pS16 = (mat->fptr + mat->step*i)[j];
				pS16++;
			}
		}
	}
}

void printHiImage(IVE_IMAGE_S *image, char *path) {
	HI_U16 height = image->u16Height;
	HI_U16 width = image->u16Width;
	HI_U8 *pU8 = image->pu8VirAddr[0];
	int i, j;
	FILE *fp = fopen(path, "w");
	for (i = 0; i < height; ++i) {
		for (j = 0; j < width; ++j) {
			fprintf(fp, "%d ", *pU8);
			pU8++;
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

void printHiImage16(IVE_IMAGE_S *image, char *path) {
	HI_U16 height = image->u16Height;
	HI_U16 width = image->u16Width;
	HI_U16 *pU16 = image->pu8VirAddr[0];
	int i, j;
	FILE *fp = fopen(path, "w");
	for (i = 0; i < height; ++i) {
		for (j = 0; j < width; ++j) {
			fprintf(fp, "%d ", *pU16);
			pU16++;
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

void readYUV(IVE_SRC_IMAGE_S *image, char *path) {
	FILE * fp= fopen(path, "rb");
	SAMPLE_COMM_IVE_ReadFile(image, fp);
}

void saveYUV(IVE_SRC_IMAGE_S *image, char *path) {
	FILE * fp= fopen(path, "wb");
	SAMPLE_COMM_IVE_WriteFile(image, fp);
}
void doSobel(IMat *mat) {
	IVE_HANDLE IveHandle;
	IVE_THRESH_CTRL_S stThrCtrl = {IVE_THRESH_MODE_BINARY, 0, 0, 0, 0, 255};
	HI_BOOL bInstant = HI_TRUE;
	HI_S32 s32Ret;

	IVE_SRC_IMAGE_S *image;
	s32Ret = SAMPLE_COMM_IVE_CreateImage(image, IVE_IMAGE_TYPE_U8C1, mat->width, mat->height);
	if (s32Ret != HI_SUCCESS) {
		SAMPLE_PRT("SAMPLE_COMM_IVE_CreateImage,Error(%#x)\n", s32Ret);
		return;
	}
	HIIMat2HIImage(mat, image);

	IVE_SOBEL_CTRL_S stSobelCtrl;
	HI_S8 as8Mask[25] = {0, 0, 0, 0, 0,
						 0, -1, 0, 1, 0,
						 0, -2, 0, 2, 0,
						 0, -1, 0, 1, 0,
						 0, 0, 0, 0, 0
	};

	IVE_SRC_IMAGE_S DstH, DstV;
	s32Ret = SAMPLE_COMM_IVE_CreateImage((&DstH), IVE_IMAGE_TYPE_S16C1, image->u16Width, image->u16Height);
	s32Ret = SAMPLE_COMM_IVE_CreateImage((&DstV), IVE_IMAGE_TYPE_S16C1, image->u16Width, image->u16Height);

	stSobelCtrl.enOutCtrl = IVE_SOBEL_OUT_CTRL_BOTH;
	memcpy(stSobelCtrl.as8Mask, as8Mask, 25);

	s32Ret = HI_MPI_SYS_MmzFlushCache(image->u32PhyAddr[0], image->pu8VirAddr[0],
									  image->u16Stride[0] * image->u16Height);
	if (s32Ret != HI_SUCCESS)
	{
		SAMPLE_PRT("HI_MPI_SYS_MmzFlushCache fail,Error(%#x)\n", s32Ret);
		return;
	}
	s32Ret = HI_MPI_IVE_Sobel(&IveHandle, &image, &DstH, &DstV, &stSobelCtrl, bInstant);

	if (s32Ret != HI_SUCCESS) {
		SAMPLE_PRT("HI_MPI_IVE_Sobel fail,Error(%#x)\n", s32Ret);
		return;
	}

	printHiImage16(&DstV, "image1.txt");
	printHiImage16(&DstH, "image2.txt");
}

