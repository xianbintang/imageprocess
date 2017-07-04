//
// Created by xianb on 2017/5/22.
//

#include <TemplateMatch.h>
#include <iostream>


void HIIMage2IMat(IVE_IMAGE_S *image, IMat *mat);
void HIIMat2HIImage(IMat *mat, IVE_IMAGE_S *image);

int main(int argc, char ** argv) {
    IPoint point;
    double degree;
//    DrawContours(dst, tpl, CvScalar(255,255,255), 1);
    TemplateMatch tpl;
    IVE_SRC_IMAGE_S imageSrc;
    SAMPLE_COMM_IVE_CheckIveMpiInit();

    SAMPLE_COMM_IVE_CreateImageByCached(&imageSrc, IVE_IMAGE_TYPE_U8C1, 135, 35);
    readYUV(imageSrc, argv[1]);

    IMat *src = ICreateMat(imageSrc.u16Height, imageSrc.u16Width, U8C1);
    HIIMage2IMat(&imageSrc, src);
    CreateGeoMatchModel(&tpl, src, 30, 170);


    IVE_SRC_IMAGE_S imageDst;
    SAMPLE_COMM_IVE_CreateImageByCached(&imageDst, IVE_IMAGE_TYPE_U8C1, 320, 240);
    readYUV(&imageDst, argv[2]);

    IMat *searchIMat= ICreateMat(imageDst.u16Height, imageDst.u16Width, U8C1);
    HIIMage2IMat(&imageDst, searchIMat);

    double score = FindGeoMatchModel(&tpl,searchIMat, 0.995, 0.99999995, &point, 30, &degree);
//    std::cout << " found at degree: " << degree << " Points: (" << point.x << "," << point.y << ")" << " score: " << score
//              << " Time used: " << tt.duration() << std::endl;
    return 0;
}
