#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "types.h"
#include "TemplateMatch.h"
#include <cmath>


//allocate memory for doubel matrix
void CreateDoubleMatrix(double ***matrix,ISize size)
{
//    matrix = new double*[size.height];
//    for(int iInd = 0; iInd < size.height; iInd++)
//        matrix[iInd] = new double[size.width];
    (*matrix) = (double **)malloc(size.height * sizeof(double *));
    int i;
    for (i = 0; i < size.height; ++i) {
        (*matrix)[i] = (double *)malloc(size.width * sizeof(double));
    }
}
void ResumeTpl(TemplateMatch *tpl)
{
    tpl->noOfCordinates = 433;
    int i;
    int x, y;
    double d;
    FILE *fp1 = fopen("cordinates.txt", "r");
    FILE *fp2 = fopen("edgeMag.txt", "r");
    FILE *fp3 = fopen("ex.txt", "r");
    FILE *fp4 = fopen("ey.txt", "r");
    tpl->modelWidth = 140;
    tpl->modelHeight= 42;

    for(i = 0; i < tpl->noOfCordinates; i++) {
        fscanf(fp1, "%d", &x);
        fscanf(fp1, "%d", &y);
        tpl->cordinates[i].x = x;
        tpl->cordinates[i].y = y;
        fscanf(fp2, "%lf", &d);
        tpl->edgeMagnitude[i] = d;
        fscanf(fp3, "%lf", &d);
        tpl->edgeDerivativeX[i] = d;
        fscanf(fp4, "%lf", &d);
        tpl->edgeDerivativeY[i] = d;

    }
    tpl->centerOfGravity.x = 20;
    tpl->centerOfGravity.y = 71;


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


int main(int argc, char ** argv)
{
    TemplateMatch tpl;

    tpl.cordinates = (IPoint *)malloc(500 * sizeof(IPoint));
    tpl.edgeMagnitude = (double *)malloc(500 * sizeof(double));		//Allocate memory for edge magnitude for selected points
    tpl.edgeDerivativeX = 	(double *)malloc(500 * sizeof(double));		//Allocate memory for edge X derivative for selected points
    tpl.edgeDerivativeY = (double *)malloc(500 * sizeof(double));		////Allocate memory for edge Y derivative for selected points

    ResumeTpl(&tpl);

    std::cout << tpl.noOfCordinates <<std::endl;
    std::cout << tpl.modelWidth<<std::endl;
    std::cout << tpl.modelHeight<<std::endl;
    std::cout << tpl.centerOfGravity.x<<std::endl;
    std::cout << tpl.centerOfGravity.y<<std::endl;



    double partialSum;
    double partialScore;
    int curX, curY;
    int i,j,m;
    const short *_Sdx;
    const short *_Sdy;
    double iTx, iTy, iSx, iSy;
    IMat *Sdx = 0, *Sdy = 0;
    int sumOfCoords = 0;

    ISize Ssize;
    Ssize.width = 320;
    Ssize.height = 240;

   double **edgeX;  //Gradient magnitude matrix
    double **edgeY;  //Gradient magnitude matrix
    CreateDoubleMatrix(&edgeX, Ssize); // create image to save gradient magnitude  values
    CreateDoubleMatrix(&edgeY, Ssize); // create image to save gradient magnitude  values

    Sdx = fICreateMat(Ssize.height, Ssize.width, S16C1); // X derivatives
    Sdy = fICreateMat(Ssize.height, Ssize.width, S16C1); // y derivatives

    double minScore = 0.95;
    double greediness = 0.999995;
    double resultScore = 0;
    IPoint resultPoint;


    resumeIMatf("Sdx.txt", Sdx);
    resumeIMatf("Sdy.txt", Sdy);

    for (i = 0; i < Ssize.height; i++) {
        _Sdx = (short *) (Sdx->fptr + Sdx->step * (i));
        _Sdy = (short *) (Sdy->fptr + Sdy->step * (i));

        for (j = 0; j < Ssize.width; j++) {

            iSx = _Sdx[j];  // X derivative of Source image
            iSy = _Sdy[j];  // Y derivative of Source image
            double vector_length = sqrt(_Sdx[j] * _Sdx[j] + _Sdy[j] * _Sdy[j]);
            edgeX[i][j] = _Sdx[j] / vector_length;
            edgeY[i][j] = _Sdy[j] / vector_length;
        }
    }



    double normMinScore = minScore / tpl.noOfCordinates; // precompute minumum score
    double normGreediness =
            ((1 - greediness * minScore) / (1 - greediness)) / tpl.noOfCordinates; // precompute greedniness

    for (i = 0; i < 240; i++) {
        for (j = 0; j < 320; j++) {
            partialSum = 0; // initilize partialSum measure
            int count = 0;
            for (m = 0; m < tpl.noOfCordinates; m++) {
                count++;
                curX = i + tpl.cordinates[m].x;    // template X coordinate
                curY = j + tpl.cordinates[m].y; // template Y coordinate
                iTx = tpl.edgeDerivativeX[m];    // template X derivative
                iTy = tpl.edgeDerivativeY[m];    // template Y derivative

                if (curX < 0 || curY < 0 || curX > Ssize.height - 1 || curY > Ssize.width - 1)
                    continue;

                _Sdx = (short *) (Sdx->fptr + Sdx->step * (curX));
                _Sdy = (short *) (Sdy->fptr + Sdy->step * (curX));

                iSx = _Sdx[curY]; // get curresponding  X derivative from source image
                iSy = _Sdy[curY];// get curresponding  Y derivative from source image

                if ((iSx != 0 || iSy != 0) && (iTx != 0 || iTy != 0)) {
                    //partial Sum  = Sum of(((Source X derivative* Template X drivative) + Source Y derivative * Template Y derivative)) / Edge magnitude of(Template)* edge magnitude of(Source))
//                        partialSum =
//                                partialSum + ((iSx * iTx) + (iSy * iTy)) * (tpl->edgeMagnitude[m] * matGradMag[curX][curY]);

                    partialSum =
                            partialSum + ((edgeX[curX][curY]* iTx) + (edgeY[curX][curY]* iTy));
                }

                sumOfCoords = m + 1;
                partialScore = partialSum / sumOfCoords;
                // check termination criteria
                // if partial score score is less than the score than needed to make the required score at that position
                // break serching at that coordinate.
                if (partialScore < (MIN((minScore - 1) + normGreediness * sumOfCoords, normMinScore * sumOfCoords)))
                    break;

            }
//                std::cout << "no of cordinates computed: " << count << std::endl;
            if (partialScore > resultScore) {
                resultScore = partialScore; //  Match score
                resultPoint.x = i;            // result coordinate X
                resultPoint.y = j;            // result coordinate Y
            }
        }
    }
    printf("score: %lf, where: (%d, %d)", resultScore, resultPoint.x, resultPoint.y);

}



