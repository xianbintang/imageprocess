#include "types.h"
#include <math.h>
#include "TemplateMatch.h"
#include <stdlib.h>
#include <cstring>

short saturate_cast(float v)
{
    int iv = round(v);
    return saturate_cast(iv);
}

int MatAt(IMat  mat, int y, int x) {
    return mat.ptr[y * mat.step + x];
}

//int RotateGrayImage(IMat matSrc, IMat *matDst, const double degree);
//int main(int argc, char **argv)
//{
//
//    IMat matSrc = Iimread(argv[1]);
//    IMat matDst;
//    RotateGrayImage(matSrc, matDst, 50);
//    Iimshow("rotated", matDst);
//
//    cvWaitKey(-1);
//}
int RotateGrayImage(IMat matSrc, IMat **matDst, const double degree, IPoint *rect)
{
//     IMat matSrc = Iimread("images//KOYO.bmp", 2 | 4);

//    if (matSrc.empty())
//        return 0;

//    const double degree = 90;
    double angle = degree * M_PI / 180.;
    double alpha = cos(angle);
    double beta = sin(angle);
    int iWidth = matSrc.width;
    int iHeight = matSrc.height;
    int iNewWidth = round(iWidth * fabs(alpha) + iHeight * fabs(beta));
    int iNewHeight = round(iHeight * fabs(alpha) + iWidth * fabs(beta));

    if(iNewHeight %2==1) {
        iNewHeight--;
    }
    if (iNewWidth%2 == 1) {
        iNewWidth--;
    }


    double m[6];
    m[0] = alpha;
    m[1] = beta;
    m[2] = (1 - alpha) * (iWidth / 2.) - beta * (iHeight / 2.);
    m[3] = -m[1];
    m[4] = m[0];
    m[5] = beta * (iWidth / 2.) + (1 - alpha) * (iHeight / 2.);

    /* 将图片移到中心, opencv的仿射变换是不包含这个的, 这个在博客的作者回答读者问中踢到怎么做的 */
    /*
     * 原始的旋转是以图像的中心来旋转的，所以可能会导致一些内容被移除图像外，但是通过下面的平移后可以将其
     * 移到图像中心。
     * 具体挪多少可能需要算，这样挪好像还是会出现超出范围的问题。
     * 最简单的办法就是先直接平移将重心移到图片的重心，这样可能不会超出范围，但是这也不确定。
     * */
    m[2] += (iNewWidth - iWidth) / 2.0f;	m[5] += (iNewHeight - iHeight) / 2.0f;


    int plx = rect[0].x, ply=rect[0].y, prx=rect[3].x, pry=rect[3].y;
    int plxb = rect[1].x, plyb=rect[1].y, prxb=rect[2].x, pryb=rect[2].y;

    rect[0].y = m[0] * plx + m[1] * ply + m[2];
    rect[0].x = m[3] * plx + m[4] * ply + m[5];

    rect[1].y = m[0] * plxb + m[1] * plyb + m[2];
    rect[1].x = m[3] * plxb + m[4] * plyb + m[5];

    rect[2].y = m[0] * prxb + m[1] * pryb + m[2];
    rect[2].x = m[3] * prxb + m[4] * pryb + m[5];

    rect[3].y = m[0] * prx + m[1] * pry + m[2];
    rect[3].x = m[3] * prx + m[4] * pry + m[5];

//    IMat M = IMat(2, 3, CV_64F, m);
/*    IMat *tmpM = ICreateMat(2, 3, U8C1);
    int i, j;
    for (i = 0; i < tmpM->height; i++) {
        for (j = 0; j < tmpM->width ; ++j) {
            (tmpM->ptr + tmpM->step*i)[j-1] = m[i * tmpM->height + j];
        }
    }
    IMat M = *tmpM;*/
//    IMat matDst1;
//            = IMat(iNewWidth, iNewHeight, U8C1);

    *matDst = (IMat *)malloc(sizeof(IMat));
    (*matDst)->width = iNewWidth;
    (*matDst)->height = iNewHeight;
    (*matDst)->step = iNewWidth;
    (*matDst)->type = U8C1;
    (*matDst)->ptr = (int *)malloc(iNewHeight * iNewWidth * sizeof(int));
    /*
     * 旋转后的图像的边缘可能会导致两条白色边缘直线
     *
     * */
    memset((*matDst)->ptr, 0, iNewHeight * iNewWidth * sizeof(int));
    /*
     *  这其实是在求逆。因为按照您文中的alpha，beta，M的公式，计算出来的M，
     *  应该去左乘源图像的坐标，得到的结果是目标图像的坐标。但是后面的程序可以看出来，
     *  其实是用M（现在已经不是公式计算出来的M了）左乘目标图像坐标，得出源图像坐标，再进行插值。
     *  所以那一段代码，就是在进行一个“求逆”，首先把M分开成A，t，然后两边同乘以A逆，就可以得到“新的M”。
     * */
    double D = m[0]*m[4] - m[1]*m[3];
    D = D != 0 ? 1./D : 0;
    double A11 = m[4]*D, A22 = m[0]*D;
    m[0] = A11; m[1] *= -D;
    m[3] *= -D; m[4] = A22;
    double b1 = -m[0]*m[2] - m[1]*m[5];
    double b2 = -m[3]*m[2] - m[4]*m[5];
    m[2] = b1; m[5] = b2;

    for (int y=0; y<iNewHeight + 5; ++y)
    {
        for (int x=0; x<iNewWidth + 5; ++x)
        {
//int tmpx = cvFloor(m[0] * x + m[1] * y + m[2]);
//int tmpy = cvFloor(m[3] * x + m[4] * y + m[5]);

            /* x 和 y是目标图像中的坐标 */
            /* sy和sx是原始图像中的坐标, 在这里用逆矩阵得到(x, y)对应的原图像中的坐标点，获得其周围的点的像素值，从而才能进行插值运算 */
            float fx = m[0] * x + m[1] * y + m[2];
            float fy = m[3] * x + m[4] * y + m[5];

            int sy  = floor(fy);
            fy -= sy;
//sy = std::min(sy, iHeight-2);
//sy = std::max(0, sy);
            if (sy < 0 || sy >= iHeight) continue;

            /* cbufy和cbufx是用来求原始图像中对应点的像素值的一个中间值，具体什么意思要看公式 */
            short cbufy[2];
            /* 这里的2048是为了将浮点数转换成整数从而方便运算，加速运算的，和后面的移位22互相抵消了 */
//            cbufy[0] = saturate_cast((1.f - fy) * 2048);
            cbufy[0] = (1.f - fy) * 2048;
            cbufy[1] = 2048 - cbufy[0];

            int sx = floor(fx);
            fx -= sx;
//if (sx < 0) {
//	fx = 0, sx = 0;
//}
//if (sx >= iWidth - 1) {
//	fx = 0, sx = iWidth - 2;
//}
            if (sx < 0 || sx >= iWidth) continue;

            short cbufx[2];
//            cbufx[0] = saturate_cast((1.f - fx) * 2048);
            cbufx[0] = (1.f - fx) * 2048;
            cbufx[1] = 2048 - cbufx[0];

            for (int k=0; k < 1; ++k)
            {
                /*
                 * 最下面的和最右边的不做处理,所以会导致两条白线
                 * */
                if (sy == iHeight - 1 || sx == iWidth - 1) {
                    continue;
                } else {
                    ((*matDst)->ptr + y * (*matDst)->step)[x] = (MatAt(matSrc,sy, sx) * cbufx[0] * cbufy[0] +
                                                             MatAt(matSrc,sy+1, sx) * cbufx[0] * cbufy[1] +
                                                             MatAt(matSrc,sy, sx+1) * cbufx[1] * cbufy[0] +
                                                             MatAt(matSrc,sy+1, sx+1) * cbufx[1] * cbufy[1]) >> 22;
                }
            }
        }
    }
//    *matDst = matDst1;
//    Iimwrite("rotate_bilinear_1.jpg", matDst1);
//    Iimshow("we", matDst1);
    return 0;
}

