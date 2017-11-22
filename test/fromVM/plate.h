#ifndef _PLATE_H_
#define _PLATE_H_

#include "opencv/cv.h"
#include "opencv/highgui.h"
#include "List.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define SCALE_MAX 4.2
#define SCALE_MIN 2.8
#define AREA_MAX 28000
#define AREA_MIN 7000
#define FILENAME_LEN 50

#define BIN_MIN 0.75

/*接下来这几个宏定义是在车牌的左右去边框的时候用到的*/
#define RESIZED_HEIGHT 200
#define AREA_MIN_1 1700
#define AREA_MAX_1 2700
#define SCALE_MAX_1 0.17
#define SCALE_MIN_1 0.1

void preprocess_car_img(IplImage * car_img);
List get_location(IplImage * img_car, IplImage * org_img_car);
void get_plate_image(IplImage * img_car, List rects);
void resize_image(IplImage * img_to_resize, IplImage * img_after_resize, float scale);
void preprocess_plate_image(IplImage * img_after_resize);

void get_contour_rect(IplImage * src_img, List  rects, CvMemStorage * storage, CvSeq * contours);
void draw_contour_rect(IplImage * src_img, List  rects);
void filter_rect_by_shape(List src_rects, List dst_rects);
void print_area_of_rect(CvRect rect);


void remove_border_ul(IplImage * img_plate);
void get_character(IplImage * img);

int character_recognizing(IplImage * img_char);

#endif