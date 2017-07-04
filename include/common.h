//
// Created by xianb on 2017/4/19.
//

#ifndef TESTOPENCV_COMMON_H
#define TESTOPENCV_COMMON_H

#include <vector>

void displayTextImage(std::vector<std::vector<int>> TextImg, int width, int height);

std::vector<std::vector<int>> GetUnions(std::string input, int width, int height);
#endif //TESTOPENCV_COMMON_H
