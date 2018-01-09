//
// Created by xianb on 2018/1/5.
//

#include <iostream>
#include "../src/create_template_pc/contour_detection.h"

int main(void)
{
    std::cout << "before" << std::endl;
    Koyo_Tool_Contour_Parameter koyo_tool_contour_parameter;
    int a;
    create_template(NULL, koyo_tool_contour_parameter, &a);
    std::cout << "after" << std::endl;
    std::cout << a << std::endl;
    return 0;
}

