//
// Created by xianb on 2017/3/31.
//

#include <opencv/cv.h>
using namespace cv;
using namespace std;
#include "../include/WeightedQuickUnionWithPC.h"
#include "../include/common.h"

std::vector<std::vector<int>> GetUnions(std::string input, int width, int height) {
    /*read graph into array*/
    std::vector<int> graph;
    std::ifstream fin(input.c_str());
    for (int i = 0; i < width * height; ++i) {
        int num = 0;
        fin >> num;
        graph.push_back(num);
    }
    WeightedQuickUnionWithPC uf(width* height, graph);
    /*Traversal the graph*/
    for (int i = 0; i < graph.size(); ++i) {
        if(i < width || i >= width * height - width|| i % width == 0 || i % width == width - 1)
            continue;
        /* is this point already connected with its adjacent points? */
        std::vector<int> adj;

        /*
         * (i - width - 1) (i - width) (i - width + 1)
         * (i - 1       ) (i        ) (i  + 1       )
         * (i + width - 1) (i + width) (i + width + 1)
         *
         * */
        adj.push_back(i - width - 1);
        adj.push_back(i - width);
        adj.push_back(i - width + 1);

        adj.push_back(i - 1);
        adj.push_back(i + 1);

        adj.push_back(i + width - 1);
        adj.push_back(i + width);
        adj.push_back(i + width + 1);

        if (graph[i] != 0) {
            /* 遍历该点的八领域，若相连，则将其进行union操作 */
            for (int a : adj) {
                /* 若领域相连则进行union操作 */
                if(graph[a] != 0) {
                    uf.unionThem(i, a);
                }
            }

        }
    }
    uf.putThemOnVectors();
    return uf.getUnions();
}


#ifdef DEBUG
int main(void)
{
    std::vector<std::vector<int>> unions = GetUnions("ImageAsText.txt", 139, 44);
    std::cout << unions.size() << " components." << std::endl;
    displayTextImage(unions, 139, 44);
    for (const auto &v : unions) {
        std::cout << "component: ";
        for (const auto &p : v) {
            std::cout << p << " ";
        }
        std::cout << std::endl;
    }

    for (int i = 0; i < 10; ++i) {
        for (int j = 0; j < 10; ++j) {
            std::cout << i * 10 + j << " ";
        }
        std::cout << std::endl;
    }
    return 0;
}
#endif
