//
// Created by xianb on 2017/3/31.
//

#include "WeightedQuickUnionWithPC.h"
#include <vector>
#include <fstream>
#include <iostream>

int main(void) {
    /*read graph into array*/
    int width = 320, height = 240;
    std::vector<int> graph;
    std::ifstream fin("E:\\ClionProjects\\CppPrimer\\image.txt");
    std::ofstream fout("E:\\ClionProjects\\CppPrimer\\out.txt");
    for (int i = 0; i < width * height; ++i) {
        int num = 0;
        fin >> num;
        graph.push_back(num);
        fout << num;
    }
    fout.close();
    WeightedQuickUnionWithPC uf(320 * 240);
    /*Traversal the graph*/
    for (int i = 0; i < graph.size(); ++i) {
        if(i < 320|| i >= 76480|| i % 320== 0 || i % 320== 319)
            continue;
        uf.getAccessedTimes();
        /* is this point already connected with its adjacent points? */
        std::vector<int> adj;
        int x = i / width;
        int y = i % width;
        /*
        int xm1 = x - 1;
        int xp1 = x + 1;
        int ym1 = y - 1;
        int yp1 = y + 1;
         */
        /*
         * (x - 1, y - 1) (x - 1, y) (x - 1, y + 1)
         * (x    , y - 1) (x    , y) (x    , y + 1)
         * (x + 1, y - 1) (x + 1, y) (x + 1, y + 1)
         *
         * */
        /*
         * (i - width - 1) (i - width) (i - width + 1)
        /* (i - 1       ) (i        ) (i  + 1       )
        /* (i + width - 1) (i + width) (i + width + 1)
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

//    uf.printId();
//    uf.printSize();
    std::cout << uf.howMany() << " components." << std::endl;
    std::cout << uf.putThemOnVectors() << " vectors." << std::endl;
    uf.SaveEachVectors();

    /*
    while(!StdIn.isEmpty()) {
        uf.setAccessedTimes();
        if(uf.connected(p, q)) {
            uf.printId();
            uf.printSize();
            StdOut.println(p + " "+ q + " array accessed: " + uf.getAccessedTimes() + " times");
            continue;
        }
        uf.unionThem(p, q);
        uf.printId();
        uf.printSize();
        StdOut.println(p + " "+ q + " array accessed: " + uf.getAccessedTimes() + " times");
        StdOut.println(p + " " + q);
    }
    StdOut.println(uf.count() + " components");
     */
}