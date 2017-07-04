//
// Created by xianb on 2017/3/31.
//

#ifndef CPPPRIMER_WEIGHTEDQUICKUNIONWITHPC_H
#define CPPPRIMER_WEIGHTEDQUICKUNIONWITHPC_H
#include <iostream>
#include <fstream>
#include <vector>
#include <map>

class WeightedQuickUnionWithPC {
private:
    std::vector<int> id;
    std::vector<int> size;
    int count = 0;
    std::map<int, std::vector<int>> unions;
    std::vector<int> graph;

    int find(int p, std::vector<int> &a) {
        while(id[p] != p) {
            a.push_back(p);
            p = id[p];
        }
        a.push_back(p);
        return p;
    }
    int find(int p) {
        while(id[p] != p) {
            p = id[p];
        }
        return p;
    }
    /*
     * 判断两个点是否连通
     * */
    bool connected(int p, int q) {
        return find(p) == find(q);
    }


public:
    WeightedQuickUnionWithPC(int N, std::vector<int> graph) {
        this->graph = graph;
        count = N;
        /*
         * 初始状态下所有的像素点各代表一个不同的ID，而每个ID上的size都是1.
         * */
        for(int i = 0; i < N; i++) {
            id.push_back(i);
            size.push_back(1);
        }
    }

    /*
     * 图像中的count计数，每连通一次就减少一个值。
     * */
    int counts() {
        return count;
    }

    /*
     * 将两个连通域进行合并，主要时将size和id进行设置。
     * 这种情况下所有在同一个连通域中的点的ID值都是一样的。
     * */
    void unionThem(int p, int q) {
        std::vector<int> sitesp;
        std::vector<int> sitesq;
        int pRoot = find(p, sitesp);
        int qRoot = find(q, sitesq);
        if(qRoot == pRoot)
            return;
        if(size[qRoot] > size[pRoot]) {
            id[pRoot] = qRoot;
            for(int i:sitesp) {
                id[i] = qRoot;
            }
            size[qRoot] = size[pRoot] + size[qRoot];
        } else {
            id[qRoot] = pRoot;
            for(int i:sitesq) {
                id[i] = pRoot;
            }
            size[pRoot] = size[qRoot] + size[pRoot];
        }
        count--;
    }

    void printId() {
        for(int i:id) {
            std::cout << i << " ";
        }
        std::cout << std::endl;
    }
    void printSize() {
        for(int i:size) {
            std::cout << i << " ";
        }
        std::cout << std::endl;
    }
    /*
     * 统计有多少个size不为1的点。
     * */
//    int howMany() {
//        int ct = 0;
//        for (int i = 0; i < size.size(); ++i) {
//            if (size[i] != 1)
//                ct++;
//        }
//        return ct;
//    }

    std::vector<std::vector<int>> getUnions(){
        std::vector<std::vector<int>> uns;
        for (const auto &u : unions) {
            std::vector<int> tmp;
            for (const auto &p : u.second) {
                tmp.push_back(p);
            }
            uns.push_back(tmp);
        }
        return uns;
    }
    /*
     * 将所有的在同一个连通域中的点放到同一个vector中, 所以不同的vector的个数就是连通域的个数。
     *
     * */

    void putThemOnVectors() {
        for (int i = 0; i < id.size(); ++i) {
            /*
             * 若一个点的id值不等于本身(表示该点属于一个union中, 这种情况下会忽略掉root节点)
             * size[i]不等于1(将union的根节点加入，不要忽略根节点,将上面被忽略的root节点加入，
             * 若只判断这个条件，那么最底层叶子的size是1，这些点会被忽略，因此前两个条件缺一不可)
             * size[i] == 1 && graph[i] != 0，一个点的size是1但是该点在图中不是0值，表示这是一个独立的单点union
             * */
            if((id[i] != i || size[i] != 1) || (size[i] == 1 && graph[i] != 0))
                unions[find(i)].push_back(i);
        }
    }

//    void SaveEachVectors() {
//        std::ofstream out("allpoints.txt");
//        for(auto &v : unions) {
//            for(auto i:v.second)
//                out << i << " ";
//        }
//        out.close();
//    }
};

#endif //CPPPRIMER_WEIGHTEDQUICKUNIONWITHPC_H

