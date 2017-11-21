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
    int accessedTimes = 0;
    std::map<int, std::vector<int>> mp;

public:
    WeightedQuickUnionWithPC(int N) {
        count = N;
        for(int i = 0; i < N; i++) {
            id.push_back(i);
            size.push_back(1);
        }
    }
    bool connected(int p, int q) {
        return find(p) == find(q);
    }

    int counts() {
        return count;
    }

    int find(int p, std::vector<int> &a) {
        while(id[p] != p) {
            a.push_back(p);
            p = id[p];
            accessedTimes += 2;
        }
        a.push_back(p);
        accessedTimes++;
        return p;
    }
    int find(int p) {
        while(id[p] != p) {
            p = id[p];
            accessedTimes += 2;
        }
        accessedTimes++;
        return p;
    }

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
        accessedTimes++;
        count--;
    }

    int getAccessedTimes() {
        return accessedTimes;
    }

    void setAccessedTimes() {
        accessedTimes = 0;
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
    int howMany() {
        int ct = 0;
        for (int i = 0; i < size.size(); ++i) {
            if (size[i] != 1)
                ct++;
        }
        return ct;
    }
    int putThemOnVectors() {
        for (int i = 0; i < id.size(); ++i) {
            if(id[i] != i)
                mp[id[i]].push_back(i);
        }
        return mp.size();
    }

    void SaveEachVectors() {
        std::ofstream out("E:\\ClionProjects\\CppPrimer\\allpoints.txt");
        for(auto &v : mp) {
            for(auto i:v.second)
                out << i << " ";
        }
        out.close();
    }
};

#endif //CPPPRIMER_WEIGHTEDQUICKUNIONWITHPC_H

