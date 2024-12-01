#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <limits>
#include <map>
#include <algorithm>
#include <iomanip>
using namespace std;
const vector<int> input_size = {100, 500, 1000, 2000, 4000};
const vector<int> input_dim  = {3, 4, 5, 6, 7, 8};

struct TestData {
    double itr;
    double n;
    double num;
};
void re(string name, map<int, map<int, vector<TestData>>>& MapC, map<int, map<int, vector<TestData>>>& MapA,
        map<int, map<int, vector<TestData>>>& MapI) {
    ifstream fileT(name);
    string line;
    int lineCount = 0;
    int setCount  = 0;
    while (getline(fileT, line)) {
        getline(fileT, line);
        TestData data;
        istringstream iss(line);
        string value;
        getline(iss, value, ':');
        getline(iss, value, ':');
        data.itr = stoi(value.substr(0, value.find_last_of(",")));
        getline(iss, value, ':');
        data.n = stoi(value.substr(0, value.find_last_of(",")));
        getline(iss, value, ':');
        data.num = stoi(value);
        int d    = input_dim[lineCount / 15];
        int s    = input_size[setCount / 3];
        if (setCount % 3 == 0) {
            vector<TestData>& t = MapC[d][s];
            t.push_back(data);
        }
        else if (setCount % 3 == 1) {
            vector<TestData>& t = MapA[d][s];
            t.push_back(data);
        }
        else if (setCount % 3 == 2) {
            vector<TestData>& t = MapI[d][s];
            t.push_back(data);
        }
        getline(fileT, line);
        lineCount++;
        setCount++;
        if (lineCount >= 15 && lineCount % 15 == 0) {
            setCount = 0;
        }
        if (lineCount >= 90 && lineCount % 90 == 0) {
            lineCount = 0;
        }
    }
    fileT.close();
}

int main() {
    int Testsize = 20;
    map<int, map<int, vector<TestData>>> MapC;
    map<int, map<int, vector<TestData>>> MapA;
    map<int, map<int, vector<TestData>>> MapI;
    string name = "CostNumTest.txt";
    re(name, MapC, MapA, MapI);
    map<int, map<int, vector<TestData>>> MapC1;
    map<int, map<int, vector<TestData>>> MapA1;
    map<int, map<int, vector<TestData>>> MapI1;
    string name1 = "CostNumTest-1.txt";
    re(name1, MapC1, MapA1, MapI1);
    map<int, map<int, vector<TestData>>> MapC12;
    map<int, map<int, vector<TestData>>> MapA12;
    map<int, map<int, vector<TestData>>> MapI12;
    string name12 = "CostNumTest-12.txt";
    re(name12, MapC12, MapA12, MapI12);
    map<int, map<int, vector<TestData>>> MapC121;
    map<int, map<int, vector<TestData>>> MapA121;
    map<int, map<int, vector<TestData>>> MapI121;
    string name121 = "CostNumTest-12-1.txt";
    re(name121, MapC121, MapA121, MapI121);
    map<int, map<int, vector<double>>> MapCt;
    map<int, map<int, vector<double>>> MapAt;
    map<int, map<int, vector<double>>> MapIt;
    for (size_t d = 0; d < input_dim.size(); d++) {
        int dim = input_dim[d];
        for (size_t tn = 0; tn < input_size.size(); tn++) {
            int len = input_size[tn];
            if (!(len == 1000 || dim == 4)) {
                continue;
            }
            double itrC = 0, nC = 0, numC = 0;
            double itrA = 0, nA = 0, numA = 0;
            double itrI = 0, nI = 0, numI = 0;
            vector<TestData>& Ct = MapC1[dim][len];
            vector<TestData>& At = MapA1[dim][len];
            vector<TestData>& It = MapI1[dim][len];
            for (size_t t = 0; t < Testsize; t++) {
                itrC += MapC[dim][len][t].itr - Ct[t].itr;
                nC += Ct[t].n - MapC[dim][len][t].n;
                numC += Ct[t].num / MapC[dim][len][t].num;
                itrA += MapA[dim][len][t].itr - At[t].itr;
                nA += At[t].n - MapA[dim][len][t].n;
                numA += At[t].num / MapA[dim][len][t].num;
                itrI += MapI[dim][len][t].itr - It[t].itr;
                nI += It[t].n - MapI[dim][len][t].n;
                numI += It[t].num / MapI[dim][len][t].num;
                // itrC += Ct[t].itr;
                // nC += Ct[t].n;
                // numC += Ct[t].num;
                // itrA += At[t].itr;
                // nA += At[t].n;
                // numA += At[t].num;
                // itrI += It[t].itr;
                // nI += It[t].n;
                // numI += It[t].num;
            }
            vector<double> nPairC(3);
            vector<double> nPairA(3);
            vector<double> nPairI(3);
            nPairC[0]       = itrC / Testsize;
            nPairC[1]       = nC / Testsize;
            nPairC[2]       = numC / Testsize;
            nPairA[0]       = itrA / Testsize;
            nPairA[1]       = nA / Testsize;
            nPairA[2]       = numA / Testsize;
            nPairI[0]       = itrI / Testsize;
            nPairI[1]       = nI / Testsize;
            nPairI[2]       = numI / Testsize;
            MapCt[dim][len] = nPairC;
            MapAt[dim][len] = nPairA;
            MapIt[dim][len] = nPairI;
        }
    }
    cout << "CORR" << endl;
    cout << "d=4" << endl;
    for (size_t i = 0; i < 3; i++)
    {
        for (auto& nPair : MapCt[4]) {
            cout << fixed << setprecision(2) << nPair.second[i] << ",";
        }
        cout << endl;
    }
    cout << "n=1000" << endl;
    for (size_t i = 0; i < 3; i++)
    {
        for (auto& itrPair : MapCt) {
            cout << fixed << setprecision(2) << itrPair.second[1000][i] << ",";
        }
        cout << endl;
    }
    cout << "========================" << endl;
    cout << "ANTI" << endl;
    cout << "d=4" << endl;
    for (size_t i = 0; i < 3; i++)
    {
        for (auto& nPair : MapAt[4]) {
            cout << fixed << setprecision(2) << nPair.second[i] << ",";
        }
        cout << endl;
    }
    cout << "n=1000" << endl;
    for (size_t i = 0; i < 3; i++)
    {
        for (auto& itrPair : MapAt) {
            cout << fixed << setprecision(2) << itrPair.second[1000][i] << ",";
        }
        cout << endl;
    }
    cout << "========================" << endl;
    cout << "INDE" << endl;
    cout << "d=4" << endl;
    for (size_t i = 0; i < 3; i++)
    {
        for (auto& nPair : MapIt[4]) {
            cout << fixed << setprecision(2) << nPair.second[i] << ",";
        }
        cout << endl;
    }
    cout << "n=1000" << endl;
    for (size_t i = 0; i < 3; i++)
    {
        for (auto& itrPair : MapIt) {
            cout << fixed << setprecision(2) << itrPair.second[1000][i] << ",";
        }
        cout << endl;
    }
    cout << "========================" << endl;
    // for (auto& nPair : MapCt[4]) {
    //     cout << nPair.first << endl;
    //     for (double t : nPair.second) {
    //         cout << t << ",";
    //     }
    //     cout << endl;
    // }
    // for (auto& itrPair : MapCt) {
    //     cout << itrPair.first << endl;
    //     vector<double> nPair = itrPair.second[1000];
    //     for (double t : nPair) {
    //         cout << t << ",";
    //     }
    //     cout << endl;
    // }
    // cout << "========================" << endl;
    // for (auto& nPair : MapAt[4]) {
    //     cout << nPair.first << endl;
    //     for (double t : nPair.second) {
    //         cout << t << ",";
    //     }
    //     cout << endl;
    // }
    // for (auto& itrPair : MapAt) {
    //     cout << itrPair.first << endl;
    //     vector<double> nPair = itrPair.second[1000];
    //     for (double t : nPair) {
    //         cout << t << ",";
    //     }
    //     cout << endl;
    // }
    // cout << "========================" << endl;
    // for (auto& nPair : MapIt[4]) {
    //     cout << nPair.first << endl;
    //     for (double t : nPair.second) {
    //         cout << t << ",";
    //     }
    //     cout << endl;
    // }
    // for (auto& itrPair : MapIt) {
    //     cout << itrPair.first << endl;
    //     vector<double> nPair = itrPair.second[1000];
    //     for (double t : nPair) {
    //         cout << t << ",";
    //     }
    //     cout << endl;
    // }
    // for (auto& itrPair : MapC1) {
    //     cout << itrPair.first << endl;
    //     for (auto& nPair : itrPair.second) {
    //         cout << nPair.first << endl;
    //         for (TestData t : nPair.second) {
    //             cout << "itr:" << t.itr << ",n:" << t.n << ",num:" << t.num << endl;
    //         }
    //     }
    // }

    return 0;
}