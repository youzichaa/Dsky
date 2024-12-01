#include <iostream>
#include <vector>
#include <random>
#include <limits>
#include <ctime>
#include "openfhe.h"
#include "cryptocontext-ser.h"
#include <map>
#include <omp.h>
#include <pthread.h>
#include <chrono>
using namespace std;
using namespace lbcrypto;
using namespace std::chrono;

vector<vector<uint64_t>> P;
uint32_t n = 1000, d = 3, dt = 0, dv = 0;
PlaintextModulus plaintextModulus = 65537;
uint32_t ringDim                  = 1 << 13;
usint batchSize                   = 16;
vector<KeyPair<DCRTPoly>> kps;
KeyPair<DCRTPoly> keys;
CryptoContext<DCRTPoly> cc;
TimeVar t;
double processingTime(0.000);
double offline_process, query_process;
vector<int32_t> indexList(1);
vector<uint32_t> NumItr;
uint32_t ts;
uint32_t ItrNum = 0, SkyNum = 0, ComNum = 0;
double Utime = 0, Ucost = 0, Stime = 0, Scost = 0;
const vector<uint32_t> input_size = {100, 500, 1000, 2000, 4000};
const vector<uint32_t> input_dim  = {3, 4, 5, 6, 7, 8};
const string relaxfile            = "./data/relax1.txt";
map<uint32_t, map<uint32_t, vector<int64_t>>> relaxData;

uint64_t* plainT(vector<vector<uint64_t>> P, uint32_t n, uint32_t d) {
    map<uint32_t, vector<uint64_t>> T;
    map<uint32_t, vector<uint64_t>> Result;
    // Process
    clock_t s21 = clock();
    for (uint32_t i = 0; i < n; i++) {
        vector<uint64_t> pt(d);
        for (uint32_t j = 0; j < d; j++) {
            pt[j] = P[j][i];
        }
        T[i] = P[i];
    }
    // Select
    while (!T.empty()) {
        n = T.size();
        vector<vector<uint64_t>> S(n, vector<uint64_t>(2));
        uint32_t i0 = 0;
        for (auto& entry : T) {
            S[i0][0]           = entry.first;
            vector<uint64_t> s = entry.second;
            S[i0][1]           = s[0];
            for (uint32_t j = 1; j < d; j++) {
                S[i0][1] += s[j];
            }
            i0++;
        }
        uint32_t pos   = S[0][0];
        uint64_t STMin = S[0][1];
        for (uint32_t i = 1; i < n; i++) {
            if (S[i][1] < STMin) {
                STMin = S[i][1];
                pos   = S[i][0];
            }
        }
        vector<uint64_t> Pmin = P[pos];
        vector<uint64_t> Tmin = T[pos];
        Result[Result.size()] = Pmin;
        for (auto it = T.begin(); it != T.end();) {
            vector<uint64_t> s = it->second;
            uint32_t sdom;
            vector<uint64_t> leq(d);
            for (uint32_t j = 0; j < d; j++) {
                leq[j] = (Tmin[j] <= s[j]) ? 1 : 0;
            }
            uint64_t LEQ = leq[0];
            for (uint32_t j = 1; j < d; j++) {
                LEQ = LEQ * leq[j];
            }
            uint64_t EQ  = 0;
            uint64_t tmp = 0;
            for (uint32_t j = 0; j < d; j++) {
                tmp = tmp + s[j];
            }
            if (STMin < tmp) {
                EQ = 1;
            }
            sdom = LEQ * EQ;
            if (sdom == 1) {
                it = T.erase(it);
            }
            else {
                ++it;
            }
        }
        T.erase(pos);
    }
    uint32_t k = Result.size();
    cout << k << endl;
    uint64_t* resD = new uint64_t[k * d];
    for (uint32_t i = 0; i < k; i++) {
        for (uint32_t j = 0; j < d; j++) {
            resD[i * d + j] = Result[i][j];
        }
    }
    clock_t s22 = clock();
    double s2   = double(s22 - s21) * 1000.0 / CLOCKS_PER_SEC;
    cout << "\nPlain time: " << s2 << "ms" << endl;
    return resD;
}

struct ElementIndex {
    uint64_t element;
    uint32_t index;

    ElementIndex(uint64_t e, uint32_t i) : element(e), index(i) {}

    bool operator<(const ElementIndex& other) const {
        return element < other.element;
    }
};

struct IndexElement {
    uint64_t element;
    uint32_t index;
    uint32_t visit;

    IndexElement() : element(0), index(0), visit(0) {}
};

void loadP(const string& filename) {
    P.resize(n, vector<uint64_t>(d, 0));
    ifstream in(filename);
    string line;
    getline(in, line);
    uint32_t i = 0;
    while (getline(in, line)) {
        istringstream iss(line);
        string token;
        uint32_t j = 0;
        while (getline(iss, token, ',')) {
            P[i][j] = stoull(token);
            ++j;
        }
        ++i;
    }
    in.close();
    ifstream inr(relaxfile);
    for (size_t d = 0; d < input_dim.size(); d++) {
        uint32_t dim = input_dim[d];
        map<uint32_t, vector<int64_t>> h1;
        for (size_t n = 0; n < input_size.size(); n++) {
            uint32_t len = input_size[n];
            vector<int64_t> ht1(3);
            getline(inr, line);
            istringstream iss(line);
            string token;
            uint32_t j = 0;
            while (getline(iss, token, ',')) {
                ht1[j] = stoull(token);
                ++j;
            }
            h1[len] = move(ht1);
        }
        relaxData[dim] = move(h1);
    }
    inr.close();
}

void Dsky(vector<vector<uint64_t>> P) {
    map<uint32_t, vector<uint64_t>> T;
    clock_t s01 = clock();
    vector<vector<IndexElement>> Ind(d, vector<IndexElement>(n));
    vector<vector<uint32_t>> IndMap(d, vector<uint32_t>(n));
    vector<uint32_t> Left(d, 0);
    for (uint32_t j = 0; j < d; ++j) {
        vector<ElementIndex> list;
        for (uint32_t i = 0; i < n; ++i) {
            list.emplace_back(P[i][j], i);
        }
        sort(list.begin(), list.end());
        for (uint32_t i = 0; i < n; ++i) {
            Ind[j][i].element = list[i].element;
            Ind[j][i].index   = list[i].index;
            Ind[j][i].visit   = 0;
        }
        for (uint32_t i = 0; i < n; ++i) {
            IndMap[j][list[i].index] = i;
        }
    }
    for (uint32_t i = 0; i < n; ++i) {
        T[i] = P[i];
    }
    // Select
    uint32_t flag = 1;
    map<uint32_t, map<uint32_t, vector<uint64_t>>> LES, LESQ;
    for (uint32_t j = 0; j < d; ++j) {
        map<uint32_t, vector<uint64_t>> h1, h2;
        for (uint32_t k = 0; k < n; ++k) {
            vector<uint64_t> ht1(n), ht2(n);
            for (uint32_t i = 0; i < n; ++i) {
                ht1[i] = (P[k][j] < P[i][j]) ? 0 : 1;
                ht2[i] = (P[k][j] <= P[i][j]) ? 0 : 1;
            }
            h1[k] = move(ht1);
            h2[k] = move(ht2);
        }
        LES[j]  = move(h1);
        LESQ[j] = move(h2);
    }
    clock_t s02     = clock();
    offline_process = double(s02 - s01) * 1000.0 / CLOCKS_PER_SEC;
    int itr = 0, num = 0;
    vector<uint32_t> posF(d, 0);
    vector<uint32_t> posI(d, 0);
    NumItr.resize(n + 1, 0);
    while (flag == 1) {
        n = T.size();
        cout << "itr:" << itr << "\tn:" << n << endl;
        itr++;
        num += n;
        NumItr[n]++;
        // cout << "LES" << endl;
        // for (uint32_t j = 0; j < LES.size(); ++j) {
        //     map<uint32_t, vector<uint64_t>> LESt = LES[j];
        //     for (auto& inner_pair : LESt) {
        //         auto& values = inner_pair.second;
        //         for (auto value : values) {
        //         cout << value << "\t";
        //         }
        //         cout << endl;
        //     }
        //     cout << endl;
        // }
        // cout << "LESQ" << endl;
        // for (uint32_t j = 0; j < d; ++j) {
        //     map<uint32_t, vector<uint64_t>> LESQt = LESQ[j];
        //     for (auto& inner_pair : LESQt) {
        //         auto& values = inner_pair.second;
        //         for (auto value : values) {
        //         cout << value << "\t";
        //         }
        //         cout << endl;
        //     }
        //     cout << endl;
        // }
        srand(time(0));
        for (uint32_t j = 0; j < d; ++j) {
            if (posF[j] == Ind[j].size()) {
                size_t randomIndex = rand() % T.size();
                auto it            = T.begin();
                next(it, randomIndex);
                posI[j] = it->first;
            }
            else {
                auto& row = Ind[j];
                for (uint32_t i = Left[j]; i < row.size(); ++i) {
                    if (row[i].visit == 0) {
                        posI[j]      = row[i].index;
                        row[i].visit = 1;
                        Left[j]      = i;
                        break;
                    }
                }
            }
        }
        // cout << "Anch:" << endl;
        // for (uint32_t j = 0; j < d; ++j) {
        //     cout << posI[j] << "\t" << T.at(posI[j])[j] << endl;
        // }
        vector<vector<uint64_t>> DomT(n, vector<uint64_t>(d, 0));
        vector<vector<vector<uint64_t>>> DomRow(d, vector<vector<uint64_t>>(n, vector<uint64_t>(d)));
        for (uint32_t j = 0; j < d; ++j) {
            for (uint32_t k = 0; k < d; ++k) {
                vector<uint64_t> innerMap = (j == k) ? LES[j][posI[k]] : LESQ[j][posI[k]];
                for (uint32_t i = 0; i < n; ++i) {
                    DomRow[j][i][k] = innerMap[i];
                }
            }
        }
        // cout << "DomRow:" << endl;
        // for (uint32_t j = 0; j < d; ++j) {
        //     for (uint32_t i = 0; i < n; ++i) {
        //         for (uint32_t k = 0; k < d; ++k) {
        //             cout << DomRow[j][i][k] << "\t";
        //         }
        //         cout << endl;
        //     }
        //     cout << endl;
        // }
        // compute DomT
        for (uint32_t i = 0; i < n; ++i) {
            for (uint32_t k = 0; k < d; ++k) {
                for (uint32_t j = 0; j < d; ++j) {
                    DomT[i][k] += DomRow[j][i][k];
                }
            }
        }
        // cout << "DomT:" << endl;
        // for (uint32_t i = 0; i < n; ++i) {
        //     for (uint32_t j = 0; j < d; ++j) {
        //         cout << DomT[i][j] << "\t";
        //     }
        //     cout << endl;
        // }
        // compute Dom
        vector<uint64_t> Dom(n);
        for (uint32_t i = 0; i < n; ++i) {
            uint64_t t = DomT[i][0];
            for (uint32_t k = 1; k < d; ++k) {
                t *= DomT[i][k];
            }
            Dom[i] = t;
        }
        // cout << "Dom:" << endl;
        // for (uint32_t i = 0; i < n; ++i) {
        //     cout << Dom[i] << "\t";
        // }
        // cout << endl;
        // index of T
        vector<uint32_t> indT(n);
        uint32_t i0 = 0;
        for (auto entry = T.begin(); entry != T.end(); ++entry) {
            indT[i0++] = entry->first;
        }
        // cout << "indT" << endl;
        // for (uint32_t i = 0; i < n; i++) {
        //     cout << indT[i] << "\t";
        // }
        // cout << endl;
        for (uint32_t i = 0; i < n; ++i) {
            if (Dom[i] == 0) {
                uint32_t keyToRemove = indT[i];
                T.erase(keyToRemove);
                for (uint32_t j = 0; j < d; ++j) {
                    auto& innerMapS = LES[j];
                    innerMapS.erase(keyToRemove);
                    auto& innerMapQ = LESQ[j];
                    innerMapQ.erase(keyToRemove);
                }
            }
        }
        // update LES and LESQ，remove 0 in Dom
        for (uint32_t i = n - 1; i >= 0 && i < n; --i) {
            if (Dom[i] == 0) {
                for (uint32_t j = 0; j < d; ++j) {
                    // auto& innerMapS = LES[j];
                    for (auto& entry : LES[j]) {
                        auto& vectorS = entry.second;
                        vectorS.erase(vectorS.begin() + i);
                    }
                    // auto& innerMapQ = LESQ[j];
                    for (auto& entry : LESQ[j]) {
                        auto& vectorQ = entry.second;
                        vectorQ.erase(vectorQ.begin() + i);
                    }
                }
            }
        }
        for (uint32_t j = 0; j < d; ++j) {
            for (uint32_t i = 0; i < n; ++i) {
                if (Dom[i] == 0)
                    // Update Ind
                    Ind[j][IndMap[j][indT[i]]].visit = 1;
            }
        }

        uint32_t posT = 0;
        for (uint32_t j = 0; j < d; ++j) {
            posF[j]         = 0;
            auto& Ind_j     = Ind[j];
            uint32_t size_j = Ind_j.size();
            for (uint32_t i = 0; i < size_j; ++i) {
                posF[j] += Ind_j[i].visit;
            }
            if (posF[j] == size_j)
                ++posT;
        }
        if (posT == d || T.empty())
            flag = 0;
        // cout << "index" << endl;
        // for (uint32_t j = 0; j < d; j++) {
        //     cout << posF[j] << "\t";
        // }
        // cout << endl;
    }
    clock_t s03   = clock();
    query_process = double(s03 - s02) * 1000.0 / CLOCKS_PER_SEC;
    cout << "num:" << num << endl;
    uint32_t k = T.size();
    cout << "k:" << k << endl;
    SkyNum = k;
    ItrNum = num;
    ComNum = itr;
    // for (auto& entry : T) {
    //     for (auto& r : P[entry.first]) {
    //         cout << r << "\t";
    //     }
    //     cout << endl;
    // }
}

void Dsky1(vector<vector<uint64_t>> P) {
    map<uint32_t, vector<uint64_t>> T;
    clock_t s01 = clock();
    vector<vector<IndexElement>> Ind(d, vector<IndexElement>(n));
    vector<vector<uint32_t>> IndMap(d, vector<uint32_t>(n));
    vector<vector<int64_t>> FlagF(d, vector<int64_t>(n, 0));
    vector<uint32_t> Left(d, 0);
    for (uint32_t j = 0; j < d; ++j) {
        vector<ElementIndex> list;
        for (uint32_t i = 0; i < n; ++i) {
            list.emplace_back(P[i][j], i);
        }
        sort(list.begin(), list.end());
        for (uint32_t i = 0; i < n; ++i) {
            Ind[j][i].element = list[i].element;
            Ind[j][i].index   = list[i].index;
            Ind[j][i].visit   = 0;
        }
        for (uint32_t i = 0; i < n; ++i) {
            IndMap[j][list[i].index] = i;
        }
    }
    for (uint32_t i = 0; i < n; ++i) {
        T[i] = P[i];
    }
    // Select
    uint32_t flag = 1;
    map<uint32_t, map<uint32_t, vector<uint64_t>>> LES, LESQ;
    for (uint32_t j = 0; j < d; ++j) {
        map<uint32_t, vector<uint64_t>> h1, h2;
        for (uint32_t k = 0; k < n; ++k) {
            vector<uint64_t> ht1(n), ht2(n);
            for (uint32_t i = 0; i < n; ++i) {
                ht1[i] = (P[k][j] < P[i][j]) ? 0 : 1;
                ht2[i] = (P[k][j] <= P[i][j]) ? 0 : 1;
            }
            h1[k] = move(ht1);
            h2[k] = move(ht2);
        }
        LES[j]  = move(h1);
        LESQ[j] = move(h2);
    }
    clock_t s02     = clock();
    offline_process = double(s02 - s01) * 1000.0 / CLOCKS_PER_SEC;
    int itr = 0, num = 0;
    vector<uint32_t> posF(d, 0);
    vector<uint32_t> posI(d, 0);
    NumItr.resize(n + 1, 0);
    usint Alln = n;
    while (flag == 1) {
        n = T.size();
        cout << "itr:" << itr << "\tn:" << n << endl;
        itr++;
        num += n;
        NumItr[n]++;
        // cout << "LES" << endl;
        // for (uint32_t j = 0; j < LES.size(); ++j) {
        //     map<uint32_t, vector<uint64_t>> LESt = LES[j];
        //     for (auto& inner_pair : LESt) {
        //         auto& values = inner_pair.second;
        //         for (auto value : values) {
        //         cout << value << "\t";
        //         }
        //         cout << endl;
        //     }
        //     cout << endl;
        // }
        // cout << "LESQ" << endl;
        // for (uint32_t j = 0; j < d; ++j) {
        //     map<uint32_t, vector<uint64_t>> LESQt = LESQ[j];
        //     for (auto& inner_pair : LESQt) {
        //         auto& values = inner_pair.second;
        //         for (auto value : values) {
        //         cout << value << "\t";
        //         }
        //         cout << endl;
        //     }
        //     cout << endl;
        // }
        srand(time(0));
        for (uint32_t j = 0; j < d; ++j) {
            if (posF[j] == Ind[j].size()) {
                size_t randomIndex = rand() % T.size();
                auto it            = T.begin();
                next(it, randomIndex);
                posI[j] = it->first;
            }
            else {
                auto& row = Ind[j];
                for (uint32_t i = Left[j]; i < row.size(); ++i) {
                    if (row[i].visit == 0) {
                        posI[j]           = row[i].index;
                        row[i].visit      = 1;
                        Left[j]           = i;
                        FlagF[j][posI[j]] = 1;
                        break;
                    }
                }
            }
        }
        // cout << "Anch:" << endl;
        // for (uint32_t j = 0; j < d; ++j) {
        //     cout << posI[j] << "\t" << T.at(posI[j])[j] << endl;
        // }
        vector<vector<uint64_t>> DomT(n, vector<uint64_t>(d, 0));
        vector<vector<vector<uint64_t>>> DomRow(d, vector<vector<uint64_t>>(n, vector<uint64_t>(d)));
        for (uint32_t j = 0; j < d; ++j) {
            for (uint32_t k = 0; k < d; ++k) {
                vector<uint64_t> innerMap = (j == k) ? LES[j][posI[k]] : LESQ[j][posI[k]];
                for (uint32_t i = 0; i < n; ++i) {
                    DomRow[j][i][k] = innerMap[i];
                }
            }
        }
        // cout << "DomRow:" << endl;
        // for (uint32_t j = 0; j < d; ++j) {
        //     for (uint32_t i = 0; i < n; ++i) {
        //         for (uint32_t k = 0; k < d; ++k) {
        //             cout << DomRow[j][i][k] << "\t";
        //         }
        //         cout << endl;
        //     }
        //     cout << endl;
        // }
        // compute DomT
        for (uint32_t i = 0; i < n; ++i) {
            for (uint32_t k = 0; k < d; ++k) {
                for (uint32_t j = 0; j < d; ++j) {
                    DomT[i][k] += DomRow[j][i][k];
                }
            }
        }
        // cout << "DomT:" << endl;
        // for (uint32_t i = 0; i < n; ++i) {
        //     for (uint32_t j = 0; j < d; ++j) {
        //         cout << DomT[i][j] << "\t";
        //     }
        //     cout << endl;
        // }
        // compute Dom
        vector<uint64_t> Dom(n);
        for (uint32_t i = 0; i < n; ++i) {
            uint64_t t = DomT[i][0];
            for (uint32_t k = 1; k < d; ++k) {
                t *= DomT[i][k];
            }
            Dom[i] = t;
        }
        // cout << "Dom:" << endl;
        // for (uint32_t i = 0; i < n; ++i) {
        //     cout << Dom[i] << "\t";
        // }
        // cout << endl;
        // index of T
        vector<uint32_t> indT(n);
        uint32_t i0 = 0;
        for (auto entry = T.begin(); entry != T.end(); ++entry) {
            indT[i0++] = entry->first;
        }
        // cout << "indT" << endl;
        // for (uint32_t i = 0; i < n; i++) {
        //     cout << indT[i] << "\t";
        // }
        // cout << endl;
        for (uint32_t i = 0; i < n; ++i) {
            if (Dom[i] == 0) {
                uint32_t keyToRemove = indT[i];
                T.erase(keyToRemove);
                for (uint32_t j = 0; j < d; ++j) {
                    auto& innerMapS = LES[j];
                    innerMapS.erase(keyToRemove);
                    auto& innerMapQ = LESQ[j];
                    innerMapQ.erase(keyToRemove);
                }
            }
        }
        // update LES and LESQ，remove 0 in Dom
        for (uint32_t i = n - 1; i >= 0 && i < n; --i) {
            if (Dom[i] == 0) {
                for (uint32_t j = 0; j < d; ++j) {
                    // auto& innerMapS = LES[j];
                    for (auto& entry : LES[j]) {
                        auto& vectorS = entry.second;
                        vectorS.erase(vectorS.begin() + i);
                    }
                    // auto& innerMapQ = LESQ[j];
                    for (auto& entry : LESQ[j]) {
                        auto& vectorQ = entry.second;
                        vectorQ.erase(vectorQ.begin() + i);
                    }
                }
            }
        }
        for (uint32_t j = 0; j < d; ++j) {
            for (uint32_t i = 0; i < n; ++i) {
                if (Dom[i] == 0)
                    // Update Ind
                    Ind[j][IndMap[j][indT[i]]].visit = 1;
            }
        }
        vector<int64_t> FlagF_ALL(FlagF[0].size(), 0);
        int64_t FlagFinal = 0;
        for (uint32_t i = 0; i < Alln; i++) {
            for (uint32_t j = 0; j < d; j++) {
                FlagF_ALL[i] += FlagF[j][i];
            }
            if (FlagF_ALL[i] == d) {
                FlagFinal = 1;
                break;
            }
        }
        if (FlagFinal >= 1)
            break;
        uint32_t posT = 0;
        for (uint32_t j = 0; j < d; ++j) {
            posF[j]         = 0;
            auto& Ind_j     = Ind[j];
            uint32_t size_j = Ind_j.size();
            for (uint32_t i = 0; i < size_j; ++i) {
                posF[j] += Ind_j[i].visit;
            }
            if (posF[j] == size_j)
                ++posT;
        }
        if (posT == d || T.empty())
            flag = 0;
        // cout << "index" << endl;
        // for (uint32_t j = 0; j < d; j++) {
        //     cout << posF[j] << "\t";
        // }
        // cout << endl;
    }
    clock_t s03   = clock();
    query_process = double(s03 - s02) * 1000.0 / CLOCKS_PER_SEC;
    cout << "num:" << num << endl;
    uint32_t k = T.size();
    cout << "k:" << k << endl;
    SkyNum = k;
    ItrNum = num;
    ComNum = itr;
    // for (auto& entry : T) {
    //     for (auto& r : P[entry.first]) {
    //         cout << r << "\t";
    //     }
    //     cout << endl;
    // }
}

void Dsky12(vector<vector<uint64_t>> P) {
    map<uint32_t, vector<uint64_t>> T;
    clock_t s01 = clock();
    vector<vector<IndexElement>> Ind(d, vector<IndexElement>(n));
    vector<vector<uint32_t>> IndMap(d, vector<uint32_t>(n));
    vector<vector<int64_t>> FlagF(d, vector<int64_t>(n, 0));
    vector<uint32_t> Left(d, 0);
    for (uint32_t j = 0; j < d; ++j) {
        vector<ElementIndex> list;
        for (uint32_t i = 0; i < n; ++i) {
            list.emplace_back(P[i][j], i);
        }
        sort(list.begin(), list.end());
        for (uint32_t i = 0; i < n; ++i) {
            Ind[j][i].element = list[i].element;
            Ind[j][i].index   = list[i].index;
            Ind[j][i].visit   = 0;
        }
        for (uint32_t i = 0; i < n; ++i) {
            IndMap[j][list[i].index] = i;
        }
    }
    for (uint32_t i = 0; i < n; ++i) {
        T[i] = P[i];
    }
    // Select
    uint32_t flag = 1;
    map<uint32_t, map<uint32_t, vector<uint64_t>>> LES, LESQ;
    for (uint32_t j = 0; j < d; ++j) {
        map<uint32_t, vector<uint64_t>> h1, h2;
        for (uint32_t k = 0; k < n; ++k) {
            vector<uint64_t> ht1(n), ht2(n);
            for (uint32_t i = 0; i < n; ++i) {
                ht1[i] = (P[k][j] < P[i][j]) ? 0 : 1;
                ht2[i] = (P[k][j] <= P[i][j]) ? 0 : 1;
            }
            h1[k] = move(ht1);
            h2[k] = move(ht2);
        }
        LES[j]  = move(h1);
        LESQ[j] = move(h2);
    }
    clock_t s02     = clock();
    offline_process = double(s02 - s01) * 1000.0 / CLOCKS_PER_SEC;
    int itr = 0, num = 0;
    vector<uint32_t> posF(d, 0);
    vector<uint32_t> posI(d, 0);
    NumItr.resize(n + 1, 0);
    usint Alln = n;
    while (flag == 1) {
        n = T.size();
        cout << "itr:" << itr << "\tn:" << n << endl;
        itr++;
        num += n;
        NumItr[n]++;
        if (NumItr[n] > ts)
            break;
        // cout << "LES" << endl;
        // for (uint32_t j = 0; j < LES.size(); ++j) {
        //     map<uint32_t, vector<uint64_t>> LESt = LES[j];
        //     for (auto& inner_pair : LESt) {
        //         auto& values = inner_pair.second;
        //         for (auto value : values) {
        //         cout << value << "\t";
        //         }
        //         cout << endl;
        //     }
        //     cout << endl;
        // }
        // cout << "LESQ" << endl;
        // for (uint32_t j = 0; j < d; ++j) {
        //     map<uint32_t, vector<uint64_t>> LESQt = LESQ[j];
        //     for (auto& inner_pair : LESQt) {
        //         auto& values = inner_pair.second;
        //         for (auto value : values) {
        //         cout << value << "\t";
        //         }
        //         cout << endl;
        //     }
        //     cout << endl;
        // }
        srand(time(0));
        for (uint32_t j = 0; j < d; ++j) {
            if (posF[j] == Ind[j].size()) {
                size_t randomIndex = rand() % T.size();
                auto it            = T.begin();
                next(it, randomIndex);
                posI[j] = it->first;
            }
            else {
                auto& row = Ind[j];
                for (uint32_t i = Left[j]; i < row.size(); ++i) {
                    if (row[i].visit == 0) {
                        posI[j]           = row[i].index;
                        row[i].visit      = 1;
                        Left[j]           = i;
                        FlagF[j][posI[j]] = 1;
                        break;
                    }
                }
            }
        }
        // cout << "Anch:" << endl;
        // for (uint32_t j = 0; j < d; ++j) {
        //     cout << posI[j] << "\t" << T.at(posI[j])[j] << endl;
        // }
        vector<vector<uint64_t>> DomT(n, vector<uint64_t>(d, 0));
        vector<vector<vector<uint64_t>>> DomRow(d, vector<vector<uint64_t>>(n, vector<uint64_t>(d)));
        for (uint32_t j = 0; j < d; ++j) {
            for (uint32_t k = 0; k < d; ++k) {
                vector<uint64_t> innerMap = (j == k) ? LES[j][posI[k]] : LESQ[j][posI[k]];
                for (uint32_t i = 0; i < n; ++i) {
                    DomRow[j][i][k] = innerMap[i];
                }
            }
        }
        // cout << "DomRow:" << endl;
        // for (uint32_t j = 0; j < d; ++j) {
        //     for (uint32_t i = 0; i < n; ++i) {
        //         for (uint32_t k = 0; k < d; ++k) {
        //             cout << DomRow[j][i][k] << "\t";
        //         }
        //         cout << endl;
        //     }
        //     cout << endl;
        // }
        // compute DomT
        for (uint32_t i = 0; i < n; ++i) {
            for (uint32_t k = 0; k < d; ++k) {
                for (uint32_t j = 0; j < d; ++j) {
                    DomT[i][k] += DomRow[j][i][k];
                }
            }
        }
        // cout << "DomT:" << endl;
        // for (uint32_t i = 0; i < n; ++i) {
        //     for (uint32_t j = 0; j < d; ++j) {
        //         cout << DomT[i][j] << "\t";
        //     }
        //     cout << endl;
        // }
        // compute Dom
        vector<uint64_t> Dom(n);
        for (uint32_t i = 0; i < n; ++i) {
            uint64_t t = DomT[i][0];
            for (uint32_t k = 1; k < d; ++k) {
                t *= DomT[i][k];
            }
            Dom[i] = t;
        }
        // cout << "Dom:" << endl;
        // for (uint32_t i = 0; i < n; ++i) {
        //     cout << Dom[i] << "\t";
        // }
        // cout << endl;
        // index of T
        vector<uint32_t> indT(n);
        uint32_t i0 = 0;
        for (auto entry = T.begin(); entry != T.end(); ++entry) {
            indT[i0++] = entry->first;
        }
        // cout << "indT" << endl;
        // for (uint32_t i = 0; i < n; i++) {
        //     cout << indT[i] << "\t";
        // }
        // cout << endl;
        for (uint32_t i = 0; i < n; ++i) {
            if (Dom[i] == 0) {
                uint32_t keyToRemove = indT[i];
                T.erase(keyToRemove);
                for (uint32_t j = 0; j < d; ++j) {
                    auto& innerMapS = LES[j];
                    innerMapS.erase(keyToRemove);
                    auto& innerMapQ = LESQ[j];
                    innerMapQ.erase(keyToRemove);
                }
            }
        }
        // update LES and LESQ，remove 0 in Dom
        for (uint32_t i = n - 1; i >= 0 && i < n; --i) {
            if (Dom[i] == 0) {
                for (uint32_t j = 0; j < d; ++j) {
                    // auto& innerMapS = LES[j];
                    for (auto& entry : LES[j]) {
                        auto& vectorS = entry.second;
                        vectorS.erase(vectorS.begin() + i);
                    }
                    // auto& innerMapQ = LESQ[j];
                    for (auto& entry : LESQ[j]) {
                        auto& vectorQ = entry.second;
                        vectorQ.erase(vectorQ.begin() + i);
                    }
                }
            }
        }
        for (uint32_t j = 0; j < d; ++j) {
            for (uint32_t i = 0; i < n; ++i) {
                if (Dom[i] == 0)
                    // Update Ind
                    Ind[j][IndMap[j][indT[i]]].visit = 1;
            }
        }
        vector<int64_t> FlagF_ALL(FlagF[0].size(), 0);
        int64_t FlagFinal = 0;
        for (uint32_t i = 0; i < Alln; i++) {
            for (uint32_t j = 0; j < d; j++) {
                FlagF_ALL[i] += FlagF[j][i];
            }
            if (FlagF_ALL[i] == d) {
                FlagFinal = 1;
                break;
            }
        }
        if (FlagFinal >= 1)
            break;
        uint32_t posT = 0;
        for (uint32_t j = 0; j < d; ++j) {
            posF[j]         = 0;
            auto& Ind_j     = Ind[j];
            uint32_t size_j = Ind_j.size();
            for (uint32_t i = 0; i < size_j; ++i) {
                posF[j] += Ind_j[i].visit;
            }
            if (posF[j] == size_j)
                ++posT;
        }
        if (posT == d || T.empty())
            flag = 0;
        // cout << "index" << endl;
        // for (uint32_t j = 0; j < d; j++) {
        //     cout << posF[j] << "\t";
        // }
        // cout << endl;
    }
    clock_t s03   = clock();
    query_process = double(s03 - s02) * 1000.0 / CLOCKS_PER_SEC;
    cout << "num:" << num << endl;
    uint32_t k = T.size();
    cout << "k:" << k << endl;
    SkyNum = k;
    ItrNum = num;
    ComNum = itr;
    // for (auto& entry : T) {
    //     for (auto& r : P[entry.first]) {
    //         cout << r << "\t";
    //     }
    //     cout << endl;
    // }
}

CryptoContext<DCRTPoly> GenerateBFVrnsContext() {
    CCParams<CryptoContextBFVRNS> parameters;
    parameters.SetSecurityLevel(SecurityLevel::HEStd_128_classic);
    // parameters.SetSecurityLevel(SecurityLevel::HEStd_NotSet);
    parameters.SetPlaintextModulus(plaintextModulus);
    vector<uint32_t> Deptht = {3, 3, 3, 4, 4, 4};
    usint multDepth         = Deptht[d - 2];
    usint digitSize         = 5;
    usint dcrtBits  = 50;
    parameters.SetMultiplicativeDepth(multDepth);
    parameters.SetSecretKeyDist(UNIFORM_TERNARY);
    parameters.SetDigitSize(digitSize);
    parameters.SetStandardDeviation(3.2);
    parameters.SetScalingModSize(dcrtBits);
    // parameters.SetRingDim(ringDim);
    // parameters.SetMaxRelinSkDeg(4);
    parameters.SetMultiplicationTechnique(HPSPOVERQLEVELED);
    // parameters.SetBatchSize(16);
    CryptoContext<DCRTPoly> cryptoContext = GenCryptoContext(parameters);
    // Enable the features that you wish to use.
    cryptoContext->Enable(PKE);
    cryptoContext->Enable(KEYSWITCH);
    cryptoContext->Enable(LEVELEDSHE);
    cryptoContext->Enable(ADVANCEDSHE);
    cryptoContext->Enable(MULTIPARTY);
    cout << "Security Level: " << parameters.GetSecurityLevel() << endl;
    cout << "Deviation: " << parameters.GetStandardDeviation() << "\t";
    cout << "Depth = " << parameters.GetMultiplicativeDepth() << "\t";
    cout << "p = " << cryptoContext->GetCryptoParameters()->GetPlaintextModulus() << "\t";
    cout << "n = " << cryptoContext->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2 << "\t";
    cout << "log2 q = "
         << log2(cryptoContext->GetCryptoParameters()->GetElementParams()->GetModulus().ConvertToDouble()) << endl;
    return cryptoContext;
}

void Init0() {
    TIC(t);
    cc             = GenerateBFVrnsContext();
    processingTime = TOC(t);
    cout << "Key generation time: " << processingTime << "ms" << endl;
    TIC(t);
    // Initialize Public Key Containers
    vector<EvalKey<DCRTPoly>> evalMultKeys(d);
    vector<EvalKey<DCRTPoly>> evalMult_i_ALL(d);
    vector<shared_ptr<map<usint, EvalKey<DCRTPoly>>>> evalSumKeys(d);
    vector<shared_ptr<map<usint, EvalKey<DCRTPoly>>>> evalAtIndexKeys(d);
    cout << "Running key generation..." << endl;
    kps.resize(d);
    kps[0] = cc->KeyGen();
    for (uint32_t i = 1; i < d; i++) {
        kps[i] = cc->MultipartyKeyGen(kps[i - 1].publicKey);
    }
    cout << "Joint public key..." << endl;
    evalMultKeys[0] = cc->KeySwitchGen(kps[0].secretKey, kps[0].secretKey);
    for (uint32_t i = 1; i < d; i++) {
        evalMultKeys[i] = cc->MultiKeySwitchGen(kps[i].secretKey, kps[i].secretKey, evalMultKeys[0]);
    }
    cout << "Joint evaluation multiplication key..." << endl;
    // vector<EvalKey<DCRTPoly>> evalMultALL_i(d-1);
    // evalMultALL_i[0] = cc->MultiAddEvalKeys(evalMultKeys[0], evalMultKeys[1], kps[1].publicKey->GetKeyTag());
    // for (uint32_t i = 2; i < d; i++) {
    //     evalMultALL_i[i-1] = cc->MultiAddEvalKeys(evalMultALL_i[i-2], evalMultKeys[i], kps[i].publicKey->GetKeyTag());
    // }
    // auto evalMultALL = evalMultALL_i[d-1];
    auto evalMultALL = cc->MultiAddEvalKeys(evalMultKeys[0], evalMultKeys[1], kps[1].publicKey->GetKeyTag());
    for (uint32_t i = 2; i < d; i++) {
        auto evalMultTmp = cc->MultiAddEvalKeys(evalMultALL, evalMultKeys[i], kps[i].publicKey->GetKeyTag());
        evalMultALL      = evalMultTmp;
    }
    for (uint32_t i = 0; i < d; i++) {
        evalMult_i_ALL[i] = cc->MultiMultEvalKey(kps[i].secretKey, evalMultALL, kps[d - 1].publicKey->GetKeyTag());
    }
    cout << "Final evaluation multiplication key..." << endl;
    // if (d == 2) {
    //     auto evalMultFinal = cc->MultiAddEvalMultKeys(evalMult_i_ALL[0], evalMult_i_ALL[1], evalMultALL->GetKeyTag());
    //     cc->InsertEvalMultKey({evalMultFinal});
    // }
    // else {
    //     auto evalMultTmp = cc->MultiAddEvalMultKeys(evalMult_i_ALL[d - 1], evalMult_i_ALL[d - 2], evalMult_i_ALL[d - 1]->GetKeyTag());
    //     for (uint32_t i = d - 3; i > 1; i--) {
    //         auto evalMultTmpt = cc->MultiAddEvalMultKeys(evalMult_i_ALL[i], evalMultTmp, evalMult_i_ALL[i]->GetKeyTag());
    //         evalMultTmp     = evalMultTmpt;
    //     }
    //     auto evalMultFinal = cc->MultiAddEvalMultKeys(evalMult_i_ALL[0], evalMultTmp, kps[d-1].publicKey->GetKeyTag());
    //     cc->InsertEvalMultKey({evalMultFinal});
    // }
    auto evalMultFinal = cc->MultiAddEvalMultKeys(evalMult_i_ALL[0], evalMult_i_ALL[1], kps[1].publicKey->GetKeyTag());
    for (uint32_t i = 2; i < d; i++) {
        auto evalMultTmp = cc->MultiAddEvalMultKeys(evalMultFinal, evalMult_i_ALL[i], kps[i].publicKey->GetKeyTag());
        evalMultFinal    = evalMultTmp;
    }
    cc->InsertEvalMultKey({evalMultFinal});
    cout << "Running evalsum key generation..." << endl;
    cc->EvalSumKeyGen(kps[0].secretKey);
    evalSumKeys[0] = make_shared<map<usint, EvalKey<DCRTPoly>>>(cc->GetEvalSumKeyMap(kps[0].secretKey->GetKeyTag()));
    for (uint32_t i = 1; i < d; i++) {
        evalSumKeys[i] = cc->MultiEvalSumKeyGen(kps[i].secretKey, evalSumKeys[0], kps[i].publicKey->GetKeyTag());
    }
    auto evalSumKeysJoin = cc->MultiAddEvalSumKeys(evalSumKeys[0], evalSumKeys[1], kps[1].publicKey->GetKeyTag());
    for (uint32_t i = 2; i < d; i++) {
        auto evalMultTmp = cc->MultiAddEvalSumKeys(evalSumKeysJoin, evalSumKeys[i], kps[i].publicKey->GetKeyTag());
        evalSumKeysJoin  = evalMultTmp;
    }
    cc->InsertEvalSumKey(evalSumKeysJoin);
    cout << "Running evalAutomorphism key generation..." << endl;
    vector<int32_t> indexList(n);
    for (uint32_t i = 1; i < n + 1; ++i) {
        indexList[i - 1] = 0 - i;
    }
    cc->EvalAtIndexKeyGen(kps[0].secretKey, indexList);
    evalAtIndexKeys[0] =
        make_shared<map<usint, EvalKey<DCRTPoly>>>(cc->GetEvalAutomorphismKeyMap(kps[0].secretKey->GetKeyTag()));
    for (uint32_t i = 1; i < d; i++) {
        evalAtIndexKeys[i] =
            cc->MultiEvalAtIndexKeyGen(kps[i].secretKey, evalAtIndexKeys[0], indexList, kps[i].publicKey->GetKeyTag());
    }
    auto evalAtIndexKeysJoin =
        cc->MultiAddEvalAutomorphismKeys(evalAtIndexKeys[0], evalAtIndexKeys[1], kps[1].publicKey->GetKeyTag());
    for (uint32_t i = 2; i < d; i++) {
        auto evalMultTmp =
            cc->MultiAddEvalAutomorphismKeys(evalAtIndexKeysJoin, evalAtIndexKeys[i], kps[i].publicKey->GetKeyTag());
        evalAtIndexKeysJoin = evalMultTmp;
    }
    cc->InsertEvalAutomorphismKey(evalAtIndexKeysJoin);
    cout << "EvalAutomorphismKey " << endl;
    keys           = kps[d - 1];
    processingTime = TOC(t);
    cout << "Joint Key generation time: " << processingTime << "ms" << endl;
    // vector<PrivateKey<DCRTPoly>> secretKeys(d);
    // for (size_t i = 0; i < d; i++)
    // {
    //     secretKeys[i] = kps[i].secretKey;
    // }
    // KeyPair<DCRTPoly> kpMultiparty = cc->MultipartyKeyGen(secretKeys);
    // cc->Decrypt(kpMultiparty.secretKey, ciphertext, &plaintext);
}

void Init1() {
    // TIC(t);
    // auto start1 = high_resolution_clock::now();
    cc = GenerateBFVrnsContext();
    // processingTime = TOC(t);
    // auto end1      = high_resolution_clock::now();
    // auto duration1 = duration_cast<milliseconds>(end1 - start1).count();
    // cout << "Key generation time: " << duration1 << " ms" << endl;
    // TIC(t);
    // auto start2 = high_resolution_clock::now();
    // Initialize Public Key Containers
    vector<EvalKey<DCRTPoly>> evalMultKeys(d);
    vector<EvalKey<DCRTPoly>> evalMult_i_ALL(d);
    vector<shared_ptr<map<usint, EvalKey<DCRTPoly>>>> evalSumKeys(d);
    // cout << "Running key generation..." << endl;
    kps.resize(d);
    kps[0] = cc->KeyGen();
    for (uint32_t i = 1; i < d; i++) {
        kps[i] = cc->MultipartyKeyGen(kps[i - 1].publicKey);
    }
    // cout << "Joint public key..." << endl;
    evalMultKeys[0] = cc->KeySwitchGen(kps[0].secretKey, kps[0].secretKey);
    for (uint32_t i = 1; i < d; i++) {
        evalMultKeys[i] = cc->MultiKeySwitchGen(kps[i].secretKey, kps[i].secretKey, evalMultKeys[0]);
    }
    // cout << "Joint evaluation multiplication key..." << endl;
    // vector<EvalKey<DCRTPoly>> evalMultALL_i(d-1);
    // evalMultALL_i[0] = cc->MultiAddEvalKeys(evalMultKeys[0], evalMultKeys[1], kps[1].publicKey->GetKeyTag());
    // for (uint32_t i = 2; i < d; i++) {
    //     evalMultALL_i[i-1] = cc->MultiAddEvalKeys(evalMultALL_i[i-2], evalMultKeys[i], kps[i].publicKey->GetKeyTag());
    // }
    // auto evalMultALL = evalMultALL_i[d-1];
    auto evalMultALL = cc->MultiAddEvalKeys(evalMultKeys[0], evalMultKeys[1], kps[1].publicKey->GetKeyTag());
    for (uint32_t i = 2; i < d; i++) {
        auto evalMultTmp = cc->MultiAddEvalKeys(evalMultALL, evalMultKeys[i], kps[i].publicKey->GetKeyTag());
        evalMultALL      = evalMultTmp;
    }
    for (uint32_t i = 0; i < d; i++) {
        evalMult_i_ALL[i] = cc->MultiMultEvalKey(kps[i].secretKey, evalMultALL, kps[d - 1].publicKey->GetKeyTag());
    }
    // cout << "Final evaluation multiplication key..." << endl;
    auto evalMultFinal = cc->MultiAddEvalMultKeys(evalMult_i_ALL[0], evalMult_i_ALL[1], kps[1].publicKey->GetKeyTag());
    for (uint32_t i = 2; i < d; i++) {
        auto evalMultTmp = cc->MultiAddEvalMultKeys(evalMultFinal, evalMult_i_ALL[i], kps[i].publicKey->GetKeyTag());
        evalMultFinal    = evalMultTmp;
    }
    cc->InsertEvalMultKey({evalMultFinal});
    // cout << "Running evalsum key generation..." << endl;
    cc->EvalSumKeyGen(kps[0].secretKey);
    evalSumKeys[0] = make_shared<map<usint, EvalKey<DCRTPoly>>>(cc->GetEvalSumKeyMap(kps[0].secretKey->GetKeyTag()));
    for (uint32_t i = 1; i < d; i++) {
        evalSumKeys[i] = cc->MultiEvalSumKeyGen(kps[i].secretKey, evalSumKeys[0], kps[i].publicKey->GetKeyTag());
    }
    auto evalSumKeysJoin = cc->MultiAddEvalSumKeys(evalSumKeys[0], evalSumKeys[1], kps[1].publicKey->GetKeyTag());
    for (uint32_t i = 2; i < d; i++) {
        auto evalMultTmp = cc->MultiAddEvalSumKeys(evalSumKeysJoin, evalSumKeys[i], kps[i].publicKey->GetKeyTag());
        evalSumKeysJoin  = evalMultTmp;
    }
    cc->InsertEvalSumKey(evalSumKeysJoin);
    keys = kps[d - 1];
    // processingTime = TOC(t);
    // auto end2      = high_resolution_clock::now();
    // auto duration2 = duration_cast<milliseconds>(end2 - start2).count();
    // cout << "Joint Key generation time: " << duration2 << " ms" << endl;
}

Plaintext Dec(Ciphertext<DCRTPoly> ciphertext, size_t n) {
    Plaintext res;
    // const shared_ptr<CryptoParametersBase<DCRTPoly>> cryptoParams = kps[0].secretKey->GetCryptoParameters();
    // const shared_ptr<typename DCRTPoly::Params> elementParams     = cryptoParams->GetElementParams();
    // cout << "\npartial decryption..." << endl;
    vector<Ciphertext<DCRTPoly>> partialCiphertextVec;
    map<usint, vector<Ciphertext<DCRTPoly>>> ciphertextPartial;
    ciphertextPartial[0] = cc->MultipartyDecryptLead({ciphertext}, kps[0].secretKey);
    partialCiphertextVec.push_back(ciphertextPartial[0][0]);
    for (uint32_t i = 1; i < d; i++) {
        ciphertextPartial[i] = cc->MultipartyDecryptMain({ciphertext}, kps[i].secretKey);
        partialCiphertextVec.push_back(ciphertextPartial[i][0]);
    }
    cc->MultipartyDecryptFusion(partialCiphertextVec, &res);
    res->SetLength(n);
    // processingTime = TOC(t);
    // cout << "Dec Completed\n";
    // cout << "Decryption time: " << processingTime << "ms" << endl;
    // cout << res << endl;
    return res;
}

vector<uint32_t> GenerateIndVec(uint32_t n) {
    vector<uint32_t> x(n);
    uint32_t mid      = n / 2;
    uint32_t current  = mid;
    uint32_t currentL = mid;
    uint32_t currentR = mid;
    bool goingDown    = true;
    for (uint32_t i = 0; i < n; ++i) {
        x[i] = current;
        if (goingDown) {
            --currentL;
            current   = currentL;
            goingDown = false;
        }
        else {
            ++currentR;
            current   = currentR;
            goingDown = true;
        }
    }
    return x;
}

void MKsky(vector<vector<uint64_t>> P) {
    offline_process = 0;
    query_process   = 0;
    map<uint32_t, vector<uint64_t>> T;
    double s01 = omp_get_wtime();
    vector<vector<IndexElement>> Ind(d, vector<IndexElement>(n));
    vector<vector<uint32_t>> IndMap(d, vector<uint32_t>(n));
    vector<uint32_t> Left(d, 0);
    vector<double> comm(d + 1), comp(d + 1);
    for (uint32_t j = 0; j < d; ++j) {
        vector<ElementIndex> list;
        for (uint32_t i = 0; i < n; ++i) {
            list.emplace_back(P[i][j], i);
        }
        sort(list.begin(), list.end());
        for (uint32_t i = 0; i < n; ++i) {
            Ind[j][i].element = list[i].element;
            Ind[j][i].index   = list[i].index;
            Ind[j][i].visit   = 0;
        }
        for (uint32_t i = 0; i < n; ++i) {
            IndMap[j][list[i].index] = i;
        }
    }
    for (uint32_t i = 0; i < n; ++i) {
        T[i] = P[i];
    }
    // Select
    uint32_t flag = 1;
    // store in column
    map<uint32_t, map<uint32_t, vector<int64_t>>> LES, LESQ;
    for (uint32_t j = 0; j < d; ++j) {
        map<uint32_t, vector<int64_t>> h1, h2;
        for (uint32_t k = 0; k < n; ++k) {
            vector<int64_t> ht1(n), ht2(n);
            for (uint32_t i = 0; i < n; ++i) {
                ht1[i] = (P[k][j] < P[i][j]) ? 0 : 1;
                ht2[i] = (P[k][j] <= P[i][j]) ? 0 : 1;
            }
            h1[k] = move(ht1);
            h2[k] = move(ht2);
        }
        LES[j]  = move(h1);
        LESQ[j] = move(h2);
    }
    unordered_map<uint32_t, Plaintext> splitVec;
    for (uint32_t i = 0; i < n; ++i) {
        vector<int64_t> atIndex(2 * n, 0);
        atIndex[n + i] = 1;
        splitVec[i]    = cc->MakePackedPlaintext(atIndex);
    }
    double s02      = omp_get_wtime();
    offline_process = s02 - s01;
    double s03      = omp_get_wtime();
    int itr = 0, num = 0;
    vector<uint32_t> posF(d, 0);
    vector<uint32_t> posI(d, 0);
    vector<Ciphertext<DCRTPoly>> IV_C(d);
    vector<Ciphertext<DCRTPoly>> DomRow(d * d);
    vector<Ciphertext<DCRTPoly>> DomT(d);
    usint Alln = n;
    while (flag == 1) {
        n = T.size();
        cout << "itr:" << itr << "\tn:" << n << endl;
        itr++;
        num += n;
        // cout << "LES" << endl;
        // for (uint32_t j = 0; j < LES.size(); ++j) {
        //     map<uint32_t, vector<int64_t>> LESt = LES[j];
        //     for (auto& inner_pair : LESt) {
        //         auto& values = inner_pair.second;
        //         for (auto value : values) {
        //         cout << value << "\t";
        //         }
        //         cout << endl;
        //     }
        //     cout << endl;
        // }
        // cout << "LESQ" << endl;
        // for (uint32_t j = 0; j < d; ++j) {
        //     map<uint32_t, vector<int64_t>> LESQt = LESQ[j];
        //     for (auto& inner_pair : LESQt) {
        //         auto& values = inner_pair.second;
        //         for (auto value : values) {
        //         cout << value << "\t";
        //         }
        //         cout << endl;
        //     }
        //     cout << endl;
        // }
        double st1 = omp_get_wtime();
        // vector<vector<int64_t>> Anch(d, vector<int64_t>(n, 0));
        // map<uint32_t, vector<int64_t>> Anch;
        map<uint32_t, Plaintext> LES_P;
        map<uint32_t, map<uint32_t, Plaintext>> LESQ_P;
        // Ciphertext<DCRTPoly> CZero = cc->Encrypt(keys.publicKey, cc->MakePackedPlaintext(vector<int64_t>(n, 0)));
        double st2                 = omp_get_wtime();
        for (uint32_t j = 0; j < d; ++j)
            comp[j] += st2 - st1;
        srand(time(0));
        for (uint32_t j = 0; j < d; ++j) {
            // compute the index "posI" about min value
            double st1 = omp_get_wtime();
            if (posF[j] == Ind[j].size()) {
                // size_t randomIndex = rand() % T.size();
                // auto it            = T.begin();
                // next(it, randomIndex);
                // posI[j] = it->first;
                posI[j] = Ind[j][0].index;
            }
            else {
                auto& row = Ind[j];
                for (uint32_t i = Left[j]; i < row.size(); ++i) {
                    if (row[i].visit == 0) {
                        posI[j]      = row[i].index;
                        row[i].visit = 1;
                        Left[j]      = i;
                        break;
                    }
                }
            }
            // compute the anchor "Anch" in the origin IV
            uint32_t closestIndex = distance(T.begin(), T.find(posI[j]));
            // Anch[j] = move(vector<int64_t>(n, 0));
            vector<int64_t> Ancht(2 * Alln, 0);
            // Anch[j][closestIndex] = 1;
            Ancht[Alln + closestIndex] = 1;
            LES_P[j]                   = cc->MakePackedPlaintext(LES[j][posI[j]]);
            // compute the encrypted "Anch" in the origin IV
            IV_C[j] = cc->Encrypt(keys.publicKey, cc->MakePackedPlaintext(Ancht));
            // cout << j << "-IV:" << Serial::SerializeToString(IV_C[j]).size() << endl;
            // comm[j] += Serial::SerializeToString(IV_C[j]).size();
            // packing the column value about LESQ
            map<uint32_t, vector<int64_t>> LESQt = LESQ[j];
            map<uint32_t, Plaintext> LESQ_Pt;
            uint32_t i = 0;
            for (auto& inner_pair : LESQt) {
                LESQ_Pt[i] = cc->MakePackedPlaintext(inner_pair.second);
                i++;
            }
            LESQ_P[j]  = move(LESQ_Pt);
            double st2 = omp_get_wtime();
            comp[j] += st2 - st1;
            size_t tmp = Serial::SerializeToString(IV_C[j]).size();
            comm[j] += tmp;
            comm[d] += tmp;
        }
        // compute the inner product between matrix and IV
        for (uint32_t j = 0; j < d; ++j) {
            double st1 = omp_get_wtime();
            for (uint32_t k = 0; k < d; ++k) {
                if (j == k) {
                    DomRow[j * d + k] = cc->Encrypt(keys.publicKey, LES_P[j]);
                }
                else {
                    // const auto& innerMap = LESQ_P[j];
                    // Ciphertext<DCRTPoly> Domi = CZero->Clone();
                    vector<Ciphertext<DCRTPoly>> Domi(n);
                    vector<Plaintext> LESQt;
                    LESQt.reserve(n);
                    for (const auto& entry : LESQ_P[j]) {
                        LESQt.push_back(entry.second);
                    }
#pragma omp parallel for
                    for (uint32_t i = 0; i < n; ++i) {
                        Ciphertext<DCRTPoly> splitInd = cc->EvalMult(IV_C[k], splitVec[i]);
                        Domi[i] = cc->EvalSum(splitInd, batchSize);
                        Domi[i] = cc->EvalMult(Domi[i], LESQt[i]);
                    }
                    // uint32_t i = 0;
                    // for (const auto& entry : innerMap) {
                    //     // // first try
                    //     // Ciphertext<DCRTPoly> splitInd = cc->EvalMult(IV_C[k], splitVec[i]);
                    //     // Ciphertext<DCRTPoly> IndexSum = cc->EvalSum(splitInd, batchSize);
                    //     // second try
                    //     // Ciphertext<DCRTPoly> splitInd   = cc->EvalInnerProduct(IV_C[k], splitVec[i], n);
                    //     // Ciphertext<DCRTPoly> Indsum     = cc->EvalMult(splitInd, splitVec[0]);
                    //     // Ciphertext<DCRTPoly> splitIndex = cc->EvalRotate(splitInd, indexList[0]);
                    //     // Ciphertext<DCRTPoly> IndexSum   = cc->EvalSum(splitIndex, Alln);
                    //     // Domi[i] = cc->EvalMult(IndexSum, entry.second);
                    //     Domi[i] = cc->EvalMult(Domi[i], entry.second);
                    //     i++;
                    // }
                    DomRow[j * d + k] = cc->EvalAddMany(Domi);
                    // DomRow[j * d + k] = Domi;
                }
                // comm[j] += Serial::SerializeToString(DomRow[j * d + k]).size();
            }
            double st2 = omp_get_wtime();
            comp[j] += st2 - st1;
            for (uint32_t k = 0; k < d; ++k) {
                comm[j] += Serial::SerializeToString(DomRow[j * d + k]).size();
            }
        }
        // cout << "Anch:" << endl;
        // for (uint32_t j = 0; j < d; ++j) {
        //     cout << posI[j] << endl;
        //     cout << Dec(IV_C[j], 2 * Alln) << endl;
        // }
        // cout << "DomRow:" << endl;
        // for (uint32_t j = 0; j < d; ++j) {
        //     for (uint32_t k = 0; k < d; ++k) {
        //         cout << Dec(DomRow[j * d + k], n) << endl;
        //     }
        //     cout << endl;
        // }
        // compute DomT
        double ss1 = omp_get_wtime();
        for (uint32_t k = 0; k < d; k++) {
            auto Domt = DomRow[k];
            for (uint32_t j = 1; j < d; j++) {
                cc->EvalAddInPlace(Domt, DomRow[j * d + k]);
            }
            DomT[k] = Domt;
        }
        // cout << "DomT:" << endl;
        // for (uint32_t k = 0; k < d; ++k) {
        //     cout << Dec(DomT[k], n) << endl;
        // }
        // compute Dom
        // Ciphertext<DCRTPoly> Dom_C = DomT[0]->Clone();
        // for (uint32_t j = 1; j < d; j++) {
        //     Dom_C = cc->EvalMult(Dom_C, DomT[j]);
        //     cout << Dec(Dom_C->Clone(), n) << endl;
        // }
        Ciphertext<DCRTPoly> Dom_C = cc->EvalMultMany(DomT);
        double ss2                 = omp_get_wtime();
        comp[d] += ss2 - ss1;
        double st3                = omp_get_wtime();
        vector<int64_t> DomValues = Dec(Dom_C, n)->GetPackedValue();
        // cout << "Dom:" << endl;
        // cout << DomValues << endl;
        // index of T
        vector<uint32_t> indT(n);
        uint32_t i0 = 0;
        for (auto entry = T.begin(); entry != T.end(); ++entry) {
            indT[i0++] = entry->first;
        }
        double st4 = omp_get_wtime();
        for (uint32_t j = 0; j < d; ++j)
            comp[j] += st4 - st3;
        for (uint32_t j = 0; j < d + 1; ++j) {
            comm[j] += Serial::SerializeToString(Dom_C).size();
        }
        for (uint32_t i = 0; i < n; ++i) {
            if (DomValues[i] == 0) {
                double ss1           = omp_get_wtime();
                uint32_t keyToRemove = indT[i];
                T.erase(keyToRemove);
                double ss2 = omp_get_wtime();
                comp[d] += ss2 - ss1;
                for (uint32_t j = 0; j < d; ++j) {
                    double st1      = omp_get_wtime();
                    auto& innerMapS = LES[j];
                    innerMapS.erase(keyToRemove);
                    auto& innerMapQ = LESQ[j];
                    innerMapQ.erase(keyToRemove);
                    double st2 = omp_get_wtime();
                    comp[j] += st2 - st1;
                }
            }
        }
        // update LES and LESQ，remove 0 in Dom
        for (uint32_t i = n - 1; i >= 0 && i < n; --i) {
            if (DomValues[i] == 0) {
                for (uint32_t j = 0; j < d; ++j) {
                    double st1 = omp_get_wtime();
                    // auto& innerMapS = LES[j];
                    for (auto& entry : LES[j]) {
                        auto& vectorS = entry.second;
                        vectorS.erase(vectorS.begin() + i);
                    }
                    // auto& innerMapQ = LESQ[j];
                    for (auto& entry : LESQ[j]) {
                        auto& vectorQ = entry.second;
                        vectorQ.erase(vectorQ.begin() + i);
                    }
                    double st2 = omp_get_wtime();
                    comp[j] += st2 - st1;
                }
            }
        }
        for (uint32_t j = 0; j < d; ++j) {
            double st1 = omp_get_wtime();
            for (uint32_t i = 0; i < n; ++i) {
                if (DomValues[i] == 0) {
                    // Update Ind
                    Ind[j][IndMap[j][indT[i]]].visit = 1;
                }
            }
            double st2 = omp_get_wtime();
            comp[j] += st2 - st1;
        }
        uint32_t posT = 0;
        for (uint32_t j = 0; j < d; ++j) {
            double st1      = omp_get_wtime();
            posF[j]         = 0;
            auto& Ind_j     = Ind[j];
            uint32_t size_j = Ind_j.size();
            for (uint32_t i = 0; i < size_j; ++i) {
                posF[j] += Ind_j[i].visit;
            }
            if (posF[j] == size_j)
                ++posT;
            double st2 = omp_get_wtime();
            comp[j] += st2 - st1;
        }
        if (posT == d || T.empty())
            flag = 0;
        // cout << "index" << endl;
        // for (uint32_t j = 0; j < d; j++) {
        //     cout << posF[j] << "\t";
        // }
        // cout << endl;
        // for (uint32_t j = 0; j < d; j++) {
        //     cout << posI[j] << endl;
        //     for (uint32_t i = 0; i < n; i++) {
        //         cout << Ind[j][IndMap[j][indT[i]]][2] << "\t";
        //     }
        //     cout << endl;
        //     for (uint32_t i = 0; i < n; i++) {
        //         cout << Anch[i][j] << "\t";
        //     }
        //     cout << endl;
        // }
    }
    double s04    = omp_get_wtime();
    query_process = s04 - s03;
    cout << "num:" << num << endl;
    uint32_t k = T.size();
    cout << "k:" << k << endl;
    for (uint32_t j = 0; j < d; ++j) {
        cout << "Time_" << j << ":" << fixed << setprecision(2) << (comp[j]) << " s" << endl;
        cout << "Cost_" << j << ":" << fixed << setprecision(2) << (comm[j] / 1024 / 1024) << " MB" << endl;
        if (Utime < comp[j])
            Utime = comp[j];
        if (Ucost < comm[j])
            Ucost = comm[j] / 1024 / 1024;
    }
    cout << "Time_S:" << fixed << setprecision(2) << (comp[d]) << " s" << endl;
    cout << "Cost_S:" << fixed << setprecision(2) << (comm[d] / 1024 / 1024) << " MB" << endl;
    SkyNum = k;
    ItrNum = num;
    ComNum = itr;
    Stime  = comp[d];
    Scost  = comm[d] / 1024 / 1024;
    // for (auto& entry : T) {
    //     for (auto& r : P[entry.first]) {
    //         cout << r << "\t";
    //     }
    //     cout << endl;
    // }
}

void MKsky1(vector<vector<uint64_t>> P) {
    offline_process = 0;
    query_process   = 0;
    map<uint32_t, vector<uint64_t>> T;
    double s01 = omp_get_wtime();
    vector<vector<IndexElement>> Ind(d, vector<IndexElement>(n));
    vector<vector<uint32_t>> IndMap(d, vector<uint32_t>(n));
    vector<vector<int64_t>> FlagF(d, vector<int64_t>(n, 0));
    vector<Plaintext> FlagF_P(d);
    vector<Ciphertext<DCRTPoly>> FlagF_C(d);
    vector<uint32_t> Left(d, 0);
    vector<double> comm(d + 1), comp(d + 1);
    for (uint32_t j = 0; j < d; ++j) {
        vector<ElementIndex> list;
        for (uint32_t i = 0; i < n; ++i) {
            list.emplace_back(P[i][j], i);
        }
        sort(list.begin(), list.end());
        for (uint32_t i = 0; i < n; ++i) {
            Ind[j][i].element = list[i].element;
            Ind[j][i].index   = list[i].index;
            Ind[j][i].visit   = 0;
        }
        for (uint32_t i = 0; i < n; ++i) {
            IndMap[j][list[i].index] = i;
        }
    }
    for (uint32_t i = 0; i < n; ++i) {
        T[i] = P[i];
    }
    // Select
    uint32_t flag = 1;
    // store in column
    map<uint32_t, map<uint32_t, vector<int64_t>>> LES, LESQ;
    for (uint32_t j = 0; j < d; ++j) {
        map<uint32_t, vector<int64_t>> h1, h2;
        for (uint32_t k = 0; k < n; ++k) {
            vector<int64_t> ht1(n), ht2(n);
            for (uint32_t i = 0; i < n; ++i) {
                ht1[i] = (P[k][j] < P[i][j]) ? 0 : 1;
                ht2[i] = (P[k][j] <= P[i][j]) ? 0 : 1;
            }
            h1[k] = move(ht1);
            h2[k] = move(ht2);
        }
        LES[j]  = move(h1);
        LESQ[j] = move(h2);
    }
    unordered_map<uint32_t, Plaintext> splitVec;
    vector<int64_t> atIndex(2 * n, 0);
    atIndex[n] = 1;
    for (uint32_t i = 0; i < n; ++i) {
        splitVec[i]        = cc->MakePackedPlaintext(atIndex);
        atIndex[n + i]     = 0;
        atIndex[n + i + 1] = 1;
    }
    double s02      = omp_get_wtime();
    offline_process = s02 - s01;
    double s03      = omp_get_wtime();
    int itr = 0, num = 0;
    vector<uint32_t> posF(d, 0);
    vector<uint32_t> posI(d, 0);
    vector<Ciphertext<DCRTPoly>> IV_C(d);
    vector<Ciphertext<DCRTPoly>> DomRow(d * d);
    vector<Ciphertext<DCRTPoly>> DomT(d);
    usint Alln = n;
    while (flag == 1) {
        n = T.size();
        cout << "itr:" << itr << "\tn:" << n << endl;
        itr++;
        num += n;
        // cout << "LES" << endl;
        // for (uint32_t j = 0; j < LES.size(); ++j) {
        //     map<uint32_t, vector<int64_t>> LESt = LES[j];
        //     for (auto& inner_pair : LESt) {
        //         auto& values = inner_pair.second;
        //         for (auto value : values) {
        //         cout << value << "\t";
        //         }
        //         cout << endl;
        //     }
        //     cout << endl;
        // }
        // cout << "LESQ" << endl;
        // for (uint32_t j = 0; j < d; ++j) {
        //     map<uint32_t, vector<int64_t>> LESQt = LESQ[j];
        //     for (auto& inner_pair : LESQt) {
        //         auto& values = inner_pair.second;
        //         for (auto value : values) {
        //         cout << value << "\t";
        //         }
        //         cout << endl;
        //     }
        //     cout << endl;
        // }
        double st1 = omp_get_wtime();
        // vector<vector<int64_t>> Anch(d, vector<int64_t>(n, 0));
        // map<uint32_t, vector<int64_t>> Anch;
        map<uint32_t, Plaintext> LES_P;
        map<uint32_t, map<uint32_t, Plaintext>> LESQ_P;
        Ciphertext<DCRTPoly> CZero = cc->Encrypt(keys.publicKey, cc->MakePackedPlaintext(vector<int64_t>(n, 0)));
        double st2                 = omp_get_wtime();
        for (uint32_t j = 0; j < d; ++j)
            comp[j] += st2 - st1;
        srand(time(0));
        for (uint32_t j = 0; j < d; ++j) {
            // compute the index "posI" about min value
            double st1 = omp_get_wtime();
            if (posF[j] == Ind[j].size()) {
                size_t randomIndex = rand() % T.size();
                auto it            = T.begin();
                next(it, randomIndex);
                posI[j] = it->first;
            }
            else {
                auto& row = Ind[j];
                for (uint32_t i = Left[j]; i < row.size(); ++i) {
                    if (row[i].visit == 0) {
                        posI[j]           = row[i].index;
                        row[i].visit      = 1;
                        Left[j]           = i;
                        FlagF[j][posI[j]] = 1;
                        break;
                    }
                }
            }
            // compute the anchor "Anch" in the origin IV
            uint32_t closestIndex = distance(T.begin(), T.find(posI[j]));
            // Anch[j] = move(vector<int64_t>(n, 0));
            vector<int64_t> Ancht(2 * Alln, 0);
            // Anch[j][closestIndex] = 1;
            Ancht[Alln + closestIndex] = 1;
            LES_P[j]                   = cc->MakePackedPlaintext(LES[j][posI[j]]);
            // compute the encrypted "Anch" in the origin IV
            IV_C[j] = cc->Encrypt(keys.publicKey, cc->MakePackedPlaintext(Ancht));
            // cout << j << "-IV:" << Serial::SerializeToString(IV_C[j]).size() << endl;
            // comm[j] += Serial::SerializeToString(IV_C[j]).size();
            // packing the column value about LESQ
            map<uint32_t, vector<int64_t>> LESQt = LESQ[j];
            map<uint32_t, Plaintext> LESQ_Pt;
            uint32_t i = 0;
            for (auto& inner_pair : LESQt) {
                LESQ_Pt[i] = cc->MakePackedPlaintext(inner_pair.second);
                i++;
            }
            LESQ_P[j]  = move(LESQ_Pt);
            double st2 = omp_get_wtime();
            comp[j] += st2 - st1;
            size_t tmp = Serial::SerializeToString(IV_C[j]).size();
            comm[j] += tmp;
            comm[d] += tmp;
        }
        // compute the inner product between matrix and IV
        for (uint32_t j = 0; j < d; ++j) {
            double st1 = omp_get_wtime();
            for (uint32_t k = 0; k < d; ++k) {
                if (j == k) {
                    DomRow[j * d + k] = cc->Encrypt(keys.publicKey, LES_P[j]);
                }
                else {
                    const auto& innerMap = LESQ_P[j];
                    // Ciphertext<DCRTPoly> Domi = CZero->Clone();
                    vector<Ciphertext<DCRTPoly>> Domi(n);
                    uint32_t i = 0;
                    for (const auto& entry : innerMap) {
                        // first try
                        Ciphertext<DCRTPoly> splitInd = cc->EvalMult(IV_C[k], splitVec[i]);
                        Ciphertext<DCRTPoly> IndexSum = cc->EvalSum(splitInd, batchSize);
                        // second try
                        // Ciphertext<DCRTPoly> splitInd   = cc->EvalInnerProduct(IV_C[k], splitVec[i], n);
                        // Ciphertext<DCRTPoly> Indsum     = cc->EvalMult(splitInd, splitVec[0]);
                        // Ciphertext<DCRTPoly> splitIndex = cc->EvalRotate(splitInd, indexList[0]);
                        // Ciphertext<DCRTPoly> IndexSum   = cc->EvalSum(splitIndex, Alln);
                        Domi[i] = cc->EvalMult(IndexSum, entry.second);
                        i++;
                    }
                    DomRow[j * d + k] = cc->EvalAddMany(Domi);
                    // DomRow[j * d + k] = Domi;
                }
                // comm[j] += Serial::SerializeToString(DomRow[j * d + k]).size();
            }
            double st2 = omp_get_wtime();
            comp[j] += st2 - st1;
            for (uint32_t k = 0; k < d; ++k) {
                comm[j] += Serial::SerializeToString(DomRow[j * d + k]).size();
            }
        }
        // cout << "Anch:" << endl;
        // for (uint32_t j = 0; j < d; ++j) {
        //     cout << posI[j] << endl;
        //     cout << Dec(IV_C[j], 2 * Alln) << endl;
        // }
        // cout << "DomRow:" << endl;
        // for (uint32_t j = 0; j < d; ++j) {
        //     for (uint32_t k = 0; k < d; ++k) {
        //         cout << Dec(DomRow[j * d + k], n) << endl;
        //     }
        //     cout << endl;
        // }
        // compute DomT
        double ss1 = omp_get_wtime();
        for (uint32_t k = 0; k < d; k++) {
            auto Domt = DomRow[k];
            for (uint32_t j = 1; j < d; j++) {
                cc->EvalAddInPlace(Domt, DomRow[j * d + k]);
            }
            DomT[k] = Domt;
        }
        // cout << "DomT:" << endl;
        // for (uint32_t k = 0; k < d; ++k) {
        //     cout << Dec(DomT[k], n) << endl;
        // }
        // compute Dom
        // Ciphertext<DCRTPoly> Dom_C = DomT[0]->Clone();
        // for (uint32_t j = 1; j < d; j++) {
        //     Dom_C = cc->EvalMult(Dom_C, DomT[j]);
        //     cout << Dec(Dom_C->Clone(), n) << endl;
        // }
        Ciphertext<DCRTPoly> Dom_C = cc->EvalMultMany(DomT);
        double ss2                 = omp_get_wtime();
        comp[d] += ss2 - ss1;
        double st3                = omp_get_wtime();
        vector<int64_t> DomValues = Dec(Dom_C, n)->GetPackedValue();
        // cout << "Dom:" << endl;
        // cout << DomValues << endl;
        // index of T
        vector<uint32_t> indT(n);
        uint32_t i0 = 0;
        for (auto entry = T.begin(); entry != T.end(); ++entry) {
            indT[i0++] = entry->first;
        }
        double st4 = omp_get_wtime();
        for (uint32_t j = 0; j < d; ++j)
            comp[j] += st4 - st3;
        for (uint32_t j = 0; j < d + 1; ++j) {
            // joint decryption requires communication
            comm[j] += Serial::SerializeToString(Dom_C).size();
        }
        for (uint32_t i = 0; i < n; ++i) {
            if (DomValues[i] == 0) {
                double ss1           = omp_get_wtime();
                uint32_t keyToRemove = indT[i];
                T.erase(keyToRemove);
                double ss2 = omp_get_wtime();
                comp[d] += ss2 - ss1;
                for (uint32_t j = 0; j < d; ++j) {
                    double st1      = omp_get_wtime();
                    auto& innerMapS = LES[j];
                    innerMapS.erase(keyToRemove);
                    auto& innerMapQ = LESQ[j];
                    innerMapQ.erase(keyToRemove);
                    double st2 = omp_get_wtime();
                    comp[j] += st2 - st1;
                }
            }
        }
        // update LES and LESQ，remove 0 in Dom
        for (uint32_t i = n - 1; i >= 0 && i < n; --i) {
            if (DomValues[i] == 0) {
                for (uint32_t j = 0; j < d; ++j) {
                    double st1 = omp_get_wtime();
                    // auto& innerMapS = LES[j];
                    for (auto& entry : LES[j]) {
                        auto& vectorS = entry.second;
                        vectorS.erase(vectorS.begin() + i);
                    }
                    // auto& innerMapQ = LESQ[j];
                    for (auto& entry : LESQ[j]) {
                        auto& vectorQ = entry.second;
                        vectorQ.erase(vectorQ.begin() + i);
                    }
                    double st2 = omp_get_wtime();
                    comp[j] += st2 - st1;
                }
            }
        }
        for (uint32_t j = 0; j < d; ++j) {
            double st1 = omp_get_wtime();
            for (uint32_t i = 0; i < n; ++i) {
                if (DomValues[i] == 0) {
                    // Update Ind
                    Ind[j][IndMap[j][indT[i]]].visit = 1;
                }
            }
            double st2 = omp_get_wtime();
            comp[j] += st2 - st1;
        }
        // double Flag1 = omp_get_wtime();
        for (uint32_t j = 0; j < d; ++j) {
            FlagF_P[j] = cc->MakePackedPlaintext(FlagF[j]);
            FlagF_C[j] = cc->Encrypt(keys.publicKey, FlagF_P[j]);
        }
        auto FlagF_ALL    = cc->EvalMultMany(FlagF_C);
        int64_t FlagFinal = Dec(cc->EvalSum(FlagF_ALL, Alln), Alln)->GetPackedValue()[0];
        // cout << "Flag:" << FlagFinal;
        // double Flag2 = omp_get_wtime();
        // cout << "\tTime_Flag:" << (Flag2 - Flag1) << "ms" << endl;
        if (FlagFinal >= 1)
            break;
        uint32_t posT = 0;
        for (uint32_t j = 0; j < d; ++j) {
            double st1      = omp_get_wtime();
            posF[j]         = 0;
            auto& Ind_j     = Ind[j];
            uint32_t size_j = Ind_j.size();
            for (uint32_t i = 0; i < size_j; ++i) {
                posF[j] += Ind_j[i].visit;
            }
            if (posF[j] == size_j)
                ++posT;
            double st2 = omp_get_wtime();
            comp[j] += st2 - st1;
        }
        if (posT == d || T.empty())
            flag = 0;
        // cout << "index" << endl;
        // for (uint32_t j = 0; j < d; j++) {
        //     cout << posF[j] << "\t";
        // }
        // cout << endl;
        // for (uint32_t j = 0; j < d; j++) {
        //     cout << posI[j] << endl;
        //     for (uint32_t i = 0; i < n; i++) {
        //         cout << Ind[j][IndMap[j][indT[i]]][2] << "\t";
        //     }
        //     cout << endl;
        //     for (uint32_t i = 0; i < n; i++) {
        //         cout << Anch[i][j] << "\t";
        //     }
        //     cout << endl;
        // }
    }
    double s04    = omp_get_wtime();
    query_process = s04 - s03;
    cout << "num:" << num << endl;
    uint32_t k = T.size();
    cout << "k:" << k << endl;
    for (uint32_t j = 0; j < d; ++j) {
        cout << "Time_" << j << ":" << fixed << setprecision(2) << (comp[j]) << " s" << endl;
        cout << "Cost_" << j << ":" << fixed << setprecision(2) << (comm[j] / 1024 / 1024) << " MB" << endl;
        if (Utime < comp[j])
            Utime = comp[j];
        if (Ucost < comm[j])
            Ucost = comm[j] / 1024 / 1024;
    }
    cout << "Time_S:" << fixed << setprecision(2) << (comp[d]) << " s" << endl;
    cout << "Cost_S:" << fixed << setprecision(2) << (comm[d] / 1024 / 1024) << " MB" << endl;
    SkyNum = k;
    ItrNum = num;
    ComNum = itr;
    Stime  = comp[d];
    Scost  = comm[d] / 1024 / 1024;
    // for (auto& entry : T) {
    //     for (auto& r : P[entry.first]) {
    //         cout << r << "\t";
    //     }
    //     cout << endl;
    // }
}

void MKsky1T(vector<vector<uint64_t>> P) {
    offline_process = 0;
    query_process   = 0;
    map<uint32_t, vector<uint64_t>> T;
    double s01 = omp_get_wtime();
    vector<vector<IndexElement>> Ind(d, vector<IndexElement>(n));
    vector<vector<uint32_t>> IndMap(d, vector<uint32_t>(n));
    vector<vector<int64_t>> FlagF(d, vector<int64_t>(n, 0));
    vector<Plaintext> FlagF_P(d);
    vector<Ciphertext<DCRTPoly>> FlagF_C(d);
    vector<uint32_t> Left(d, 0);
    vector<double> comm(d + 1), comp(d + 1);
    for (uint32_t j = 0; j < d; ++j) {
        vector<ElementIndex> list;
        for (uint32_t i = 0; i < n; ++i) {
            list.emplace_back(P[i][j], i);
        }
        sort(list.begin(), list.end());
        for (uint32_t i = 0; i < n; ++i) {
            Ind[j][i].element = list[i].element;
            Ind[j][i].index   = list[i].index;
            Ind[j][i].visit   = 0;
        }
        for (uint32_t i = 0; i < n; ++i) {
            IndMap[j][list[i].index] = i;
        }
    }
    for (uint32_t i = 0; i < n; ++i) {
        T[i] = P[i];
    }
    // Select
    uint32_t flag = 1;
    // store in column
    map<uint32_t, map<uint32_t, vector<int64_t>>> LES, LESQ;
    for (uint32_t j = 0; j < d; ++j) {
        map<uint32_t, vector<int64_t>> h1, h2;
        for (uint32_t k = 0; k < n; ++k) {
            vector<int64_t> ht1(n), ht2(n);
            for (uint32_t i = 0; i < n; ++i) {
                ht1[i] = (P[k][j] < P[i][j]) ? 0 : 1;
                ht2[i] = (P[k][j] <= P[i][j]) ? 0 : 1;
            }
            h1[k] = move(ht1);
            h2[k] = move(ht2);
        }
        LES[j]  = move(h1);
        LESQ[j] = move(h2);
    }
    unordered_map<uint32_t, Plaintext> splitVec;
    for (uint32_t i = 0; i < n; ++i) {
        vector<int64_t> atIndex(2 * n, 0);
        atIndex[n + i] = 1;
        splitVec[i]    = cc->MakePackedPlaintext(atIndex);
    }
    double s02      = omp_get_wtime();
    offline_process = s02 - s01;
    double s03      = omp_get_wtime();
    int itr = 0, num = 0;
    vector<uint32_t> posF(d, 0);
    vector<uint32_t> posI(d, 0);
    vector<Ciphertext<DCRTPoly>> IV_C(d);
    vector<Ciphertext<DCRTPoly>> DomRow(d * d);
    vector<Ciphertext<DCRTPoly>> DomT(d);
    usint Alln = n;
    while (flag == 1) {
        n = T.size();
        cout << "itr:" << itr << "\tn:" << n << endl;
        itr++;
        num += n;
        // cout << "LES" << endl;
        // for (uint32_t j = 0; j < LES.size(); ++j) {
        //     map<uint32_t, vector<int64_t>> LESt = LES[j];
        //     for (auto& inner_pair : LESt) {
        //         auto& values = inner_pair.second;
        //         for (auto value : values) {
        //         cout << value << "\t";
        //         }
        //         cout << endl;
        //     }
        //     cout << endl;
        // }
        // cout << "LESQ" << endl;
        // for (uint32_t j = 0; j < d; ++j) {
        //     map<uint32_t, vector<int64_t>> LESQt = LESQ[j];
        //     for (auto& inner_pair : LESQt) {
        //         auto& values = inner_pair.second;
        //         for (auto value : values) {
        //         cout << value << "\t";
        //         }
        //         cout << endl;
        //     }
        //     cout << endl;
        // }
        double st1 = omp_get_wtime();
        // vector<vector<int64_t>> Anch(d, vector<int64_t>(n, 0));
        // map<uint32_t, vector<int64_t>> Anch;
        map<uint32_t, Plaintext> LES_P;
        map<uint32_t, map<uint32_t, Plaintext>> LESQ_P;
        // Ciphertext<DCRTPoly> CZero = cc->Encrypt(keys.publicKey, cc->MakePackedPlaintext(vector<int64_t>(n, 0)));
        double st2                 = omp_get_wtime();
        for (uint32_t j = 0; j < d; ++j)
            comp[j] += st2 - st1;
        srand(time(0));
        for (uint32_t j = 0; j < d; ++j) {
            // compute the index "posI" about min value
            double st1 = omp_get_wtime();
            if (posF[j] == Ind[j].size()) {
                size_t randomIndex = rand() % T.size();
                auto it            = T.begin();
                next(it, randomIndex);
                posI[j] = it->first;
            }
            else {
                auto& row = Ind[j];
                for (uint32_t i = Left[j]; i < row.size(); ++i) {
                    if (row[i].visit == 0) {
                        posI[j]           = row[i].index;
                        row[i].visit      = 1;
                        Left[j]           = i;
                        FlagF[j][posI[j]] = 1;
                        break;
                    }
                }
            }
            // compute the anchor "Anch" in the origin IV
            uint32_t closestIndex = distance(T.begin(), T.find(posI[j]));
            // Anch[j] = move(vector<int64_t>(n, 0));
            vector<int64_t> Ancht(2 * Alln, 0);
            // Anch[j][closestIndex] = 1;
            Ancht[Alln + closestIndex] = 1;
            LES_P[j]                   = cc->MakePackedPlaintext(LES[j][posI[j]]);
            // compute the encrypted "Anch" in the origin IV
            IV_C[j] = cc->Encrypt(keys.publicKey, cc->MakePackedPlaintext(Ancht));
            // cout << j << "-IV:" << Serial::SerializeToString(IV_C[j]).size() << endl;
            // comm[j] += Serial::SerializeToString(IV_C[j]).size();
            // packing the column value about LESQ
            map<uint32_t, vector<int64_t>> LESQt = LESQ[j];
            map<uint32_t, Plaintext> LESQ_Pt;
            uint32_t i = 0;
            for (auto& inner_pair : LESQt) {
                LESQ_Pt[i] = cc->MakePackedPlaintext(inner_pair.second);
                i++;
            }
            LESQ_P[j]  = move(LESQ_Pt);
            double st2 = omp_get_wtime();
            comp[j] += st2 - st1;
            size_t tmp = Serial::SerializeToString(IV_C[j]).size();
            comm[j] += tmp;
            comm[d] += tmp;
        }
        // compute the inner product between matrix and IV
        for (uint32_t j = 0; j < d; ++j) {
            double st1 = omp_get_wtime();
            for (uint32_t k = 0; k < d; ++k) {
                if (j == k) {
                    DomRow[j * d + k] = cc->Encrypt(keys.publicKey, LES_P[j]);
                }
                else {
                    // const auto& innerMap = LESQ_P[j];
                    // Ciphertext<DCRTPoly> Domi = CZero->Clone();
                    vector<Ciphertext<DCRTPoly>> Domi(n);
                    vector<Plaintext> LESQt;
                    LESQt.reserve(n);
                    for (const auto& entry : LESQ_P[j]) {
                        LESQt.push_back(entry.second);
                    }
#pragma omp parallel for
                    for (uint32_t i = 0; i < n; ++i) {
                        Ciphertext<DCRTPoly> splitInd = cc->EvalMult(IV_C[k], splitVec[i]);
                        Domi[i] = cc->EvalSum(splitInd, batchSize);
                        Domi[i] = cc->EvalMult(Domi[i], LESQt[i]);
                    }
                    // uint32_t i = 0;
                    // for (const auto& entry : innerMap) {
                    //     // // first try
                    //     // Ciphertext<DCRTPoly> splitInd = cc->EvalMult(IV_C[k], splitVec[i]);
                    //     // Ciphertext<DCRTPoly> IndexSum = cc->EvalSum(splitInd, batchSize);
                    //     // second try
                    //     // Ciphertext<DCRTPoly> splitInd   = cc->EvalInnerProduct(IV_C[k], splitVec[i], n);
                    //     // Ciphertext<DCRTPoly> Indsum     = cc->EvalMult(splitInd, splitVec[0]);
                    //     // Ciphertext<DCRTPoly> splitIndex = cc->EvalRotate(splitInd, indexList[0]);
                    //     // Ciphertext<DCRTPoly> IndexSum   = cc->EvalSum(splitIndex, Alln);
                    //     // Domi[i] = cc->EvalMult(IndexSum, entry.second);
                    //     Domi[i] = cc->EvalMult(Domi[i], entry.second);
                    //     i++;
                    // }
                    DomRow[j * d + k] = cc->EvalAddMany(Domi);
                    // DomRow[j * d + k] = Domi;
                }
                // comm[j] += Serial::SerializeToString(DomRow[j * d + k]).size();
            }
            double st2 = omp_get_wtime();
            comp[j] += st2 - st1;
            for (uint32_t k = 0; k < d; ++k) {
                comm[j] += Serial::SerializeToString(DomRow[j * d + k]).size();
            }
        }
        // cout << "Anch:" << endl;
        // for (uint32_t j = 0; j < d; ++j) {
        //     cout << posI[j] << endl;
        //     cout << Dec(IV_C[j], 2 * Alln) << endl;
        // }
        // cout << "DomRow:" << endl;
        // for (uint32_t j = 0; j < d; ++j) {
        //     for (uint32_t k = 0; k < d; ++k) {
        //         cout << Dec(DomRow[j * d + k], n) << endl;
        //     }
        //     cout << endl;
        // }
        // compute DomT
        double ss1 = omp_get_wtime();
        for (uint32_t k = 0; k < d; k++) {
            auto Domt = DomRow[k];
            for (uint32_t j = 1; j < d; j++) {
                cc->EvalAddInPlace(Domt, DomRow[j * d + k]);
            }
            DomT[k] = Domt;
        }
        // cout << "DomT:" << endl;
        // for (uint32_t k = 0; k < d; ++k) {
        //     cout << Dec(DomT[k], n) << endl;
        // }
        // compute Dom
        // Ciphertext<DCRTPoly> Dom_C = DomT[0]->Clone();
        // for (uint32_t j = 1; j < d; j++) {
        //     Dom_C = cc->EvalMult(Dom_C, DomT[j]);
        //     cout << Dec(Dom_C->Clone(), n) << endl;
        // }
        Ciphertext<DCRTPoly> Dom_C = cc->EvalMultMany(DomT);
        double ss2                 = omp_get_wtime();
        comp[d] += ss2 - ss1;
        double st3                = omp_get_wtime();
        vector<int64_t> DomValues = Dec(Dom_C, n)->GetPackedValue();
        // cout << "Dom:" << endl;
        // cout << DomValues << endl;
        // index of T
        vector<uint32_t> indT(n);
        uint32_t i0 = 0;
        for (auto entry = T.begin(); entry != T.end(); ++entry) {
            indT[i0++] = entry->first;
        }
        double st4 = omp_get_wtime();
        for (uint32_t j = 0; j < d; ++j)
            comp[j] += st4 - st3;
        for (uint32_t j = 0; j < d + 1; ++j) {
            // joint decryption requires communication
            comm[j] += Serial::SerializeToString(Dom_C).size();
        }
        for (uint32_t i = 0; i < n; ++i) {
            if (DomValues[i] == 0) {
                double ss1           = omp_get_wtime();
                uint32_t keyToRemove = indT[i];
                T.erase(keyToRemove);
                double ss2 = omp_get_wtime();
                comp[d] += ss2 - ss1;
                for (uint32_t j = 0; j < d; ++j) {
                    double st1      = omp_get_wtime();
                    auto& innerMapS = LES[j];
                    innerMapS.erase(keyToRemove);
                    auto& innerMapQ = LESQ[j];
                    innerMapQ.erase(keyToRemove);
                    double st2 = omp_get_wtime();
                    comp[j] += st2 - st1;
                }
            }
        }
        // update LES and LESQ，remove 0 in Dom
        for (uint32_t i = n - 1; i >= 0 && i < n; --i) {
            if (DomValues[i] == 0) {
                for (uint32_t j = 0; j < d; ++j) {
                    double st1 = omp_get_wtime();
                    // auto& innerMapS = LES[j];
                    for (auto& entry : LES[j]) {
                        auto& vectorS = entry.second;
                        vectorS.erase(vectorS.begin() + i);
                    }
                    // auto& innerMapQ = LESQ[j];
                    for (auto& entry : LESQ[j]) {
                        auto& vectorQ = entry.second;
                        vectorQ.erase(vectorQ.begin() + i);
                    }
                    double st2 = omp_get_wtime();
                    comp[j] += st2 - st1;
                }
            }
        }
        for (uint32_t j = 0; j < d; ++j) {
            double st1 = omp_get_wtime();
            for (uint32_t i = 0; i < n; ++i) {
                if (DomValues[i] == 0) {
                    // Update Ind
                    Ind[j][IndMap[j][indT[i]]].visit = 1;
                }
            }
            double st2 = omp_get_wtime();
            comp[j] += st2 - st1;
        }
        double Flag1 = omp_get_wtime();
        for (uint32_t j = 0; j < d; ++j) {
            FlagF_P[j] = cc->MakePackedPlaintext(FlagF[j]);
            FlagF_C[j] = cc->Encrypt(keys.publicKey, FlagF_P[j]);
        }
        auto FlagF_ALL    = cc->EvalMultMany(FlagF_C);
        int64_t FlagFinal = Dec(cc->EvalSum(FlagF_ALL, Alln), Alln)->GetPackedValue()[0];
        // cout << "Flag:" << FlagFinal;
        double Flag2 = omp_get_wtime();
        for (uint32_t j = 0; j < d; ++j)
            comp[j] += Flag2 - Flag1;
        for (uint32_t j = 0; j < d + 1; ++j) {
            comm[j] += Serial::SerializeToString(FlagF_ALL).size();
        }
        // cout << "\tTime_Flag:" << (Flag2 - Flag1) << "ms" << endl;
        if (FlagFinal >= 1)
            break;
        uint32_t posT = 0;
        for (uint32_t j = 0; j < d; ++j) {
            double st1      = omp_get_wtime();
            posF[j]         = 0;
            auto& Ind_j     = Ind[j];
            uint32_t size_j = Ind_j.size();
            for (uint32_t i = 0; i < size_j; ++i) {
                posF[j] += Ind_j[i].visit;
            }
            if (posF[j] == size_j)
                ++posT;
            double st2 = omp_get_wtime();
            comp[j] += st2 - st1;
        }
        if (posT == d || T.empty())
            flag = 0;
        // cout << "index" << endl;
        // for (uint32_t j = 0; j < d; j++) {
        //     cout << posF[j] << "\t";
        // }
        // cout << endl;
        // for (uint32_t j = 0; j < d; j++) {
        //     cout << posI[j] << endl;
        //     for (uint32_t i = 0; i < n; i++) {
        //         cout << Ind[j][IndMap[j][indT[i]]][2] << "\t";
        //     }
        //     cout << endl;
        //     for (uint32_t i = 0; i < n; i++) {
        //         cout << Anch[i][j] << "\t";
        //     }
        //     cout << endl;
        // }
    }
    double s04    = omp_get_wtime();
    query_process = s04 - s03;
    cout << "num:" << num << endl;
    uint32_t k = T.size();
    cout << "k:" << k << endl;
    for (uint32_t j = 0; j < d; ++j) {
        cout << "Time_" << j << ":" << fixed << setprecision(2) << (comp[j]) << " s" << endl;
        cout << "Cost_" << j << ":" << fixed << setprecision(2) << (comm[j] / 1024 / 1024) << " MB" << endl;
        if (Utime < comp[j])
            Utime = comp[j];
        if (Ucost < comm[j])
            Ucost = comm[j] / 1024 / 1024;
    }
    cout << "Time_S:" << fixed << setprecision(2) << (comp[d]) << " s" << endl;
    cout << "Cost_S:" << fixed << setprecision(2) << (comm[d] / 1024 / 1024) << " MB" << endl;
    SkyNum = k;
    ItrNum = num;
    ComNum = itr;
    Stime  = comp[d];
    Scost  = comm[d] / 1024 / 1024;
    // for (auto& entry : T) {
    //     for (auto& r : P[entry.first]) {
    //         cout << r << "\t";
    //     }
    //     cout << endl;
    // }
}

void MKsky12(vector<vector<uint64_t>> P) {
    offline_process = 0;
    query_process   = 0;
    map<uint32_t, vector<uint64_t>> T;
    double s01 = omp_get_wtime();
    vector<vector<IndexElement>> Ind(d, vector<IndexElement>(n));
    vector<vector<uint32_t>> IndMap(d, vector<uint32_t>(n));
    vector<vector<int64_t>> FlagF(d, vector<int64_t>(n, 0));
    vector<Plaintext> FlagF_P(d);
    vector<Ciphertext<DCRTPoly>> FlagF_C(d);
    vector<uint32_t> Left(d, 0);
    vector<double> comm(d + 1), comp(d + 1);
    for (uint32_t j = 0; j < d; ++j) {
        vector<ElementIndex> list;
        for (uint32_t i = 0; i < n; ++i) {
            list.emplace_back(P[i][j], i);
        }
        sort(list.begin(), list.end());
        for (uint32_t i = 0; i < n; ++i) {
            Ind[j][i].element = list[i].element;
            Ind[j][i].index   = list[i].index;
            Ind[j][i].visit   = 0;
        }
        for (uint32_t i = 0; i < n; ++i) {
            IndMap[j][list[i].index] = i;
        }
    }
    for (uint32_t i = 0; i < n; ++i) {
        T[i] = P[i];
    }
    // Select
    uint32_t flag = 1;
    // store in column
    map<uint32_t, map<uint32_t, vector<int64_t>>> LES, LESQ;
    for (uint32_t j = 0; j < d; ++j) {
        map<uint32_t, vector<int64_t>> h1, h2;
        for (uint32_t k = 0; k < n; ++k) {
            vector<int64_t> ht1(n), ht2(n);
            for (uint32_t i = 0; i < n; ++i) {
                ht1[i] = (P[k][j] < P[i][j]) ? 0 : 1;
                ht2[i] = (P[k][j] <= P[i][j]) ? 0 : 1;
            }
            h1[k] = move(ht1);
            h2[k] = move(ht2);
        }
        LES[j]  = move(h1);
        LESQ[j] = move(h2);
    }
    unordered_map<uint32_t, Plaintext> splitVec;
    vector<int64_t> atIndex(2 * n, 0);
    atIndex[n] = 1;
    for (uint32_t i = 0; i < n; ++i) {
        splitVec[i]        = cc->MakePackedPlaintext(atIndex);
        atIndex[n + i]     = 0;
        atIndex[n + i + 1] = 1;
    }
    double s02      = omp_get_wtime();
    offline_process = s02 - s01;
    // cout << "offline_process time: " << fixed << setprecision(2) << offline_process << " s" << endl;
    double s03 = omp_get_wtime();
    int itr = 0, num = 0;
    vector<uint32_t> posF(d, 0);
    vector<uint32_t> posI(d, 0);
    vector<Ciphertext<DCRTPoly>> IV_C(d);
    vector<Ciphertext<DCRTPoly>> DomRow(d * d);
    vector<Ciphertext<DCRTPoly>> DomT(d);
    usint Alln = n;
    NumItr.resize(n + 1, 0);
    while (flag == 1) {
        n = T.size();
        cout << "itr:" << itr << "\tn:" << n << endl;
        itr++;
        num += n;
        NumItr[n]++;
        if (NumItr[n] > ts)
            break;
        // cout << "LES" << endl;
        // for (uint32_t j = 0; j < LES.size(); ++j) {
        //     map<uint32_t, vector<int64_t>> LESt = LES[j];
        //     for (auto& inner_pair : LESt) {
        //         auto& values = inner_pair.second;
        //         for (auto value : values) {
        //         cout << value << "\t";
        //         }
        //         cout << endl;
        //     }
        //     cout << endl;
        // }
        // cout << "LESQ" << endl;
        // for (uint32_t j = 0; j < d; ++j) {
        //     map<uint32_t, vector<int64_t>> LESQt = LESQ[j];
        //     for (auto& inner_pair : LESQt) {
        //         auto& values = inner_pair.second;
        //         for (auto value : values) {
        //         cout << value << "\t";
        //         }
        //         cout << endl;
        //     }
        //     cout << endl;
        // }
        double st1 = omp_get_wtime();
        // vector<vector<int64_t>> Anch(d, vector<int64_t>(n, 0));
        // map<uint32_t, vector<int64_t>> Anch;
        map<uint32_t, Plaintext> LES_P;
        map<uint32_t, map<uint32_t, Plaintext>> LESQ_P;
        Ciphertext<DCRTPoly> CZero = cc->Encrypt(keys.publicKey, cc->MakePackedPlaintext(vector<int64_t>(n, 0)));
        double st2                 = omp_get_wtime();
        for (uint32_t j = 0; j < d; ++j)
            comp[j] += st2 - st1;
        srand(time(0));
        for (uint32_t j = 0; j < d; ++j) {
            // compute the index "posI" about min value
            double st1 = omp_get_wtime();
            if (posF[j] == Ind[j].size()) {
                size_t randomIndex = rand() % T.size();
                auto it            = T.begin();
                next(it, randomIndex);
                posI[j] = it->first;
            }
            else {
                auto& row = Ind[j];
                for (uint32_t i = Left[j]; i < row.size(); ++i) {
                    if (row[i].visit == 0) {
                        posI[j]           = row[i].index;
                        row[i].visit      = 1;
                        Left[j]           = i;
                        FlagF[j][posI[j]] = 1;
                        break;
                    }
                }
            }
            // compute the anchor "Anch" in the origin IV
            uint32_t closestIndex = distance(T.begin(), T.find(posI[j]));
            // Anch[j] = move(vector<int64_t>(n, 0));
            vector<int64_t> Ancht(2 * Alln, 0);
            // Anch[j][closestIndex] = 1;
            Ancht[Alln + closestIndex] = 1;
            LES_P[j]                   = cc->MakePackedPlaintext(LES[j][posI[j]]);
            // compute the encrypted "Anch" in the origin IV
            IV_C[j] = cc->Encrypt(keys.publicKey, cc->MakePackedPlaintext(Ancht));
            // cout << j << "-IV:" << Serial::SerializeToString(IV_C[j]).size() << endl;
            // comm[j] += Serial::SerializeToString(IV_C[j]).size();
            // packing the column value about LESQ
            map<uint32_t, vector<int64_t>> LESQt = LESQ[j];
            map<uint32_t, Plaintext> LESQ_Pt;
            uint32_t i = 0;
            for (auto& inner_pair : LESQt) {
                LESQ_Pt[i] = cc->MakePackedPlaintext(inner_pair.second);
                i++;
            }
            LESQ_P[j]  = move(LESQ_Pt);
            double st2 = omp_get_wtime();
            comp[j] += st2 - st1;
            size_t tmp = Serial::SerializeToString(IV_C[j]).size();
            comm[j] += tmp;
            comm[d] += tmp;
        }
        // compute the inner product between matrix and IV
        for (uint32_t j = 0; j < d; ++j) {
            double st1 = omp_get_wtime();
            for (uint32_t k = 0; k < d; ++k) {
                if (j == k) {
                    DomRow[j * d + k] = cc->Encrypt(keys.publicKey, LES_P[j]);
                }
                else {
                    const auto& innerMap = LESQ_P[j];
                    // Ciphertext<DCRTPoly> Domi = CZero->Clone();
                    vector<Ciphertext<DCRTPoly>> Domi(n);
                    uint32_t i = 0;
                    for (const auto& entry : innerMap) {
                        // first try
                        Ciphertext<DCRTPoly> splitInd = cc->EvalMult(IV_C[k], splitVec[i]);
                        Ciphertext<DCRTPoly> IndexSum = cc->EvalSum(splitInd, batchSize);
                        // second try
                        // Ciphertext<DCRTPoly> splitInd   = cc->EvalInnerProduct(IV_C[k], splitVec[i], n);
                        // Ciphertext<DCRTPoly> Indsum     = cc->EvalMult(splitInd, splitVec[0]);
                        // Ciphertext<DCRTPoly> splitIndex = cc->EvalRotate(splitInd, indexList[0]);
                        // Ciphertext<DCRTPoly> IndexSum   = cc->EvalSum(splitIndex, Alln);
                        // cout << j << "," << k << "-1:" << Dec(splitInd, 2 * Alln) << endl;
                        // cout << j << "," << k << "-2:" << Dec(IndexSum, 2 * Alln) << endl;
                        Domi[i] = cc->EvalMult(IndexSum, entry.second);
                        i++;
                    }
                    DomRow[j * d + k] = cc->EvalAddMany(Domi);
                    // DomRow[j * d + k] = Domi;
                }
                // comm[j] += Serial::SerializeToString(DomRow[j * d + k]).size();
            }
            double st2 = omp_get_wtime();
            comp[j] += st2 - st1;
            for (uint32_t k = 0; k < d; ++k) {
                comm[j] += Serial::SerializeToString(DomRow[j * d + k]).size();
            }
        }
        // cout << "Anch:" << endl;
        // for (uint32_t j = 0; j < d; ++j) {
        //     cout << posI[j] << endl;
        //     cout << Dec(IV_C[j], 2 * Alln) << endl;
        // }
        // cout << "DomRow:" << endl;
        // for (uint32_t j = 0; j < d; ++j) {
        //     for (uint32_t k = 0; k < d; ++k) {
        //         cout << Dec(DomRow[j * d + k], n) << endl;
        //     }
        //     cout << endl;
        // }
        // compute DomT
        double ss1 = omp_get_wtime();
        for (uint32_t k = 0; k < d; k++) {
            // DomT[k] = DomRow[k];
            // for (uint32_t j = 1; j < d; j++) {
            //     cc->EvalAddInPlace(DomT[k], DomRow[j * d + k]);
            // }
            auto Domt = DomRow[k];
            for (uint32_t j = 1; j < d; j++) {
                cc->EvalAddInPlace(Domt, DomRow[j * d + k]);
            }
            DomT[k] = Domt;
        }
        // cout << "DomT:" << endl;
        // for (uint32_t k = 0; k < d; ++k) {
        //     cout << Dec(DomT[k], n) << endl;
        // }
        // compute Dom
        // Ciphertext<DCRTPoly> Dom_C = DomT[0]->Clone();
        // for (uint32_t j = 1; j < d; j++) {
        //     Dom_C = cc->EvalMult(Dom_C, DomT[j]);
        //     cout << Dec(Dom_C->Clone(), n) << endl;
        // }
        Ciphertext<DCRTPoly> Dom_C = cc->EvalMultMany(DomT);
        double ss2                 = omp_get_wtime();
        comp[d] += ss2 - ss1;
        double st3                = omp_get_wtime();
        vector<int64_t> DomValues = Dec(Dom_C, n)->GetPackedValue();
        // cout << "Dom:" << endl;
        // cout << DomValues << endl;
        // cout << endl;
        // index of T
        vector<uint32_t> indT(n);
        uint32_t i0 = 0;
        for (auto entry = T.begin(); entry != T.end(); ++entry) {
            indT[i0++] = entry->first;
        }
        double st4 = omp_get_wtime();
        for (uint32_t j = 0; j < d; ++j)
            comp[j] += st4 - st3;
        for (uint32_t j = 0; j < d + 1; ++j) {
            comm[j] += Serial::SerializeToString(Dom_C).size();
        }
        for (uint32_t i = 0; i < n; ++i) {
            if (DomValues[i] == 0) {
                double ss1           = omp_get_wtime();
                uint32_t keyToRemove = indT[i];
                T.erase(keyToRemove);
                double ss2 = omp_get_wtime();
                comp[d] += ss2 - ss1;
                for (uint32_t j = 0; j < d; ++j) {
                    double st1      = omp_get_wtime();
                    auto& innerMapS = LES[j];
                    innerMapS.erase(keyToRemove);
                    auto& innerMapQ = LESQ[j];
                    innerMapQ.erase(keyToRemove);
                    double st2 = omp_get_wtime();
                    comp[j] += st2 - st1;
                }
            }
        }
        // update LES and LESQ，remove 0 in Dom
        for (uint32_t i = n - 1; i >= 0 && i < n; --i) {
            if (DomValues[i] == 0) {
                for (uint32_t j = 0; j < d; ++j) {
                    double st1 = omp_get_wtime();
                    // auto& innerMapS = LES[j];
                    for (auto& entry : LES[j]) {
                        auto& vectorS = entry.second;
                        vectorS.erase(vectorS.begin() + i);
                    }
                    // auto& innerMapQ = LESQ[j];
                    for (auto& entry : LESQ[j]) {
                        auto& vectorQ = entry.second;
                        vectorQ.erase(vectorQ.begin() + i);
                    }
                    double st2 = omp_get_wtime();
                    comp[j] += st2 - st1;
                }
            }
        }
        for (uint32_t j = 0; j < d; ++j) {
            double st1 = omp_get_wtime();
            for (uint32_t i = 0; i < n; ++i) {
                if (DomValues[i] == 0) {
                    // Update Ind
                    Ind[j][IndMap[j][indT[i]]].visit = 1;
                }
            }
            double st2 = omp_get_wtime();
            comp[j] += st2 - st1;
        }
        // double Flag1 = omp_get_wtime();
        for (uint32_t j = 0; j < d; ++j) {
            FlagF_P[j] = cc->MakePackedPlaintext(FlagF[j]);
            FlagF_C[j] = cc->Encrypt(keys.publicKey, FlagF_P[j]);
        }
        auto FlagF_ALL    = cc->EvalMultMany(FlagF_C);
        int64_t FlagFinal = Dec(cc->EvalSum(FlagF_ALL, Alln), Alln)->GetPackedValue()[0];
        // cout << "Flag:" << FlagFinal;
        // double Flag2 = omp_get_wtime();
        // cout << "\tTime_Flag:" << (Flag2 - Flag1) << "ms" << endl;
        if (FlagFinal >= 1)
            break;
        uint32_t posT = 0;
        for (uint32_t j = 0; j < d; ++j) {
            double st1      = omp_get_wtime();
            posF[j]         = 0;
            auto& Ind_j     = Ind[j];
            uint32_t size_j = Ind_j.size();
            for (uint32_t i = 0; i < size_j; ++i) {
                posF[j] += Ind_j[i].visit;
            }
            if (posF[j] == size_j)
                ++posT;
            double st2 = omp_get_wtime();
            comp[j] += st2 - st1;
        }
        if (posT == d || T.empty())
            flag = 0;
        // cout << "index" << endl;
        // for (uint32_t j = 0; j < d; j++) {
        //     cout << posF[j] << "\t";
        // }
        // cout << endl;
        // for (uint32_t j = 0; j < d; j++) {
        //     cout << posI[j] << endl;
        //     for (uint32_t i = 0; i < n; i++) {
        //         cout << Ind[j][IndMap[j][indT[i]]][2] << "\t";
        //     }
        //     cout << endl;
        //     for (uint32_t i = 0; i < n; i++) {
        //         cout << Anch[i][j] << "\t";
        //     }
        //     cout << endl;
        // }
    }
    double s04    = omp_get_wtime();
    query_process = s04 - s03;
    cout << "num:" << num << endl;
    uint32_t k = T.size();
    cout << "k:" << k << endl;
    for (uint32_t j = 0; j < d; ++j) {
        cout << "Time_" << j << ":" << fixed << setprecision(2) << (comp[j]) << " s" << endl;
        cout << "Cost_" << j << ":" << fixed << setprecision(2) << (comm[j] / 1024 / 1024) << " MB" << endl;
        if (Utime < comp[j])
            Utime = comp[j];
        if (Ucost < comm[j])
            Ucost = comm[j] / 1024 / 1024;
    }
    cout << "Time_S:" << fixed << setprecision(2) << (comp[d]) << " s" << endl;
    cout << "Cost_S:" << fixed << setprecision(2) << (comm[d] / 1024 / 1024) << " MB" << endl;
    SkyNum = k;
    ItrNum = num;
    ComNum = itr;
    Stime  = comp[d];
    Scost  = comm[d] / 1024 / 1024;
    // for (auto& entry : T) {
    //     for (auto& r : P[entry.first]) {
    //         cout << r << "\t";
    //     }
    //     cout << endl;
    // }
}

void MKsky12T(vector<vector<uint64_t>> P) {
    offline_process = 0;
    query_process   = 0;
    map<uint32_t, vector<uint64_t>> T;
    double s01 = omp_get_wtime();
    vector<vector<IndexElement>> Ind(d, vector<IndexElement>(n));
    vector<vector<uint32_t>> IndMap(d, vector<uint32_t>(n));
    vector<vector<int64_t>> FlagF(d, vector<int64_t>(n, 0));
    vector<Plaintext> FlagF_P(d);
    vector<Ciphertext<DCRTPoly>> FlagF_C(d);
    vector<uint32_t> Left(d, 0);
    vector<double> comm(d + 1), comp(d + 1);
    for (uint32_t j = 0; j < d; ++j) {
        vector<ElementIndex> list;
        for (uint32_t i = 0; i < n; ++i) {
            list.emplace_back(P[i][j], i);
        }
        sort(list.begin(), list.end());
        for (uint32_t i = 0; i < n; ++i) {
            Ind[j][i].element = list[i].element;
            Ind[j][i].index   = list[i].index;
            Ind[j][i].visit   = 0;
        }
        for (uint32_t i = 0; i < n; ++i) {
            IndMap[j][list[i].index] = i;
        }
    }
    for (uint32_t i = 0; i < n; ++i) {
        T[i] = P[i];
    }
    // Select
    uint32_t flag = 1;
    // store in column
    map<uint32_t, map<uint32_t, vector<int64_t>>> LES, LESQ;
    for (uint32_t j = 0; j < d; ++j) {
        map<uint32_t, vector<int64_t>> h1, h2;
        for (uint32_t k = 0; k < n; ++k) {
            vector<int64_t> ht1(n), ht2(n);
            for (uint32_t i = 0; i < n; ++i) {
                ht1[i] = (P[k][j] < P[i][j]) ? 0 : 1;
                ht2[i] = (P[k][j] <= P[i][j]) ? 0 : 1;
            }
            h1[k] = move(ht1);
            h2[k] = move(ht2);
        }
        LES[j]  = move(h1);
        LESQ[j] = move(h2);
    }
    unordered_map<uint32_t, unordered_map<uint32_t, Plaintext>> splitVec;
    for (uint32_t it = 1; it < log2(ringDim); ++it) {
        unordered_map<uint32_t, Plaintext> splitVecT;
        uint32_t len = (uint32_t)pow(2, it);
        for (uint32_t i = 0; i < len; ++i) {
            vector<int64_t> atIndex(2 * len, 0);
            atIndex[len + i] = 1;
            splitVecT[i]     = cc->MakePackedPlaintext(atIndex);
        }
        splitVec[it] = move(splitVecT);
    }
    double s02      = omp_get_wtime();
    offline_process = s02 - s01;
    // cout << "offline_process time: " << fixed << setprecision(2) << offline_process << " s" << endl;
    double s03 = omp_get_wtime();
    int itr = 0, num = 0;
    vector<uint32_t> posF(d, 0);
    vector<uint32_t> posI(d, 0);
    vector<Plaintext> LES_P(d);
    vector<Ciphertext<DCRTPoly>> IV_C(d);
    vector<vector<Ciphertext<DCRTPoly>>> DomRow(d, vector<Ciphertext<DCRTPoly>>(d));
    vector<Ciphertext<DCRTPoly>> DomT(d);
    usint Alln = n;
    NumItr.resize(n + 1, 0);
    while (flag == 1) {
        n = T.size();
        cout << "itr:" << itr << "\tn:" << n << endl;
        itr++;
        num += n;
        NumItr[n]++;
        if (NumItr[n] > ts || n == 1)
            break;
        // cout << "LES" << endl;
        // for (uint32_t j = 0; j < LES.size(); ++j) {
        //     map<uint32_t, vector<int64_t>> LESt = LES[j];
        //     for (auto& inner_pair : LESt) {
        //         auto& values = inner_pair.second;
        //         for (auto value : values) {
        //         cout << value << "\t";
        //         }
        //         cout << endl;
        //     }
        //     cout << endl;
        // }
        // cout << "LESQ" << endl;
        // for (uint32_t j = 0; j < d; ++j) {
        //     map<uint32_t, vector<int64_t>> LESQt = LESQ[j];
        //     for (auto& inner_pair : LESQt) {
        //         auto& values = inner_pair.second;
        //         for (auto value : values) {
        //         cout << value << "\t";
        //         }
        //         cout << endl;
        //     }
        //     cout << endl;
        // }
        double st1 = omp_get_wtime();
        // vector<vector<int64_t>> Anch(d, vector<int64_t>(n, 0));
        // map<uint32_t, vector<int64_t>> Anch;
        // map<uint32_t, Plaintext> LES_P;
        map<uint32_t, map<uint32_t, Plaintext>> LESQ_P;
        uint32_t batchI = static_cast<size_t>(ceil(log2(n)));
        batchSize = 1 << batchI;
        // cout << batchI << "," << batchSize << endl;
        unordered_map<uint32_t, Plaintext> &splitVecT = splitVec.at(batchI);
        // Ciphertext<DCRTPoly> CZero = cc->Encrypt(keys.publicKey, cc->MakePackedPlaintext(vector<int64_t>(n, 0)));
        double st2 = omp_get_wtime();
        for (uint32_t j = 0; j < d; ++j)
            comp[j] += st2 - st1;
        srand(time(0));
        for (uint32_t j = 0; j < d; ++j) {
            // compute the index "posI" about min value
            double st1 = omp_get_wtime();
            if (posF[j] == Ind[j].size()) {
                size_t randomIndex = rand() % T.size();
                auto it            = T.begin();
                next(it, randomIndex);
                posI[j] = it->first;
            }
            else {
                auto& row = Ind[j];
                for (uint32_t i = Left[j]; i < row.size(); ++i) {
                    if (row[i].visit == 0) {
                        posI[j]           = row[i].index;
                        row[i].visit      = 1;
                        Left[j]           = i;
                        FlagF[j][posI[j]] = 1;
                        break;
                    }
                }
            }
            // compute the anchor "Anch" in the origin IV
            uint32_t closestIndex = distance(T.begin(), T.find(posI[j]));
            // Anch[j] = move(vector<int64_t>(n, 0));
            vector<int64_t> Ancht(2 * batchSize, 0);
            // Anch[j][closestIndex] = 1;
            Ancht[batchSize + closestIndex] = 1;
            LES_P[j]                   = cc->MakePackedPlaintext(LES[j][posI[j]]);
            // compute the encrypted "Anch" in the origin IV
            IV_C[j] = cc->Encrypt(keys.publicKey, cc->MakePackedPlaintext(Ancht));
            // cout << j << "-IV:" << Serial::SerializeToString(IV_C[j]).size() << endl;
            // comm[j] += Serial::SerializeToString(IV_C[j]).size();
            // packing the column value about LESQ
            // map<uint32_t, vector<int64_t>> LESQt = LESQ[j];
            map<uint32_t, Plaintext> LESQ_Pt;
            uint32_t i = 0;
            for (auto& inner_pair : LESQ[j]) {
                LESQ_Pt[i] = cc->MakePackedPlaintext(inner_pair.second);
                i++;
            }
            LESQ_P[j]  = move(LESQ_Pt);
            double st2 = omp_get_wtime();
            comp[j] += st2 - st1;
            size_t tmp = Serial::SerializeToString(IV_C[j]).size();
            comm[j] += tmp;
            comm[d] += tmp;
        }
        // compute the inner product between matrix and IV
        for (uint32_t j = 0; j < d; ++j) {
            double st1 = omp_get_wtime();
            for (uint32_t k = 0; k < d; ++k) {
                if (j == k) {
                    DomRow[k][j] = cc->Encrypt(keys.publicKey, LES_P[j]);
                }
                else {
                    // const auto& innerMap = LESQ_P[j];
                    // Ciphertext<DCRTPoly> Domi = CZero->Clone();
                    vector<Ciphertext<DCRTPoly>> Domi(n);
                    vector<Plaintext> LESQt;
                    LESQt.reserve(n);
                    for (const auto& entry : LESQ_P[j]) {
                        LESQt.push_back(entry.second);
                    }
#pragma omp parallel for
                    for (uint32_t i = 0; i < n; ++i) {
                        Ciphertext<DCRTPoly> splitInd = cc->EvalMult(IV_C[k], splitVecT[i]);
                        Domi[i]                       = cc->EvalSum(splitInd, 2 * batchSize);
                        Domi[i]                       = cc->EvalMult(Domi[i], LESQt[i]);
                    }
                    // uint32_t i = 0;
                    // for (const auto& entry : innerMap) {
                    //     // // first try
                    //     // Ciphertext<DCRTPoly> splitInd = cc->EvalMult(IV_C[k], splitVec[i]);
                    //     // Ciphertext<DCRTPoly> IndexSum = cc->EvalSum(splitInd, batchSize);
                    //     // second try
                    //     // Ciphertext<DCRTPoly> splitInd   = cc->EvalInnerProduct(IV_C[k], splitVec[i], n);
                    //     // Ciphertext<DCRTPoly> Indsum     = cc->EvalMult(splitInd, splitVec[0]);
                    //     // Ciphertext<DCRTPoly> splitIndex = cc->EvalRotate(splitInd, indexList[0]);
                    //     // Ciphertext<DCRTPoly> IndexSum   = cc->EvalSum(splitIndex, Alln);
                    //     // Domi[i] = cc->EvalMult(IndexSum, entry.second);
                    //     Domi[i] = cc->EvalMult(Domi[i], entry.second);
                    //     i++;
                    // }
                    DomRow[k][j] = cc->EvalAddMany(Domi);
                    // DomRow[j * d + k] = Domi;
                }
                // comm[j] += Serial::SerializeToString(DomRow[j * d + k]).size();
            }
            double st2 = omp_get_wtime();
            comp[j] += st2 - st1;
            for (uint32_t k = 0; k < d; ++k) {
                comm[j] += Serial::SerializeToString(DomRow[k][j]).size();
            }
        }
        // cout << "Anch:" << endl;
        // for (uint32_t j = 0; j < d; ++j) {
        //     cout << posI[j] << endl;
        //     cout << Dec(IV_C[j], 2 * Alln) << endl;
        // }
        // cout << "DomRow:" << endl;
        // for (uint32_t j = 0; j < d; ++j) {
        //     for (uint32_t k = 0; k < d; ++k) {
        //         cout << Dec(DomRow[j * d + k], n) << endl;
        //     }
        //     cout << endl;
        // }
        // compute DomT
        double ss1 = omp_get_wtime();
        for (uint32_t k = 0; k < d; k++) {
            DomT[k] = cc->EvalAddMany(DomRow[k]);
        }
        // cout << "DomT:" << endl;
        // for (uint32_t k = 0; k < d; ++k) {
        //     cout << Dec(DomT[k], n) << endl;
        // }
        // compute Dom
        // Ciphertext<DCRTPoly> Dom_C = DomT[0]->Clone();
        // for (uint32_t j = 1; j < d; j++) {
        //     Dom_C = cc->EvalMult(Dom_C, DomT[j]);
        //     cout << Dec(Dom_C->Clone(), n) << endl;
        // }
        Ciphertext<DCRTPoly> Dom_C = cc->EvalMultMany(DomT);
        double ss2                 = omp_get_wtime();
        comp[d] += ss2 - ss1;
        double st3                = omp_get_wtime();
        vector<int64_t> DomValues = Dec(Dom_C, n)->GetPackedValue();
        // cout << "Dom:" << endl;
        // cout << DomValues << endl;
        // cout << endl;
        // index of T
        vector<uint32_t> indT(n);
        uint32_t i0 = 0;
        for (auto entry = T.begin(); entry != T.end(); ++entry) {
            indT[i0] = entry->first;
            i0++;
        }
        double st4 = omp_get_wtime();
        for (uint32_t j = 0; j < d; ++j)
            comp[j] += st4 - st3;
        for (uint32_t j = 0; j < d + 1; ++j) {
            comm[j] += Serial::SerializeToString(Dom_C).size();
        }
        for (uint32_t i = 0; i < n; ++i) {
            if (DomValues[i] == 0) {
                double ss1           = omp_get_wtime();
                uint32_t keyToRemove = indT[i];
                T.erase(keyToRemove);
                double ss2 = omp_get_wtime();
                comp[d] += ss2 - ss1;
                for (uint32_t j = 0; j < d; ++j) {
                    double st1      = omp_get_wtime();
                    auto& innerMapS = LES[j];
                    innerMapS.erase(keyToRemove);
                    auto& innerMapQ = LESQ[j];
                    innerMapQ.erase(keyToRemove);
                    double st2 = omp_get_wtime();
                    comp[j] += st2 - st1;
                }
            }
        }
        // update LES and LESQ，remove 0 in Dom
        for (uint32_t i = n - 1; i >= 0 && i < n; --i) {
            if (DomValues[i] == 0) {
                for (uint32_t j = 0; j < d; ++j) {
                    double st1 = omp_get_wtime();
                    // auto& innerMapS = LES[j];
                    for (auto& entry : LES[j]) {
                        auto& vectorS = entry.second;
                        vectorS.erase(vectorS.begin() + i);
                    }
                    // auto& innerMapQ = LESQ[j];
                    for (auto& entry : LESQ[j]) {
                        auto& vectorQ = entry.second;
                        vectorQ.erase(vectorQ.begin() + i);
                    }
                    double st2 = omp_get_wtime();
                    comp[j] += st2 - st1;
                }
            }
        }
        for (uint32_t j = 0; j < d; ++j) {
            double st1 = omp_get_wtime();
            for (uint32_t i = 0; i < n; ++i) {
                if (DomValues[i] == 0) {
                    // Update Ind
                    Ind[j][IndMap[j][indT[i]]].visit = 1;
                }
            }
            double st2 = omp_get_wtime();
            comp[j] += st2 - st1;
        }
        for (uint32_t j = 0; j < d; ++j) {
            double st1 = omp_get_wtime();
            FlagF_P[j] = cc->MakePackedPlaintext(FlagF[j]);
            FlagF_C[j] = cc->Encrypt(keys.publicKey, FlagF_P[j]);
            double st2 = omp_get_wtime();
            comp[j] += st2 - st1;
        }
        double Flag1      = omp_get_wtime();
        auto FlagF_ALL    = cc->EvalMultMany(FlagF_C);
        int64_t FlagFinal = Dec(cc->EvalSum(FlagF_ALL, Alln), Alln)->GetPackedValue()[0];
        // cout << "Flag:" << FlagFinal;
        double Flag2 = omp_get_wtime();
        for (uint32_t j = 0; j < d; ++j)
            comp[j] += Flag2 - Flag1;
        for (uint32_t j = 0; j < d + 1; ++j) {
            comm[j] += Serial::SerializeToString(FlagF_ALL).size();
        }
        // cout << "\tTime_Flag:" << (Flag2 - Flag1) << "ms" << endl;
        if (FlagFinal >= 1)
            break;
        uint32_t posT = 0;
        for (uint32_t j = 0; j < d; ++j) {
            double st1      = omp_get_wtime();
            posF[j]         = 0;
            auto& Ind_j     = Ind[j];
            uint32_t size_j = Ind_j.size();
            for (uint32_t i = 0; i < size_j; ++i) {
                posF[j] += Ind_j[i].visit;
            }
            if (posF[j] == size_j)
                ++posT;
            double st2 = omp_get_wtime();
            comp[j] += st2 - st1;
        }
        if (posT == d || T.empty())
            flag = 0;
        // cout << "index" << endl;
        // for (uint32_t j = 0; j < d; j++) {
        //     cout << posF[j] << "\t";
        // }
        // cout << endl;
        // for (uint32_t j = 0; j < d; j++) {
        //     cout << posI[j] << endl;
        //     for (uint32_t i = 0; i < n; i++) {
        //         cout << Ind[j][IndMap[j][indT[i]]][2] << "\t";
        //     }
        //     cout << endl;
        //     for (uint32_t i = 0; i < n; i++) {
        //         cout << Anch[i][j] << "\t";
        //     }
        //     cout << endl;
        // }
    }
    double s04    = omp_get_wtime();
    query_process = s04 - s03;
    cout << "num:" << num << endl;
    uint32_t k = T.size();
    cout << "k:" << k << endl;
    for (uint32_t j = 0; j < d; ++j) {
        cout << "Time_" << j << ":" << fixed << setprecision(2) << (comp[j]) << " s" << endl;
        cout << "Cost_" << j << ":" << fixed << setprecision(2) << (comm[j] / 1024 / 1024) << " MB" << endl;
        if (Utime < comp[j])
            Utime = comp[j];
        if (Ucost < comm[j])
            Ucost = comm[j] / 1024 / 1024;
    }
    cout << "Time_S:" << fixed << setprecision(2) << (comp[d]) << " s" << endl;
    cout << "Cost_S:" << fixed << setprecision(2) << (comm[d] / 1024 / 1024) << " MB" << endl;
    SkyNum = k;
    ItrNum = num;
    ComNum = itr;
    Stime  = comp[d];
    Scost  = comm[d] / 1024 / 1024;
    // for (auto& entry : T) {
    //     for (auto& r : P[entry.first]) {
    //         cout << r << "\t";
    //     }
    //     cout << endl;
    // }
}

void MKsky12T_2(vector<vector<uint64_t>> P) {
    offline_process = 0;
    query_process   = 0;
    map<uint32_t, vector<uint64_t>> T;
    double s01 = omp_get_wtime();
    vector<vector<IndexElement>> Ind(d, vector<IndexElement>(n));
    vector<vector<uint32_t>> IndMap(d, vector<uint32_t>(n));
    vector<vector<int64_t>> FlagF(d, vector<int64_t>(n, 0));
    vector<Plaintext> FlagF_P(d);
    vector<Ciphertext<DCRTPoly>> FlagF_C(d);
    vector<uint32_t> Left(d, 0);
    vector<double> comm(d + 1), comp(d + 1);
    for (uint32_t j = 0; j < d; ++j) {
        vector<ElementIndex> list;
        for (uint32_t i = 0; i < n; ++i) {
            list.emplace_back(P[i][j], i);
        }
        sort(list.begin(), list.end());
        for (uint32_t i = 0; i < n; ++i) {
            Ind[j][i].element = list[i].element;
            Ind[j][i].index   = list[i].index;
            Ind[j][i].visit   = 0;
        }
        for (uint32_t i = 0; i < n; ++i) {
            IndMap[j][list[i].index] = i;
        }
    }
    for (uint32_t i = 0; i < n; ++i) {
        T[i] = P[i];
    }
    // Select
    uint32_t flag = 1;
    // store in column
    map<uint32_t, map<uint32_t, vector<int64_t>>> LES, LESQ;
    for (uint32_t j = 0; j < d; ++j) {
        map<uint32_t, vector<int64_t>> h1, h2;
        for (uint32_t k = 0; k < n; ++k) {
            vector<int64_t> ht1(n), ht2(n);
            for (uint32_t i = 0; i < n; ++i) {
                ht1[i] = (P[k][j] < P[i][j]) ? 0 : 1;
                ht2[i] = (P[k][j] <= P[i][j]) ? 0 : 1;
            }
            h1[k] = move(ht1);
            h2[k] = move(ht2);
        }
        LES[j]  = move(h1);
        LESQ[j] = move(h2);
    }
    unordered_map<uint32_t, unordered_map<uint32_t, Plaintext>> splitVec;
    for (uint32_t it = 1; it < log2(ringDim); ++it) {
        unordered_map<uint32_t, Plaintext> splitVecT;
        uint32_t len = (uint32_t)pow(2, it);
        for (uint32_t i = 0; i < len; ++i) {
            vector<int64_t> atIndex(2 * len, 0);
            atIndex[len + i] = 1;
            splitVecT[i]     = cc->MakePackedPlaintext(atIndex);
        }
        splitVec[it] = move(splitVecT);
    }
    double s02      = omp_get_wtime();
    offline_process = s02 - s01;
    // cout << "offline_process time: " << fixed << setprecision(2) << offline_process << " s" << endl;
    double s03 = omp_get_wtime();
    int itr = 0, num = 0;
    vector<uint32_t> posF(d, 0);
    vector<uint32_t> posI(d, 0);
    vector<Plaintext> LES_P(d);
    vector<Ciphertext<DCRTPoly>> IV_C(d);
    vector<vector<Ciphertext<DCRTPoly>>> DomRow(d, vector<Ciphertext<DCRTPoly>>(d));
    vector<Ciphertext<DCRTPoly>> DomT(d);
    usint Alln = n;
    NumItr.resize(n + 1, 0);
    while (flag == 1) {
        n = T.size();
        cout << "itr:" << itr << "\tn:" << n << endl;
        itr++;
        num += n;
        NumItr[n]++;
        if (NumItr[n] > ts || n == 1)
            break;
        // cout << "LES" << endl;
        // for (uint32_t j = 0; j < LES.size(); ++j) {
        //     map<uint32_t, vector<int64_t>> LESt = LES[j];
        //     for (auto& inner_pair : LESt) {
        //         auto& values = inner_pair.second;
        //         for (auto value : values) {
        //         cout << value << "\t";
        //         }
        //         cout << endl;
        //     }
        //     cout << endl;
        // }
        // cout << "LESQ" << endl;
        // for (uint32_t j = 0; j < d; ++j) {
        //     map<uint32_t, vector<int64_t>> LESQt = LESQ[j];
        //     for (auto& inner_pair : LESQt) {
        //         auto& values = inner_pair.second;
        //         for (auto value : values) {
        //         cout << value << "\t";
        //         }
        //         cout << endl;
        //     }
        //     cout << endl;
        // }
        double st1 = omp_get_wtime();
        // vector<vector<int64_t>> Anch(d, vector<int64_t>(n, 0));
        // map<uint32_t, vector<int64_t>> Anch;
        // map<uint32_t, Plaintext> LES_P;
        map<uint32_t, vector<Plaintext>> LESQ_P;
        uint32_t batchI = static_cast<size_t>(ceil(log2(n)));
        batchSize = 1 << batchI;
        // cout << batchI << "," << batchSize << endl;
        unordered_map<uint32_t, Plaintext> &splitVecT = splitVec.at(batchI);
        // Ciphertext<DCRTPoly> CZero = cc->Encrypt(keys.publicKey, cc->MakePackedPlaintext(vector<int64_t>(n, 0)));
        double st2 = omp_get_wtime();
        for (uint32_t j = 0; j < d; ++j)
            comp[j] += st2 - st1;
        srand(time(0));
        for (uint32_t j = 0; j < d; ++j) {
            // compute the index "posI" about min value
            double st1 = omp_get_wtime();
            if (posF[j] == Ind[j].size()) {
                size_t randomIndex = rand() % T.size();
                auto it            = T.begin();
                next(it, randomIndex);
                posI[j] = it->first;
            }
            else {
                auto& row = Ind[j];
                for (uint32_t i = Left[j]; i < row.size(); ++i) {
                    if (row[i].visit == 0) {
                        posI[j]           = row[i].index;
                        row[i].visit      = 1;
                        Left[j]           = i;
                        FlagF[j][posI[j]] = 1;
                        break;
                    }
                }
            }
            // compute the anchor "Anch" in the origin IV
            uint32_t closestIndex = distance(T.begin(), T.find(posI[j]));
            // Anch[j] = move(vector<int64_t>(n, 0));
            vector<int64_t> Ancht(2 * batchSize, 0);
            // Anch[j][closestIndex] = 1;
            Ancht[batchSize + closestIndex] = 1;
            LES_P[j]                   = cc->MakePackedPlaintext(LES[j][posI[j]]);
            // compute the encrypted "Anch" in the origin IV
            IV_C[j] = cc->Encrypt(keys.publicKey, cc->MakePackedPlaintext(Ancht));
            // cout << j << "-IV:" << Serial::SerializeToString(IV_C[j]).size() << endl;
            // comm[j] += Serial::SerializeToString(IV_C[j]).size();
            // packing the column value about LESQ
            // map<uint32_t, vector<int64_t>> LESQt = LESQ[j];
            vector<Plaintext> LESQ_Pt;
            uint32_t i = 0;
            for (auto& inner_pair : LESQ[j]) {
                LESQ_Pt[i] = cc->MakePackedPlaintext(inner_pair.second);
                i++;
            }
            LESQ_P[j]  = move(LESQ_Pt);
            double st2 = omp_get_wtime();
            comp[j] += st2 - st1;
            size_t tmp = Serial::SerializeToString(IV_C[j]).size();
            comm[j] += tmp;
            comm[d] += tmp;
        }
        // compute the inner product between matrix and IV
        for (uint32_t j = 0; j < d; ++j) {
            double st1 = omp_get_wtime();
            if(j<2){
            for (uint32_t k = 0; k < d; ++k) {
                if (j == k) {
                    DomRow[k][j] = cc->Encrypt(keys.publicKey, LES_P[j]);
                }
                else {
                    // const auto& innerMap = LESQ_P[j];
                    // Ciphertext<DCRTPoly> Domi = CZero->Clone();
                    vector<Ciphertext<DCRTPoly>> Domi(n);
                    vector<Plaintext> LESQt = LESQ_P[j];
#pragma omp parallel for
                    for (uint32_t i = 0; i < n; ++i) {
                        Ciphertext<DCRTPoly> splitInd = cc->EvalMult(IV_C[k], splitVecT[i]);
                        Domi[i]                       = cc->EvalSum(splitInd, 2 * batchSize);
                        Domi[i]                       = cc->EvalMult(Domi[i], LESQt[i]);
                    }
                    // uint32_t i = 0;
                    // for (const auto& entry : innerMap) {
                    //     // // first try
                    //     // Ciphertext<DCRTPoly> splitInd = cc->EvalMult(IV_C[k], splitVec[i]);
                    //     // Ciphertext<DCRTPoly> IndexSum = cc->EvalSum(splitInd, batchSize);
                    //     // second try
                    //     // Ciphertext<DCRTPoly> splitInd   = cc->EvalInnerProduct(IV_C[k], splitVec[i], n);
                    //     // Ciphertext<DCRTPoly> Indsum     = cc->EvalMult(splitInd, splitVec[0]);
                    //     // Ciphertext<DCRTPoly> splitIndex = cc->EvalRotate(splitInd, indexList[0]);
                    //     // Ciphertext<DCRTPoly> IndexSum   = cc->EvalSum(splitIndex, Alln);
                    //     // Domi[i] = cc->EvalMult(IndexSum, entry.second);
                    //     Domi[i] = cc->EvalMult(Domi[i], entry.second);
                    //     i++;
                    // }
                    DomRow[k][j] = cc->EvalAddMany(Domi);
                    // DomRow[j * d + k] = Domi;
                }
                // comm[j] += Serial::SerializeToString(DomRow[j * d + k]).size();
            }
            }else {
                for (uint32_t k = 0; k < d; ++k) {
                    if (j == k) {
                        DomRow[k][j] = cc->Encrypt(keys.publicKey, LES_P[j]);
                    }else {
                        DomRow[k][j] = cc->Encrypt(keys.publicKey, cc->MakePackedPlaintext(LESQ[j][posI[k]]));
                    }
                }
            }
            double st2 = omp_get_wtime();
            comp[j] += st2 - st1;
            for (uint32_t k = 0; k < d; ++k) {
                comm[j] += Serial::SerializeToString(DomRow[k][j]).size();
            }
        }
        // cout << "Anch:" << endl;
        // for (uint32_t j = 0; j < d; ++j) {
        //     cout << posI[j] << endl;
        //     cout << Dec(IV_C[j], 2 * Alln) << endl;
        // }
        // cout << "DomRow:" << endl;
        // for (uint32_t j = 0; j < d; ++j) {
        //     for (uint32_t k = 0; k < d; ++k) {
        //         cout << Dec(DomRow[j * d + k], n) << endl;
        //     }
        //     cout << endl;
        // }
        // compute DomT
        double ss1 = omp_get_wtime();
        for (uint32_t k = 0; k < d; k++) {
            DomT[k] = cc->EvalAddMany(DomRow[k]);
        }
        // cout << "DomT:" << endl;
        // for (uint32_t k = 0; k < d; ++k) {
        //     cout << Dec(DomT[k], n) << endl;
        // }
        // compute Dom
        // Ciphertext<DCRTPoly> Dom_C = DomT[0]->Clone();
        // for (uint32_t j = 1; j < d; j++) {
        //     Dom_C = cc->EvalMult(Dom_C, DomT[j]);
        //     cout << Dec(Dom_C->Clone(), n) << endl;
        // }
        Ciphertext<DCRTPoly> Dom_C = cc->EvalMultMany(DomT);
        double ss2                 = omp_get_wtime();
        comp[d] += ss2 - ss1;
        double st3                = omp_get_wtime();
        vector<int64_t> DomValues = Dec(Dom_C, n)->GetPackedValue();
        // cout << "Dom:" << endl;
        // cout << DomValues << endl;
        // cout << endl;
        // index of T
        vector<uint32_t> indT(n);
        uint32_t i0 = 0;
        for (auto entry = T.begin(); entry != T.end(); ++entry) {
            indT[i0] = entry->first;
            i0++;
        }
        double st4 = omp_get_wtime();
        for (uint32_t j = 0; j < d; ++j)
            comp[j] += st4 - st3;
        for (uint32_t j = 0; j < d + 1; ++j) {
            comm[j] += Serial::SerializeToString(Dom_C).size();
        }
        for (uint32_t i = 0; i < n; ++i) {
            if (DomValues[i] == 0) {
                double ss1           = omp_get_wtime();
                uint32_t keyToRemove = indT[i];
                T.erase(keyToRemove);
                double ss2 = omp_get_wtime();
                comp[d] += ss2 - ss1;
                for (uint32_t j = 0; j < d; ++j) {
                    double st1      = omp_get_wtime();
                    auto& innerMapS = LES[j];
                    innerMapS.erase(keyToRemove);
                    auto& innerMapQ = LESQ[j];
                    innerMapQ.erase(keyToRemove);
                    double st2 = omp_get_wtime();
                    comp[j] += st2 - st1;
                }
            }
        }
        // update LES and LESQ，remove 0 in Dom
        for (uint32_t i = n - 1; i >= 0 && i < n; --i) {
            if (DomValues[i] == 0) {
                for (uint32_t j = 0; j < d; ++j) {
                    double st1 = omp_get_wtime();
                    // auto& innerMapS = LES[j];
                    for (auto& entry : LES[j]) {
                        auto& vectorS = entry.second;
                        vectorS.erase(vectorS.begin() + i);
                    }
                    // auto& innerMapQ = LESQ[j];
                    for (auto& entry : LESQ[j]) {
                        auto& vectorQ = entry.second;
                        vectorQ.erase(vectorQ.begin() + i);
                    }
                    double st2 = omp_get_wtime();
                    comp[j] += st2 - st1;
                }
            }
        }
        for (uint32_t j = 0; j < d; ++j) {
            double st1 = omp_get_wtime();
            for (uint32_t i = 0; i < n; ++i) {
                if (DomValues[i] == 0) {
                    // Update Ind
                    Ind[j][IndMap[j][indT[i]]].visit = 1;
                }
            }
            double st2 = omp_get_wtime();
            comp[j] += st2 - st1;
        }
        for (uint32_t j = 0; j < d; ++j) {
            double st1 = omp_get_wtime();
            FlagF_P[j] = cc->MakePackedPlaintext(FlagF[j]);
            FlagF_C[j] = cc->Encrypt(keys.publicKey, FlagF_P[j]);
            double st2 = omp_get_wtime();
            comp[j] += st2 - st1;
        }
        double Flag1      = omp_get_wtime();
        auto FlagF_ALL    = cc->EvalMultMany(FlagF_C);
        int64_t FlagFinal = Dec(cc->EvalSum(FlagF_ALL, Alln), Alln)->GetPackedValue()[0];
        // cout << "Flag:" << FlagFinal;
        double Flag2 = omp_get_wtime();
        for (uint32_t j = 0; j < d; ++j)
            comp[j] += Flag2 - Flag1;
        for (uint32_t j = 0; j < d + 1; ++j) {
            comm[j] += Serial::SerializeToString(FlagF_ALL).size();
        }
        // cout << "\tTime_Flag:" << (Flag2 - Flag1) << "ms" << endl;
        if (FlagFinal >= 1)
            break;
        uint32_t posT = 0;
        for (uint32_t j = 0; j < d; ++j) {
            double st1      = omp_get_wtime();
            posF[j]         = 0;
            auto& Ind_j     = Ind[j];
            uint32_t size_j = Ind_j.size();
            for (uint32_t i = 0; i < size_j; ++i) {
                posF[j] += Ind_j[i].visit;
            }
            if (posF[j] == size_j)
                ++posT;
            double st2 = omp_get_wtime();
            comp[j] += st2 - st1;
        }
        if (posT == d || T.empty())
            flag = 0;
        // cout << "index" << endl;
        // for (uint32_t j = 0; j < d; j++) {
        //     cout << posF[j] << "\t";
        // }
        // cout << endl;
        // for (uint32_t j = 0; j < d; j++) {
        //     cout << posI[j] << endl;
        //     for (uint32_t i = 0; i < n; i++) {
        //         cout << Ind[j][IndMap[j][indT[i]]][2] << "\t";
        //     }
        //     cout << endl;
        //     for (uint32_t i = 0; i < n; i++) {
        //         cout << Anch[i][j] << "\t";
        //     }
        //     cout << endl;
        // }
    }
    double s04    = omp_get_wtime();
    query_process = s04 - s03;
    cout << "num:" << num << endl;
    uint32_t k = T.size();
    cout << "k:" << k << endl;
    for (uint32_t j = 0; j < d; ++j) {
        cout << "Time_" << j << ":" << fixed << setprecision(2) << (comp[j]) << " s" << endl;
        cout << "Cost_" << j << ":" << fixed << setprecision(2) << (comm[j] / 1024 / 1024) << " MB" << endl;
        if (Utime < comp[j])
            Utime = comp[j];
        if (Ucost < comm[j])
            Ucost = comm[j] / 1024 / 1024;
    }
    cout << "Time_S:" << fixed << setprecision(2) << (comp[d]) << " s" << endl;
    cout << "Cost_S:" << fixed << setprecision(2) << (comm[d] / 1024 / 1024) << " MB" << endl;
    SkyNum = k;
    ItrNum = num;
    ComNum = itr;
    Stime  = comp[d];
    Scost  = comm[d] / 1024 / 1024;
    // for (auto& entry : T) {
    //     for (auto& r : P[entry.first]) {
    //         cout << r << "\t";
    //     }
    //     cout << endl;
    // }
}

int main(int argc, char* argv[]) {
    for (int i = 1; i < argc; ++i) {
        string param = argv[i];
        size_t pos   = param.find('=');
        if (pos != string::npos) {
            string key   = param.substr(0, pos);
            string value = param.substr(pos + 1);
            if (key == "n") {
                n = stoi(value);
            }
            else if (key == "d") {
                d = stoi(value);
            }
            else if (key == "t") {
                dt = stoi(value);
            }
        }
    }
    // n = 1000;
    // d = 8;
    // cout << "n:" << n << ", d:" << d << ", dataset:" << dt << endl;
    string filename = "./data/data" + to_string(d) + "/size=" + to_string(n);
    // string filename = "./data/data"+to_string(d)+"/size="+to_string(n)+"/small_anti_corr.txt";
    // string filename = "./data/data"+to_string(d)+"/size="+to_string(n)+"/small_uniformly.txt";
    if (dt == 0) {
        filename += "/small_correlate.txt";
    }
    else if (dt == 1) {
        filename += "/small_anti_corr.txt";
    }
    else if (dt == 2) {
        filename += "/small_uniformly.txt";
    }
    else if (dt == 3) {
        filename += "/nba.txt";
    }
    cout << filename << endl;
    loadP(filename);
    while (batchSize < (2 * n)) {
        batchSize <<= 1;
    }
    ts               = relaxData[d][n][dt%3];
    double TimeInit1 = omp_get_wtime();
    Init1();
    double TimeInit2 = omp_get_wtime();
    cout << "Init time: " << fixed << setprecision(2) << (TimeInit2 - TimeInit1) << " s" << endl;
    double query1 = omp_get_wtime();
    MKsky12T(P);
    double query2 = omp_get_wtime();
    cout << "Overall time: " << fixed << setprecision(2) << (query2 - query1) << " s" << endl;
    cout << "offline_process: " << fixed << setprecision(2) << offline_process << " s" << endl;
    cout << "query_process: " << fixed << setprecision(2) << query_process << " s" << endl;
    ofstream outfile;
    // outfile.open("./dataTest/CostNumTrain.txt", ios::app | ios::in);
    outfile.open("./dataTest/Time-12-1T.txt", ios::app | ios::in);
    outfile << filename << endl;
    outfile << "itr:" << ComNum << "\tn:" << SkyNum << "\tnum:" << ItrNum << endl;
    outfile << "offline: " << fixed << setprecision(2) << (TimeInit2 - TimeInit1 + offline_process) << " s\t";
    outfile << "Utime: " << fixed << setprecision(2) << Utime << " s\t";
    outfile << "Ucost: " << fixed << setprecision(2) << Ucost << " MB\t";
    outfile << "Stime: " << fixed << setprecision(2) << Stime << " s\t";
    outfile << "Scost: " << fixed << setprecision(2) << Scost << " MB" << endl;
    outfile.close();
    return 0;
}

int mainExtend(int argc, char* argv[]) {
    for (int i = 1; i < argc; ++i) {
        string param = argv[i];
        size_t pos   = param.find('=');
        if (pos != string::npos) {
            string key   = param.substr(0, pos);
            string value = param.substr(pos + 1);
            if (key == "n") {
                n = stoi(value);
            }
            else if (key == "d") {
                d = stoi(value);
            }
            else if (key == "t") {
                dt = stoi(value);
            }
            else if (key == "v") {
                dv = stoi(value);
            }
        }
    }
    uint32_t nn = n;
    // n = 1000;
    // d = 8;
    cout << "v:" << dv << ",n:" << n << ", d:" << d << ", dataset:" << dt << endl;
    string filename = "./dataTest/data" + to_string(dv) + "/data" + to_string(d) + "/size=" + to_string(n);
    // string filename = "./data/data"+to_string(d)+"/size="+to_string(n)+"/small_anti_corr.txt";
    // string filename = "./data/data"+to_string(d)+"/size="+to_string(n)+"/small_uniformly.txt";
    if (dt == 0) {
        filename += "/small_correlate.txt";
    }
    else if (dt == 1) {
        filename += "/small_anti_corr.txt";
    }
    else {
        filename += "/small_uniformly.txt";
    }
    cout << filename << endl;
    loadP(filename);
    while (batchSize < (2 * n)) {
        batchSize <<= 1;
    }
    Dsky1(P);
    ofstream outfile;
    int paral = dv / 4;
    // outfile.open("./dataTest/CostNum.txt", ios::app | ios::in);
    outfile.open("./dataTest/" + to_string(paral) + "CostNum-" + to_string(d) + "-" + to_string(nn) + ".txt",
                 ios::app | ios::in);
    if (dt == 0) {
        outfile << "CORR-";
    }
    else if (dt == 1) {
        outfile << "ANTI-";
    }
    else {
        outfile << "INDE-";
    }
    uint32_t maxnum = 0;
    uint32_t maxi   = 0;
    for (size_t i = 0; i < NumItr.size(); i++) {
        if (NumItr[i] > maxnum) {
            maxnum = NumItr[i];
            maxi   = i;
        }
    }
    outfile << "max:" << maxi << ":" << maxnum << "@";
    for (size_t i = 0; i < NumItr.size(); i++) {
        if (NumItr[i] > 1)
            outfile << i << ":" << NumItr[i] << ";";
    }
    outfile << endl;
    outfile.close();
    return 0;
}

int mainNum(int argc, char* argv[]) {
    for (int i = 1; i < argc; ++i) {
        string param = argv[i];
        size_t pos   = param.find('=');
        if (pos != string::npos) {
            string key   = param.substr(0, pos);
            string value = param.substr(pos + 1);
            if (key == "n") {
                n = stoi(value);
            }
            else if (key == "d") {
                d = stoi(value);
            }
            else if (key == "t") {
                dt = stoi(value);
            }
            else if (key == "v") {
                dv = stoi(value);
            }
        }
    }
    // n = 1000;
    // d = 8;
    cout << "v:" << dv << ",n:" << n << ", d:" << d << ", dataset:" << dt << endl;
    string filename = "./dataTest/data" + to_string(dv) + "/data" + to_string(d) + "/size=" + to_string(n);
    // string filename = "./data/data"+to_string(d)+"/size="+to_string(n)+"/small_anti_corr.txt";
    // string filename = "./data/data"+to_string(d)+"/size="+to_string(n)+"/small_uniformly.txt";
    if (dt == 0) {
        filename += "/small_correlate.txt";
    }
    else if (dt == 1) {
        filename += "/small_anti_corr.txt";
    }
    else {
        filename += "/small_uniformly.txt";
    }
    cout << filename << endl;
    loadP(filename);
    while (batchSize < (2 * n)) {
        batchSize <<= 1;
    }
    // for (size_t d = 0; d < input_dim.size(); d++) {
    //     uint32_t dim                          = input_dim[d];
    //     map<uint32_t, vector<int64_t>> relaxt = relaxData[dim];
    //     for (auto& inner_pair : relaxt) {
    //         auto& values = inner_pair.second;
    //         for (auto value : values) {
    //             cout << value << "\t";
    //         }
    //         cout << endl;
    //     }
    // }
    ts = relaxData[d][n][dt];
    Dsky(P);
    ofstream outfile;
    // outfile.open("./dataTest/CostNumTrain.txt", ios::app | ios::in);
    outfile.open("./dataTest/CostNumTest.txt", ios::app | ios::in);
    outfile << filename << endl;
    outfile << "itr:" << ComNum << "\tn:" << SkyNum << "\tnum:" << ItrNum << endl;
    outfile << "offline: " << fixed << setprecision(2) << offline_process << " s\t";
    outfile << "query: " << fixed << setprecision(2) << query_process << " s" << endl;
    outfile.close();
    return 0;
}
