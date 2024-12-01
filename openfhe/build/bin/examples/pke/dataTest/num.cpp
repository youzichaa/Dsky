#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <limits>
#include <algorithm>
using namespace std;

struct Data {
    int y;
    vector<pair<int, int>> values;
};

int main() {
    const vector<string> input_size = {"100", "500", "1000", "2000", "4000"};
    const vector<int> input_dim     = {3, 4, 5, 6, 7, 8};
    // const vector<string> input_size = {"100"};
    // const vector<int> input_dim     = {3};

    int Testsize = 80;
    for (size_t d = 0; d < input_dim.size(); d++) {
        int dim = input_dim[d];
        for (size_t n = 0; n < input_size.size(); n++) {
            string len      = input_size[n];
            string filename = "CostNum-" + to_string(dim) + "-" + len + ".txt";
            // cout << filename << "\t";
            ifstream fileT(filename);
            vector<int> corr(81), anti(81), inde(81), corr1(81), anti1(81), inde1(81);
            string line;
            int lineCount = 0;
            int yt        = std::numeric_limits<int>::max();
            int bt        = std::numeric_limits<int>::min();
            std::vector<int> ys;
            while (getline(fileT, line)) {
                istringstream iss(line);
                string token;
                int y, yi, b = 0, b1 = 0;
                // Find max:x:y
                getline(iss, token, '@');
                yi = stoi(token.substr(token.find_first_of(":") + 1, token.find_last_of(":")));
                y  = stoi(token.substr(token.find_last_of(":")).substr(1));
                // Find @a:b
                getline(iss, token, '@');
                if (!token.empty()) {
                    token = token.substr(token.find_first_of(";") + 1);
                    if (!token.empty()) {
                        // bi = stoi(token.substr(0, token.find_first_of(":")));
                        b = stoi(token.substr(token.find_first_of(":") + 1, token.find_first_of(";")));
                    }
                    token = token.substr(token.find_first_of(";") + 1);
                    if (!token.empty()) {
                        // bi = stoi(token.substr(0, token.find_first_of(":")));
                        b1 = stoi(token.substr(token.find_first_of(":") + 1, token.find_first_of(";")));
                    }
                    // bi = stoi(token.substr(0,token.find_first_of(":")));
                    // b  = stoi(token.substr(token.find_first_of(":") + 1, subToken.find_last_of(";")));
                }
                if (lineCount % 3 == 0) {
                    corr[lineCount / 3]  = b;
                    corr1[lineCount / 3] = b1;
                }
                else if (lineCount % 3 == 1) {
                    anti[lineCount / 3]  = b;
                    anti1[lineCount / 3] = b1;
                }
                else if (lineCount % 3 == 2) {
                    inde[lineCount / 3]  = b;
                    inde1[lineCount / 3] = b1;
                }
                lineCount++;
            }
            // extra 1 to avoid filter ts=1
            corr[Testsize]  = 1;
            corr1[Testsize] = 1;
            anti[Testsize]  = 1;
            anti1[Testsize] = 1;
            inde[Testsize]  = 1;
            inde1[Testsize] = 1;
            // for (int number : anti) {
            //     cout << number << " ";
            // }
            // cout << endl;
            sort(corr.begin(), corr.end());
            sort(anti.begin(), anti.end());
            sort(inde.begin(), inde.end());
            cout << corr[corr.size()-1] << "," << anti[anti.size()-1] << ","
            << inde[inde.size()-1] << endl;

            // sort(corr1.begin(), corr1.end());
            // sort(anti1.begin(), anti1.end());
            // sort(inde1.begin(), inde1.end());
            // cout << corr1[corr1.size() - 1] << "," << anti1[anti1.size() - 1] << "," << inde1[inde1.size() - 1] << endl;
            fileT.close();
        }
    }
    return 0;
}