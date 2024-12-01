#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

int main() {
    const vector<string> input_size = {"100", "500", "1000", "2000", "4000"};
    const vector<int> input_dim     = {3, 4, 5, 6, 7, 8};

    for (size_t d = 0; d < input_dim.size(); d++) {
        int dim = input_dim[d];
        for (size_t n = 0; n < input_size.size(); n++) {
            string len         = input_size[n];
            string outfilename = "CostNum-" + to_string(dim) + "-" + len + ".txt";
            ofstream fileCost(outfilename);
            for (size_t i = 0; i <= 20; i++) {
                string filename = to_string(i) + "CostNum-" + to_string(dim) + "-" + len + ".txt";
                cout << filename << endl;
                ifstream file1Cost(filename);
                string line;
                while (getline(file1Cost, line)) {
                    fileCost << line << endl;
                }
                file1Cost.close();
            }
            fileCost.close();
        }
    }

    return 0;
}