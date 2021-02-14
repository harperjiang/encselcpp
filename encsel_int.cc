//
// Created by harper on 2/13/21.
//

#include <chrono>
#include <iostream>
#include <fstream>
#include "feature.h"
#include <vector>
#include <sstream>
#include <string>

using namespace std;
using namespace std::chrono;
using namespace lqf::encsel;

int main(int argc, char **argv) {

    std::ifstream infile(argv[1]);

    std::vector<unique_ptr<Feature>> features;
    features.push_back(unique_ptr<Feature>(new Sparsity()));
    features.push_back(unique_ptr<Feature>(new Entropy()));
    features.push_back(unique_ptr<Feature>(new Length()));
    features.push_back(unique_ptr<Feature>(new Distinct()));
    features.push_back(unique_ptr<Feature>(new Sortness(50)));
    features.push_back(unique_ptr<Feature>(new Sortness(100)));
    features.push_back(unique_ptr<Feature>(new Sortness(200)));

    MemFeatureRecorder recorder;

    auto start = high_resolution_clock::now();

    std::string line;
    while (std::getline(infile, line)) {
        for (auto &f: features) {
            f->Add(line);
        }
    }

    for (auto &f: features) {
        f->Close(recorder);
    }

    // Get ending timepoint
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "Time taken: " << duration.count() << " microseconds" << endl;
}