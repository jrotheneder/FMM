#ifndef DEBUGGING_H
#define DEBUGGING_H

#include <iostream>
#include <string>
#include <sstream>
#include <iomanip> 
#include <limits> 
#include <vector> 
#include <fstream> 

using namespace std;

string time_proc_stamp(size_t iteration, size_t procId) {

    stringstream s;
    s << "[" << setw(4) << setfill('0') << iteration << "][p" <<
        procId << "]: ";
    return s.str();
}

template <typename T>
void print_iterable(T obj) {
    cout << "[ ";
    for(auto& i : obj) {
        cout << i << ", "; 
    }
    cout << "]\n";
}

template <typename T>
void print_vec(T * vec, size_t len) {
    cout << "[ ";
    for(size_t i = 0; i < len; i++) {
        cout << vec[i] << ", "; 
    }
    cout << "]\n";
}

template <typename T>
void VectorToFile(std::vector<T> vec, string filename, string sep = "\n") {

    ofstream ofile; 
    ofile.open(filename);
    ofile << std::setprecision(std::numeric_limits<double>::digits10) << std::scientific;

    for(size_t i = 0; i < vec.size(); ++i) {
        ofile << vec[i] << sep;
    }

    ofile.close();
}


template <typename T>
void VectorToFile(T * vec, size_t len, string filename, string sep = "\n") {

    ofstream ofile; 
    ofile.open(filename);
    ofile << std::setprecision(std::numeric_limits<double>::digits10) << std::scientific;

    for(size_t i = 0; i < len; ++i) {
        ofile << vec[i] << sep;
    }

    ofile.close();
}

#endif 
