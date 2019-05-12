#ifndef DEBUGGING_H
#define DEBUGGING_H

#include <iostream>
#include <string>
#include <sstream>
#include <iomanip> 
#include <limits> 
#include <vector> 
#include <fstream> 
#include <chrono> 

template<typename ChronoDuration>
double chrono_duration(ChronoDuration d) {
    return std::chrono::duration_cast<std::chrono::microseconds>(d).count()/1E6;  
}

std::string time_proc_stamp(size_t iteration, size_t procId) {

    std::stringstream s;
    s << "[" << std::setw(4) << std::setfill('0') << iteration << "][p" <<
        procId << "]: ";
    return s.str();
}

template <typename T>
void printIterable(T obj, std::string sep = ",") {
    std::cout << "[ ";
    for(auto& i : obj) {
        std::cout << i << " " << sep; 
    }
    std::cout << "]\n";
}

template <typename T>
void printVec(T& vec, size_t len, std::string sep = ",") {
    std::cout << "[ ";
    for(size_t i = 0; i < len; i++) {
        std::cout << vec[i] << " " << sep; 
    }
    std::cout << "]\n";
}

template <typename T>
void iterableToFile(T& t, std::string filename, std::string sep = "\n") {

    std::ofstream ofile; 
    ofile.open(filename);
    ofile << std::setprecision(std::numeric_limits<double>::digits10) << std::scientific;

    for(auto el : t) {
        ofile << el << sep;
    }

    ofile.close();
}


template <typename T>
void vecToFile(T& vec, size_t len, std::string filename, std::string sep = "\n") {

    std::ofstream ofile; 
    ofile.open(filename);
    ofile << std::setprecision(std::numeric_limits<double>::digits10) 
        << std::scientific;

    for(size_t i = 0; i < len; ++i) {
        ofile << vec[i] << sep;
    }

    ofile.close();
}

#endif 
