#include <iostream>
#include <string>
#include <sstream>
#include <iomanip> 
#include <vector> 

using namespace std;

string time_proc_stamp(size_t iteration, size_t procId) {

    stringstream s;
    s << "[" << setw(4) << setfill('0') << iteration << "][p" <<
        procId << "]: ";
    return s.str();
}

template <class T>
void print_vec(vector<T> vec) {
    cout << "[ ";
    for(auto i : vec) {
        cout << i << " "; 
    }
    cout << "]\n";
}

template <class T>
void print_vec(T * vec, size_t len) {
    cout << "[ ";
    for(size_t i = 0; i < len; i++) {
        cout << vec[i] << " "; 
    }
    cout << "]\n";
}
