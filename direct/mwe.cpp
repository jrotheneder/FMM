#include <functional> 
#include <iostream> 

struct Test {
    const std::function <void()>& testFunction;
    Test(const std::function <void()>& testFunction): testFunction(testFunction) {}
    void call() { testFunction(); }
};

int main() {

    Test test1([]()  { std::cout << "test1!\n";});
    Test test2([]() { std::cout << "test2!\n";});

    test1.call(); // "test2!" 
    test2.call(); // "test2!"
}

