#include "MyClass.hpp"
#include <iostream>

MyClass::MyClass(const std::string& name) : name(name) {}

void MyClass::greet() {
    std::cout << "Hello from " << name << std::endl;
}