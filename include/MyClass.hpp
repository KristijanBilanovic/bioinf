#ifndef MYCLASS_HPP
#define MYCLASS_HPP

#include <string>

class MyClass {
private:
    std::string name;
    
public:
    MyClass(const std::string& name);
    void greet();
};

#endif