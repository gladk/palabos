#include <iostream>
#include <cxxabi.h>

using namespace std;

const char* demangle(const char* name) {
    char buf[1024];
    size_t size=1024;
    int status;
    char* res = abi::__cxa_demangle(name,buf, &size, &status);
    return res;
}

int main(int argc, char* argv[0])
{
    if (argc != 2) {
        cout << "Syntax: " << argv[0] << " mangled_name" << endl;
    }
    cout << demangle(argv[1]) << endl;
}
