/*
 * core/plbInit.i
 */

%{
#include "PALABOS_ROOT/src/core/plbInit.h"
%}


namespace plb {

void plbInit();

namespace global {
class MainArgv {
public:
    MainArgv(std::string argument_, int whichArg_);
    template<typename T> void read(T& variable);
    %template(readInt)    read<int>;
    %template(readDouble) read<double>;
    %template(readString) read<std::string>;
private:
    MainArgv();
};

class MainArgs {
public:
    int argc() const;
    MainArgv argv(int whichArg) const;
    void setArgs(int argcValue_, char*** argvPointer_);
private:
    MainArgs();
};


MainArgs& mainArguments();
int argc();
MainArgv argv(int whichArg);

}  // namespace global

}  // namespace plb
