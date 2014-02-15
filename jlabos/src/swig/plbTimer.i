namespace plb {
    template<typename T> class MultiNTensorField3D;
}
/*
 * core/plbTimer.i
 */

%{
#include "PALABOS_ROOT/src/core/plbTimer.h"
%}

namespace plb {

namespace global {

class PlbTimer {
public:
    PlbTimer();
    void start();
    void restart();
    double stop();
    void reset();
    double getTime() const;
};

PlbTimer& timer(std::string nameOfTimer);

}  // namespace global

}  // namespace plb
