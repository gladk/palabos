namespace plb {
    template<typename T> class MultiNTensorField3D;
}
/*
 * core/blockStatistics.i
 */

%{
#include "PALABOS_ROOT/src/core/blockStatistics.h"
%}

namespace plb {

class StatSubscriber {
public:
    virtual ~StatSubscriber() { }
    virtual plint subscribeAverage() =0;
    virtual plint subscribeSum() =0;
    virtual plint subscribeMax() =0;
    virtual plint subscribeIntSum() =0;
};

class BlockStatistics {
public:
    BlockStatistics();
    void gatherAverage(plint whichAverage, double value);
    void gatherSum(plint whichSum, double value);
    void gatherMax(plint whichMax, double value);
    void gatherIntSum(plint whichSum, plint value);
    void incrementStats();
    pluint const& getNumCells() const;

    double getAverage(plint whichAverage) const;
    double getSum(plint whichSum) const;
    double getMax(plint whichMax) const;
    plint getIntSum(plint whichSum) const;

    plint subscribeAverage();
    plint subscribeSum();
    plint subscribeMax();
    plint subscribeIntSum();
};

}  // namespace pbl
