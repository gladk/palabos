namespace plb {
    template<typename T> class MultiNTensorField3D;
}
/*
 * core/globalDefs.i
 */

%{
#include "PALABOS_ROOT/src/core/globalDefs.h"
%}

namespace plb {

/*typedef ptrdiff_t plint;*/
typedef long plint;
typedef size_t pluint;

namespace IndexOrdering {
    enum OrderingT {forward, backward, memorySaving};
}

namespace global {

class IOpolicyClass {
public:
    void setIndexOrderingForStreams(IndexOrdering::OrderingT streamOrdering_);
    IndexOrdering::OrderingT getIndexOrderingForStreams() const;
    void setEndianSwitchOnBase64out(bool doSwitch);
    bool getEndianSwitchOnBase64out();
    void setEndianSwitchOnBase64in(bool doSwitch);
    bool getEndianSwitchOnBase64in();
private:
    IOpolicyClass();
};

class Directories {
public:
    void setOutputDir(std::string outputDir);
    void setLogOutDir(std::string logOutDir_);
    void setImageOutDir(std::string imageOutDir_);
    void setVtkOutDir(std::string inputDir_);
    void setInputDir(std::string inputDir);

    std::string getLogOutDir() const;
    std::string getImageOutDir() const;
    std::string getVtkOutDir() const;
    std::string getInputDir() const;
private:
    Directories();
};

Directories& directories();
IOpolicyClass& IOpolicy();

}   // namespace global

}   // namespace plb
