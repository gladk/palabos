#ifndef KBCTEMPLATES_H
#define KBCTEMPLATES_H

#include "core/globalDefs.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include <cmath>

namespace plb {

template<typename T, template<typename U> class Descriptor>
struct kbcTemplates
{

};

}

#include "latticeBoltzmann/kbcTemplates2D.h"

#endif // KBCTEMPLATES_H
