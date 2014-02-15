namespace plb {
    template<typename T> class MultiNTensorField3D;
}
/*
 * core/runTimeDiagnostics.i
 */

%{
#include "PALABOS_ROOT/src/core/runTimeDiagnostics.h"
%}

namespace plb {

class PlbException : public std::exception
{ };

class PlbMemoryException : public PlbException {
public:
    PlbMemoryException(std::string message_) throw();
    virtual ~PlbMemoryException() throw() { }
    virtual const char* what() const throw();
};

class PlbIOException : public PlbException {
public:
    PlbIOException(std::string message_) throw();
    virtual ~PlbIOException() throw() { }
    virtual const char* what() const throw();
};

class PlbNetworkException : public PlbException {
public:
    PlbNetworkException(std::string message_) throw();
    virtual ~PlbNetworkException() throw() { }
    virtual const char* what() const throw();
};

class PlbLogicException : public PlbException {
public:
    PlbLogicException(std::string message_) throw();
    virtual ~PlbLogicException() throw() { }
    virtual const char* what() const throw();
};

class PlbOutOfRangeException : public PlbException {
public:
    PlbOutOfRangeException(std::string message_) throw();
    virtual ~PlbOutOfRangeException() throw() { }
    virtual const char* what() const throw();
};

}  // namespace plb
