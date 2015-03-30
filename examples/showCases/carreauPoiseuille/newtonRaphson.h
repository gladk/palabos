#include <cmath>

#include "functions.h"

#ifndef NEWTON_RAPHSON_H
#define NEWTON_RAPHSON_H


/** This class is a Newton-Raphson root finding solver.
 * It find y such that F(y)=0. It takes the function type as template
 * parameter and F(y), F'(y), the tolerance and the maximum iterations
 * the user allows in order for the NR solver to find the root.
 */

template <typename T>
class NewtonRaphson : public Function<T>
{
public :
	NewtonRaphson(FunctionAndDerivative<T> *function_, T tol_, int maxIter_)
		: function(function_), tol(tol_), maxIter(maxIter_)
	{   }
	
	virtual T operator()(T y) const 
	{
		T root0 = T();
		T root = T();
// 		std::cout << root0 << std::endl;
		for (int iPop = 0; iPop < maxIter; ++iPop)
		{
			root = root0 - (*(function))(y, root0) / function->derivative(y, root0);
// 			std::cout << std::fabs((root - root0)) << std::endl;
			if ((std::fabs((root - root0)/root0) < tol) || (root-root0) == T())
			{
// 				std::cout << root << std::endl;
				return root;
			}
			root0 = root;
		}
		
		std::cout << "Error Newton-Raphson never converged." << std::endl;
		exit(1);
        
		return root;
    }
	
private :
	FunctionAndDerivative<T> *function;
	T tol;
	int maxIter;
};

#endif
