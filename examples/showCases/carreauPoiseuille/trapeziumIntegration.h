#include <cmath>

#include "functions.h"

#ifndef TRAPEZIUM_INTEGRATION_H
#define TRAPEZIUM_INTEGRATION_H


/** This class is a Newton-Raphson root finding solver.
 * It find y such that F(y)=0. It takes the function type as template
 * parameter and F(y), F'(y), the tolerance and the maximum iterations
 * the user allows in order for the NR solver to find the root.
 */

template <typename T>
class TrapeziumIntegration
{
public :
	TrapeziumIntegration(Function<T> *function_, T y0_, int numberSteps_)
		: function(function_), y0(y0_), numberSteps(numberSteps_)
	{   }
	
	T operator()(T y) const 
	{	
		T dy = (y-y0)/(T)numberSteps;
		
		if (dy <= T() || y==y0)
		{
			return T();
		}
		T integral = ((*(function))(y0)+(*(function))(y))/(T)2;
// 		std::cout << y << ", " << integral << std::endl;
		for (int iS = 1; iS < numberSteps; ++iS)
		{
// 			std::cout << integral << std::endl;
			T ty = y0+(T)iS*dy;
			integral += (*(function))(ty);
		}
		integral *= dy;
        
// 		std::cout << y << ", " << integral << std::endl;
		return integral;
    }
	
private :
	Function<T> *function;
	T y0;
	int numberSteps;
};

#endif
