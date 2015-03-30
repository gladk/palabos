#include <cmath>
#include <vector>

#ifndef FUNCTIONS_2_H
#define FUNCTIONS_2_H


template <typename T>
T factorial(int num)
{
	T result = (T)1;
	for (int iPop = 1; iPop <= num; ++iPop)
	{
		result = result *= (T)iPop;
	}
	return result;
}



/** This class computes the hypergeometric function
 * _pF_q(a;b;z)=sum_k=0^\infty z^k/k!*(Product_i=0^{p-1} (a[i])_k)/(Product_i=0^{q-1} (b[i])_k)
 * where (c)_k=c*(c+1)*...*(c+k-1).
 */

template <typename T>
class Hypergeom
{
public :
	Hypergeom(std::vector<T> a_, std::vector<T> b_, T tol_, int maxIter_)
		: a(a_), b(b_), tol(tol_), maxIter(maxIter_)
	{ }
	
	T operator()(T z) const 
	{
		T hypergeom = T();
		T hypergeomRef = 100;
		for (int iK = 0; iK < maxIter; ++iK)
		{
			T aPoch = (T)1;
			for (unsigned iPop = 0; iPop < a.size(); ++iPop)
			{
				aPoch *= pochHammer(a[iPop],iK);
			}
			
			T bPoch = (T)1;
			for (unsigned iPop = 0; iPop < b.size(); ++iPop)
			{
				bPoch *= pochHammer(b[iPop],iK);
			}
			
			hypergeom += std::pow(z,(T)iK) / factorial<T>(iK) * aPoch/bPoch;
			
// 			std::cout << "Factorial = " << factorial<T>(iK) << std::endl;
// 			std::cout << "aPoch = " << aPoch << std::endl;
// 			std::cout << "bPoch = " << bPoch << std::endl;
// 			std::cout << "power = " << std::pow(z,(T)iK) << std::endl;
// 			std::cout << "z = " << z << std::endl;
// 			
			
			T diff = T();
			if (hypergeomRef != T())
			{
				diff = std::fabs((hypergeom-hypergeomRef)/hypergeomRef);
			}
			else 
			{
				diff = std::fabs(hypergeom-hypergeomRef);
			}
// 			std::cout << hypergeom << " " << diff << std::endl;
			if ( diff < tol)
			{
// 				std::cout << "Hypergeom converged" << std::endl;
				return hypergeom;
			}
			hypergeomRef = hypergeom;
		}
        
		
		std::cout << "Error Hypergeom never converged." << std::endl;
		exit(1);
		
		return hypergeom;
    }
	
private :
	T pochHammer(T start, int iK) const
	{
		T poch = (T)1;
		for (int iPop = 0; iPop < iK; ++iPop)
		{
			poch *= (start + (T)iPop);
		}
		return poch;
	}
	
private :
	std::vector<T> a, b;
	
	T tol;
	int maxIter;
};

template <typename T>
class FunctionAndDerivative
{
public :
	virtual ~FunctionAndDerivative() { }
	
	virtual T operator()(T y, T z) const = 0;
	
	virtual T derivative(T y, T z) const = 0;
};

template <typename T>
struct Function
{
public :
	virtual ~Function() { }
	
	virtual T operator()(T y) const = 0;
};


/** This class computes the the function related with the derivative of the velocity
 * that one needs to compute u(y) for the case of the planar poiseille 
 * with a carreau-line viscosity.
 */

template <typename T>
class CarreauFunction : public FunctionAndDerivative<T>
{
public :
	CarreauFunction(T nu0_, T lambda_, T n_, T kappa_, T ly_, T tol_, int maxIter_)
		: nu0(nu0_), lambda(lambda_), n(n_), kappa(kappa_), ly(ly_), 
		  tol(tol_), maxIter(maxIter_)
	{  }
	
	virtual T operator()(T y, T z) const 
	{
		std::vector<T> a, b; 
		a.push_back((T)1/(T)2);
		a.push_back(((T)1-n)/(T)2);
		b.push_back((T)3/(T)2);
		Hypergeom<T> hyp(a,b,tol,maxIter);
		
		return nu0 * z * hyp(-lambda*lambda*z*z) + kappa*(ly/(T)2-y);
    }
	
	virtual T derivative(T y, T z) const
	{
		std::vector<T> a, b; 
		a.push_back((T)1/(T)2);
		a.push_back(((T)1-n)/(T)2);
		b.push_back((T)3/(T)2);
		
		std::vector<T> c, d; 
		a.push_back((T)3/(T)2);
		a.push_back(((T)3-n)/(T)2);
		b.push_back((T)5/(T)2);
		Hypergeom<T> hyp(a,b,tol,maxIter), hyp2(c,d,tol,maxIter);
		
		return nu0 * hyp(-lambda*lambda*z*z) 
				-((T)1 - n)/(T)3*lambda*lambda*z*z*nu0*hyp2(-lambda*lambda*z*z);
	}
		
private :
	T nu0, lambda, n, kappa, ly;
	
	T tol;
	int maxIter;
};

template <typename T>
class CarreauFunction_2 : public FunctionAndDerivative<T>
{
public :
	CarreauFunction_2(T nu0_, T lambda_, T n_, T kappa_, T ly_)
	: nu0(nu0_), lambda(lambda_), n(n_), kappa(kappa_), ly(ly_)
	  {  }
	
	  virtual T operator()(T y, T z) const 
	  {
		  return (ly-2*y)*kappa+2*std::pow((T)1+lambda*lambda*z*z,(n-(T)1)/(T)2)*nu0*z;
	  }

	  virtual T derivative(T y, T z) const
	  {
		  T lz2 = lambda*lambda*z*z;
		  
		  return (T)4*std::pow((T)1+lz2,(n-(T)1)/(T)2)*((n-(T)1)/(T)2)*lz2*nu0/((T)1+lz2)
				  + (T)2*std::pow((T)1+lz2,(n-(T)1)/(T)2)*nu0;
	  }
	
private :
	T nu0, lambda, n, kappa, ly;
};



#endif
