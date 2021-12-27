//TODO: Constraints
// BZZMATH: Release 7.0

//	==================< BzzNumericalDifferentiation.hpp >=======================
//	* Class BzzNumericalDifferentiation for numerical differentiation				*
// * Description:																					*
//	*					Metodi Numerici e Software in C++ (Capitolo 18)					*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
// * Examples: c:\BzzMath\examples\BzzAdvanced\										 	*
// *				NumericalDifferentiation\NumericalDifferentiation.cpp				*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	05-1993	Date Written in double precision.
//	02-1994	Added the Bzz prefix to the names of the classes.

//	============================================================================
//	****** Constructors for BzzNumericalDifferentiation:								*
//	* BzzNumericalDifferentiation D(Func);	// Func name of function				*
//	****************************************************************************
//	***** Access functions :																	*
//	* double der = D(a);																			*
//	* double error = D.ErrorEstimation();													*
//	****************************************************************************

#ifndef NUMERICAL_DIFFERENTIATION_HPP
#define NUMERICAL_DIFFERENTIATION_HPP

//	============================================================================
//	===================< class BzzNumericalDifferentiation >====================
//	============================================================================
class BzzNumericalDifferentiation
{
private:
	static const char* const BZZ_ERROR;
	static const double ETA;
	static const double ZERO_DER;
	static const double SAFETY;
	static const int MAX_ITERATIONS;

	double (*ptrFunc)(double x);
	double errorEstimation;
	char stop;

public:
	//	============================================================================
	//	******************************< constructors >******************************
	//	============================================================================
		// default
	BzzNumericalDifferentiation(void)
	{
		BzzError("%s%s", BzzNumericalDifferentiation::BZZ_ERROR,
			BZZ_ERR_IMPLEMENTATION);
	}
	BzzNumericalDifferentiation(double (*ptr)(double x));

	//	============================================================================
	//	*******************************< destructor >*******************************
	//	============================================================================
	~BzzNumericalDifferentiation(void) {};

	//	============================================================================
	//	*******************************< functions >********************************
	//	============================================================================
	double ErrorEstimation(void) { return errorEstimation; }
	double operator () (double a);
};

#endif // NUMERICAL_DIFFERENTIATION_HPP