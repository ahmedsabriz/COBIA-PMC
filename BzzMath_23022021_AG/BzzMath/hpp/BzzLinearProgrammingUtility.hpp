// BZZMATH: Release 7.0

//	=====================< BzzLinearProgrammingUtility.hpp >================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	05-2013	Date Written.

#ifndef BZZ_LINEAR_PROGRAMMING_UTILITY
#define BZZ_MATRIX_DOUBLE_HPP

//	============================================================================
//	=================< class BzzLinearProgrammingUtility >======================
//	============================================================================
class BzzLinearProgrammingUtility
{
private:
	int numVariables;
	int numEquality;
	int numInequality;
	int numActiveInequality, numActiveBound;
	BzzVector xLower, xUpper, x, xx, dxLU, d, a, aa, aaa, columnNull;
	BzzVectorInt activeBound, iActiveBound,
		activeInequality, iActiveInequality;
	double eLower, eUpper;
	double dLower, dUpper;
	BzzMatrix A;
	BzzFactorizedLQ LQ;

public:
	//	============================================================================
	//	**************************< constructors >**********************************
	//	============================================================================
		// default
	BzzLinearProgrammingUtility(void) {}

	//	============================================================================
	//	*****************************< destructor >*********************************
	//	============================================================================
	~BzzLinearProgrammingUtility(void) {}

	//	============================================================================
	//	=========================< Non-modifying functions >========================
	//	============================================================================

	//	*****************************< BzzPrint >***********************************
	void FindBaricenter(BzzVector& xL, BzzVector& xU,
		BzzMatrix& E, BzzVector& e,
		BzzVector& dL, BzzMatrix& D, BzzVector& dU,
		BzzVector* xB);
};

#endif // BZZ_LINEAR_PROGRAMMING_UTILITY