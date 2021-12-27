// BZZMATH: Release 7.0

// ===========================< BzzMonoAttic.hpp >================================
// * Class BzzMonoAttic for attic linear programming										*
// ============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	06-2007	Date Written.
//	10-2009	Revision.
//	05-2013	Revision.

// ****************************************************************************

#ifndef BZZ_MONO_ATTIC_HPP
#define BZZ_MONO_ATTIC_HPP

// ============================================================================
// =========================< class BzzMonoAttic >=============================
// ============================================================================

class BzzMonoAttic : public BzzBaseClass
{
private:
	static const char* const BZZ_ERROR;
	static const double BOUND_ABSOLUTE_ERROR;
	static const double BOUND_RELATIVE_ERROR;
	static const double T_ERROR;
	static const double LAMBDA_ABSOLUTE_ERROR;
	static const double BZZ_BIG_ATTIC;
	static int count; // for whoAmI
	static int countInScope;
	int whoAmI;
	int	numVariables,
		unfeasible,
		solutionFounded, // 0 not yet 1 OK
//			restart,
vertex,
workingVariable,
numVertex,
numPivot;

	double	boundAbsoluteError,
		boundRelativeError,
		tError,
		lambdaAbsoluteError,
		bigAttic;

	BzzVector	x,
		lambda,
		lowerBounds,
		upperBounds,
		eq, // coefficients of equation
		s, // coefficients of function
		dxUL,
		dxL,
		dxU,
		d,
		sorted;

	BzzVectorInt	iSorted,
		activeBound,
		satisfiedBound;

public:
	void SetBoundAbsoluteError(double ba)
	{
		boundAbsoluteError = ba;
	}
	void SetBoundRelativeError(double br)
	{
		boundRelativeError = br;
	}
	void SetLambdaAbsoluteError(double la)
	{
		lambdaAbsoluteError = la;
	}
	void SetBigAttic(double big)
	{
		bigAttic = fabs(big);
	}
	double	f,
		fSolution,
		bBest,
		fBest,
		bEq,
		bLower,
		bUpper,
		bFloor,
		bRoof,
		fLower,
		fUpper,
		fFloor,
		fRoof;

	BzzVector	xSolution,
		xBest,
		xLower,
		xUpper,
		xFloor,
		xRoof;

	int WhoAmI(void) { return whoAmI; }
	BzzMonoAttic(void);
	void operator()(BzzVector* cc, BzzVector* ss,
		BzzVector* xxL, BzzVector* xxU);
	int operator()(double b);
	virtual void ObjectBzzPrint(void);
};

#endif // BZZ_MONO_ATTIC_HPP