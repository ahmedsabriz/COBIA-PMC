// BZZMATH: Release 7.0

//	=========================< BzzFactorizedSVD.HPP >===============================
//	* Class BzzFactorizedSVD derived from BzzFactorized								*
// * Description:																					*
//	*					Metodi Numerici e Software in C++ (Capitoli 8, 9)				*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
// * Examples: BzzMath\examples\BzzMathBasic\LinearSystems\							*
// *				FactorizedSVD\FactorizedSVD.cpp								*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	08-1991	Date Written
//	03-1994	Added BzzBaseClass for BzzPrint and BzzMessage.
//	11-1994	Conversion to double precision done.
//	02-1995	Added the Bzz prefix to the names of the classes.
// 02-1995	Added GetSortedBzzVectorDAndBzzMatricesVU function.

////////////////// Release 4.0
//	09-2001	Added TransposeSolve.

//	============================================================================
//	* BzzFactorizedSVD permits linear system solution and solutions to related	*
// * problems when numRows >= numColumns													*
//	* The solution minimises the residuals norm											*
//	****************************************************************************
//	* Other functions specifics for BzzFactorizedSVD:										*
//	* void Solve(A,b,&x,numD);																	*
//	* void Solve(&A,b,&x,numD);																*
//	* float GetMaxD(void);																		*
//	* float GetMinD(void);																		*
//	* void GetBzzVectorD(BzzVector *d);												*
//	* void GetBzzVectorDAndBzzMatrixV(BzzVector *d,BzzMatrix *V);	*
//	* void GetBzzVectorDAndBzzMatricesVU(BzzVector *d,							*
// *		BzzMatrix *V,BzzMatrix *U);										*
//	* void GetSortedBzzVectorDAndBzzMatricesVU(BzzVector *d,					*
// *		BzzMatrix *V,BzzMatrix *U);										*
//	* void GetBzzMatrixN(BzzMatrix *N);												*
//	* void GetRange(BzzMatrix *G);													*
//	* int PseudoInverse(BzzMatrix *inv);											*
//	* int PseudoInverse(BzzMatrix *inv,int numD);								*
//	* float ConditionNumber(void);															*
//	============================================================================

#ifndef BZZ_FACTORIZEDSVD_DOUBLE_HPP
#define BZZ_FACTORIZEDSVD_DOUBLE_HPP

//	============================================================================
//	====================< prototype functions for SVD >=========================
//	============================================================================

//#include "factordd.hpp"
//#include "factqrld.hpp"

char SVDFactorization
(int required, int m, int n, double** a, double* d, double** v);

void SVDEigenvalues(int required, int m, int n, double** a,
	double* d, double* ds, double** v);

void SVDSolution
(double dmax, int m, int n, double** u, double* d, double** v,
	double* b, double* x);

void SVDTransposeSolution
(double dmax, int m, int n, double** u, double* d, double** v,
	double* b, double* x);

//	============================================================================
//	=====================< class BzzFactorizedSVD >=========================
//	============================================================================
class BzzVector;
class BzzMatrix;
class BzzFactorizedSVD : public BzzFactorized
{
	friend void Householder(int j, int k, double* aux, double* d);
	friend void HouseholderApplyLeft(int ri, int rf, int ci, int cf, double** a,
		double* aux);
	friend void HouseholderApplyRight(int ri, int rf, int ci, int cf, double** a,
		double* aux);
	friend void Givens(double* x1, double* x2, double* c, double* s, double* w);
	friend void GivensApply(double* y1, double* y2, double c, double s, double w);
	friend void Hessenberg(int m, int n, double** a);
	friend void HessenbergP(int m, int n, double** a, double** p);
	friend void Bidiagonal(int m, int n, double** a, double* d, double* ds);
	friend void BidiagonalR(int m, int n, double** a, double* d,
		double* ds, double** an);

protected:

	double* dsvd, * bsvd, ** vBzzMatrix;
	double dMax, dMin;
	void FurtherInit(void);
	virtual void SpecificInitialize(void);
	virtual void SpecificDeinitialize(void);
	BzzFactorizedSVD(char ch, int rows, int columns)
		: BzzFactorized(ch, rows, columns) {
		FurtherInit();
	}

	//	============================================================================
	//	*********************< Functions for linear algebra >***********************
	//	============================================================================
	friend void Solve(BzzFactorizedSVD& A, const BzzVector& b,
		BzzVector* x, int numD);
	friend void Solve(BzzFactorizedSVD* A, const BzzVector& b,
		BzzVector* x, int numD);
	friend void TransposeSolve(BzzFactorizedSVD& A,
		const BzzVector& b, BzzVector* x, int numD);
	friend void TransposeSolve(BzzFactorizedSVD* A,
		const BzzVector& b, BzzVector* x, int numD);

	virtual void Factorization(void);
	virtual void Solution(BzzVector* bx);
	virtual void Solution(const BzzVector& b, BzzVector* x);
	void Solution(const BzzVector& b, BzzVector* x, double dSmall);
	virtual void TransposeSolution(BzzVector* bx);
	virtual void TransposeSolution(const BzzVector& b, BzzVector* x);
	void TransposeSolution(const BzzVector& b,
		BzzVector* x, double dSmall);
	virtual void Solution(BzzMatrix* BX);
	virtual void Solution(const BzzMatrix& B, BzzMatrix* X);
	virtual void TransposeSolution(const BzzMatrix& B,
		BzzMatrix* X);
	virtual void TransposeSolution(BzzMatrix* BX);
	virtual double	Condition(void) { return 1.; };
	virtual double	DeterminantEvaluation(void);

public:
	//	============================================================================
	//	****************************< constructors >********************************
	//	============================================================================
		// default constructor
	BzzFactorizedSVD(void)
		: BzzFactorized() {
		FurtherInit();
	}

	// copy constructor
	BzzFactorizedSVD(const BzzFactorized& rval)
		: BzzFactorized(rval) {
		FurtherInit();
	}

	// constructor A(3,5);
	BzzFactorizedSVD(int rows, int columns)
		: BzzFactorized(rows, columns) {
		FurtherInit();
	}

	// constructor A(3,2,1.,2.,3.,4.,5.,6.);
	BzzFactorizedSVD(int rows, int columns, double a11, ...);

	// constructor A(5,3,w);
	BzzFactorizedSVD(int rows, int columns, double* initvalues)
		: BzzFactorized(rows, columns, initvalues) {
		FurtherInit();
	}

	// constructor from BzzMatrix
	BzzFactorizedSVD(const BzzMatrix& rval)
		: BzzFactorized(rval) {
		FurtherInit();
	}

	// make a submatrix with rows,columns
	BzzFactorizedSVD(int rows, int columns, const BzzMatrix& rval)
		: BzzFactorized(rows, columns, rval) {
		FurtherInit();
	}

	// as above, commencing from irow,jcol
	BzzFactorizedSVD(int rows, int columns,
		int irow, int jcol, const BzzMatrix& rval)
		: BzzFactorized(rows, columns, irow, jcol, rval) {
		FurtherInit();
	}

	//	============================================================================
	//	******************************< destructor >********************************
	//	============================================================================
	~BzzFactorizedSVD(void);

	//	============================================================================
	//	**************************< assignment operators >**************************
	//	============================================================================
	void operator = (const BzzFactorizedSVD& rval);
	BzzFactorizedSVD& operator = (const BzzMatrix& rval);

	//	============================================================================
	//	============================< Functions >===================================
	//	============================================================================
	double GetMaxD(void);
	double GetMinD(void);
	void GetBzzVectorD(BzzVector* d); // A = UDVT
	void GetBzzVectorDAndBzzMatrixV(BzzVector* d, BzzMatrix* V);
	void GetBzzVectorDAndBzzMatrixU(BzzVector* d, BzzMatrix* U);
	void GetBzzVectorDAndBzzMatricesVU(BzzVector* d,
		BzzMatrix* V, BzzMatrix* U);
	void GetSortedBzzVectorDAndBzzMatricesVU
	(BzzVector* d, BzzMatrix* V, BzzMatrix* U);
	void GetBzzMatrixN(BzzMatrix* N);
	void GetRange(BzzMatrix* R);
	int PseudoInverse(BzzMatrix* inv);
	int PseudoInverse(BzzMatrix* inv, int numD);
	double ConditionNumber(void);
};

#endif // BZZ_FACTORIZEDSVD_DOUBLE_HPP