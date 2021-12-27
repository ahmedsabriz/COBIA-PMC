// TODO: NTProduct, NProduct.
// TODO: NTSNProduct.

// BZZMATH: Release 7.0

//	=============================< FACTQRLD.HPP >===============================
//	* Class BzzFactorizedQRLQ derived from BzzFactorized					*
//	* Class BzzFactorizedQR derived from BzzFactorizedQRLQ					*
//	* Class BzzFactorizedLQ derived from BzzFactorizedQRLQ					*
// * Description:																					*
// *					Dal Fortan al C++	(Capitolo 15)										*
// *					by G. Buzzi-Ferraris														*
// *					Addison Wesley(1991)														*
// *							and																	*
// *					Scientific C++	(Chapter 15)											*
// *					by G. Buzzi-Ferraris														*
// *					Addison Wesley(1991)														*
// * Examples: BzzMath\Examples\BzzMathBasic\LinearSystems\FactorizedQR\	*
// *				FactorizedQR.cpp															*
// * Examples: BzzMath\Examples\BzzMathBasic\LinearSystems\FactorizedLQ\	*
// *				FactorizedLQ.cpp															*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	06-1991	Date Written
//	11-1992	English version
//	03-1994	Added BzzBaseClass for BzzPrint and BzzMessage.
// 04-1994	Added functions TransposeSolve, GetBzzMatrixN, GetBzzVectorD
//	11-1994	Conversion to double precision done.
//	02-1995	Added the Bzz prefix to the names of the classes.
// 02-1995	Modified GetBzzMatrixQ, P, L, N in BzzFactorizedLQ.
// 08-1995	Added SolveSquare in BzzFactorizedLQ.
// 07-1996	Added LinearLeastSquaresWithLinearEqualityConstraints.
// 07-1996	Added NonNegativeLinearLeastSquares.
// 09-1996	Added BoundedLinearLeastSquares.
// 09-1996	Added LinearWeightedLeastSquares.
// 09-1996	Added BoundedLinearWeightedLeastSquares.
// 10-1996	Added LeastDistanceProgramming.
// 03-1999	Added GetLinearCombinations and PrintLinearCombinations.
// 03-1999	Added GetTheNumberOfIndependentEquations.

////////////////// Release 4.0
//	07-2000	Added SetPrecision.
//	08-2000	Added AnalyzeDependencies.
//	08-2000	Added StartLinProMat,SolveLinProMat.
// 10-2001	Added function GetBzzVectorD in QR.
// 10-2001	Added GetLinearCombinations and PrintLinearCombinations in QR.
// 10-2001	Added new version of GetLinearCombinations in LQ e QR.
// 10-2001	Added GetTheNumberOfIndependentColumns in QR.
// 10-2001	Added GetNormColumns QR.
// 10-2001	Added GetNormRows LQ.
// 03-2002	Added NTProduct, NProduct.
////////////////// Release 5.0
// 01-2007	Added IntrinsicConditionNumber function.
// 01-2007	Added BzzLinearReconciliation function.
// 01-2007	Added BzzLinearReconciliationWithConstraints function.

//	============================================================================
//	* BzzFactorizedQR permits linear system solution and solutions to		*
// * related problems when numRows >= numColumns										*
//	* The solution minimises the residuals norm											*
//	****************************************************************************
//	* Other functions specifics for BzzFactorizedQR:								*
//	* void GetBzzMatrixQ(BzzMatrix *Q);												*
//	*	 calculates the matrix Q of the factorisation A = QR							*
//	* void GetBzzMatrixP(BzzMatrix *P);												*
//	*	 calculates the matrix P of the factorisation PA = R							*
//	* void GetBzzMatrixR(BzzMatrixRight *R);										*
//	*	 calculates the matrix R of the factorisation A = QR							*
//	* void GetResiduals(BzzVector *residuals);									*
//	*	 provides the residuals (only if numRows > numColumns)						*
//	* void GetBzzVectorD(BzzVector *d);														*
//	*	 calculates the diagonal vector														*
//	============================================================================
//	* BzzFactorizedLQ permits linear system solution and solutions to		*
//	* related problems when numRows <= numColumns										*
//	* The solution minimises the solution norm											*
//	****************************************************************************
//	* Other specific functions for BzzFactorizedLQ:								*
//	* void GetBzzMatrixQ(BzzMatrix *Q);												*
//	*	 calculates the matrix Q of the factorisation A = LQ							*
//	* void GetBzzMatrixP(BzzMatrix *P);												*
//	*	 calculates the matrix P of the factorisation AP = L							*
//	* void GetBzzMatrixL(BzzMatrixLeft *L);										*
//	*	 calculates the matrix L of the factorisation A = LQ							*
//	* void GetBzzMatrixN(BzzMatrix *N);												*
//	*	 calculates the matrix N null space of A											*
//	* void GetBzzVectorD(BzzVector *d);												*
//	*	 calculates the diagonal vector														*
//	****************************************************************************
//	* Other special functions:																	*
//	* **** SolveSquare																			*
//	* void SolveSquare(const BzzMatrix &A,const BzzVector &b,		*
//	* 				Bzzmatrix *X);																	*
//	*	 Given a matrix A(m,n) (m >= n) and vector b(n) solve							*
//	*	 m - n + 1 systems with the first n - 1 rows are constant				 	*
//	* **** LinearLeastSquaresWithLinearEqualityConstraints							*
//	* BzzVector LinearLeastSquaresWithLinearEqualityConstraints				*
//	*		(const BzzMatrix &A,const BzzVector &b,						*
//	*		const BzzMatrix &C,const BzzVector &d);						*
//	*	 Minimize the norm Ax - b and solve the linear system Cx = d				*
//	* **** NonNegativeLinearLeastSquares													*
//	* BzzVector NonNegativeLinearLeastSquares										*
//	* 		(BzzMatrix &A,BzzVector &b);										*
//	* BzzVector NonNegativeLinearLeastSquares										*
//	* 		(BzzMatrix *A,BzzVector *b);										*
//	*	 Minimize the norm Ax - b subject to x >= 0.										*
//	* **** BoundedLinearLeastSquares															*
//	* BzzVector BoundedLinearLeastSquares(BzzMatrix &A,				*
//	*		BzzVector &b,BzzVector &xMin,BzzVector &xMax);		*
//	*	 Minimize the norm Ax - b subject to xMin <= x <= xMax						*
//	* **** LinearWeightedLeastSquares														*
//	* BzzVector LinearWeightedLeastSquares(BzzMatrix &A,				*
//	*		BzzVector &b,BzzMatrixDiagonal W);								*
//	*	 Minimize the function (Ax - b)WTW(Ax - b)										*
//	* **** BoundedLinearWeightedLeastSquares												*
//	* BzzVector BoundedLinearWeightedLeastSquares(BzzMatrix &A,				*
//	*		BzzVector &b,BzzMatrixDiagonal W,								*
//	*		BzzVector &xMin,BzzVector &xMax);								*
//	*	 Minimize the function (Ax - b)WTW(Ax - b) subject to xMin <= x <= xMax	*
//	* **** LeastDistanceProgramming															*
//	* BzzVector LeastDistanceProgramming(BzzMatrix &A,BzzVector &b)*
//	*	 Minimize the x norm subject to Ax >= b											*
//	* **** LeastDistanceProgramming															*
//	* BzzVector LeastDistanceProgramming(BzzVector &xMin,						*
//	*		BzzVector &xMax,BzzMatrix &A,BzzVector &b)			*
//	*	 Minimize the x norm subject to xMin <= x <= xMax and Ax >= b				*
//	****************************************************************************

#ifndef BZZ_FACTORIZEDQRLQ_DOUBLE_HPP
#define BZZ_FACTORIZEDQRLQ_DOUBLE_HPP

//	============================================================================
//	======================< prototype functions for QR and LQ >=================
//	============================================================================
char QRFactorization
(int m, int n, double** a, double* d, double sT, double* norm);
void QRSolution
(int m, int n, double** a, double* d, double* b, double* x);
void QRTransposeSolution
(int m, int n, double** a, double* d, double* b, double* x);
double QRCondition
(int n, double** a, double* d);

char LQFactorization
(int m, int n, double** a, double* d, double sT, double* norm, double* intr);
void LQSolution
(int m, int n, double** a, double* d, double* b, double* x);
void LQTransposeSolution
(int m, int n, double** a, double* d, double* b, double* x);
double LQCondition
(int n, double** a, double* d);
void SolveSquare(BzzMatrix& A, BzzVector& b, BzzMatrix* X);
BzzVector LinearLeastSquaresWithLinearEqualityConstraints
(BzzMatrix& A, BzzVector& b,
	BzzMatrix& C, BzzVector& d);
BzzVector LinearLeastSquaresWithLinearEqualityConstraints
(BzzMatrix* A, BzzVector* b,
	BzzMatrix* C, BzzVector* d);
BzzVector NonNegativeLinearLeastSquares
(BzzMatrix& A, BzzVector& b);
BzzVector NonNegativeLinearLeastSquares
(BzzMatrix* A, BzzVector* b);
BzzVector BoundedLinearLeastSquares(BzzMatrix& A,
	BzzVector& b, BzzVector& xMin, BzzVector& xMax);
BzzVector LinearWeightedLeastSquares(BzzMatrix& A,
	BzzVector& b, BzzMatrixDiagonal W);
BzzVector BoundedLinearWeightedLeastSquares(BzzMatrix& A,
	BzzVector& b, BzzMatrixDiagonal W, BzzVector& xMin,
	BzzVector& xMax);
BzzVector LeastDistanceProgramming(BzzMatrix& A,
	BzzVector& b);
BzzVector LeastDistanceProgramming
(BzzVector& xMin, BzzVector& xMax, BzzMatrix& A,
	BzzVector& b);
BzzVector BzzLinearReconciliation(BzzMatrix* A, BzzVector* b,
	BzzVectorInt* iv, BzzVector* y);
BzzVector BzzLinearReconciliationWithConstraints(BzzMatrix* A, BzzVector* b,
	BzzVectorInt* iv, BzzVector* y, BzzVector* xMin, BzzVector* xMax);
BzzVector BzzWeightedLinearReconciliationWithConstraints(BzzMatrix* EE, BzzVector* ee,
	BzzVectorInt* iv, BzzVector* y, BzzVector* w, BzzVector* xxMin, BzzVector* xxMax);

//	============================================================================
//	=====================< class BzzFactorizedQRLQ >========================
//	============================================================================
class BzzVector;
class BzzMatrix;
class BzzFactorizedQRLQ : public BzzFactorized
{
	friend void Householder
	(int j, int k, double* aux, double* d);
	friend void HouseholderApplyLeft
	(int ri, int rf, int ci, int cf, double** a, double* aux);
	friend void HouseholderApplyRight
	(int ri, int rf, int ci, int cf, double** a, double* aux);
	friend void Givens
	(double* x1, double* x2, double* c, double* s, double* w);
	friend void GivensApply
	(double* y1, double* y2, double c, double s, double w);
	friend void Hessenberg
	(int m, int n, double** a);
	friend void HessenbergP
	(int m, int n, double** a, double** p);
	friend void Bidiagonal
	(int m, int n, double** a, double* d, double* ds);
	friend void BidiagonalR
	(int m, int n, double** a, double* d,
		double* ds, double** an);

protected:
	static const double BZZ_SIN_TINY;
	int* indx;
	double* dqr;
	double* bqr;
	double sinTiny;

	// initialisation of constructors
	void FurtherInit(void);

	//	initialisation of special vectors
	virtual void SpecificInitialize(void);

	//	deinitialisation of special vectors
	virtual void SpecificDeinitialize(void);

	// constructor A('*',3,5);
	BzzFactorizedQRLQ(char ch, int rows, int columns)
		: BzzFactorized(ch, rows, columns) {
		FurtherInit();
	}
	int matrixForLinProMat, startingNumRows;

	//	============================================================================
	//	*********************< Functions for linear algebra >***********************
	//	============================================================================
	virtual void	 Factorization(void) = 0;
	virtual void	 Solution(BzzVector* bx) = 0;
	virtual void	 Solution(const BzzVector& b, BzzVector* x) = 0;
	virtual void	 TransposeSolution(BzzVector* bx) = 0;
	virtual void	 TransposeSolution
	(const BzzVector& b, BzzVector* x) = 0;
	virtual void	 Solution(BzzMatrix* BX) = 0;
	virtual void	 Solution(const BzzMatrix& B, BzzMatrix* X) = 0;
	virtual void	 TransposeSolution(BzzMatrix* BX) = 0;
	virtual void	 TransposeSolution(const BzzMatrix& B,
		BzzMatrix* X) = 0;

	virtual double	Condition(void) = 0;
	virtual double	DeterminantEvaluation(void);

public:
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
		// default constructor
	BzzFactorizedQRLQ(void)
		: BzzFactorized() {
		FurtherInit();
	}

	// copy constructor
	BzzFactorizedQRLQ(const BzzFactorized& rval)
		: BzzFactorized(rval) {
		FurtherInit();
	}

	// constructor A(3,5);
	BzzFactorizedQRLQ(int rows, int columns)
		: BzzFactorized(rows, columns) {
		FurtherInit();
	}

	// constructor A(3,5,w);
	BzzFactorizedQRLQ(int rows, int columns, double* initvalues)
		: BzzFactorized(rows, columns, initvalues) {
		FurtherInit();
	}

	// constructor from BzzMatrix
	BzzFactorizedQRLQ(const BzzMatrix& rval)
		: BzzFactorized(rval) {
		FurtherInit();
	}

	// make a submatrix with rows,columns
	BzzFactorizedQRLQ(int rows, int columns, const BzzMatrix& rval)
		: BzzFactorized(rows, columns, rval) {
		FurtherInit();
	}

	// as above, commencing from irow,jcol
	BzzFactorizedQRLQ(int rows, int columns,
		int irow, int jcol, const BzzMatrix& rval)
		: BzzFactorized(rows, columns, irow, jcol, rval) {
		FurtherInit();
	}

	//	============================================================================
	//	*****************************< destructor >*********************************
	//	============================================================================
	~BzzFactorizedQRLQ(void);

	//	*****************************< precision >*********************************
	void SetPrecision(double precision);
	void SetSinTiny(double sinT) { sinTiny = sinT; }
};

//	============================================================================
//	======================< class BzzFactorizedQR >=========================
//	============================================================================
class BzzFactorizedQR : public BzzFactorizedQRLQ
{
	friend void Householder
	(int j, int k, double* aux, double* d);
	friend void HouseholderApplyLeft
	(int ri, int rf, int ci, int cf, double** a, double* aux);
	friend void HouseholderApplyRight
	(int ri, int rf, int ci, int cf, double** a, double* aux);
	friend void Givens
	(double* x1, double* x2, double* c, double* s, double* w);
	friend void GivensApply
	(double* y1, double* y2, double c, double s, double w);
	friend void Hessenberg
	(int m, int n, double** a);
	friend void HessenbergP
	(int m, int n, double** a, double** p);
	friend void Bidiagonal
	(int m, int n, double** a, double* d, double* ds);
	friend void BidiagonalR
	(int m, int n, double** a, double* d,
		double* ds, double** an);

private:
	//	============================================================================
	//	*********************< Functions for linear algebra >***********************
	//	============================================================================
	virtual void	 Factorization(void);
	virtual void	 Solution(BzzVector* bx);
	virtual void	 Solution(const BzzVector& b, BzzVector* x);
	virtual void	 TransposeSolution(BzzVector* bx);
	virtual void	 TransposeSolution
	(const BzzVector& b, BzzVector* x);
	virtual void	 Solution(BzzMatrix* BX);
	virtual void	 Solution(const BzzMatrix& B, BzzMatrix* X);
	virtual void	 TransposeSolution(BzzMatrix* BX);
	virtual void	 TransposeSolution(const BzzMatrix& B,
		BzzMatrix* X);

	virtual double	Condition(void);
	BzzVector normColumn;

public:
	//	============================================================================
	//	***************************< constructors >*********************************
	//	============================================================================
		// default constructor type BzzFactorized A;
	BzzFactorizedQR(void)
		: BzzFactorizedQRLQ() {}

	// copy constructor
	BzzFactorizedQR(const BzzFactorized& rval)
		: BzzFactorizedQRLQ(rval) {}

	// constructor A(3,5);
	BzzFactorizedQR(int rows, int columns)
		: BzzFactorizedQRLQ(rows, columns) {}

	// sized and initialised
	//	A(2,3,1.,2.,3.,4.,5.6.);
	BzzFactorizedQR(int rows, int columns, double a11, ...);

	// constructor A(3,5,w);
	BzzFactorizedQR(int rows, int columns, double* initvalues)
		: BzzFactorizedQRLQ(rows, columns, initvalues) {}

	// constructor from BzzMatrix
	BzzFactorizedQR(const BzzMatrix& rval)
		: BzzFactorizedQRLQ(rval) {}

	// make a submatrix with rows,columns
	BzzFactorizedQR(int rows, int columns, const BzzMatrix& rval)
		: BzzFactorizedQRLQ(rows, columns, rval) {}

	//as above, commencing from irow,jcol
	BzzFactorizedQR(int rows, int columns,
		int irow, int jcol, const BzzMatrix& rval)
		: BzzFactorizedQRLQ(rows, columns, irow, jcol, rval) {}

	//	============================================================================
	//	******************************< destructor >********************************
	//	============================================================================
	~BzzFactorizedQR(void);

	//	============================================================================
	//	**************************< assignment operators >**************************
	//	============================================================================
	void operator = (const BzzFactorizedQR& rval);
	BzzFactorizedQR& operator = (const BzzMatrix& rval);

	//	============================================================================
	//	================< Functions which do not modify the matrix >================
	//	============================================================================
	void GetBzzMatrixQ(BzzMatrix* Q); // A = QR
	void GetBzzMatrixP(BzzMatrix* P);	// PA = R
	void GetBzzMatrixR(BzzMatrixRight* R);
	void GetBzzVectorD(BzzVector* d);
	void GetNormColumns(BzzVector* normC);
	BzzVectorInt GetLinearCombinations(void);
	void GetLinearCombinations(BzzVectorInt* eq);
	void PrintLinearCombinations(void);
	int GetTheNumberOfIndependentColumns(void);
	void GetResiduals(BzzVector* residuals); // nRows > nColumns
};

//	============================================================================
//	===================< class BzzFactorizedLQ >============================
//	============================================================================
class BzzFactorizedLQ : public BzzFactorizedQRLQ
{
	friend void Householder
	(int j, int k, double* aux, double* d);
	friend void HouseholderApplyLeft
	(int ri, int rf, int ci, int cf, double** a, double* aux);
	friend void HouseholderApplyRight
	(int ri, int rf, int ci, int cf, double** a, double* aux);
	friend void Givens
	(double* x1, double* x2, double* c, double* s, double* w);
	friend void GivensApply
	(double* y1, double* y2, double c, double s, double w);
	friend void Hessenberg
	(int m, int n, double** a);
	friend void HessenbergP
	(int m, int n, double** a, double** p);
	friend void Bidiagonal
	(int m, int n, double** a, double* d, double* ds);
	friend void BidiagonalR
	(int m, int n, double** a, double* d,
		double* ds, double** an);
	friend char AnalyzeDependencies(BzzFactorizedLQ* E, BzzVector* e,
		BzzMatrix& D, BzzVector& d, BzzVector& L, BzzVector& U,
		BzzVectorInt* ie, BzzVectorInt* iv, BzzVectorInt* id);

private:
	//	============================================================================
	//	*******************< Functions for linear algebra >*************************
	//	============================================================================
	virtual void	 Factorization(void);
	virtual void	 Solution(BzzVector* bx);
	virtual void	 Solution(const BzzVector& b, BzzVector* x);
	virtual void	 TransposeSolution(BzzVector* bx);
	virtual void	 TransposeSolution
	(const BzzVector& b, BzzVector* x);
	virtual void	 Solution(BzzMatrix* BX);
	virtual void	 Solution(const BzzMatrix& B, BzzMatrix* X);
	virtual void	 TransposeSolution(BzzMatrix* BX);
	virtual void	 TransposeSolution(const BzzMatrix& B,
		BzzMatrix* X);
	virtual double	Condition(void);
	BzzVector normRow;
	int iIntrinsicConditionNumber;
	double intrinsicConditionNumber;

public:
	char StartingSystem(BzzMatrix& W);
	char StartingSystem(BzzMatrixSparseLockedByRows& W);
	char SolveSystem(BzzVector* bx);
	char SolveSystem(BzzVector* newRow, BzzVector* bx);
	void ResetSystem(void);
	void DeleteLastRow(void);

	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
		// default constructor type BzzFactorized A;
	BzzFactorizedLQ(void)
		: BzzFactorizedQRLQ() {}

	// copy constructor
	BzzFactorizedLQ(const BzzFactorized& rval)
		: BzzFactorizedQRLQ(rval) {}

	// constructor A(3,5);
	BzzFactorizedLQ(int rows, int columns)
		: BzzFactorizedQRLQ(rows, columns) {}

	// sized and initialised
	//	A(2,3,1.,2.,3.,4.,5.6.);
	BzzFactorizedLQ(int rows, int columns, double a11, ...);

	// constructor A(3,5,w);
	BzzFactorizedLQ(int rows, int columns, double* initvalues)
		: BzzFactorizedQRLQ(rows, columns, initvalues) {}

	// constructor from BzzMatrix
	BzzFactorizedLQ(const BzzMatrix& rval)
		: BzzFactorizedQRLQ(rval) {}

	// makes a submatrix with rows,columns
	BzzFactorizedLQ(int rows, int columns, const BzzMatrix& rval)
		: BzzFactorizedQRLQ(rows, columns, rval) {}

	//as above, commencing from irow,jcol
	BzzFactorizedLQ(int rows, int columns,
		int irow, int jcol, const BzzMatrix& rval)
		: BzzFactorizedQRLQ(rows, columns, irow, jcol, rval) {}

	//	============================================================================
	//	*****************************< destructor >*********************************
	//	============================================================================
	~BzzFactorizedLQ(void);

	//	============================================================================
	//	************************< assignment operators >****************************
	//	============================================================================
	void operator = (const BzzFactorizedLQ& rval);
	BzzFactorizedLQ& operator = (const BzzMatrix& rval);

	//	============================================================================
	//	==============< Functions which do not modify the matrix >==================
	//	============================================================================
	void GetBzzMatrixQ(BzzMatrix* Q); // A = LQ
	void GetBzzMatrixP(BzzMatrix* P);	// AP = L
	void GetBzzMatrixL(BzzMatrixLeft* L);
	void GetBzzMatrixN(BzzMatrix* N); // Null Space
	void GetBzzVectorD(BzzVector* d);
	void GetNormRows(BzzVector* normR);
	BzzVectorInt GetLinearCombinations(void);
	void GetLinearCombinations(BzzVectorInt* eq);
	void PrintLinearCombinations(void);
	int GetTheNumberOfIndependentEquations(void);
	void NTProduct(BzzVector* x);
	void NTProduct(BzzVector& x, BzzVector* y);
	void NProduct(BzzVector* x);
	void NProduct(BzzVector& x, BzzVector* y);
	//	void NTSNProduct(BzzMatrixSymmetric &S,BzzMatrixSymmetric *W);
	double IntrinsicConditionNumber(void)
	{
		if (iIntrinsicConditionNumber != 0)
			return intrinsicConditionNumber;
		else
			return 0.;
	}
};
#endif // BZZ_FACTORIZEDQRLQ_DOUBLE_HPP