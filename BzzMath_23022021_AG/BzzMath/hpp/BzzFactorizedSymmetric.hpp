// BxxMath release 3.1 version alfa
// TODO: BzzFactorizedSymmetric::ConditionNumber
// TODO: BzzFactorizedSymmetric::Determinant
// TODO: BzzFactorizedSymmetric::ChangeFactorization
// TODO: BzzFactorizedSymmetric::GetBzzVectorDAndBzzMatrixLeftL

//	======================< BzzFactorizedSymmetric.hpp >====================
//	* BzzFactorizedSymmetricBase: class for symmetrical matrices				*
// *               factorization in double precision									*
// * Description:																					*
//	*					Metodi Numerici e Software in C++ (Capitoli 6, 7)				*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
// *																									*
// * Examples: BzzMath\examples\BzzMathBasic\LinearSystems\							*
// *				FactorizedSymmetric\FactorizedSymmetric.cpp				*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	06-1994	date written.
//	11-1994	Conversion to double precision done.
//	02-1994	Added the Bzz prefix to the names of the classes.
//	06-1996	New version for non definite positive matrix factorization
//	01-1998	Added BzzFactorizedSymmetricPositive class
//	01-1998	New version for non definite positive matrix factorization
//	03-1999	Added GetBzzVectorDAndBzzMatrixLeftL

//	============================================================================
//	*						factorizationStatus		matrixStatus							*
//	* initially				 UNFACTORIZED			 AVAILABLE									*
//	*																									*
//	* OnlySolve					FACTORIZED				DESTROYED								*
//	* Solve						FACTORIZED				AVAILABLE								*
//	* TRansposeSolve			FACTORIZED				AVAILABLE								*
//	* HyperSolve				FACTORIZED				AVAILABLE								*
//	****************************************************************************
//	* Base class for solving linear systems with symmetrical defite            *
// * positive matrix																				*
//	* Valid functions for all factorisations:												*
//	* Solve(A,b,&x);																				*
//	* Solve(A,&bx);	else	Solve(A,b,&b);													*
//	* Solve(&A,b,&x);																				*
//	* Solve(&A,&bx); else	Solve(&A,b,&b);												*
//	* Solve(A,B,&X);																				*
//	* Solve(A,&BX);	else	Solve(A,B,&B);													*
//	* Solve(&A,B,&X);																				*
//	* Solve(&A,&BX); else	Solve(&A,B,&B);												*
//	* HyperSolve(A,b,&x);																		*
//	* HyperSolve(A,&bx); or HyperSolve(A,b,&b);											*
//	* double det = A.Determinant(void)														*
//	* double cond = A.ConditionNumber(void)												*
//	* char s = A.Singular(void)																*
//	* char s = A.Positive(void)																*
//	* BzzMatrix inv; A.Inverse(&inv)													*
//	* BzzMatrix inv = A.Inverse(void)												*
//	* Delete(&A);																					*
//	* ChangeDimensions(newr,newc,&A);														*
//	* SumRankOne(&A,alfa,&z);																	*
//	* SumRankTwo(&A,alfa,&z,beta,&s);														*
//	****************************************************************************
//	* Can be defined as objects of the classes:											*
//	* BzzFactorizedSymmetric, BzzFactorizedCholesky, and BzzFactorizedGillMurray		*
//	****************************************************************************
//	* The matrix A is DESTROYED if you use													*
//	* Solve(&A,...)																				*
//	* because it is overwritten with factorization										*
//	****************************************************************************
//	* The coefficients of the matrix are accessible if									*
//	* matrixStatus != DESTROYED																*
//	* float val = A(i,j);																		*
//	****************************************************************************
//	* The matrix can be modified only if matrixStatus != DESTROYED					*
//	* A(i,j) = val;																				*
//	****************************************************************************
//	* Other functions for objects of BzzFactorizedGillMurray class:					*
//	* char positive = A.Positive();															*
//	* positive = 1 if definite positive														*
//	* positive = 0 if not																		*
//	* char positive = A.NegativeCurvature(&x);											*
//	****************************************************************************

#ifndef BZZ_FACTORIZEDSYMM_DOUBLE_HPP
#define BZZ_FACTORIZEDSYMM_DOUBLE_HPP

//	============================================================================
//	====================< class BzzFactorizedSymmetricBase >=========================
//	============================================================================

class BzzFactorizedSymmetricBase : public BzzBaseClass
{
	//	============================================================================
	//	======================< Functions for linear algebra >======================
	//	============================================================================

	friend void Solve
	(BzzFactorizedSymmetricBase& A, const BzzVector& b, BzzVector* x);
	friend void Solve
	(BzzFactorizedSymmetricBase& A, BzzVector* bx);
	friend void Solve
	(BzzFactorizedSymmetricBase* A, const BzzVector& b, BzzVector* x);
	friend void Solve
	(BzzFactorizedSymmetricBase* A, BzzVector* bx);
	friend void Solve
	(BzzFactorizedSymmetricBase& A, const BzzMatrix& B, BzzMatrix* X);
	friend void Solve
	(BzzFactorizedSymmetricBase& A, BzzMatrix* BX);
	friend void Solve
	(BzzFactorizedSymmetricBase* A, const BzzMatrix& B, BzzMatrix* X);
	friend void Solve
	(BzzFactorizedSymmetricBase* A, BzzMatrix* BX);

	friend void HyperSolve
	(BzzFactorizedSymmetricBase& A, const BzzVector& b, BzzVector* x);
	friend void HyperSolve
	(BzzFactorizedSymmetricBase& A, BzzVector* bx);

	friend void CholeskyFactorization
	(int n, double** a);
	friend void CholeskySolution
	(int n, double** a, double* b);
	friend double CholeskyNormInvEst
	(int n, double** a);

	friend void SymmetricFactorization
	(int n, double** a);
	friend void SymmetricSolution
	(int n, double** a, double* b);
	friend double SymmetricNormInvEst
	(int n, double** a);

	friend void SymmetricPositiveFactorization
	(int n, double** a);

	friend char GillMurrayFactorization
	(int n, double** a, BzzVector& cor);
	friend void GillMurraySolution
	(int n, double** a, double* b);
	friend double GillMurrayNormInvEst
	(int n, double** a);

protected:
	enum SymmFactorizationStatus
	{
		UNFACTORIZED,
		FACTORIZED
	}factorizationStatus;
	enum SymmBzzMatrixStatus
	{
		AVAILABLE,
		DESTROYED,
		MODIFIED
	}matrixStatus;

	static const char* const BZZ_ERROR;
	static const char* const BZZ_ERR_DESTROYED;
	static int count; // for whoAmI

	double** matrix;
	double** factorized;
	int numRows, numColumns;
	int size;
	int whoAmI;
	double norm;
	double tollerance;
	char positive;
	// eigenvalues only for BzzFactorizedSymmetric
	// 0 not yet evaluated
	// 1 definte positive								0 0 1
	// 2 definite negative								0 1 0
	// 3 some positive some negative					0 1 1
	// 4 all eigenvalues are null						1 0 0
	// 5 some positive some null						1 0 1
	// 6 some negative some null						1 1 0
	// 7 some positive some negative aome null	1 1 1
	int eigenvalues;

	// initialising constructors
	void Initialize(int mrows, int ncolumns);

	// initialise factorized if requested
	virtual void BzzFactorizedInitialize(void);

	// deinitialisation
	void Deinitialize(void);

	// deinitialise factorized if requested
	void BzzFactorizedDeinitialize(void);

	// preparing assignments
	void PrepCopy(int rRows, int rColumns);

	// protected constructor type BzzFactorized A('*',3,5);
	BzzFactorizedSymmetricBase(char, int rows, int columns);

	// preparing the OnlySolve
	void PrepOnlySolve(void);

	// preparing the Solve
	void PrepSolve(void);

	// calculating the norm
	virtual void Norm(void);

	//	============================================================================
	//	======================< Functions for linear algebra >======================
	//	============================================================================
		// For each factorisation
	virtual void Factorization(void) = 0;
	virtual void Solution(int n, double** a, double* b) = 0;
	virtual void Solution(BzzVector* bx);
	void Solution(const BzzMatrix& B, BzzMatrix* X);
	void Solution(BzzMatrix* BX);
	virtual void ChangeFactorization
	(int n, double** a, double alfa, double* z) = 0;

public:
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
		// default constructor BzzFactorizedSymmetricBase A;
	BzzFactorizedSymmetricBase(void);

	// copy constructor
	BzzFactorizedSymmetricBase(const BzzFactorizedSymmetricBase& rval);

	// constructor A(3,3);
	BzzFactorizedSymmetricBase(int rows, int columns);

	// constructor A(3,3,w);
	BzzFactorizedSymmetricBase(int rows, int columns, double* initvalues);

	// constructor from BzzMatrixSymmetric
	BzzFactorizedSymmetricBase(const BzzMatrixSymmetric& rval);

	//	============================================================================
	//	******************************< destructor >********************************
	//	============================================================================
	virtual ~BzzFactorizedSymmetricBase(void);

	//	============================================================================
	//	**********************< Non-modifying access functions >********************
	//	============================================================================
		// number of rows
	int Rows(void) const { return numRows; }

	// number of columns
	int Columns(void) const { return numColumns; }

	int WhoAmI(void) const { return whoAmI; }

	//	============================================================================
	//	***********************< Modifying access functions >***********************
	//	============================================================================
		// assigning and receiving values with control
	double& operator () (int row, int col);

	//	============================================================================
	//	**************************< assignment operators >**************************
	//	============================================================================
		// transforms a BzzMatrixSymmetric in BzzFactorizedSymmetricBase
		// creates an independent matrix
	virtual void CopyFromBzzMatrixSymm
	(const BzzMatrixSymmetric& rval);

	// transforms a BzzMatrixSymmetric in BzzFactorizedSymmetricBase
	// BzzMatrixSymmetric is destroyed
	friend void ReplaceBzzMatrixWithBzzFactorized
	(BzzMatrixSymmetric* lval, BzzFactorizedSymmetricBase* rval);

	//	============================================================================
	//	======================< Non-modifying functions >===========================
	//	============================================================================
	//	*********************************< BzzPrint >*******************************
	virtual void ObjectBzzPrint(void);

	//	============================================================================
	//	=========================< Modifying functions >============================
	//	============================================================================
	friend void Delete(BzzFactorizedSymmetricBase* result);
	friend void ChangeDimensions(int rows, int columns,
		BzzFactorizedSymmetricBase* result, char);
	double Determinant(void);
	double ConditionNumber(void);
	char Positive(void);
	void Inverse(BzzMatrix* inv);
	BzzMatrix Inverse(void);
	friend void SumRankOne(BzzFactorizedSymmetricBase* A,
		double a, BzzVector* z);
	friend void SumRankTwo(BzzFactorizedSymmetricBase* A,
		double a, BzzVector* z, double b, BzzVector* s);
};

//	============================================================================
//	=======================< class BzzFactorizedCholesky >==================
//	============================================================================

class BzzFactorizedCholesky : public BzzFactorizedSymmetricBase
{
public:
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
		// default constructor
	BzzFactorizedCholesky(void)
		: BzzFactorizedSymmetricBase() {
		tollerance = 0.;
	}

	// copy constructor
	BzzFactorizedCholesky(const BzzFactorizedSymmetricBase& rval)
		: BzzFactorizedSymmetricBase(rval) {
		tollerance = 0.;
	}

	// constructor A(3,5);
	BzzFactorizedCholesky(int rows, int columns)
		: BzzFactorizedSymmetricBase(rows, columns) {
		tollerance = 0.;
	}

	// constructor A(2,2,1.,2.,3.,4.);
	BzzFactorizedCholesky(int rows, int columns, double a11, ...);

	// constructor A(3,5,w);
	BzzFactorizedCholesky(int rows, int columns, double* initvalues)
		: BzzFactorizedSymmetricBase(rows, columns, initvalues) {
		tollerance = 0.;
	}

	// constructor from BzzMatrixSymmetric
	BzzFactorizedCholesky(const BzzMatrixSymmetric& rval)
		: BzzFactorizedSymmetricBase(rval) {
		tollerance = 0.;
	}

	//	============================================================================
	//	******************************< destructor >********************************
	//	============================================================================
	~BzzFactorizedCholesky(void);

	//	============================================================================
	//	**************************< assignment operators >**************************
	//	============================================================================
	void operator = (const BzzFactorizedCholesky& rval);
	BzzFactorizedCholesky& operator =
		(const BzzMatrixSymmetric& rval);

	//	============================================================================
	//	======================< Functions for linear algebra >======================
	//	============================================================================
	virtual void Factorization(void);
#if BZZ_COMPILER >= 100
	using BzzFactorizedSymmetricBase::Solution;
#endif
	virtual void Solution(int n, double** a, double* b);
	virtual void ChangeFactorization
	(int n, double** a, double alfa, double* z);
	void GetBzzVectorDAndBzzMatrixLeftL(BzzVector* d,
		BzzMatrixLeft* L);
};

//	============================================================================
//	================< class BzzFactorizedSymmetricPositive >================
//	============================================================================

class BzzFactorizedSymmetricPositive : public BzzFactorizedCholesky
{
public:
	//	============================================================================
	//	******************************< constructors >******************************
	//	============================================================================
		// default constructor
	BzzFactorizedSymmetricPositive(void)
		: BzzFactorizedCholesky() {
		tollerance = 0.;
	}

	// copy constructor
	BzzFactorizedSymmetricPositive(const BzzFactorizedSymmetricBase& rval)
		: BzzFactorizedCholesky(rval) {
		tollerance = 0.;
	}

	// constructor A(3,5);
	BzzFactorizedSymmetricPositive(int rows, int columns)
		: BzzFactorizedCholesky(rows, columns) {
		tollerance = 0.;
	}

	// constructor A(2,2,1.,2.,3.,4.);
	BzzFactorizedSymmetricPositive(int rows, int columns, double a11, ...);

	// constructor A(3,5,w);
	BzzFactorizedSymmetricPositive(int rows, int columns, double* initvalues)
		: BzzFactorizedCholesky(rows, columns, initvalues) {
		tollerance = 0.;
	}

	// constructor from BzzMatrixSymmetric
	BzzFactorizedSymmetricPositive(const BzzMatrixSymmetric& rval)
		: BzzFactorizedCholesky(rval) {
		tollerance = 0.;
	}

	//	============================================================================
	//	**************************< assignment operators >**************************
	//	============================================================================
	void operator = (const BzzFactorizedSymmetricPositive& rval);
	BzzFactorizedSymmetricPositive& operator =
		(const BzzMatrixSymmetric& rval);

	//	============================================================================
	//	=======================< Functions for linear algebra >=====================
	//	============================================================================
	virtual void Factorization(void);
};

//	============================================================================
//	====================< class BzzFactorizedSymmetric >====================
//	============================================================================

class BzzFactorizedSymmetric : public BzzFactorizedCholesky
{
private:
	BzzVectorInt kMemo, pivot;
	BzzVector lMemo, rMemo;
	char pos, neg, nul;

public:
	//	============================================================================
	//	******************************< constructors >******************************
	//	============================================================================
		// default constructor
	BzzFactorizedSymmetric(void)
		: BzzFactorizedCholesky() {
		tollerance = 0.;
	}

	// copy constructor
	BzzFactorizedSymmetric(const BzzFactorizedSymmetricBase& rval)
		: BzzFactorizedCholesky(rval) {
		tollerance = 0.;
	}

	// constructor A(3,5);
	BzzFactorizedSymmetric(int rows, int columns)
		: BzzFactorizedCholesky(rows, columns) {
		tollerance = 0.;
	}

	// constructor A(2,2,1.,2.,3.,4.);
	BzzFactorizedSymmetric(int rows, int columns, double a11, ...);

	// constructor A(3,5,w);
	BzzFactorizedSymmetric(int rows, int columns, double* initvalues)
		: BzzFactorizedCholesky(rows, columns, initvalues) {
		tollerance = 0.;
	}

	// constructor from BzzMatrixSymmetric
	BzzFactorizedSymmetric(const BzzMatrixSymmetric& rval)
		: BzzFactorizedCholesky(rval) {
		tollerance = 0.;
	}

	//	============================================================================
	//	**************************< assignment operators >**************************
	//	============================================================================
	void operator = (const BzzFactorizedSymmetric& rval);
	BzzFactorizedSymmetric& operator =
		(const BzzMatrixSymmetric& rval);
	int Eigenvalues(void) { return eigenvalues; }

	//	============================================================================
	//	=======================< Functions for linear algebra >=====================
	//	============================================================================
	virtual void Factorization(void);
#if BZZ_COMPILER >= 100
	using BzzFactorizedSymmetricBase::Solution;
#endif
	virtual void Solution(int n, double** a, double* b);
	virtual void ChangeFactorization
	(int n, double** a, double alfa, double* z);
	double ConditionNumber(void);
	double Determinant(void);
	void GetBzzVectorDAndBzzMatrixLeftL(BzzVector* d,
		BzzMatrixLeft* L)
	{
		BzzError("%s Function not available", BZZ_ERROR);
	}
};

//	============================================================================
//	==================< class BzzFactorizedGillMurray >=====================
//	============================================================================

class BzzFactorizedGillMurray : public BzzFactorizedSymmetricBase
{
private:
	BzzVector correction;
public:
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
		// default constructor
	BzzFactorizedGillMurray(void)
		: BzzFactorizedSymmetricBase() {
		tollerance = sqrt(BZZ_TINY_FLOAT);
	}

	// copy constructor
	BzzFactorizedGillMurray(const BzzFactorizedSymmetricBase& rval)
		: BzzFactorizedSymmetricBase(rval) {
		tollerance = sqrt(BZZ_TINY_FLOAT);
	}

	// constructor A(3,3);
	BzzFactorizedGillMurray(int rows, int columns)
		: BzzFactorizedSymmetricBase(rows, columns)
	{
		tollerance = sqrt(BZZ_TINY_FLOAT);
	}

	// constructor A(2,2,1.,2.,3.,4.);
	BzzFactorizedGillMurray(int rows, int columns, double a11, ...);

	// constructor A(3,5,w);
	BzzFactorizedGillMurray
	(int rows, int columns, double* initvalues)
		: BzzFactorizedSymmetricBase(rows, columns, initvalues)
	{
		tollerance = sqrt(BZZ_TINY_FLOAT);
	}

	// constructor from BzzMatrixSymmetric
	BzzFactorizedGillMurray(const BzzMatrixSymmetric& rval)
		: BzzFactorizedSymmetricBase(rval)
	{
		tollerance = sqrt(BZZ_TINY_FLOAT);
	}

	//	============================================================================
	//	******************************< destructor >********************************
	//	============================================================================
	~BzzFactorizedGillMurray(void);

	//	============================================================================
	//	**************************< assignment operators >**************************
	//	============================================================================
	void operator = (const BzzFactorizedGillMurray& rval);
	BzzFactorizedGillMurray& operator =
		(const BzzMatrixSymmetric& rval);

	//	============================================================================
	//	======================< Functions for linear algebra >======================
	//	============================================================================
	virtual void Factorization(void);
#if BZZ_COMPILER >= 100
	using BzzFactorizedSymmetricBase::Solution;
#endif
	virtual void Solution(int n, double** a, double* b);
	virtual void ChangeFactorization
	(int n, double** a, double alfa, double* z);
	char NegativeCurvature(BzzVector* x);
	void GetDiagonalCorrections(BzzVector* cor)
	{
		*cor = correction;
	}
};

//	============================================================================
//	================< class BzzFactorizedSymmetricBand >================
//	============================================================================

class BzzFactorizedSymmetricBand : public BzzFactorizedSymmetricPositive
{
private:
	int band;
	void InitializeBand(int r, int c, int b);
	friend void Delete(BzzFactorizedSymmetricBand* S);
	//	void Norm(void){norm = 1.;BzzError();}
public:
	//	============================================================================
	//	******************************< constructors >******************************
	//	============================================================================
		// default constructor
	BzzFactorizedSymmetricBand(void);

	// copy constructor
	BzzFactorizedSymmetricBand(const BzzFactorizedSymmetricBase& rval);

	// constructor A(5,5,2);
	BzzFactorizedSymmetricBand(int rows, int columns, int b);

	// constructor A(1,5,5,1.,2.,3.,4.);
	BzzFactorizedSymmetricBand(int rows, int columns, int b, double a11, ...);

	// constructor from BzzMatrixSymmetric
//	BzzFactorizedSymmetricBand(BzzMatrixSymmetricBand &rval);

//	============================================================================
//	**************************< assignment operators >**************************
//	============================================================================
	double& operator ()(int row, int col);
	virtual void ObjectBzzPrint(void);
	int Band(void) { return band; }
	//	void operator = (BzzFactorizedSymmetricBand &rval);
	BzzFactorizedSymmetricBand& operator =
		(BzzMatrixSymmetricBand& rval);
	void CopyFromBzzMatrixSymmetricBand(BzzMatrixSymmetricBand& rval);
	friend void ReplaceBzzMatrixWithBzzFactorized(BzzMatrixSymmetricBand* M,
		BzzFactorizedSymmetricBand* F);
	void SetDiagonal(int j, BzzVector& rval);
	void SetDiagonal(int j, double xf);

	//	============================================================================
	//	=======================< Functions for linear algebra >=====================
	//	============================================================================
	virtual void Norm(void) { norm = 1.; }
	virtual void Factorization(void);
#if BZZ_COMPILER >= 100
	using BzzFactorizedSymmetricBase::Solution;
#endif
	virtual void Solution(int n, double** a, double* b);
	virtual void BzzFactorizedInitialize(void);
};

// Friend functions with default arguments
void ChangeDimensions(int rows, int columns,
	BzzFactorizedSymmetricBase* result, char zero = 0);


#endif // BZZ_FACTORIZEDSYMM_DOUBLE_HPP