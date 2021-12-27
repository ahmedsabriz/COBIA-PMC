// BZZMATH: Release 7.0

//	============================< BzzFactorized.hpp >=======================
//	* Class BzzFactorized: Base class for factorisations							*
// * in double precision																		*
// * Description:																					*
// *					Dal Fortan al C++	(Capitolo 15)										*
// *					by G. Buzzi-Ferraris														*
// *					Addison Wesley(1991)														*
// *							and																	*
// *					Scientific C++	(Chapter 15)											*
// *					by G. Buzzi-Ferraris														*
// *					Addison Wesley(1991)														*
// *							and																	*
//	*					Metodi Numerici e Software in C++ (Capitolo 4, 5)				*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	06-1991	Date Written
//	11-1992	English version
//	03-1994	Added BzzBaseClass for BzzPrint and BzzMessage.
//	11-1994	Conversion to double precision done.
//	02-1995	Added the Bzz prefix to the names of the classes.

//	============================================================================
//	*						factorizationStatus		matrixStatus							*
//	* initially					UNFACTORIZED				AVAILABLE							*
//	*																									*
//	* OnlySolve					FACTORIZED					DESTROYED							*
//	* Solve						FACTORIZED					AVAILABLE							*
//	* TRansposeSolve			FACTORIZED					AVAILABLE							*
//	* HyperSolve				FACTORIZED					AVAILABLE							*
//	****************************************************************************
//	* These classes are used for solving linear systems and connected problems	*
//	* Valid functions for all factorisations:												*
//	* Solve(A,b,&x);																				*
//	* Solve(A,&bx);	else	Solve(A,b,&b);													*
//	* Solve(&A,b,&x);																				*
//	* Solve(&A,&bx); else	Solve(&A,b,&b);												*
//	* Solve(A,B,&X);																				*
//	* Solve(A,&BX);	else	Solve(A,B,&B);													*
//	* Solve(&A,B,&X);																				*
//	* Solve(&A,&BX); else	Solve(&A,B,&B);												*
//	* TransposeSolve(A,b,&x);																	*
//	* TransposeSolve(A,&bx);	or TransposeSolve(A,b,&b);								*
//	* TransposeSolve(&A,b,&x);																	*
//	* TransposeSolve(&A,&bx); or TransposeSolve(&A,b,&b);								*
//	* TransposeSolve(A,B,&X);																	*
//	* TransposeSolve(A,&BX);	else	TransposeSolve(A,B,&B);							*
//	* TransposeSolve(&A,B,&X);																	*
//	* TransposeSolve(&A,&BX); else	TransposeSolve(&A,B,&B);						*
//	* HyperSolve(A,b,&x);																		*
//	* HyperSolve(A,&bx); or HyperSolve(A,b,&b);											*
//	* double det = A.Determinant(void)														*
//	* double cond = A.ConditionNumber(void)													*
//	* char s = A.Singular(void)																*
//	* BzzMatrix inv; A.Inverse(&inv)													*
//	* BzzMatrix inv = A.Inverse(void)												*
//	* Delete(&A);																					*
//	* ChangeDimensions(newr,newc,&A);														*
//	****************************************************************************
//	* Can be defined as objects of the classes:											*
//	* BzzFactorizedGauss, BzzFactorizedCrout									*
//	* BzzFactorizedQR, BzzFactorizedLQ											*
//	* See c:\bzzmath\exampled\dxgauss.cpp for initialising							*
//	****************************************************************************
//	* The matrix A is DESTROYED if you use													*
//	* Solve(&A,...) or TransposeSolve(&A,...)												*
//	* because it is overwritten with factorization										*
//	****************************************************************************
//	* The coefficients of the matrix are accessible if									*
// * matrixStatus != DESTROYED																*
//	* double val = A(i,j);																		*
//	* val =	A.GetValue(i,j);																	*
//	* BzzVector v = A.GetRow(i);														*
//	* v =	A.GetColumn(j);																		*
//	* v = GetDiagonal(j); // j = 0 principal positive at right						*
//	****************************************************************************
//	* The matrix can be modified oly if														*
//	* matrixStatus != DESTROYED																*
//	* A.SetValue(i,j,val);																		*
//	* A(i,j) = val;																				*
//	* A.SetRow(i,v);																				*
//	* A.SetRow(i,xf);																				*
//	* A.SetColumn(j,v);																			*
//	* A.SetColumn(j,xf);																			*
//	* A.SetDiagonal(j,xf); // j = 0 principal positive at right						*
//	* A.SetDiagonal(j,v);  // j = 0 principal positive at right						*
//	****************************************************************************

#ifndef BZZ_FACTORIZED_DOUBLE_HPP
#define BZZ_FACTORIZED_DOUBLE_HPP

//	============================================================================
//	=============================< class BzzFactorized >====================
//	============================================================================

class BzzFactorized : public BzzBaseClass
{
	friend class BzzSave;
	friend class BzzLoad;

	//	============================================================================
	//	======================< Functions for linear algebra >======================
	//	============================================================================
	friend void SolveWithOrdering
	(BzzFactorizedGauss* A, BzzVector* bx);

	friend void Solve
	(BzzFactorized& A, const BzzVector& b, BzzVector* x);
	friend void Solve
	(BzzFactorized& A, BzzVector* bx);
	friend void Solve
	(BzzFactorized* A, const BzzVector& b, BzzVector* x);
	friend void Solve
	(BzzFactorized* A, BzzVector* bx);
	friend void Solve
	(BzzFactorized& A, const BzzMatrix& B, BzzMatrix* X);
	friend void Solve
	(BzzFactorized& A, BzzMatrix* BX);
	friend void Solve
	(BzzFactorized* A, const BzzMatrix& B, BzzMatrix* X);
	friend void Solve
	(BzzFactorized* A, BzzMatrix* BX);

	friend void Solve(BzzFactorizedSVD& A, const BzzVector& b,
		BzzVector* x, int numD);
	friend void Solve(BzzFactorizedSVD* A, const BzzVector& b,
		BzzVector* x, int numD);
	friend void TransposeSolve(BzzFactorizedSVD& A,
		const BzzVector& b, BzzVector* x, int numD);
	friend void TransposeSolve(BzzFactorizedSVD* A,
		const BzzVector& b, BzzVector* x, int numD);

	friend void TransposeSolve
	(BzzFactorized& A, const BzzVector& b, BzzVector* x);
	friend void TransposeSolve
	(BzzFactorized& A, BzzVector* bx);
	friend void TransposeSolve
	(BzzFactorized* A, const BzzVector& b, BzzVector* x);
	friend void TransposeSolve
	(BzzFactorized* A, BzzVector* bx);
	friend void TransposeSolve
	(BzzFactorized& A, const BzzMatrix& B, BzzMatrix* X);
	friend void TransposeSolve
	(BzzFactorized& A, BzzMatrix* BX);
	friend void TransposeSolve
	(BzzFactorized* A, const BzzMatrix& B, BzzMatrix* X);
	friend void TransposeSolve
	(BzzFactorized* A, BzzMatrix* BX);

	friend void HyperSolve
	(BzzFactorized& A, BzzVector& b, BzzVector* x);
	friend void HyperSolve
	(BzzFactorized& A, BzzVector* bx);

protected:
	enum FactorizationStatus
	{
		UNFACTORIZED,
		FACTORIZED
	}factorizationStatus;
	enum BzzMatrixStatus
	{
		AVAILABLE,
		DESTROYED,
		MODIFIED,
		RESETTED
	}matrixStatus;

	static const char* const BZZ_ERROR;
	static const char* const BZZ_ERR_DESTROYED;
	static int count; // for whoAmI

	double** matrix;
	double** factorized;
	int numRows, numColumns;
	int size;
	int whoAmI;
	char singular;	// singular = 1 when singular
	double norm;

	// initialising constructors
	void Initialize(int mrows, int ncolumns);

	// initialise factorized if requested
	void BzzFactorizedInitialize(void);

	// deinitialisation
	void Deinitialize(void);

	// deinitialise factorized if requested
	void BzzFactorizedDeinitialize(void);

	// preparing assignments
	void PrepCopy(int rRows, int rColumns);

	// protected constructor type BzzFactorized A('*',3,5);
	BzzFactorized(char, int rows, int columns);

	// preparing the OnlySolve
	void PrepOnlySolve(void);

	// preparing the Solve
	void PrepSolve(void);

	// calculating the norm
	void Norm(void);

	int iSystemBalancing;
	BzzVector	rowCoefficients,
		columnCoefficients,
		rhsVector,
		diffCoeff;
	BzzMatrix originalMatrix;
	//	int RhsControl(BzzVector &bx);
	//	void BzzBalance(void);

	//	============================================================================
	//	===============< Functions for linear algebra >=============================
	//	============================================================================
		// For each PLR factorisation, Gauss, Qr, SVD
	virtual void Factorization(void) = 0;
	virtual void SpecificInitialize(void) = 0;
	virtual void SpecificDeinitialize(void) = 0;
	virtual void Solution(const BzzVector& b, BzzVector* x) = 0;
	virtual void Solution(BzzVector* bx) = 0;
	virtual void TransposeSolution
	(const BzzVector& b, BzzVector* x) = 0;
	virtual void TransposeSolution(BzzVector* bx) = 0;
	virtual void Solution(const BzzMatrix& B, BzzMatrix* X) = 0;
	virtual void Solution(BzzMatrix* BX) = 0;
	virtual void TransposeSolution(const BzzMatrix& B, BzzMatrix* X) = 0;
	virtual void TransposeSolution(BzzMatrix* BX) = 0;

	virtual double Condition(void) = 0;
	virtual double DeterminantEvaluation(void) = 0;

public:
	//	============================================================================
	//	*************************< constructors >***********************************
	//	============================================================================
		// default constructor type BzzFactorized A;
	BzzFactorized(void);

	// copy constructor
	BzzFactorized(const BzzFactorized& rval);

	// constructor A(3,5);
	BzzFactorized(int rows, int columns);

	// constructor A(3,5,w);
	BzzFactorized(int rows, int columns, double* initvalues);

	// constructor from BzzMatrix
	BzzFactorized(const BzzMatrix& rval);

	// make a submatrix with rows,columns
	BzzFactorized(int rows, int columns, const BzzMatrix& rval);

	// as above, starting from irow,jcol
	BzzFactorized(int rows, int columns,
		int irow, int jcol, const BzzMatrix& rval);

	//	============================================================================
	//	*****************************< destructor >*********************************
	//	============================================================================
	virtual ~BzzFactorized(void);

	//	============================================================================
	//	*********************< Non-modifying access functions >*********************
	//	============================================================================
		// number of rows
	int Rows(void) const { return numRows; }

	// number of columns
	int Columns(void) const { return numColumns; }

	int WhoAmI(void) const { return whoAmI; }

	// receives values with control
	double GetValue(int row, int col) const;
	BzzVector GetRow(int i) const;	// row i of matrix
	BzzVector GetColumn(int j) const; // column j of matrix

	// takes the diagonal i
	// i = 0 principal
	// i > 0 right diagonal
	// i < 0 left diagonal
	BzzVector GetDiagonal(int j) const;

	//	============================================================================
	//	********************< Modifying access functions >**************************
	//	============================================================================
		//assigning values with control
	void SetValue(int row, int column, double val);

	// assigning and receiving values with control
	double& operator () (int row, int col);

	// assigns and receives vector values without control
	double* operator [] (int r)
	{
		return matrix[r];
	}

	//substitutes row i with xf
	void SetRow(int i, const double& xf);

	// row i = rval
	void SetRow(int i, const BzzVector& rval);

	// substitutes column j with xf
	void SetColumn(int j, const double& xf);

	// column j = rval
	void SetColumn(int j, const BzzVector& rval);

	// substitutes diagonal j with xf
	void SetDiagonal(int j, const double& xf);

	// substitutes diagonal j with Vector rval
	void SetDiagonal(int j, BzzVector& rval);

	void SetToZeroAndReset(void);

	//	============================================================================
	//	***********************< assignment operators >*****************************
	//	============================================================================
		// transforms a BzzMatrix in BzzFactorized
		// creates an independent matrix
	virtual void CopyFromBzzMatrix(const BzzMatrix& rval);

	// transforms a BzzMatrix in BzzFactorized
	// BzzMatrix is destroyed
	friend void ReplaceBzzMatrixWithBzzFactorized
	(BzzMatrix* lval, BzzFactorized* rval);

	//	============================================================================
	//	===================< Non-modifying functions >==============================
	//	============================================================================

	//	*******************************< BzzPrint >*********************************
	virtual void ObjectBzzPrint(void);

	//	============================================================================
	//	=======================< Modifying functions >==============================
	//	============================================================================
	friend void Delete(BzzFactorized* result);
	friend void ChangeDimensions(int rows, int columns,
		BzzFactorized* result, char);
	double Determinant(void);
	double ConditionNumber(void);
	char Singular(void); // 1 if singular 0 if not singular
	void Inverse(BzzMatrix* inv);
	BzzMatrix Inverse(void);
	void SolveWithBalancing(BzzMatrix* A, BzzVector* b);
	void SolveWithBalancing(BzzMatrix& A, BzzVector* b);
	void SolveWithBalancing(BzzVector* b);
};

// Friend functions with default arguments
void ChangeDimensions(int rows, int columns,
	BzzFactorized* result, char zero = 0);


#endif // BZZ_FACTORIZED_DOUBLE_HPP