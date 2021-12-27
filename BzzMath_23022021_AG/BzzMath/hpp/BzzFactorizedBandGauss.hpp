// BZZMATH: Release 7.0

//	=====================< BzzFactorizedBandGauss.hpp >=====================
//	* BzzFactorizedBandGauss: Class for solution of linear banded systems	*
// * Description:																					*
//	*					Metodi Numerici e Software in C++ (Capitolo 6, 7)				*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
// * Examples: BzzMath\Examples\BzzMathBasic\LinearSystems\							*
//					FactorizedBandGauss\FactorizedBandGauss.cpp				*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	10-1995	Date Written.
//	01-1998	Added SetDiagonal.
//	11-1999	Bug fixed in TransposeSolve function.

////////////////// Release 6.0
//	07-2009	BzzFactorizedBandGauss parallelized wth OpenMP.

//	============================================================================
//	******* Functions for solving linear sparse systems:								*
//	* Solve(&A,&bx);																				*
//	* TransposeSolve(&A,&bx);																	*
//	****************************************************************************

#ifndef BZZ_FACTORIZED_DOUBLE_BAND_HPP
#define BZZ_FACTORIZED_DOUBLE_BAND_HPP

class BzzMatrixSparse;
class BzzVectorInt;
class BzzMatrixIntBand;

//	============================================================================
//	========================< class BzzFactorizedBandGauss >===========================
//	============================================================================

class BzzFactorizedBandGauss : public BzzBaseClass
{
	//	============================================================================
	//	======================< Functions for linear algebra >======================
	//	============================================================================
	friend void Solve
	(BzzFactorizedBandGauss* A, BzzVector* bx);
	friend void Solve
	(BzzFactorizedBandGauss* A, BzzVector& b, BzzVector* x);

	friend void Solve
	(BzzFactorizedBandGauss* A, BzzMatrix* BX);
	friend void Solve
	(BzzFactorizedBandGauss* A, BzzMatrix& B, BzzMatrix* X);

	friend void TransposeSolve
	(BzzFactorizedBandGauss* A, BzzVector* bx);
	friend void TransposeSolve
	(BzzFactorizedBandGauss* A, BzzVector& b, BzzVector* x);

	friend void TransposeSolve
	(BzzFactorizedBandGauss* A, BzzMatrix* BX);
	friend void TransposeSolve
	(BzzFactorizedBandGauss* A, BzzMatrix& B, BzzMatrix* X);

private:
	enum BandFactorizationStatus
	{
		UNFACTORIZED,
		FACTORIZED
	}factorizationStatus;

	enum BandTransposeStatus
	{
		REGULAR,
		TRANSPOSE
	}transposeStatus;

	static const char* const BZZ_ERROR;
	static const char* const BZZ_ERR_DESTROYED;
	static int count; // per whoAmI
	int numRows, numColumns;
	int size;
	int	lowerBand,
		upperBand,
		totalBand,
		originalUpperBand,
		swapBand;

	double** matrix;

	BzzVectorInt	lowerRow,
		upperRow,
		lowerColumn,
		upperColumn,
		originalUpperColumn,
		indx,
		indxy;

	BzzMatrixIntBand K;

	int whoAmI;
	char singular;	// singular = 1 when singular
	char pivoting;
	int signd;
	double norm;

	// initialising constructors
	void Initialize(int m, int n, int low, int up);

	// re-initialising
//	void ReInitialize(int m,int n,int low,int up);

	// deinitialisation
	void Deinitialize(void);

	// preparing assignments
	void PrepCopy(int rRows, int rColumns);

	// preparing the OnlySolve
	void PrepOnlySolve(void);

	// preparing the Solve
	void PrepSolve(void);

	// calculating the norm
	void Norm(void);

	void Copy(BzzMatrixSparse& rval);

	//	void PrepCopy(int rRows,int rColumns);
	//	BzzFactorizedBandGauss(char,int rows,int columns);
	void Factorization(void);
	void Solution(BzzVector* bx);
	void TransposeSolution(BzzVector* bx);

public:
	double factorizationTime, solutionTime;
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
		// default constructor
	BzzFactorizedBandGauss(void);

	// copy constructor
	BzzFactorizedBandGauss(const BzzFactorizedBandGauss& rval);

	// constructor from BzzMatrixBand
	BzzFactorizedBandGauss(BzzMatrixBand& rval);

	// constructor from BzzMatrixSparse
	BzzFactorizedBandGauss(BzzMatrixSparse& rval);

	// sizing constructor
	BzzFactorizedBandGauss(int rows, int columns, int low, int up);

	BzzFactorizedBandGauss(char* file);
	void operator()(char* file);

	//	============================================================================
	//	********************************< destructor >******************************
	//	============================================================================
	~BzzFactorizedBandGauss(void);

	//	============================================================================
	//	*****************************< Access functions >***************************
	//	============================================================================
		// number of rows
	int Rows(void) const
	{
		return numRows;
	}

	// number of columns
	int Columns(void) const
	{
		return numColumns;
	}

	int WhoAmI(void) const { return whoAmI; }

	// assigns and receives vector values with control
	double& operator ()
		(int row, int column);

	// assigns and receives vector values without control
	double* operator [] (int r)
	{
		return matrix[r];
	}

	//	============================================================================
	//	***************************< assignment operators >*************************
	//	============================================================================
	BzzFactorizedBandGauss& operator =
		(BzzFactorizedBandGauss& rval);

	BzzFactorizedBandGauss& operator =
		(BzzMatrixBand& rval);

	BzzFactorizedBandGauss& operator =
		(BzzMatrixSparse& rval);

	// transforms a BzzMatrixBand in BzzFactorizedBandGauss
	// BzzMatrixSparse is destroyed
	friend void ReplaceBzzMatrixWithBzzFactorized
	(BzzMatrixBand* lval, BzzFactorizedBandGauss* rval);

	// transforms a BzzMatrixSparse in BzzFactorizedBandGauss
	// BzzMatrixSparse is destroyed
	friend void ReplaceBzzMatrixWithBzzFactorized
	(BzzMatrixSparse* lval, BzzFactorizedBandGauss* rval);

	// re-initialising
	void ReInitialize(int m, int n, int low, int up);
	void operator()(int m, int n, int low, int up)
	{
		ReInitialize(m, n, low, up);
	}

	//	============================================================================
	//	========================< Non-modifying functions	=========================
	//	============================================================================

	//	*********************************< BzzPrint >*******************************
	virtual void ObjectBzzPrint(void);

	//	*****************************< Max and Min >********************************
		//to have the position Max(im,jc)
	double Max(int* imax = 0, int* jmax = 0);
	double MaxAbs(int* imax = 0, int* jmax = 0);

	double Min(int* imin = 0, int* jmin = 0);
	double MinAbs(int* imin = 0, int* jmin = 0);

	//	============================================================================
	//	=====================< Modifying Functions >================================
	//	============================================================================
	//	void RemoveElement(int row,int column);
	//	void RemoveAllElementsInRow(int row);
	//	void RemoveAllElementsInColumn(int column);
	friend void Delete(BzzFactorizedBandGauss* A);
	void SetToZeroAndReset(void);// All the matrix elements are zeroed
//	void CleanMatrix(double eps); // eliminates those <= eps
	double Determinant(void);
	double ConditionNumber(void);
	char Singular(void) { return singular; } // 1 if singular 0 if not
//	int CountElements(void);
	void SetDiagonal(int j, BzzVector& rval);
	void SetDiagonal(int j, float xf);
	void SetDiagonal(int j, double xf);

	//	============================================================================
	//	===================< System solution functions >============================
	//	============================================================================
	//	void SolveRight(BzzVector *bx);
	//	void SolveLeft(BzzVector *bx);
};

#endif // BZZ_FACTORIZED_DOUBLE_BAND_HPP