// BZZMATH: Release 7.0
/*
Linear systems with rectangular band structure
WITHOUT PIVOTING!!!!!
Example:
*xxx.........
x*xxxx.......
xx*xxx.......
..x*xxxxx....
....*xxxxxx..
.xxxx*xxxxx..
......*xxxxx.
......x*xxxx.
......xx*xxxx
..xxxxxxx*xxx
.......xxx*xx
.........xx*x
......xxxxxx*
*/
///////////////////////////////////////
// TODO:
// ConditionNumber
///////////////////////////////////////

//	================< BzzFactorizedStaircaseGauss.hpp >==================
//	* BzzFactorizedStaircaseGauss: Class for solution of square linear	*
// * systems with rectangle band structure												*
// * Examples: c:\bzzmath\examples\exfactRectangle.cpp									*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	08-2002	Date Written.

////////////////// Release 6.0
// 11-2008	Added TransposeSolve function.

//	============================================================================
//	******* Functions for solving square linear sparse systems:						*
//	* Solve(&A,b,&x);																				*
//	* Solve(&A,&bx);																				*
//	* TransposeSolve(&A,b,&x);																	*
//	* TransposeSolve(&A,&bx);																	*
//	****************************************************************************

#ifndef BZZ_FACTORIZED_DOUBLE_STAIRCASE_GAUSS_HPP
#define BZZ_FACTORIZED_DOUBLE_STAIRCASE_GAUSS_HPP

//	============================================================================
//	=====================< class BzzFactorizedStaircaseGauss >==========================
//	============================================================================
class BzzFactorizedStaircaseGauss : public BzzBaseClass
{
	//	============================================================================
	//	======================< Functions for linear algebra >======================
	//	============================================================================
	friend void Solve
	(BzzFactorizedStaircaseGauss* A, BzzVector* bx);
	friend void Solve
	(BzzFactorizedStaircaseGauss* A, BzzVector& bx, BzzVector* x);

	friend void TransposeSolve
	(BzzFactorizedStaircaseGauss* A, BzzVector* bx);
	friend void TransposeSolve
	(BzzFactorizedStaircaseGauss* A, BzzVector& bx, BzzVector* x);

protected:
	enum OdebFactorizationStatus
	{
		MATRIX_NON_AVAILABLE,
		UNFACTORIZED,
		FACTORIZED
	}factorizationStatus;

	static const char* const BZZ_ERROR;
	static const char* const BZZ_ERR_DESTROYED;
	static int count; // for whoAmI
	static int countInScope;

	int	numRows,
		numColumns,
		size;
	char pivoting;

	double** matrix;

	BzzVectorInt lowerBand, upperBand;
	char singular;	// singular = 1 when singular
	int whoAmI;

	void Factorization(void);
	void Solution(BzzVector* bx);
	void TransposeSolution(BzzVector* bx);

public:
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
		// default constructor
	BzzFactorizedStaircaseGauss(void);

	// copy constructor
//	BzzFactorizedStaircaseGauss(const BzzFactorizedStaircaseGauss &rval);

	BzzFactorizedStaircaseGauss(BzzMatrixStaircase& rval);

	BzzFactorizedStaircaseGauss(BzzVectorInt& low, BzzVectorInt& up);
	void operator()(BzzVectorInt& low, BzzVectorInt& up);
	// assigns and receives vector values with control
	double& operator()(int i, int j);
	// assigns and receives vector values without control
	double* operator [] (int i)
	{
		return matrix[i];
	}

	// FILE ASCII
//	BzzFactorizedStaircaseGauss(char *filematrix);

//	============================================================================
//	***************************< assignment operators >*************************
//	============================================================================
	BzzFactorizedStaircaseGauss& operator =
		(BzzMatrixStaircase& rval);
	// transforms a BzzMatrixStaircase in BzzFactorizedStaircaseGauss
	// BzzMatrixSparse is destroyed
	friend void ReplaceBzzMatrixWithBzzFactorized
	(BzzMatrixStaircase* lval, BzzFactorizedStaircaseGauss* rval);

	//	============================================================================
	//	********************************< destructor >******************************
	//	============================================================================
	~BzzFactorizedStaircaseGauss(void) { Delete(this); }
	friend void Delete(BzzFactorizedStaircaseGauss* rval);
	void RemoveAllElementsInMatrix(void);

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

	//	============================================================================
	//	========================< Non-modifying functions	=========================
	//	============================================================================
	BzzVectorInt GetLowerElementInRows(void);
	void GetLowerElementInRows(BzzVectorInt* low);
	BzzVectorInt GetUpperElementInRows(void);
	void GetUpperElementInRows(BzzVectorInt* up);

	//	*********************************< BzzPrint >*******************************
	virtual void ObjectBzzPrint(void);
	void BzzPrintStructure(void);
	int CountElements(void)
	{
		return size - 1;
	}
	//	double Condition(void);

	//	============================================================================
	//	==========================< Modifying functions	==========================
	//	============================================================================
	void SetToZeroAndReset(void);
	void NonPivoting(void)
	{
		if (factorizationStatus == UNFACTORIZED)
			pivoting = 0;
	}
	void Pivoting(void)
	{
		BzzError("TODO Pivoting in BzzFactorizedStaircaseGauss");
		if (factorizationStatus == UNFACTORIZED)
			pivoting = 1;
	}
};

/*
//	============================================================================
//	=====================< class BzzFactorizedOdeb >==========================
//	============================================================================
// Galerkin and Orthogonal collocation multi elements
// Example:
//	int numElements = 4;
//	int numInteriorPoints = 2;
//	int numVariables = 3;
//	numRows = numColumns = numVariables * numElements *
//		(numInteriorPoints + 1) + numVariables = 39;

//*xxxxxxxxxxx...........................
//x*xxxxxxxxxx...........................
//xx*xxxxxxxxx...........................
//xxx*xxxxxxxx...........................
//xxxx*xxxxxxx...........................
//xxxxx*xxxxxx...........................
//xxxxxx*xxxxx...........................
//xxxxxxx*xxxx...........................
//xxxxxxxx*xxx...........................
//xxxxxxxxx*xxxxxxxxxxx..................
//xxxxxxxxxx*xxxxxxxxxx..................
//xxxxxxxxxxx*xxxxxxxxx..................
//.........xxx*xxxxxxxx..................
//.........xxxx*xxxxxxx..................
//.........xxxxx*xxxxxx..................
//.........xxxxxx*xxxxx..................
//.........xxxxxxx*xxxx..................
//.........xxxxxxxx*xxx..................
//.........xxxxxxxxx*xxxxxxxxxxx.........
//.........xxxxxxxxxx*xxxxxxxxxx.........
//.........xxxxxxxxxxx*xxxxxxxxx.........
//..................xxx*xxxxxxxx.........
//..................xxxx*xxxxxxx.........
//..................xxxxx*xxxxxx.........
//..................xxxxxx*xxxxx.........
//..................xxxxxxx*xxxx.........
//..................xxxxxxxx*xxx.........
//..................xxxxxxxxx*xxxxxxxxxxx
//..................xxxxxxxxxx*xxxxxxxxxx
//..................xxxxxxxxxxx*xxxxxxxxx
//...........................xxx*xxxxxxxx
//...........................xxxx*xxxxxxx
//...........................xxxxx*xxxxxx
//...........................xxxxxx*xxxxx
//...........................xxxxxxx*xxxx
//...........................xxxxxxxx*xxx
//...........................xxxxxxxxx*xx
//...........................xxxxxxxxxx*x
//...........................xxxxxxxxxxx*

class BzzFactorizedOdeb : public BzzFactorizedStaircaseGauss
	{
//friend void Solve
//		 (BzzFactorizedOdeb *A,BzzVector *bx);
//friend void Solve
//		 (BzzFactorizedOdeb *A,BzzVector &bx,BzzVector *x);
private:
	int	numElements,
			numInteriorPoints,
			numVariables,
			numColumnsSubMatrix,
			numColumnsSubMatrixInLongRow;
public:
	BzzFactorizedOdeb(void) : BzzFactorizedStaircaseGauss(){};

	// copy constructor
//	BzzFactorizedOdeb(const BzzFactorizedOdeb &rval);

	//
	BzzFactorizedOdeb(int numE,int numI,int numV);
	void operator()(int numE,int numI,int numV);
	double &operator()(int i,int j);
	double &operator()(int e,int ie,int je,int iv,int jv);
	};

//	============================================================================
//	=================< class BzzFactorizedFiniteDifference >======================
//	============================================================================
// Generalized Finite Difference
// Example:
//	int numElements = 4;
//	int numInteriorPoints = 2;
//	int numVariables = 3;
// int numSupportPoints = (numElements * numInteriorPoints + 2) = 10;
//	numRows = numColumns = numSupportPoints * numVariables = 30;

class BzzFactorizedFiniteDifference : public BzzFactorizedStaircaseGauss
	{
private:
	int	numElements,
			numInteriorPoints,
			numVariables,
			numSupportPoints;
public:
	BzzFactorizedFiniteDifference(void) : BzzFactorizedStaircaseGauss(){};

	// copy constructor
//	BzzFactorizedFiniteDifference(const BzzFactorizedFiniteDifference &rval);

	BzzFactorizedFiniteDifference(int numE,int numI,int numV);
	void operator()(int numE,int numI,int numV);
	double &operator()(int i,int j);
	double &operator()(int e,int ie,int je,int iv,int jv);
	};
*/
#endif // BZZ_FACTORIZED_DOUBLE_STAIRCASE_GAUSS_HPP