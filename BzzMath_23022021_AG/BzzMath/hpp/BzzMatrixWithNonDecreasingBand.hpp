// BZZMATH: Release 7.0
/*
Matrices with rectangular band structure
Example:
*xxx..........
x*xxxx........
xx*xxx........
..x*xxxxx.....
....*xxxxxx...
....x*xxxxx...
......*xxxxxx.
......x*xxxxx.
......xx*xxxx.
.......xx*xxx.
.......xxx*xx.
.........xx*x.
.........xxx*.
*/

//	==================< BzzMatrixWithNonDecreasingBand.hpp >====================
//	* BzzMatrixWithNonDecreasingBand: matrices with rectangle band structure	*
// * Examples: BzzMath\examples\BzzMathBasic\LinearAlgebra\							*
// *				MatrixWithNonDecreasingBand\MatrixWithNonDecreasingBand.cpp		*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	08-2002	Date Written.

//	============================================================================

#ifndef BZZ_MATRIX_DOUBLE_WITH_NON_DECREASING_BAND_HPP
#define BZZ_MATRIX_DOUBLE_WITH_NON_DECREASING_BAND_HPP

//	============================================================================
//	=====================< class BzzMatrixWithNonDecreasingBand >==========================
//	============================================================================
class BzzMatrixWithNonDecreasingBand : public BzzBaseClass
{
	friend class BzzFactorizedWithNonDecreasingBandGauss;
	friend void Product(BzzMatrixWithNonDecreasingBand& A, BzzVector& x, BzzVector* y);
private:
	static const char* const BZZ_ERROR;
	static const char* const BZZ_ERR_DESTROYED;
	static int count; // for whoAmI
	static int countInScope;

	int	numRows,
		numColumns,
		size;
	int	numSupportPoints,
		numElements,
		numInternalPoints,
		numComponents;

	double** matrix;

	BzzVectorInt lowerBand, upperBand;
	char singular;	// singular = 1 when singular
	int whoAmI;

public:
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
		// default constructor
	BzzMatrixWithNonDecreasingBand(void);

	// copy constructor
//	BzzMatrixWithNonDecreasingBand(const BzzMatrixWithNonDecreasingBand &rval);

	//
	BzzMatrixWithNonDecreasingBand(BzzVectorInt& low, BzzVectorInt& up);
	void operator()(BzzVectorInt& low, BzzVectorInt& up);
	// assigns and receives vector values with control
	double& operator()(int i, int j);
	// assigns and receives vector values without control
	double* operator [] (int i)
	{
		return matrix[i];
	}
	// FILE ASCII
//	BzzMatrixWithNonDecreasingBand(char *filematrix);

//	============================================================================
//	********************************< destructor >******************************
//	============================================================================
	~BzzMatrixWithNonDecreasingBand(void);
	friend void Delete(BzzMatrixWithNonDecreasingBand* rval);

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
	void RemoveAllElementsInMatrix(void);
	void SetMatrixForFiniteDifference
	(int numE, int numI, int numC, BzzVectorInt* low, BzzVectorInt* up);
	void SetMatrixForFiniteElements
	(int numE, int numI, int numC, BzzVectorInt* low, BzzVectorInt* up);
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

class BzzFactorizedOdeb : public BzzMatrixWithNonDecreasingBand
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
	BzzFactorizedOdeb(void) : BzzMatrixWithNonDecreasingBand(){};

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

class BzzFactorizedFiniteDifference : public BzzMatrixWithNonDecreasingBand
	{
private:
	int	numElements,
			numInteriorPoints,
			numVariables,
			numSupportPoints;
public:
	BzzFactorizedFiniteDifference(void) : BzzMatrixWithNonDecreasingBand(){};

	// copy constructor
//	BzzFactorizedFiniteDifference(const BzzFactorizedFiniteDifference &rval);

	BzzFactorizedFiniteDifference(int numE,int numI,int numV);
	void operator()(int numE,int numI,int numV);
	double &operator()(int i,int j);
	double &operator()(int e,int ie,int je,int iv,int jv);
	};
*/
#endif // BZZ_MATRIX_DOUBLE_WITHNONDECREASING_BAND_HPP