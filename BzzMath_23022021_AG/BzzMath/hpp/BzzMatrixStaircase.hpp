// BZZMATH: Release 7.0
/*
Linear systems with Staircase structure
  Example:
*xxx............
x*xxxxxx........
xx*xxxxx........
xxx*xxxx........
xxxx*xxx........
....x*xx........
....xx*xxxxx....
....xxx*xxxx....
....xxxx*xxx....
........x*xxxxxx
........xx*xxxxx
........xxx*xxxx
........xxxx*xxx
............x*xx
............xx*x
*/

//	====================< BzzMatrixStaircase.hpp >==============================
//	* BzzMatrixStaircase: matrices with rectangle band structure					*
// * Examples: BzzMath\examples\BzzMathBasic\LinearAlgebra\							*
// *				MatrixStaircase\MatrixStaircase.cpp										*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	05-2005	Date Written.
// 02-2012	Added BzzMatrixVariableBand class
//	============================================================================

#ifndef BZZ_MATRIX_DOUBLE_STAIRCASE_HPP
#define BZZ_MATRIX_DOUBLE_STAIRCASE_HPP

//	============================================================================
//	======================< class BzzMatrixStaircase >==========================
//	============================================================================
class BzzMatrixStaircase : public BzzBaseClass
{
	friend class BzzFactorizedStaircaseGauss;
	friend void Product(BzzMatrixStaircase& A, BzzVector& x, BzzVector* y);
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
	BzzVectorInt lowerColumnsBand, upperColumnsBand;
	char singular;	// singular = 1 when singular
	int whoAmI;

public:
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
		// default constructor
	BzzMatrixStaircase(void);

	// copy constructor
//	BzzMatrixStaircase(const BzzMatrixStaircase &rval);

	//
	BzzMatrixStaircase(BzzVectorInt& low, BzzVectorInt& up);
	void operator()(BzzVectorInt& low, BzzVectorInt& up);
	// assigns and receives vector values with control
	double& operator()(int i, int j);
	// assigns and receives vector values without control
	double* operator [] (int i)
	{
		return matrix[i];
	}
	void SetToZeroAndReset(void);
	// FILE ASCII
//	BzzMatrixStaircase(char *filematrix);

//	============================================================================
//	********************************< destructor >******************************
//	============================================================================
	~BzzMatrixStaircase(void);
	friend void Delete(BzzMatrixStaircase* rval);

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
	int GetLowerElementInRow(int row) { if (row > 0 && row <= numRows)return lowerBand[row]; }
	int GetUpperElementInRow(int row) { if (row > 0 && row <= numRows)return upperBand[row]; }

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
	void RemoveAllElementsInMatrix(void);
};

// BZZMATH: Release 7.0
/*
Linear systems with Staircase structure
  Example:
*xxx............
x*xxxxxx........
xx*xxxxx........
xxx*xxxx........
xxxx*xxx........
....x*xx........
....xx*xxxxx....
....xxx*xxxx....
....xxxx*xxx....
........x*xxxxxx
........xx*xxxxx
........xxx*xxxx
........xxxx*xxx
............x*xx
............xx*x
*/

//	====================< BzzMatrixVariableBand.hpp >==============================
//	* BzzMatrixVariableBand: matrices with rectangle band structure					*
// * Examples: BzzMath\examples\BzzMathBasic\LinearAlgebra\							*
// *				MatrixVariableBand\MatrixVariableBand.cpp										*
//	============================================================================

//	============================================================================
//	======================< class BzzMatrixVariableBand >==========================
//	============================================================================
class BzzMatrixVariableBand : public BzzBaseClass
{
	friend class BzzFactorizedStaircaseGauss;
	friend void Product(BzzMatrixVariableBand& A, BzzVector& x, BzzVector* y);
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
	BzzVectorInt lowerColumnsBand, upperColumnsBand;
	char singular;	// singular = 1 when singular
	int whoAmI;

public:
	//	============================================================================
	//	*****************************< constructors >*******************************
	//	============================================================================
		// default constructor
	BzzMatrixVariableBand(void);

	// copy constructor
//	BzzMatrixVariableBand(const BzzMatrixVariableBand &rval);

	//
	BzzMatrixVariableBand(BzzVectorInt& low, BzzVectorInt& up);
	void operator()(BzzVectorInt& low, BzzVectorInt& up);
	// assigns and receives vector values with control
	double& operator()(int i, int j);
	// assigns and receives vector values without control
	double* operator [] (int i)
	{
		return matrix[i];
	}
	void SetToZeroAndReset(void);
	// FILE ASCII
//	BzzMatrixVariableBand(char *filematrix);

//	============================================================================
//	********************************< destructor >******************************
//	============================================================================
	~BzzMatrixVariableBand(void);
	friend void Delete(BzzMatrixVariableBand* rval);

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
	int GetLowerElementInRow(int row) { if (row > 0 && row <= numRows)return lowerBand[row]; }
	int GetUpperElementInRow(int row) { if (row > 0 && row <= numRows)return upperBand[row]; }

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
	void RemoveAllElementsInMatrix(void);
};

#endif // BZZ_MATRIX_DOUBLE_STAIRCASE_HPP