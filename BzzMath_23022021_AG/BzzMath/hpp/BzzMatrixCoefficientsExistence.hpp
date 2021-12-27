// BZZMATH: Release 7.0

//	=================< MatrixCoefficientsExistence.hpp >========================
//	* BzzMatrixCoefficientsExistence Class													*
//	* for giving the non null coefficients in a sparse matrix						*
// * Description:																					*
//	*					Metodi Numerici e Software in C++ (Capitolo 3, 26)				*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
//	*					Numerical Methods and Software in C++ (Chapter 3)				*
//	*					by G. Buzzi-Ferraris														*
//	*					Addison Wesley Longman (1999)											*
// * Examples: c:\bzzmath\examples\exexist.cpp											*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	07-1994	Date Written
//	01-1995	Modified FindDependence.
//	02-1995	Added the Bzz prefix to the names of the classes.
//	05-1996	Modified FindDependence.
//	03-1997	Added Analysis and PrintAnalysis.
//	03-1997	Added GaussFillingUp.
//	03-1997	MatrixAfterGaussFactorization
//	04-1997	Added OrderingQT.
//	03-1998	Added MatrixAfterLQFactorization.
//	03-1998	Added OrderingLQ.
//	04-1999	Added SetDiagonal.
//	10-1999	Added GetRowsExistence.

/////////// Release 4.0
//	01-2000	Added BeginScanning.
// 10-2000	Added GetStartingElementInRow.
//	11-2000	Added GetNumEquationsForEachVariable.
//	11-2000	Added GetNumVariablesForEachEquation.
//	11-2000	Added GetNumVariablesForEachEquationAndNumEquationsForEachVariable.

/////////// Release 5.0
//	09-2003	Added Transpose.
//	09-2003	Added MatrixAfterQRFactorization.
//	09-2003	Added CountElementsOutsideLowerBand.
//	09-2003	Added CountElementsOutsideUpperBand.
//	09-2003	Added CountElementsOutsideBands.
//	11-2003	Added GetLowerAndUpperElementInRows.
//	11-2003	Added GetLowerAndUpperElementInColumns.
//	01-2005	Added costructor from BzzMatrixSparse.
//	01-2005	Added costructor from BzzMatrixSparse.
//	01-2005	Added assignment from BzzMatrixSparse.
//	01-2005	Added assignment from BzzMatrixSparse.
//	03-2005	Added SwapRows and SwapColumns functions.
//	03-2005	Added SwapElements function.
//	03-2005	Added Transpose function.
//	12-2006	Modified Analysis and PrintAnalysis functions.

//	============================================================================
//	******	BzzMatrixCoefficientsExistence constructors:								*
//	* BzzMatrixCoefficientsExistence A; // default										*
//	* BzzMatrixCoefficientsExistence A = B; // copy-initializer						*
//	* BzzMatrixCoefficientsExistence A(m,n); // sizes									*
//	* BzzMatrixCoefficientsExistence A("EXIST.DAT"); //	file						*
//	* BzzMatrixCoefficientsExistence A = B; // from BzzMatrixSparse				*
//	* BzzMatrixCoefficientsExistence A = B; // from BzzMatrixDoubeSparse			*
//	****************************************************************************
//	***** Access functions:																		*
//	* i = A.Rows(); // numRows																	*
//	* i = A.Columns(); // numColumns															*
//	* who = A.WhoAmI();																			*
//	* A(i,j);																						*
//	* char e = A.Existing(i,j);																*
//	* char scan = A.Scanning(&i,&j);															*
//	* A.BeginScanning();																			*
//	* int countElements = A.CountElements();												*
//	* int countElementsOutsideLowerBand = A.CountElementsOutsideLowerBand(10);	*
//	* int countElementsOutsideUpperrBand =													*
//	*			A.CountElementsOutsideUpperrBand(15);										*
//	* A.CountElementsOutsideBands(10,15,&low,&up);										*
//	* BzzVectorInt numVar = A.CountElementsInRows();									*
//	* BzzVectorInt numEq = A.CountEquationsForEachVaraible();						*
//	* int countObject = BzzMatrixCoefficientsExistence::ObjectCount();			*
//	* int countInScope =	BzzMatrixCoefficientsExistence::ObjectCountInScope();	*
//	****************************************************************************
//	***** Assignments:																			*
//	* A = B;	// BzzMatrixCoefficientsExistence											*
//	* A = B; // from BzzMatrixSparse															*
//	* A = B; // from BzzMatrixDoubeSparse													*
//	****************************************************************************
//	***** BzzMessage and BzzPrint																*
//	* A.BzzPrint("Comment"); // always														*
//	* A.BzzMessage("Comment"); // only if bzzMessageYesNo == 'Y'					*
//	* A.PrintAnalysis();																			*
//	****************************************************************************
//	***** Save and Load																			*
//	* A.Save("exist.dat");	// formatted													*
//	* Load(&A,"exist.dat");																		*
//	****************************************************************************
//	***** Delete, Clean, ChangeDimensions and Swap										*
//	* A.RemoveElement(3,151);																	*
//	* A.RemoveAllElementsInRow(row);															*
//	* A.RemoveAllElementsInColumn(column);													*
//	* Delete(&A);																					*
//	* Transpose(&A);																				*
//	* A.RemoveAllElementsInMatrix();															*
//	* ChangeDimensions(newr,newc,&A);														*
//	* Swap(&A,&B);																					*
//	* A.SwapRows(123,435);																		*
//	* A.SwapColumns(215,7);																		*
//	* A.SwapElements(215,7,415,212);															*
//	****************************************************************************
//	***** SetDiagonal																				*
//	* A.SetDiagonal(0);	// Principal														*
//	* A.SetDiagonal(-5);	// Lower																*
//	* A.SetDiagonal(8);	// Upper																*
//	****************************************************************************
//	***** Functions for non linear systems													*
//	* see Chapter 25 Numerical Software in C++											*
//	* A.FindDependence();																		*
//	* A.PrintDependence();																		*
//	****************************************************************************

#ifndef BZZ_EXIST_HPP
#define BZZ_EXIST_HPP

//	============================================================================
//	==================< class BzzMatrixCoefficientsExistence >==================
//	============================================================================

struct ElementBzzMatrixCoefficientsExistence
{
	int column;
	ElementBzzMatrixCoefficientsExistence* next;
};

class BzzMatrixCoefficientsExistence : public BzzBaseClass
{
	friend class BzzOdeSparseStiff;
private:
	static const char* const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;

	ElementBzzMatrixCoefficientsExistence** elementRow;
	ElementBzzMatrixCoefficientsExistence* elemScanning;

	int rowScanning;

	int	numRows,
		numColumns,
		numElements,
		maxElementsInRow, iMaxRow,
		minElementsInRow, iMinRow,
		maxElementsInColumn, iMaxColumn,
		minElementsInColumn, iMinColumn,
		numElementsNullOnDiagonal,
		whoAmI,
		* ptrDer;

	int	lowerBand,
		upperBand;

public:
	int** dependence,
		* numEqWithSuchVariable,
		** variableInGroup,
		numGroup,
		* numVariablesInGroup;

	int	numElementsUpPrincipalDiagonal,
		numElementsDownPrincipalDiagonal;

private:
	long int numTheoreticalElements;

	char	dependenceAvailable;

	void InsertElement
	(ElementBzzMatrixCoefficientsExistence* elem,
		int row, int column, int first);
	void Initialize(int rows, int columns);
	void Copy(const BzzMatrixCoefficientsExistence& rval);
	void DeleteDependence(void);

public:
	//	============================================================================
	//	******************************< constructors >******************************
	//	============================================================================
		// default constructor
	BzzMatrixCoefficientsExistence(void);

	// copy constructor
	BzzMatrixCoefficientsExistence
	(const BzzMatrixCoefficientsExistence& rval);

	// sizing constructor
	BzzMatrixCoefficientsExistence(int rows, int columns);

	// initialises from formatted File
	// BzzMatrixCoefficientsExistence A("EXIST.DAT");
	BzzMatrixCoefficientsExistence(char* filematrix);
	BzzMatrixCoefficientsExistence(char* filematrix, float val);
	BzzMatrixCoefficientsExistence(char* filematrix, double val);

	// from BzzMatrixSparse
	BzzMatrixCoefficientsExistence(BzzMatrixSparse& rval);

	//	============================================================================
	//	***************************< destructor >***********************************
	//	============================================================================
	~BzzMatrixCoefficientsExistence(void);

	//	============================================================================
	//	*************************< Access functions >*******************************
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
	static int ObjectCount(void) { return count; }
	static int ObjectCountInScope(void) { return countInScope; }

	// assigns existence
	void operator ()
		(int row, int column);

	char Existing(int i, int j);
	int CountElements(void);
	int CountElementsOutsideLowerBand(int lowerB);
	int CountElementsOutsideUpperBand(int upB);
	void CountElementsOutsideBands(int lowerB, int upB, int* low, int* up);
	void CountElementsOutsideBands(int* l, int* u, BzzVectorInt* low, BzzVectorInt* up);
	BzzVectorInt CountElementsInRows(void);
	BzzVectorInt CountEquationsForEachVariable(void);
	char Scanning(int* i, int* j);
	ElementBzzMatrixCoefficientsExistence* GetStartingElementInRow(int i);
	void BeginScanning(void);
	void GetRowsExistence(BzzVectorInt* iRows);
	void GetNumEquationsForEachVariable
	(BzzVectorInt* numEquationsForEachVariable);
	void GetNumVariablesForEachEquation
	(BzzVectorInt* numVariablesForEachEquation);
	void GetNumVariablesForEachEquationAndNumEquationsForEachVariable
	(BzzVectorInt* numVariablesForEachEquation,
		BzzVectorInt* numEquationsForEachVariable);

	//	============================================================================
	//	*************************< assignment operators >***************************
	//	============================================================================
	BzzMatrixCoefficientsExistence& operator =
		(const BzzMatrixCoefficientsExistence& rval);
	BzzMatrixCoefficientsExistence& operator =(BzzMatrixSparse& rval);

	//	============================================================================
	//	========================< Non-modifying functions >=========================
	//	============================================================================

	//	*********************************< BzzPrint >*******************************
	virtual void ObjectBzzPrint(void);
	void	Save(char* filematrix); // formatted

	void FindDependence(void);
	void PrintDependence(void);

	void Analysis(void);
	void PrintAnalysis(void);

	void FindBands(int* low, int* up);
	void GetLowerAndUpperElementInRows(BzzVectorInt* low, BzzVectorInt* up);
	void GetLowerAndUpperElementInColumns(BzzVectorInt* low, BzzVectorInt* up);
	int GaussFillingUp(void);
	int LQFillingUp(void);
	void MatrixAfterGaussFactorization(void);
	void MatrixAfterLQFactorization(void);
	void MatrixAfterQRFactorization(void);
	void OrderingQT(BzzVectorInt* ordRows, BzzVectorInt* ordColumns);
	void OrderingLQ(BzzVectorInt* ordRows, BzzVectorInt* ordColumns);

	//	============================================================================
	//	==========================< Modifying Functions >===========================
	//	============================================================================
	void SetDiagonal(int i);
	void RemoveElement(int row, int column);
	void RemoveAllElementsInRow(int row);
	void RemoveAllElementsInColumn(int column);
	friend void Delete(BzzMatrixCoefficientsExistence* A);
	friend void Transpose(BzzMatrixCoefficientsExistence* A);
	void RemoveAllElementsInMatrix(void);// eliminates all elements
	friend void ChangeDimensions(int nr, int nc,
		BzzMatrixCoefficientsExistence* m);
	friend void Swap(BzzMatrixCoefficientsExistence* lval,
		BzzMatrixCoefficientsExistence* rval);
	void SwapRows(int i, int j);
	void SwapColumns(int i, int j);
	void SwapElements(int i1, int j1, int i2, int j2);
	friend void Load(BzzMatrixCoefficientsExistence* A, char* filevector);
	friend void Transpose(BzzMatrixCoefficientsExistence* A);
};

#endif // BZZ_EXIST_HPP