// BZZMATH: Release 7.0

//	====================< BzzMatricesExistence.hpp >============================
//	* BzzMatricesExistence Class																*
//	* for giving the non null matrices in a block matrix								*
// * Description:																					*
//	*					Metodi Numerici e Software in C++ (Capitolo 26)					*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
// * Examples: c:\bzzmath\examples\exexblo.cpp											*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	05-1997	Date Written

/////////// Release 4.0
//	01-2000	Added BeginScanning.

/////////// Release 5.0

//	============================================================================
//	******	BzzMatricesExistence constructors:											*
//	* BzzMatricesExistence A; // default													*
//	* BzzMatricesExistence A = B; // copy-initializer									*
//	* BzzMatricesExistence A(m); // BzzVectorInt m(n);									*
//	* BzzMatricesExistence A("EXISTBLO.DAT"); //	file									*
//	****************************************************************************
//	***** Access functions:																		*
//	* i = A.Rows(); // num matrices in row													*
//	* i = A.Columns(); // num matrices in column											*
//	* who = A.WhoAmI();																			*
//	* A(i,j);																						*
//	* char e = A.Existing(i,j);																*
//	* char scan = A.Scanning(&i,&j);															*
//	* A.BeginScanning();																			*
//	* int countElements = A.CountElements();												*
//	* BzzVectorInt numVar = A.CountElementsInRows();									*
//	* BzzVectorInt numEq = A.CountElementsInColumns();									*
//	* int countObject = BzzMatricesExistence::ObjectCount();							*
//	* int countInScope =	BzzMatricesExistence::ObjectCountInScope();				*
//	****************************************************************************
//	***** Assignments:																			*
//	* A = B;	// BzzMatricesExistence															*
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
//	* A.RemoveAllElementsInMatrix();															*
//	* ChangeDimensions(newr,newc,&A);														*
//	* Swap(&A,&B);																					*
//	****************************************************************************
//	***** Functions for non linear systems													*
//	* see Chapter 25 Numerical Software in C++											*
//	* A.FindDependence();																		*
//	* A.PrintDependence();																		*
//	****************************************************************************

#ifndef BZZ_EXISTBLO_HPP
#define BZZ_EXISTBLO_HPP

//	============================================================================
//	=======================< class BzzMatricesExistence >=======================
//	============================================================================

struct ElementBzzMatricesExistence
{
	int column;
	ElementBzzMatricesExistence* next;
};

class BzzMatricesExistence : public BzzBaseClass
{
private:
	static const char* const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;

	ElementBzzMatricesExistence** elementRow;
	ElementBzzMatricesExistence* elemScanning;

	int rowScanning;

	int	numRows,
		numColumns,
		numVariables,
		numElements,
		maxElementsInRow,
		minElementsInRow,
		numElementsNullOnDiagonal,
		whoAmI,
		* ptrDer;

	int	lowerBand,
		upperBand,
		lowerMeanBand,
		upperMeanBand;

	///////////////////////////////////////////////////////////////////////////////
public:
	int** dependence,
		* numEqWithSuchVariable,
		** variableInGroup,
		numGroup,
		* numVariablesInGroup;

	BzzVectorInt diagonal;

	///////////////////////////////////////////////////////////////////////////////
private:
	long int numTheoreticalElements;

	char	dependenceAvailable;

	void InsertElement
	(ElementBzzMatricesExistence* elem,
		int row, int column, int first);
	void Initialize(int rows, int columns);
	void Copy(const BzzMatricesExistence& rval);
	void DeleteDependence(void);

public:
	//	============================================================================
	//	******************************< constructors >******************************
	//	============================================================================
		// default constructor
	BzzMatricesExistence(void);

	// copy constructor
	BzzMatricesExistence
	(const BzzMatricesExistence& rval);

	// sizing constructor
	BzzMatricesExistence(BzzVectorInt& diag);

	//	============================================================================
	//	***************************< destructor >***********************************
	//	============================================================================
	~BzzMatricesExistence(void);

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

	// change dimensions
	void operator ()
		(BzzVectorInt& diag);

	char Existing(int i, int j);
	int CountElements(void);
	BzzVectorInt CountElementsInRows(void);
	BzzVectorInt CountElementsInColumns(void);
	char Scanning(int* i, int* j);
	void BeginScanning(void);
	ElementBzzMatricesExistence* GetStartingElementInRow(int i);

	//	============================================================================
	//	*************************< assignment operators >***************************
	//	============================================================================
	BzzMatricesExistence& operator =
		(const BzzMatricesExistence& rval);

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
	int GaussFillingUp(void);
	void MatrixAfterGaussFactorization(void);
	void OrderingQT(BzzVectorInt* ordRows, BzzVectorInt* ordColumns);
	void OrderingGauss(BzzMatricesExistence* A,
		BzzVectorInt* ordRows, BzzVectorInt* ordColumns);

	//	============================================================================
	//	==========================< Modifying Functions >===========================
	//	============================================================================
	void RemoveElement(int row, int column);
	void RemoveAllElementsInRow(int row);
	void RemoveAllElementsInColumn(int column);
	friend void Delete(BzzMatricesExistence* A);
	void RemoveAllElementsInMatrix(void);// eliminates all elements
	friend void ChangeDimensions(BzzVectorInt& diag,
		BzzMatricesExistence* m);
	friend void Swap(BzzMatricesExistence* lval,
		BzzMatricesExistence* rval);
	friend void Load(BzzMatricesExistence* A, char* filevector);
};

#endif // BZZ_EXISTBLO_HPP