// BZZMATH: Release 7.0

//	===========================< EXISTSYM.HPP >=================================
//	* BzzMatrixSymmetricCoefficientsExistence Class										*
//	* for giving the non null coefficients in a sparse symmetric matrix			*
// * Description:																					*
//	*					Metodi Numerici e Software in C++ (Capitolo 3)					*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
// * Examples: c:\bzzmath\examples\exexsym.cpp											*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	06-1996	Date Written

/////////// Release 4.0
//	01-2000	Added BeginScanning.

//	============================================================================
//	******	BzzMatrixSymmetricCoefficientsExistence constructors:					*
//	* BzzMatrixSymmetricCoefficientsExistence A; // default							*
//	* BzzMatrixSymmetricCoefficientsExistence A = B; // copy-initializer			*
//	* BzzMatrixSymmetricCoefficientsExistence A(n,n); // sizes						*
//	* BzzMatrixSymmetricCoefficientsExistence A("EXISTSYM.DAT"); //	file		*
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
//	* int countObject =																			*
//	* BzzMatrixSymmetricCoefficientsExistence::ObjectCount();						*
//	* int countInScope =																			*
//	* BzzMatrixSymmetricCoefficientsExistence::ObjectCountInScope();				*
//	****************************************************************************
//	***** Assignments:																			*
//	* A = B;	// BzzMatrixSymmetricCoefficientsExistence								*
//	****************************************************************************
//	***** BzzMessage and BzzPrint																*
//	* A.BzzPrint("Comment"); // always														*
//	* A.BzzMessage("Comment"); // only if bzzMessageYesNo == 'Y'					*
//	****************************************************************************
//	***** Save and Load																			*
//	* A.Save("existsym.dat");	// formatted												*
//	* Load(&A,"existsym.dat");																	*
//	****************************************************************************
//	***** Delete, Clean, ChangeDimensions and Swap										*
//	* A.RemoveElement(3,151);																	*
//	* Delete(&A);																					*
//	* A.RemoveAllElementsInMatrix();															*
//	* ChangeDimensions(newr,newc,&A);														*
//	* Swap(&A,&B);																					*
//	****************************************************************************

#ifndef BZZ_EXIST_SYMM_HPP
#define BZZ_EXIST_SYMM_HPP

//	============================================================================
//	============< class BzzMatrixSymmetricCoefficientsExistence >===============
//	============================================================================
class BzzMatrixSymmetricCoefficientsExistence : public BzzBaseClass
{
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
		maxElementsInRow,
		minElementsInRow,
		numElementsNullOnDiagonal,
		whoAmI,
		* ptrDer;

	int	symBand;
	void InsertElement
	(ElementBzzMatrixCoefficientsExistence* elem,
		int row, int column, int first);
	void Initialize(int rows, int columns);
	void Copy(const BzzMatrixSymmetricCoefficientsExistence& rval);
	void DeleteDependence(void);

public:
	//	============================================================================
	//	******************************< constructors >******************************
	//	============================================================================
		// default constructor
	BzzMatrixSymmetricCoefficientsExistence(void);

	// copy constructor
	BzzMatrixSymmetricCoefficientsExistence
	(const BzzMatrixSymmetricCoefficientsExistence& rval);

	// sizing constructor
	BzzMatrixSymmetricCoefficientsExistence(int rows, int columns);

	// initialises from formatted File
	// BzzMatrixSymmetricCoefficientsExistence A("EXIST.DAT");
	BzzMatrixSymmetricCoefficientsExistence(char* filematrix);
	BzzMatrixSymmetricCoefficientsExistence(char* filematrix, float val);
	BzzMatrixSymmetricCoefficientsExistence(char* filematrix, double val);

	//	============================================================================
	//	***************************< destructor >***********************************
	//	============================================================================
	~BzzMatrixSymmetricCoefficientsExistence(void);

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
	BzzVectorInt CountElementsInRows(void);
	char Scanning(int* i, int* j);
	void BeginScanning(void);

	//	============================================================================
	//	*************************< assignment operators >***************************
	//	============================================================================
	BzzMatrixSymmetricCoefficientsExistence& operator =
		(const BzzMatrixSymmetricCoefficientsExistence& rval);

	//	============================================================================
	//	========================< Non-modifying functions >=========================
	//	============================================================================

	//	*********************************< BzzPrint >*******************************
	virtual void ObjectBzzPrint(void);
	void	Save(char* filematrix); // formatted

	BzzVectorInt Ordering(void);
	void ReOrdering(BzzVectorInt& newOrder);

	void FindBand(int* band);
	void FindBand(int* band, int* envelope, BzzVectorInt* beta);

	//	============================================================================
	//	==========================< Modifying Functions >===========================
	//	============================================================================
	void RemoveElement(int row, int column);
	friend void Delete(BzzMatrixSymmetricCoefficientsExistence* A);
	void RemoveAllElementsInMatrix(void);// eliminates all elements
	friend void ChangeDimensions(int nr, int nc,
		BzzMatrixSymmetricCoefficientsExistence* m);
	friend void Swap(BzzMatrixSymmetricCoefficientsExistence* lval,
		BzzMatrixSymmetricCoefficientsExistence* rval);
	friend void Load(BzzMatrixSymmetricCoefficientsExistence* A, char* filevector);
};

#endif // BZZ_EXIST_SYMM_HPP