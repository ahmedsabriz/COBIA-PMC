// BZZMATH: Release 7.0

//	=============================< BzzMatrixBlocks.hpp >===============================
//	* class BzzMatrixBlocks																*
//	* This class is used for matrices with block structure							*
// * Description:																					*
//	*					Metodi Numerici e Software in C++ (Capitolo 6, 7)				*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
// * Examples: c:\bzzmath\examples\dxmatblo.cpp											*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	05-1997	Date Written
// 09-1997	Added Product(c,&A);

//	============================================================================
//	******	constructor:																		*
//	* BzzMatrixBlocks A(numVariables,numDiagonalMatrices);						*
//	****************************************************************************
//	****** Insert SubMatrices																	*
//	* A.Insert(i,j,&aij);																		*
//	****************************************************************************
//	****** Functions:																				*
//	* Product(A,x,&y);																			*
//	* Product(4.,&A);																				*
//	* TProduct(A,x,&y);																			*
//	****************************************************************************
//	* Other functions:																			*
//	* A.BzzPrint("Comments");																	*
//	* A.BzzMessage("Comments");																*
//	****************************************************************************

#ifndef MATBLOCD_HPP
#define MATBLOCD_HPP

class ElementBzzMatrixBlocks
{
public:
	int column;
	BzzMatrix alk;
	ElementBzzMatrixBlocks* next;
};

//	============================================================================
//	=========================< class BzzMatrixBlocks >=====================
//	============================================================================

class BzzMatrix;
class BzzVector;
class BzzVectorInt;

class BzzMatrixBlocks : public BzzBaseClass
{
	friend class ElementBzzMatrixBlocks;
	friend class BzzFactorizedMatrixBlocksGauss;
	friend class BzzFactorizedLQMatrixBlocks;
private:

	static const char* const BZZ_ERROR;
	static int count; // for whoAmI

	int whoAmI;
	int	numDiagonalMatrices,
		numVariables;

	BzzVector b;
	BzzVectorInt nb;
	ElementBzzMatrixBlocks** elementRow;

	void Initialize(int nvar, int numMat);	// initialising constructors
	BzzMatrix* GetSubBzzMatrix(int l, int k);
	void DeleteSubBzzMatrix(int l, int k);
	void DeleteBzzMatrix(void);

public:
	//	============================================================================
	//	*******************************< constructors >*****************************
	//	============================================================================

	BzzMatrixBlocks(void);
	BzzMatrixBlocks(int numVar, int numMat);

	//	============================================================================
	//	******************************< destructor >********************************
	//	============================================================================
	~BzzMatrixBlocks(void);

	//	============================================================================
	//	================================< Functions >===============================
	//	============================================================================
	friend void Product(BzzMatrixBlocks& A, BzzVector& x,
		BzzVector* y);
	friend void TProduct(BzzMatrixBlocks& A, BzzVector& x,
		BzzVector* y);
	friend void Product(double c, BzzMatrixBlocks* A);

	// assigns and receives vector values with control
	double& operator () (int row, int col);

	virtual void ObjectBzzPrint(void);
	int CountSubMatrices(void);
	void Insert(int l, int k, BzzMatrix* alk);
	int GetSubMatrixExistence(int row, int column);
	void Reinitialize(int nvar, int numMat);
	BzzMatrixBlocks& operator = (BzzMatrixBlocks& m);
};

#endif // MATBLOCD_HPP