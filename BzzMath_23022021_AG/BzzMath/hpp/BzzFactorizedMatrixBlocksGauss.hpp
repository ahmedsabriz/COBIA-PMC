// BZZMATH: Release 3.1

//	================< BzzFactorizedMatrixBlocksGauss.hpp >======================
//	* class BzzFactorizedMatrixBlocksGauss													*
//	* This class is used for solving linear systems with block structure			*
// * Description:																					*
//	*					Metodi Numerici e Software in C++ (Capitolo 6,7)				*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
// * Examples: c:\bzzmath\examples\exfatblo.cpp											*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	03-1995	Date Written
//	04-1997	Added Reinitialize.
//	05-1997	Added ConditionNumber.

//	============================================================================
//	******	constructor:																		*
//	* BzzFactorizedMatrixBlocksGauss A(numVariables,numDiagonalMatrices);			*
//	****************************************************************************
//	****** Insert SubMatrices																	*
//	* A.Insert(i,j,&aij);																		*
//	****************************************************************************
//	****** Functions for solving systems:													*
//	* Solve(&A,&bx);																				*
//	* A.SolveAndDEstroy(&A,&bx);																*
//	****************************************************************************
//	* Other functions:																			*
//	* A.BzzPrint("Comments");																	*
//	* A.BzzMessage("Comments");																*
//	****************************************************************************

#ifndef FACTBLOCD_HPP
#define FACTBLOCD_HPP

class ElementBzzFactorizedMatrixBlocksGauss
{
	friend class BzzFactorizedMatrixBlocksGauss;
	friend void ProductForDae(double ch, double c,
		BzzVectorInt& id, BzzFactorizedMatrixBlocksGauss* A);

private:
	int column;
	BzzMatrix alk;
	ElementBzzFactorizedMatrixBlocksGauss* next;
};

//	============================================================================
//	=================< class BzzFactorizedMatrixBlocksGauss >=====================
//	============================================================================

class BzzMatrix;
class BzzVector;
class BzzVectorInt;
class BzzFactorized;
class BzzFactorizedPLR;
class BzzFactorizedGauss;

class BzzFactorizedMatrixBlocksGauss : public BzzBaseClass
{
	friend class ElementBzzFactorizedMatrixBlocksGauss;
private:
	enum BzzFactorizedMatrixBlocksGaussStatus
	{
		AVAILABLE,
		DESTROYED,
		FACTORIZED
	}matrixBlockStatus;

	static const char* const BZZ_ERROR;
	static int count; // for whoAmI

	int whoAmI;
	int	numDiagonalMatrices,
		numVariables;

	BzzVector b;
	BzzVectorInt nb;
	BzzFactorizedGauss* All;
	ElementBzzFactorizedMatrixBlocksGauss** elementRow;

	void Initialize(int nvar, int numMat);	// initialising constructors
	BzzMatrix* GetSubBzzMatrix(int l, int k);
	void DeleteSubBzzMatrix(int l, int k);
	void DeleteBzzMatrix(void);

public:
	//	============================================================================
	//	*******************************< constructors >*****************************
	//	============================================================================

	BzzFactorizedMatrixBlocksGauss(void);
	BzzFactorizedMatrixBlocksGauss(int numVar, int numMat);

	//	============================================================================
	//	******************************< destructor >********************************
	//	============================================================================
	~BzzFactorizedMatrixBlocksGauss(void);

	//	============================================================================
	//	================================< Functions >===============================
	//	============================================================================
	friend void SolveAndDestroy(BzzFactorizedMatrixBlocksGauss* ab,
		BzzVector* bx);
	friend void Solve(BzzFactorizedMatrixBlocksGauss* ab, BzzVector* bx);
	double ConditionNumber(void);

	virtual void ObjectBzzPrint(void);
	int CountSubMatrices(void);
	void Insert(int l, int k, BzzMatrix* alk);
	int GetSubMatrixExistence(int row, int column);
	void Reinitialize(int nvar, int numMat);
	BzzFactorizedMatrixBlocksGauss& operator = (BzzMatrixBlocks& m);
	void ReplaceForOde(double c, BzzMatrixBlocks& rval);
	friend void ProductForDae(double ch, double c,
		BzzVectorInt& id, BzzFactorizedMatrixBlocksGauss* A);
};

#endif // FACTBLOCD_HPP