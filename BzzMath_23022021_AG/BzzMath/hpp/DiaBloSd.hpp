// BZZMATH: Release 3.1

//	==============================< DIABLOSD.HPP >==============================
//	* class BzzFactorizedDiagonalBlocksAndSparse: 									*
//	* This class is used for solving linear systems with block diagonal 			*
//	* structure and other sparse element													*
// * Description:																					*
//	*					Metodi Numerici e Software in C++ (Capitolo 6, 7)				*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
// * Examples: c:\bzzmath\examples\exdiablo.cpp											*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	10-1995	Date Written

//	============================================================================
//	******	constructor:																		*
//	* BzzFactorizedDiagonalBlocksAndSparse A(numVariables,numDiagonalMatrices);	*
//	****************************************************************************
//	****** Insert a double																		*
//	****** double f;																				*
//	****** BzzMatrixSparse S(1,numVariables);												*
//	****** Insert Diagonal double and Sparse SubMatrix									*
//	* A.Insert(j,f,&S);																			*
//	****************************************************************************
//	****** Insert a BzzMatrix																	*
//	****** BzzMatrix Z(n,n);																	*
//	****** BzzMatrixSparse S(n,numVariables);												*
//	****** Insert Diagonal SubMatrix and Sparse SubMatrix								*
//	* A.Insert(i,&Z,&S);																			*
//	****************************************************************************
//	****** Insert a BzzMatrixLeft																*
//	****** BzzMatrix L(n,n);																	*
//	****** BzzMatrixSparse S(n,numVariables);												*
//	****** Insert Diagonal SubMatrixLeft and Sparse SubMatrix						*
//	* A.Insert(i,&L,&S);																			*
//	****************************************************************************
//	****** Insert a BzzMatrixRight															*
//	****** BzzMatrixRight R(n,n);																*
//	****** BzzMatrixSparse S(n,numVariables);												*
//	****** Insert Diagonal SubMatrixRight and Sparse SubMatrix						*
//	* A.Insert(i,&R,&S);																			*
//	****************************************************************************
//	****** Insert a BzzMatrixSymmetric														*
//	****** BzzMatrixSymmetric W(n,n);														*
//	****** BzzMatrixSparse S(n,numVariables);												*
//	****** Insert Diagonal SubMatrixSparse and Sparse SubMatrix						*
//	* A.Insert(i,&W,&S);																			*
//	****************************************************************************
//	****** Insert a BzzFactorizedBandGauss															*
//	****** BzzFactorizedBandGauss B(m,m,l,u);														*
//	****** BzzMatrixSparse S(m,numVariables);												*
//	****** Insert Diagonal Band SubMatrix and Sparse SubMatrix						*
//	* A.Insert(k,&B,&S);																			*
//	****************************************************************************
//	****** Insert a BzzMatrixSparse															*
//	****** BzzMatrixSparse P(k,k);															*
//	****** BzzMatrixSparse S(k,numVariables);												*
//	****** Insert Diagonal Sparse SubMatrix and Sparse SubMatrix					*
//	* A.Insert(l,&P,&S);																			*
//	****************************************************************************
//	****** Function for solving systems:													*
//	* Solve(&A,&b);																				*
//	****************************************************************************
//	****** It is possible to change the diagonal elements with the functions:	*
//	* A.Insert(j,f,&S);																			*
//	* A.Insert(i,&Z,&S);																			*
//	* A.Insert(i,&L,&S);																			*
//	* A.Insert(i,&R,&S);																			*
//	* A.Insert(i,&W,&S);																			*
//	* A.Insert(k,&B,&S);																			*
//	* A.Insert(l,&P,&S);																			*
//	****** Note: You can not change the type of element.								*
//	****** Note: You can not change the extra elements (in the sparse matrix).	*
//	****** Note: It is possible to change a matrix already factorized.			*
//	****************************************************************************
//	****** Other functions:																		*
//	* iterations = A.GetIterations();														*
//	* SetStartX(&x0);																				*
//	* SetMaxIterations(numIterations);														*
//	* SetToleranceAbs(tollA);																	*
//	* SetToleranceRel(tollR);																	*
//	* A.BzzPrint("Comments");																	*
//	* A.BzzMessage("Comments");																*
//	****************************************************************************

#ifndef FACTORED_DOUBLE_DIAGONAL_BLOCK_AND_SPARSE_HPP
#define FACTORED_DOUBLE_DIAGONAL_BLOCK_AND_SPARSE_HPP

//	============================================================================
//	=====================< class BzzFactorizedDiagonalBlocksAndSparse >=======================
//	============================================================================

class BzzVector;
class BzzMatrix;
class BzzMatrixLeft;
class BzzMatrixRight;
class BzzMatrixSymmetric;
class BzzMatrixSparse;
class BzzFactorized;
class BzzFactorizedPLR;
class BzzFactorizedGauss;
class BzzFactorizedSymmetric;

class BzzFactorizedDiagonalBlocksAndSparse : public BzzBaseClass
{
private:
	static const char* const BZZ_ERROR;
	static int count; // for whoAmI

	int	whoAmI,
		numDiagonalMatrices,
		numVariables,
		numIterations,
		maxIterations;
	double	tollAbs,
		tollRel;

	BzzVector x0;

	int* diagonalBlockAndSparseType;
	double* fi;
	BzzFactorizedGauss* Ai;
	BzzFactorizedSymmetric* Wi;
	BzzFactorizedBandGauss* Bi;
	BzzFactorizedSparseGauss* Pi;
	BzzMatrixLeft* Li;
	BzzMatrixRight* Ri;
	BzzMatrixSparse* Si;

	void Initialize(int numvar, int nummat);	// initialising constructors

public:
	//	============================================================================
	//	*******************************< constructors >*****************************
	//	============================================================================

	BzzFactorizedDiagonalBlocksAndSparse(void);
	BzzFactorizedDiagonalBlocksAndSparse(int numvar, int nummat);

	//	============================================================================
	//	******************************< destructor >********************************
	//	============================================================================
	~BzzFactorizedDiagonalBlocksAndSparse(void);

	//	============================================================================
	//	================================< Functions >===============================
	//	============================================================================
	friend void Solve(BzzFactorizedDiagonalBlocksAndSparse* Q, BzzVector* b);

	virtual void ObjectBzzPrint(void);

	void SetStartX(BzzVector* xx);
	void SetMaxIterations(int nit) { numIterations = nit; }
	void SetToleranceAbs(double tolA) { tollAbs = tolA; }
	void SetToleranceRel(double tolR) { tollRel = tolR; }
	int GetIterations(void) { return numIterations; }

	void Insert(int i, double f, BzzMatrixSparse* si);
	void Insert(int i, BzzMatrix* ai, BzzMatrixSparse* si);
	void Insert(int i, BzzMatrixLeft* li, BzzMatrixSparse* si);
	void Insert(int i, BzzMatrixRight* ri, BzzMatrixSparse* si);
	void Insert(int i, BzzFactorizedBandGauss* bi, BzzMatrixSparse* si);
	void Insert(int i, BzzMatrixBand* bi, BzzMatrixSparse* si);
	void Insert(int i, BzzMatrixSymmetric* wi, BzzMatrixSparse* si);
	void Insert(int i, BzzMatrixSparse* pi, BzzMatrixSparse* si);

	void Insert(int i, double f);
	void Insert(int i, BzzMatrix* ai);
	void Insert(int i, BzzMatrixLeft* li);
	void Insert(int i, BzzMatrixRight* ri);
	void Insert(int i, BzzMatrixSymmetric* wi);
	void Insert(int i, BzzFactorizedBandGauss* bi);
	void Insert(int i, BzzMatrixBand* bi);
	void Insert(int i, BzzMatrixSparse* pi);
	friend void Delete(BzzFactorizedDiagonalBlocksAndSparse* A);
	friend void ChangeDimensions(int numv, int mat,
		BzzFactorizedDiagonalBlocksAndSparse* result);
};

#endif // FACTORED_DOUBLE_DIAGONAL_BLOCK_AND_SPARSE_HPP