// BZZMATH: Release 7.0

//	==================< BzzFactorizedFourBlocksGauss.hpp >==================
//	* class BzzFactorizedFourBlocksGauss:	 											*
//	* This class is used for solving linear systems with four matrices 			*
//	* The matrix A11 can be structured														*
//	* The matrices A12 and A21 can be sparse												*
//	* The matrix A22 is a normal matrix														*
// * Description:																					*
//	*					Metodi Numerici e Software in C++ (Capitolo 6, 7)				*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
//	*																									*
// * Examples: c:\bzzmath\examples\exfourbl.cpp											*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	04-1997	Date Written first version

//	============================================================================
//	******	constructor:																		*
//	* BzzFactorizedFourBlocksGauss(&A11,&A12,&A21,&A22);										*
//	****************************************************************************
//	****** A11 and A22 are square matrices													*
//	****************************************************************************
//	****** A11 can be:																			*
//	BzzMatrix																						*
//	BzzMatrixBand																					*
//	BzzMatrixLeft																					*
//	BzzMatrixRight																					*
//	BzzMatrixBand																					*
//	BzzMatrixSparse																				*
//	BzzMatrixSymmetric																			*
//	Bzz TODO.......																				*
//	****************************************************************************
//	****** A12 and A21																			*
//	BzzMatrix																						*
//	BzzMatrixSparse																				*
//	****************************************************************************
//	****** A12 and A21																			*
//	BzzMatrix																						*
//	****************************************************************************
//	****** Function for solving systems:													*
//	* Solve(&A,&b);																				*
//	****************************************************************************
//	****** Other functions:																		*
//	* A.BzzPrint("Comments");																	*
//	* A.BzzMessage("Comments");																*
//	****************************************************************************

#ifndef FACTORIZED_DOUBLE_FOUR_BLOCKS_HPP
#define FACTORIZED_DOUBLE_FOUR_BLOCKS_HPP

//	============================================================================
//	================< class BzzFactorizedFourBlocksGauss >==================
//	============================================================================

class BzzFactorizedFourBlocksGauss : public BzzBaseClass
{
private:
	enum FactorizedFourBlocksStatus
	{
		UNFACTORIZED,
		FACTORIZED
	}factorizedStatus;
	enum BzzFourBlocksType
	{
		A11_FACTORIZED_GAUSS_A12_DENSE_A21_DENSE_A22_DENSE,
		A11_FACTORIZED_GAUSS_A12_SPARSE_A21_DENSE_A22_DENSE,
		A11_FACTORIZED_GAUSS_A12_DENSE_A21_SPARSE_A22_DENSE,
		A11_FACTORIZED_GAUSS_A12_SPARSE_A21_SPARSE_A22_DENSE,

		A11_FACTORIZED_BAND_GAUSS_A12_DENSE_A21_DENSE_A22_DENSE,
		A11_FACTORIZED_BAND_GAUSS_A12_SPARSE_A21_DENSE_A22_DENSE,
		A11_FACTORIZED_BAND_GAUSS_A12_DENSE_A21_SPARSE_A22_DENSE,
		A11_FACTORIZED_BAND_GAUSS_A12_SPARSE_A21_SPARSE_A22_DENSE,

		A11_FACTORIZED_TRIDIAGONALBLOCK_GAUSS_A12_DENSE_A21_DENSE_A22_DENSE,
		A11_FACTORIZED_TRIDIAGONALBLOCK_GAUSS_A12_SPARSE_A21_DENSE_A22_DENSE,
		A11_FACTORIZED_TRIDIAGONALBLOCK_GAUSS_A12_DENSE_A21_SPARSE_A22_DENSE,
		A11_FACTORIZED_TRIDIAGONALBLOCK_GAUSS_A12_SPARSE_A21_SPARSE_A22_DENSE,

		A11_FACTORIZED_SPARSE_GAUSS_A12_SPARSE_A21_SPARSE_A22_DENSE,

		A11_DIAGONALBLOCK_A12_DIAGONALBLOCK_A21_DIAGONALBLOCK_A22_DIAGONALBLOCK,
		A11_DIAGONALBLOCK_A12_DIAGONALBLOCK_A21_DIAGONALBLOCK_A22_TRIDIAGONALBLOCK,
		A11_DIAGONALBLOCK_A12_DIAGONALBLOCK_A21_DIAGONALBLOCK_A22_BAND,

		A11_DIAGONAL_A12_SPARSE_A21_NORMAL,
		A11_DIAGONAL_A12_NORMAL_A21_SPARSE,
		A11_DIAGONAL_A12_SPARSE_A21_SPARSE,

		SYMMETRIC
	}fourBlockType;

	static const char* const BZZ_ERROR;
	static int count; // for whoAmI

	int	whoAmI,
		n11, n22,
		numVariables,
		numBlocks, dim11Block, dim22Block;

	BzzFactorizedGauss* F11;
	BzzFactorizedBandGauss* FB11;
	BzzFactorizedTridiagonalBlocksGauss* FTDB11;
	BzzFactorizedSparseGauss* FS11;
	BzzFactorizedDiagonalBlocksGauss* FDB11;

	BzzMatrix* M12, * M21, * M22;
	BzzMatrixSparse* MS12, * MS21;
	BzzMatrixDiagonalBlocks* MDB12, * MDB21, * MDB22;
	BzzMatrixTridiagonalBlocks* MTDB22;
	BzzMatrixBand* MB22;

	BzzFactorizedGauss F22;
	BzzFactorizedDiagonalBlocksGauss FDB22;
	BzzFactorizedTridiagonalBlocksGauss FTDB22;
	BzzFactorizedBandGauss FB22;

	void Initialize(int numv);
public:
	//	============================================================================
	//	*******************************< constructors >*****************************
	//	============================================================================

	BzzFactorizedFourBlocksGauss(void);

	// dense dense dense dense
	BzzFactorizedFourBlocksGauss(BzzFactorizedGauss* a11, BzzMatrix* a12, BzzMatrix* a21,
		BzzMatrix* a22);
	void operator()(BzzFactorizedGauss* a11, BzzMatrix* a12, BzzMatrix* a21,
		BzzMatrix* a22);

	// dense sparse dense dense
	BzzFactorizedFourBlocksGauss(BzzFactorizedGauss* a11, BzzMatrixSparse* a12, BzzMatrix* a21,
		BzzMatrix* a22);
	void operator()(BzzFactorizedGauss* a11, BzzMatrixSparse* a12, BzzMatrix* a21,
		BzzMatrix* a22);

	// dense dense sparse dense
	BzzFactorizedFourBlocksGauss(BzzFactorizedGauss* a11, BzzMatrix* a12, BzzMatrixSparse* a21,
		BzzMatrix* a22);
	void operator()(BzzFactorizedGauss* a11, BzzMatrix* a12, BzzMatrixSparse* a21,
		BzzMatrix* a22);

	// dense sparse sparse dense
	BzzFactorizedFourBlocksGauss(BzzFactorizedGauss* a11, BzzMatrixSparse* a12,
		BzzMatrixSparse* a21, BzzMatrix* a22);
	void operator()(BzzFactorizedGauss* a11, BzzMatrixSparse* a12,
		BzzMatrixSparse* a21, BzzMatrix* a22);

	// band dense dense dense
	BzzFactorizedFourBlocksGauss(BzzFactorizedBandGauss* a11, BzzMatrix* a12, BzzMatrix* a21,
		BzzMatrix* a22);
	void operator()(BzzFactorizedBandGauss* a11, BzzMatrix* a12, BzzMatrix* a21,
		BzzMatrix* a22);

	// band sparse dense dense
	BzzFactorizedFourBlocksGauss(BzzFactorizedBandGauss* a11, BzzMatrixSparse* a12, BzzMatrix* a21,
		BzzMatrix* a22);
	void operator()(BzzFactorizedBandGauss* a11, BzzMatrixSparse* a12, BzzMatrix* a21,
		BzzMatrix* a22);

	// band dense sparse dense
	BzzFactorizedFourBlocksGauss(BzzFactorizedBandGauss* a11, BzzMatrix* a12, BzzMatrixSparse* a21,
		BzzMatrix* a22);
	void operator()(BzzFactorizedBandGauss* a11, BzzMatrix* a12, BzzMatrixSparse* a21,
		BzzMatrix* a22);

	// band sparse sparse dense
	BzzFactorizedFourBlocksGauss(BzzFactorizedBandGauss* a11, BzzMatrixSparse* a12,
		BzzMatrixSparse* a21, BzzMatrix* a22);
	void operator()(BzzFactorizedBandGauss* a11, BzzMatrixSparse* a12,
		BzzMatrixSparse* a21, BzzMatrix* a22);

	// tridia dense dense dense
	BzzFactorizedFourBlocksGauss(BzzFactorizedTridiagonalBlocksGauss* a11, BzzMatrix* a12, BzzMatrix* a21,
		BzzMatrix* a22);
	void operator()(BzzFactorizedTridiagonalBlocksGauss* a11, BzzMatrix* a12, BzzMatrix* a21,
		BzzMatrix* a22);

	// tridia sparse dense dense
	BzzFactorizedFourBlocksGauss(BzzFactorizedTridiagonalBlocksGauss* a11, BzzMatrixSparse* a12, BzzMatrix* a21,
		BzzMatrix* a22);
	void operator()(BzzFactorizedTridiagonalBlocksGauss* a11, BzzMatrixSparse* a12, BzzMatrix* a21,
		BzzMatrix* a22);

	// tridia dense sparse dense
	BzzFactorizedFourBlocksGauss(BzzFactorizedTridiagonalBlocksGauss* a11, BzzMatrix* a12, BzzMatrixSparse* a21,
		BzzMatrix* a22);
	void operator()(BzzFactorizedTridiagonalBlocksGauss* a11, BzzMatrix* a12, BzzMatrixSparse* a21,
		BzzMatrix* a22);

	// tridia sparse sparse dense
	BzzFactorizedFourBlocksGauss(BzzFactorizedTridiagonalBlocksGauss* a11, BzzMatrixSparse* a12,
		BzzMatrixSparse* a21, BzzMatrix* a22);
	void operator()(BzzFactorizedTridiagonalBlocksGauss* a11, BzzMatrixSparse* a12,
		BzzMatrixSparse* a21, BzzMatrix* a22);

	BzzFactorizedFourBlocksGauss(BzzFactorizedSparseGauss* a11, BzzMatrixSparse* a12,
		BzzMatrixSparse* a21, BzzMatrix* a22);
	void operator()(BzzFactorizedSparseGauss* a11, BzzMatrixSparse* a12,
		BzzMatrixSparse* a21, BzzMatrix* a22);

	BzzFactorizedFourBlocksGauss(BzzFactorizedDiagonalBlocksGauss* a11,
		BzzMatrixDiagonalBlocks* a12, BzzMatrixDiagonalBlocks* a21,
		BzzMatrixDiagonalBlocks* a22);
	void operator()(BzzFactorizedDiagonalBlocksGauss* a11,
		BzzMatrixDiagonalBlocks* a12, BzzMatrixDiagonalBlocks* a21,
		BzzMatrixDiagonalBlocks* a22);

	BzzFactorizedFourBlocksGauss(BzzFactorizedDiagonalBlocksGauss* a11,
		BzzMatrixDiagonalBlocks* a12, BzzMatrixDiagonalBlocks* a21,
		BzzMatrixTridiagonalBlocks* a22);
	void operator()(BzzFactorizedDiagonalBlocksGauss* a11,
		BzzMatrixDiagonalBlocks* a12, BzzMatrixDiagonalBlocks* a21,
		BzzMatrixTridiagonalBlocks* a22);

	BzzFactorizedFourBlocksGauss(BzzFactorizedDiagonalBlocksGauss* a11,
		BzzMatrixDiagonalBlocks* a12, BzzMatrixDiagonalBlocks* a21,
		BzzMatrixBand* a22);
	void operator()(BzzFactorizedDiagonalBlocksGauss* a11,
		BzzMatrixDiagonalBlocks* a12, BzzMatrixDiagonalBlocks* a21,
		BzzMatrixBand* a22);

	//	============================================================================
	//	******************************< destructor >********************************
	//	============================================================================
	~BzzFactorizedFourBlocksGauss(void);

	//	============================================================================
	//	================================< Functions >===============================
	//	============================================================================
	friend void Solve(BzzFactorizedFourBlocksGauss* Q, BzzVector* b);
	virtual void ObjectBzzPrint(void);
};

#endif // FACTORIZED_FOUR_BLOCKS_HPP