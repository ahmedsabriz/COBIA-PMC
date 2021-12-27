// BZZMATH: Release 7.0
// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	07-2005	Date Written.

#ifndef BZZ_SUDOKU_HPP
#define BZZ_SUDOKU_HPP
void Sudoku(BzzMatrixInt& A);

//	============================================================================
//	============================< class BzzSudoku >=============================
//	============================================================================
class BzzSudoku
{
private:
public:
	//	============================================================================
	//	***************************< constructors >*********************************
	//	============================================================================
		// default
	BzzSudoku(void);
	~BzzSudoku(void);

	int operator()(BzzMatrixInt* B, int* minSize, int* i, int* j, BzzVectorInt* ch);
};

#endif // BZZ_MATRIX_INT_HPP