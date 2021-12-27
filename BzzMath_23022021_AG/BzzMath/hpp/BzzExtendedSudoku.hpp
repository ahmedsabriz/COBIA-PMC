extern BzzMatrixInt S;
extern BzzVectorInt s;
extern int n1, n2, n3, n4;
extern BzzVectorInt vRCM, iCross, jCross, ir, ic, im, r, c;
extern BzzMatrixInt B, IM;
extern BzzVectorIntArray V;
extern BzzVectorInt v;

// BZZMATH: Release 7.0

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	02-2006	Date Written.

#ifndef BZZ_EXTENDED_SUDOKU_HPP
#define BZZ_EXTENDED_SUDOKU_HPP

// ============================================================================
// =======================< class BzzExtendedSudoku >==========================
// ============================================================================

class BzzExtendedSudoku
{
	// ============================================================================
	// ===================< Private Data and Functions >===========================
	// ============================================================================
private:
	char printAll;
	BzzMatrixInt A;
	void InitializeSudoku(void);
	int ControlSudoku(void);
	void PrintSudoku(void);
	int RecursionSudoku(void);
	int CancelSudoku(void);
	int CrossCancelSudoku(void);
	int StopSudoku(void);
	int PrintSudokuDeletion(void);
	void FindExistence(void);
	int maxLevel;

	// ============================================================================
	// ===================< Public Data and Functions >============================
	// ============================================================================
public:

	// ============================================================================
	// ***************************< constructors >*********************************
	// ============================================================================
	BzzExtendedSudoku(char* file);
	BzzExtendedSudoku(BzzMatrixInt& AA);

	// ============================================================================
	// ==============================< Functions >=================================
	// ============================================================================
	int operator () (void);
	void PrintAll(void)
	{
		printAll = 1;
	}
};

#endif // BZZ_EXTENDED_SUDOKU_HPP