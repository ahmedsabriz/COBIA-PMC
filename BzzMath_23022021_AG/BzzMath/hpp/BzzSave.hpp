// BZZMATH: Release 7.0

// ==============================< BzzSave.hpp >===============================
// * Class BzzSave for saving BzzMath objects on file									*
// * Description:																					*
//	*					Metodi Numerici e Software in C++ (Capitolo 3)					*
//	*					by G. Buzzi-Ferraris														*
// *					Addison Wesley Longman(1998)											*
// *																									*
// * Examples: c:\bzzmath\examples\exsave.cpp										 	*
// ============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	12-1996	Date Written.
//	10-1997	Added BzzFactordeGauss, BzzFactorizedGauss.
//	11-1999	Added ASCII file.

////////////////// Release 5.0
//	04-2005	Added BzzVectorArray, BzzVectorArray, BzzVectorIntArray.
//	08-2005	Added BzzString.

// ============================================================================
// ****** Constructors for BzzSave:															*
// * BzzSave bzz("BZZ.DAT");																	*
// ****************************************************************************
// ***** Functions :																				*
// * BzzSave bzz("BZZ.DAT");		// formatted ASCII file								*
// * BzzSave bzz('*',"BZZ.BIN");	// binary file											*
// * BzzVector v(10);								 											*
// * bzz << v;																						*
// * bzz.End();																					*
// * BzzLoad load("BZZ.DAT");			// formatted ASCII file							*
// * BzzLoad load('*',"BZZ.BIN");	// binary file										*
// * load >> v;																					*
// ****************************************************************************

#ifndef BZZ_SAVE_HPP
#define BZZ_SAVE_HPP

// ============================================================================
// =============================< class BzzSave >==============================
// ============================================================================

class BzzSave
{
private:
	static const char* const BZZ_ERROR;
	FILE* fileSave;
	char type;
public:

	// ============================================================================
	// ***************************< constructors >*********************************
	// ============================================================================
		// default constructor;
	BzzSave(void);

	// copy-initializer
	BzzSave(BzzSave& rval) { BzzError("No copy constructor");rval.fileSave = 0; }

	// file constructor;
	BzzSave(char* filevector);
	BzzSave(const char* filevector);
	BzzSave(char t, char* filevector);
	BzzSave(const char t, const char* filevector);

	// ============================================================================
	// *****************************< destructor >*********************************
	// ============================================================================
	~BzzSave(void);

	// ============================================================================
	// **************************< Functions >*************************************
	// ============================================================================
	void operator()(const char* file);
	void operator()(const char* file, char);
	void operator()(char t, const char* file);
	void operator()(char t, const char* file, char);
	void End(void);

	BzzSave& operator << (int i);
	BzzSave& operator << (float x);
	BzzSave& operator << (double d);
	BzzSave& operator << (BzzVectorInt& iv);
	BzzSave& operator << (BzzVector& v);
	BzzSave& operator << (BzzMatrixInt& A);
	BzzSave& operator << (BzzMatrix& A);
	BzzSave& operator << (BzzMatrixLeft& A);
	BzzSave& operator << (BzzMatrixRight& A);
	BzzSave& operator << (BzzMatrixSymmetric& A);
	BzzSave& operator << (BzzMatrixDiagonal& A);
	BzzSave& operator << (BzzMatrixBand& A);
	BzzSave& operator << (BzzMatrixSparse& A);
	BzzSave& operator << (BzzFactorizedGauss& A);
	BzzSave& operator << (BzzVectorIntArray& iv);
	BzzSave& operator << (BzzVectorArray& iv);
	BzzSave& operator << (BzzString& str);
};

#endif // BZZ_SAVE_HPP