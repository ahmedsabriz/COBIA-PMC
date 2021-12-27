// BZZMATH: Release 7.0

// ==============================< BzzLoad.hpp >===============================
// * Class BzzLoad for saving BzzMath objects on file									*
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
// ****** Constructors for BzzLoad:															*
// * BzzLoad load("BZZ.DAT");																	*
// ****************************************************************************
// ***** Functions :																				*
// * BzzSave bzz("BZZ.DAT");		// formatted ASCII file								*
// * BzzSave bzz("'*',BZZ.BIN");	// binary file											*
// * BzzVector v(10);								 											*
// * bzz << v;																						*
// * bzz.End();																					*
// * BzzLoad load("BZZ.DAT");			// formatted ASCII file							*
// * BzzLoad load('*',"BZZ.BIN");	// binary file										*
// * load >> v;																					*
// * load.End();																					*
// ****************************************************************************
#ifndef BZZ_LOAD_HPP
#define BZZ_LOAD_HPP

// ============================================================================
// =============================< class BzzLoad >==============================
// ============================================================================

class BzzLoad
{
private:
	static const char* const BZZ_ERROR;
	FILE* fileLoad;
	char type;
public:

	// ============================================================================
	// ***************************< constructors >*********************************
	// ============================================================================
		// default constructor;
	BzzLoad(void);

	// copy-initializer
	BzzLoad(BzzLoad& rval) { BzzError("No copy constructor");rval.fileLoad = 0; }

	// file constructor;
	BzzLoad(const char* file);
	BzzLoad(char t, const char* file);

	// ============================================================================
	// *****************************< destructor >*********************************
	// ============================================================================
	~BzzLoad(void) { fclose(fileLoad); }

	// ============================================================================
	// **************************< Functions >*************************************
	// ============================================================================
	void operator()(const char* file);
	void operator()(char t, const char* file);
	void End(void) { fclose(fileLoad); }

	BzzLoad& operator >> (int& i);
	BzzLoad& operator >> (float& x);
	BzzLoad& operator >> (double& d);
	BzzLoad& operator >> (BzzVectorInt& iv);
	BzzLoad& operator >> (BzzVector& v);
	BzzLoad& operator >> (BzzMatrixInt& A);
	BzzLoad& operator >> (BzzMatrix& A);
	BzzLoad& operator >> (BzzMatrixLeft& A);
	BzzLoad& operator >> (BzzMatrixRight& A);
	BzzLoad& operator >> (BzzMatrixSymmetric& A);
	BzzLoad& operator >> (BzzMatrixBand& A);
	BzzLoad& operator >> (BzzMatrixDiagonal& A);
	BzzLoad& operator >> (BzzMatrixSparse& A);
	BzzLoad& operator >> (BzzFactorizedGauss& A);
	BzzLoad& operator >> (BzzVectorIntArray& iv);
	BzzLoad& operator >> (BzzVectorArray& dv);
	BzzLoad& operator >> (BzzString& str);
};

#endif // BZZ_LOAD_HPP