// BZZMATH: Release 7.0

// =============================< BzzString.hpp >==============================
// * Class BzzString 																			*
// ============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	08-1995	Date Written.

////////////////// Release 5.0
//	08-2005 Added GetLength, BzzPrintString functions.
//	09-2005 Added IsEmpty, Empty, MakeUpper, MakeLower functions.

#ifndef BZZ_STRING_HPP
#define BZZ_STRING_HPP

// ============================================================================
// ============================< class BzzString >=============================
// ============================================================================
class BzzString : public BzzBaseClass
{
private:
	static const char* const BZZ_ERROR;
	static int count; // for whoAmI
	static int countInScope;

	char* vector;
	int dimensions;
	int whoAmI;
	char shadow;

	// initialise constructors
	void Initialize(int nc);

	// private constructor BzzString('*',nc)
	BzzString(char, int nc);

public:

	// ============================================================================
	// ***************************< constructors >*********************************
	// ============================================================================
		// default constructor
		// BzzString s1;
	BzzString(void);

	// copy-initializer
	// BzzString s2 = sz;
	BzzString(BzzString& rval);

	// sizes and initialises at blanks
	BzzString(int nc);

	// BzzString s4 = "Milan"; // or BzzString s4("Milan");
	BzzString(char* initvalues);

	// ============================================================================
	// *****************************< destructor >*********************************
	// ============================================================================
	~BzzString(void);

	// ============================================================================
	// *******************< Non-modifying access functions >***********************
	// ============================================================================
	int Size(void) const { return dimensions; } // dimensions
	int GetLength(void) { return dimensions; }
	int WhoAmI(void) const { return whoAmI; }
	static int ObjectCount(void) { return count; }
	static int ObjectCountInScope(void) { return countInScope; }
	int IsEmpty(void);
	// receives the value of the vector with control
	char GetValue(int i) const;

	// ============================================================================
	// **********************< Modifying access functions >************************
	// ============================================================================
		// assigns and receives vector values with control
	char& operator () (int i);

	// assigns and receives values without control
	char& operator []
	(int i) {
		return vector[i];
	}
	void MakeUpper(void);
	void MakeLower(void);
	void Empty(void);

	// ============================================================================
	// *************************< assignment operators >***************************
	// ============================================================================
	BzzString& operator = (BzzString& rval);
	BzzString& operator = (char* rval);

	void operator = (char c);

	// ============================================================================
	// *************************< operators for tests >****************************
	// ============================================================================
	friend char operator ==
		(const BzzString& lval, const BzzString& rval);

	friend char operator !=
		(const BzzString& lval, const BzzString& rval);

	// ============================================================================
	// ====================< Non-modifying functions >=============================
	// ============================================================================

	// ******************************< BzzPrint >**********************************
	virtual void ObjectBzzPrint(void);
	void BzzPrintString(void);

	// ============================================================================
	// ======================< Modifying Functions >===============================
	// ============================================================================
	friend void Delete(BzzString* result); // eliminates BzzString
	friend void ChangeDimensions(int dim, BzzString* result,
		char);

	friend void Reverse(BzzString* result);	//inverts the vector

	friend void Sort(BzzString* result); // orders the vector

	// swaps the contents of two vectors
	friend void Swap(BzzString* lval, BzzString* rval);

	// ============================================================================
	// =========================< Other functions >================================
	// ============================================================================

	void FirstInLastOut(char f);
	void LastInFirstOut(char f);
	void Append(char f);
	void Append(BzzString& v);
	void Insert(int i, char f);
	void Insert(int i, BzzString& v);
	void DeleteElement(int i);
	void DeleteElements(BzzVectorInt& iv);
	void SwapElements(int i, int j);
	int Locate(char c);
	void GetString(char* str);

	// ============================================================================
	// ===========================< Save and Load >================================
	// ============================================================================

	void Save(char* filevector); // formatted
	void Save(char, char* filevector); // unformatted
	friend void Load(BzzString* result, char* filevector);
	friend void Load(BzzString* result, char, char* filevector);
};

// Friend functions with default arguments
void ChangeDimensions(int dim, BzzString* result,
	char zero = 0);


#endif // BZZ_STRING_HPP