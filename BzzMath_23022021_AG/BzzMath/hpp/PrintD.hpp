// BZZMATH: Release 3.1

//	=============================< PrintData.hpp >==============================
//	* Class BzzPrintData			*
// * Description:																					*
// *																									*
// * Examples: 							*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	07-1998	Date Written.

//	============================================================================
//	******	BzzPrintData constructors:															*

#ifndef BZZ_PRINTDATA_HPP
#define BZZ_PRINTDATA_HPP

//	============================================================================
//	==============================< class BzzData >=============================
//	============================================================================
class BzzData //: public BzzBaseClass
{
private:
	int numRows;
	int size;

public:
	CStringArray	strArray;
	BzzData(void);
	BzzData(int sz);
	void operator ()(CString str);
	int GetRows(void) { return numRows; }
	void RemoveAll(void)
	{
		if (GetRows() != 0)
			strArray.RemoveAll();
		numRows = 0;
		strArray.SetSize(size);
	}
};

//	============================================================================
//	============================< class BzzPrintData >=============================
//	============================================================================
class BzzPrintData //: public BzzBaseClass
{
private:
	CDC* pDC;
	CRect rectO;
	COLORREF colorText;
	COLORREF colorBackground;
	BOOL bFirst;

	int	minX,
		minY,
		maxX,
		maxY,
		yHeight,
		rectHeight,
		rectWidth,
		numRows,
		iRow,
		numArray;

	int pointsFont;
	//	void Initialize(CDC *pDC,CRect &rct);
	void Initialize(void);

public:
	//	============================================================================
	//	***************************< constructors >*********************************
	//	============================================================================
		// default // BzzPrintData plot;
	BzzPrintData(void);

	// copy-initializer // BzzPrintData PlotA = PlotB;
//	BzzPrintData(BzzPrintData &rval);

	// Standard
//	BzzPrintData(CDC *pDC,CRect &rect);

//	============================================================================
//	****************************< destructor >**********************************
//	============================================================================
	~BzzPrintData(void) {};

	//	============================================================================
	//	********************< Non-modifying access functions >**********************
	//	============================================================================

	//	============================================================================
	//	*********************< Modifying access functions >*************************
	//	============================================================================

	//	============================================================================
	//	**********************< assignment operators >******************************
	//	============================================================================

	//	============================================================================
	//	=======================< Non-modifying functions >==========================
	//	============================================================================

	//	******************************< BzzPrint >**********************************
	//	virtual void ObjectBzzPrint(void);

	//	*****************************< operator () >********************************
		// printData(str);
	void operator ()(CDC* pDC, CRect& rct, BzzData& data, int init = 0);

	//	============================================================================
	//	===========================< Modifying Functions >==========================
	//	============================================================================
	void SetPointsFont(int points) { pointsFont = points; }
	void SetTextColor(COLORREF colorT) { colorText = colorT; }
	void SetBkColor(COLORREF colorB) { colorBackground = colorB; };
	void Modify(void) { bFirst = TRUE; }
};

#endif // BZZ_PRINTDATA_HPP