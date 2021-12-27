// BZZMATH: Release 7.0

//	=============================< Plot.hpp >=================================
//	* Class BzzPlot			*
// * Description:																					*
// *																									*
// * Examples: 							*
//	============================================================================

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	07-1998	Date Written.

//	============================================================================
//	******	BzzPlot constructors:															*

#ifndef BZZ_PLOT_HPP
#define BZZ_PLOT_HPP

const int BZZ_SOLID = 0;
const int BZZ_DOT = 100;
const int BZZ_DASH = 200;
const int BZZ_DASHDOT = 300;
const int BZZ_DASHDOTDOT = 400;
const int BZZ_NOLINE = 500;

const int BZZ_BLACK = 0;
const int BZZ_LIGHTRED = 1;
const int BZZ_LIGHTBLUE = 2;
const int BZZ_LIGHTMAGENTA = 3;
const int BZZ_BROWN = 4;
const int BZZ_CYAN = 5;
const int BZZ_GREEN = 6;
const int BZZ_MAGENTA = 7;
const int BZZ_RED = 8;
const int BZZ_GRAY = 9;
const int BZZ_BLUE = 10;
const int BZZ_WHITE = 11;
const int BZZ_LIGHTGREEN = 12;
const int BZZ_LIGHTCYAN = 13;
const int BZZ_YELLOW = 14;
const int BZZ_BRIGHTWHITE = 15;

const int BZZ_NULL = 0;
const int BZZ_CIRCLE = 1000;
const int BZZ_SQUARE = 2000;
const int BZZ_DIAMOND = 3000;
const int BZZ_PLUS = 4000;
const int BZZ_X = 5000;
const int BZZ_TRIANGLE = 6000;

const int BZZ_NORMAL_SCALE = 1;
const int BZZ_LOG10_SCALE = 2;
const int BZZ_INVERSE_SCALE = 3;

//	============================================================================
//	============================< class BzzPlot >=============================
//	============================================================================
class BzzPlot //: public BzzBaseClass
{
private:
	CString	strTitle,
		strOrdinate,
		strAbscissa;

	CDC* pDC;
	CRect rect;

	double xMin, xMax;
	double yMin, yMax;

	BOOL	bGrid,
		bColorGrid,
		bColorBackground,
		bFirst;
	char type;

	COLORREF	colorGridRGB,
		colorBackgroundRGB;

	int	colorGrid,
		colorBackground,
		pointsFontTitle,
		pointsFontOrdinate,
		pointsFontAbscissa,
		pointsFontNumbers,
		scaleX,
		scaleY,
		thickness,
		opaque,
		pointDimensions;

	int leftD, rightD;

	//	void Initialize(CDC *pDC,CRect &rct);
	void Initialize(void);
	BzzVectorInt iL; //,iP;
	BOOL biL; //,biP;

public:
	//	============================================================================
	//	***************************< constructors >*********************************
	//	============================================================================
		// default // BzzPlot plot;
	BzzPlot(void);

	// copy-initializer // BzzPlot PlotA = PlotB;
//	BzzPlot(BzzPlot &rval);

	// Standard
//	BzzPlot(CDC *pDC,CRect &rect);

//	============================================================================
//	****************************< destructor >**********************************
//	============================================================================
	~BzzPlot(void) {};

	//	============================================================================
	//	********************< Non-modifying access functions >**********************
	//	============================================================================
	void GetDefault(void);

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
	///////////////////////////////////////////////////////////////////////////////
		// plot(&dc,rect,x,F); or plot(&dc,rect,x,F,xMin,xMax);
	void operator ()(CDC* pDC, CRect& rct,
		BzzVector& x, BzzMatrix& F,
		double xMi = BZZ_BIG, double xMa = BZZ_BIG); //,double yMin = 0.,double yMax = 0.);
	// plot(&dc,rect,x,F,j); or plot(&dc,rect,x,F,j,xMin,xMax);
	void operator ()(CDC* pDC, CRect& rct,
		BzzVector& x, BzzMatrix& F,
		BzzVectorInt& j,
		double xMi = BZZ_BIG, double xMa = BZZ_BIG); //,double yMin = 0.,double yMax = 0.);
///////////////////////////////////////////////////////////////////////////////
	void operator ()(CDC* pDC, CRect& rct,
		BzzMatrixSparse& X, BzzMatrixSparse& Y,
		BzzVectorInt& j);
	///////////////////////////////////////////////////////////////////////////////

	//	============================================================================
	//	===========================< Modifying Functions >==========================
	//	============================================================================
	void SetPointDimensions(int pd)
	{
		pointDimensions = pd;
	}

	void SetTitle(CString& str)
	{
		if (bFirst == TRUE)
			strTitle = str;
	}
	void SetOrdinate(CString& str)
	{
		if (bFirst == TRUE)
			strOrdinate = str;
	}
	void SetAbscissa(CString& str)
	{
		if (bFirst == TRUE)
			strAbscissa = str;
	}
	void SetTitle(char* str)
	{
		if (bFirst == TRUE)
			strTitle = str;
	}
	void SetOrdinate(char* str)
	{
		if (bFirst == TRUE)
			strOrdinate = str;
	}
	void SetAbscissa(char* str)
	{
		if (bFirst == TRUE)
			strAbscissa = str;
	}

	void SetSolidThickness(int th)
	{
		if (bFirst == TRUE && th >= 1 && th <= 10)
			thickness = th;
	}

	void SetOpaque(void)
	{
		opaque = 1;
	}
	void SetTransparent(void)
	{
		opaque = 0;
	}
	void SetXMin(double xMi) { xMin = xMi; }
	void SetXMax(double xMa) { xMax = xMa; }
	void SetYMin(double yMi) { yMin = yMi; }
	void SetYMax(double yMa) { yMax = yMa; }
	void SetGrid(BOOL bG) { bGrid = bG; }
	void SetGridColor(int clrG) { colorGrid = clrG; }
	void SetBackgroundColor(int clrB) { colorBackground = clrB; }
	void SetGridColor(COLORREF clrG) { colorGridRGB = clrG;bColorGrid = TRUE; }
	void SetBackgroundColor(COLORREF clrB)
	{
		colorBackgroundRGB = clrB;bColorBackground = TRUE;
	}
	void SetPointsFontTitle(int pointsT) { pointsFontTitle = pointsT; }
	void SetPointsFontOrdinate(int pointsO) { pointsFontOrdinate = pointsO; }
	void SetPointsFontAbscissa(int pointsA) { pointsFontAbscissa = pointsA; }
	void SetPointsFontNumbers(int pointsN) { pointsFontNumbers = pointsN; }
	void SetXScale(int s = BZZ_NORMAL_SCALE) { scaleX = s; }
	void SetYScale(int s = BZZ_NORMAL_SCALE) { scaleY = s; }
	void Modify(void) { bFirst = TRUE; }
	void GetLeftRight(int& lD, int& rD)
	{
		lD = leftD;rD = rightD;
	}
};

#endif // BZZ_PLOT_HPP