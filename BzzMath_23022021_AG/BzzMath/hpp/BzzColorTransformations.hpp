// BZZMATH: Release 7.0

// Revision History (MM-YYYY)
//	Author	Guido Buzzi-Ferraris
//	09-2005	Date Written.

#ifndef BZZ_COLOR_TRANSFORMATIONS
#define BZZ_COLOR_TRANSFORMATIONS

class BzzColorTransformations : public BzzBaseClass
{
private:
	char** colorName;
public:
	BzzColorTransformations(void);
	~BzzColorTransformations(void);
	void GetHLSFromRGB(int r, int g, int b, int* h, int* l, int* s);
	void GetRGBFromHLS(int h, int l, int s, int* r, int* g, int* b);
	double GetHue(double v1, double v2, double vH);
	void GetHLSFromSelectedColors(int i, char* name, int* h, int* l, int* s);
	void GetRGBFromSelectedColors(int i, char* name, int* r, int* g, int* b);
	void GetRGBHLSFromSelectedColors(int i, char* name, int* r, int* g, int* b, int* h, int* l, int* s);
	void GetSelectedColorName(int i, char* name);
	virtual void ObjectBzzPrint(void);
	void PrintSelectedColors(void);
	void PrintSelectedColors(int i);
	int FindBestHLSCheckFromSelectedColors(int h, int l, int s);
	double CheckHLS(int h1, int l1, int s1, int h2, int l2, int s2);
};

#endif // BZZ_COLOR_TRANSFORMATIONS