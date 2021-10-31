#pragma once
#include <COBIA.h>

using namespace COBIA;

class RealParameter :
	public CapeOpenObject<RealParameter>,
	public CAPEOPEN_1_2::CapeIdentificationAdapter<RealParameter>,
	public CAPEOPEN_1_2::CapeParameterAdapter<RealParameter>,
	public CAPEOPEN_1_2::CapeRealParameterAdapter<RealParameter>,
	public CAPEOPEN_1_2::CapeParameterSpecificationAdapter<RealParameter>,
	public CAPEOPEN_1_2::CapeRealParameterSpecificationAdapter<RealParameter> {

	// Members
	CapeStringImpl& unitName;
	CapeStringImpl paramName;
	CapeReal value, defaultValue, upperBound, lowerBound;
	CapeReal dimensionality[9] = { 0, 0, 0, 0, 0, 0, 0, 0 ,0 }; // TODO: Use built-in CapeEnumeration
	CAPEOPEN_1_2::CapeParamMode paramMode;
	CAPEOPEN_1_2::CapeValidationStatus paramValidationStatus;
	CAPEOPEN_1_2::CapeValidationStatus& unitValidationStatus;
	CapeBoolean& dirty;

public:
	const CapeStringImpl getDescriptionForErrorSource() {
		return paramName + COBIATEXT("Parameter of ") + unitName;
	}

	RealParameter(CapeStringImpl& _unitName, const COBIACHAR* _paramName,
		CapeReal _defaultValue, CapeReal _lowerBound, CapeReal upperBound,
		CAPEOPEN_1_2::CapeParamMode _paramMode,
		CAPEOPEN_1_2::CapeValidationStatus& _unitValidationStatus, CapeBoolean& _dirty) :
		unitName(_unitName), paramName(_paramName), value (_defaultValue),
		defaultValue(_defaultValue), lowerBound(_lowerBound), upperBound(upperBound),
		paramMode(_paramMode), unitValidationStatus(_unitValidationStatus), dirty(_dirty),
		paramValidationStatus(CAPEOPEN_1_2::CAPE_NOT_VALIDATED) {
	}

	~RealParameter() {
	}

	//CAPEOPEN_1_2::ICapeIdentification

	void getComponentName(/*out*/ CapeString name) {
		name = this->paramName;
	}

	void putComponentName(/*in*/ CapeString name) {
		throw cape_open_error(COBIAERR_Denied);
	}

	void getComponentDescription(/*out*/ CapeString desc) {
		desc = L"Real Parameter";
	}

	void putComponentDescription(/*in*/ CapeString desc) {
		throw cape_open_error(COBIAERR_Denied);
	}

	//CAPEOPEN_1_2::ICapeParameter
	
	CAPEOPEN_1_2::CapeValidationStatus getValStatus() {
		return paramValidationStatus;
	}
	
	CAPEOPEN_1_2::CapeParamMode getMode() {
		return paramMode;
	}
	
	CAPEOPEN_1_2::CapeParamType getType() {
		return CAPEOPEN_1_2::CAPE_PARAMETER_REAL;
	}

	CapeBoolean Validate(/*out*/ CapeString message) {
		CapeBoolean val = true;

		if (getType() != CAPEOPEN_1_2::CAPE_PARAMETER_REAL || getMode() != paramMode) {
			val = false;
			message = COBIATEXT("Parameter does not meet specifications");
			throw cape_open_error(COBIAERR_UnknownError);
		}
		return val;
	}

	void Reset() {
		value = this->defaultValue;
	}
	
	//CAPEOPEN_1_2::ICapeRealParameter
	
	CapeReal getValue() {
		return value;
	}
	
	void putValue(/*in*/ CapeReal value) {
		paramValidationStatus = CAPEOPEN_1_2::CAPE_NOT_VALIDATED;
		unitValidationStatus = CAPEOPEN_1_2::CAPE_NOT_VALIDATED;
		this->value = value;
		dirty = true;
	}
	
	CapeReal getDefaultValue() {
		return defaultValue;
	}
	
	CapeReal getLowerBound() {
		return lowerBound;
	}
	
	CapeReal getUpperBound() {
		return upperBound;
	}

	// This method is not part of the code generator.
	void putDimensionality(/*in*/ CapeInteger dim, /*in*/ CapeReal power) {
		dimensionality[dim] = power;
	}

	void getDimensionality(/*out*/ CapeArrayReal _dimensionality) {
		_dimensionality.resize(9);
		_dimensionality[0] = dimensionality[0];	// CAPE_METER
		_dimensionality[1] = dimensionality[1];	// CAPE_KILOGRAM
		_dimensionality[2] = dimensionality[2];	// CAPE_SECOND
		_dimensionality[3] = dimensionality[3];	// CAPE_AMPERE
		_dimensionality[4] = dimensionality[4];	// CAPE_KELVIN
		_dimensionality[5] = dimensionality[5];	// CAPE_MOLE
		_dimensionality[6] = dimensionality[6];	// CAPE_CANDELA
		_dimensionality[7] = dimensionality[7];	// CAPE_RADIAN
		_dimensionality[8] = dimensionality[8];	// CAPE_DIFFERENCE_FLAG
	}
	
	CapeBoolean Validate(/*in*/ CapeReal value,/*out*/ CapeString message) {
		CapeBoolean val = true;
		if (value < lowerBound || value > upperBound) {
			val = false;
			message = COBIATEXT("Parameter Value is out of bound");
			throw cape_open_error(COBIAERR_InvalidArgument);
		}
		paramValidationStatus = (val) ? CAPEOPEN_1_2::CAPE_VALID : CAPEOPEN_1_2::CAPE_INVALID;
		return val;
	}
};

using RealParameterPtr = CapeOpenObjectSmartPointer<RealParameter>;
