#pragma once
#include <COBIA.h>
#include "Parameter.h"

using namespace COBIA;

class ParameterReal :
	public CapeOpenObject<ParameterReal>,
	public CAPEOPEN_1_2::CapeIdentificationAdapter<ParameterReal>,
	public CAPEOPEN_1_2::CapeParameterAdapter<ParameterReal>,
	public CAPEOPEN_1_2::CapeRealParameterAdapter<ParameterReal>,
	public CAPEOPEN_1_2::CapeParameterSpecificationAdapter<ParameterReal>,
	public CAPEOPEN_1_2::CapeRealParameterSpecificationAdapter<ParameterReal>,
	public Parameter {

	// Members
	CapeReal value, defaultValue, upperBound, lowerBound;
	CapeArrayReal& dimensionality;

public:

	ParameterReal(CapeStringImpl& _unitName, CAPEOPEN_1_2::CapeValidationStatus& _unitValidationStatus,
		CapeBoolean& _dirty, const COBIACHAR* _paramName, CAPEOPEN_1_2::CapeParamMode _paramMode,
		CapeReal _defaultValue, CapeReal _lowerBound, CapeReal _upperBound, CapeArrayReal& _dimensionality) :
		Parameter(_unitName, _unitValidationStatus, _dirty, _paramName, _paramMode),
		defaultValue(_defaultValue), lowerBound(_lowerBound), upperBound(_upperBound),
		dimensionality(_dimensionality) {
		paramValidationStatus = CAPEOPEN_1_2::CAPE_NOT_VALIDATED;
		value = defaultValue;
	}

	~ParameterReal() {
	}

	//CAPEOPEN_1_2::ICapeParameter
	CAPEOPEN_1_2::CapeParamType getType() {
		return CAPEOPEN_1_2::CAPE_PARAMETER_REAL;
	}
	CapeBoolean Validate(/*out*/ CapeString message) {
		// First validation is for type and mode
		if (getType() != CAPEOPEN_1_2::CAPE_PARAMETER_REAL || getMode() != paramMode) {
			message = COBIATEXT("Parameter does not meet specifications");
			throw cape_open_error(COBIAERR_UnknownError);
		}
		return true;
	}
	void Reset() {
		this->value = defaultValue;
		dirty = true;
	}
	
	//CAPEOPEN_1_2::ICapeRealParameter
	CapeReal getValue() {
		return this->value;
	}
	void putValue(/*in*/ CapeReal value) {
		this->value = value;
	}
	CapeReal getDefaultValue() {
		return this->defaultValue;
	}
	CapeReal getLowerBound() {
		return this->lowerBound;
	}
	CapeReal getUpperBound() {
		return this->upperBound;
	}
	void getDimensionality(/*out*/ CapeArrayReal dimensionality) {
		dimensionality.resize(9);
		dimensionality[0] = this->dimensionality[0];	// CAPE_METER
		dimensionality[1] = this->dimensionality[1];	// CAPE_KILOGRAM
		dimensionality[2] = this->dimensionality[2];	// CAPE_SECOND
		dimensionality[3] = this->dimensionality[3];	// CAPE_AMPERE
		dimensionality[4] = this->dimensionality[4];	// CAPE_KELVIN
		dimensionality[5] = this->dimensionality[5];	// CAPE_MOLE
		dimensionality[6] = this->dimensionality[6];	// CAPE_CANDELA
		dimensionality[7] = this->dimensionality[7];	// CAPE_RADIAN
		dimensionality[8] = this->dimensionality[8];	// CAPE_DIFFERENCE_FLAG
	}
	CapeBoolean Validate(/*in*/ CapeReal value,/*out*/ CapeString message) {
		// Second validation is for value
		CapeBoolean val = true;
		if (value < lowerBound || value > upperBound) {
			val = false;
			message = COBIATEXT("Parameter Value is out of bound");
			throw cape_open_error(COBIAERR_UnknownError);
		}
		paramValidationStatus = (val) ? CAPEOPEN_1_2::CAPE_VALID : CAPEOPEN_1_2::CAPE_INVALID;
		return val;
	}
};

using ParameterRealPtr = CapeOpenObjectSmartPointer<ParameterReal>;
