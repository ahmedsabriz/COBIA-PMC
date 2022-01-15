#pragma once
#include <COBIA.h>
#include "Parameter.h"

using namespace COBIA;

class ParameterOption :
	public CapeOpenObject<ParameterOption>,
	public CAPEOPEN_1_2::CapeIdentificationAdapter<ParameterOption>,
	public CAPEOPEN_1_2::CapeParameterAdapter<ParameterOption>,
	public CAPEOPEN_1_2::CapeParameterSpecificationAdapter<ParameterOption>,
	public CAPEOPEN_1_2::CapeStringParameterAdapter<ParameterOption>,
	public CAPEOPEN_1_2::CapeStringParameterSpecificationAdapter<ParameterOption>,
	public Parameter {

	// Members
	CapeArrayStringImpl& optionNames;
	CapeStringImpl defaultValue, value;
	
public:

	ParameterOption(CapeStringImpl& _unitName, CAPEOPEN_1_2::CapeValidationStatus& _unitValidationStatus,
		CapeBoolean& _dirty, const COBIACHAR* _paramName, CapeArrayStringImpl& _optionNames) :
		Parameter(_unitName, _unitValidationStatus, _dirty, _paramName, CAPEOPEN_1_2::CAPE_INPUT),
		optionNames(_optionNames) {
		paramValidationStatus = CAPEOPEN_1_2::CAPE_NOT_VALIDATED;
		defaultValue = optionNames[0];
		value = defaultValue;
	}

	~ParameterOption() {
	}

	//CAPEOPEN_1_2::ICapeParameter
	CAPEOPEN_1_2::CapeParamType getType() {
		return CAPEOPEN_1_2::CAPE_PARAMETER_STRING;
	}
	CapeBoolean Validate(/*out*/ CapeString message) {
		// First validation is for type and mode
		if (getType() != CAPEOPEN_1_2::CAPE_PARAMETER_STRING || getMode() != CAPEOPEN_1_2::CAPE_INPUT) {
			message = COBIATEXT("Parameter does not meet specifications");
			throw cape_open_error(COBIAERR_UnknownError);
		}
		return true;
	}
	void Reset() {
		this->value = defaultValue;
		dirty = true;
	}
	
	//CAPEOPEN_1_2::ICapeStringParameter
	void getValue(/*out*/ CapeString value) {
		value = this->value;
	}
	void putValue(/*in*/ CapeString value) {
		this->value = value;
	}
	void getDefaultValue(/*out*/ CapeString defaultValue) {
		defaultValue = this->defaultValue;
	}
	void getOptionList(/*out*/ CapeArrayString optionNames) {
		optionNames.resize(this->optionNames.size());
		for (size_t i = 0, length = this->optionNames.size(); i < length; i++)
		{
			optionNames[i] = this->optionNames[i];
		}
	}
	CapeBoolean getRestrictedToList() {
		return true;
	}
	CapeBoolean Validate(/*in*/ CapeString value,/*out*/ CapeString message) {
		// Second validation for value
		for (CapeString option : CapeArrayString(optionNames)) {
			if (value == option) {
				paramValidationStatus = CAPEOPEN_1_2::CAPE_VALID;
				return true;
			}
		}
		message = COBIATEXT("Option is not allowed"); // TODO: indicate position
		throw cape_open_error(COBIAERR_UnknownError);
		paramValidationStatus = CAPEOPEN_1_2::CAPE_INVALID;
		return false;
	}
	
};

using ParameterOptionPtr = CapeOpenObjectSmartPointer<ParameterOption>;
