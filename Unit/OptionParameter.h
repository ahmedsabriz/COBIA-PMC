#pragma once
#include <COBIA.h>

using namespace COBIA;

class OptionParameter :
	public CapeOpenObject<OptionParameter>,
	public CAPEOPEN_1_2::CapeIdentificationAdapter<OptionParameter>,
	public CAPEOPEN_1_2::CapeParameterAdapter<OptionParameter>,
	public CAPEOPEN_1_2::CapeStringParameterAdapter<OptionParameter>,
	public CAPEOPEN_1_2::CapeParameterSpecificationAdapter<OptionParameter>,
	public CAPEOPEN_1_2::CapeStringParameterSpecificationAdapter<OptionParameter> {

	// Members
	CapeStringImpl& unitName;
	CAPEOPEN_1_2::CapeValidationStatus& unitValidationStatus;
	CapeBoolean& dirty;
	CapeArrayStringImpl& optionNames;

	CapeStringImpl paramName;
	CAPEOPEN_1_2::CapeValidationStatus paramValidationStatus;

	CapeStringImpl defaultValue, value;
	
public:

	const CapeStringImpl getDescriptionForErrorSource() {
		return paramName + COBIATEXT("Parameter of ") + unitName;
	}

	OptionParameter(CapeStringImpl& _unitName, CAPEOPEN_1_2::CapeValidationStatus& _unitValidationStatus,
		CapeBoolean& _dirty, CapeArrayStringImpl& _optionNames,
		const COBIACHAR* _paramName) :
		unitName(_unitName), unitValidationStatus(_unitValidationStatus),
		dirty(_dirty), optionNames(_optionNames),
		paramName(_paramName) {
		paramValidationStatus = CAPEOPEN_1_2::CAPE_NOT_VALIDATED;
		defaultValue = optionNames[0];
		value = defaultValue;
	}

	~OptionParameter() {
	}

	//CAPEOPEN_1_2::ICapeIdentification
	
	void getComponentName(/*out*/ CapeString name) {
		name = this->paramName;
	}
	
	void putComponentName(/*in*/ CapeString name) {
		throw cape_open_error(COBIAERR_Denied);
	}
	
	void getComponentDescription(/*out*/ CapeString desc) {
		desc = COBIATEXT("Option Parameter Array");
	}
	
	void putComponentDescription(/*in*/ CapeString desc) {
		throw cape_open_error(COBIAERR_Denied);
	}
	
	//CAPEOPEN_1_2::ICapeParameter
	
	CAPEOPEN_1_2::CapeValidationStatus getValStatus() {
		return paramValidationStatus;
	}
	
	CAPEOPEN_1_2::CapeParamMode getMode() {
		return CAPEOPEN_1_2::CAPE_INPUT;
	}
	
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

using OptionParameterPtr = CapeOpenObjectSmartPointer<OptionParameter>;
