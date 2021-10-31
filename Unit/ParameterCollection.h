/*
Implementation of Collection Common Interface for the developed unit:
"The ICapeCollection interface provides a means of collecting together
lists of CAPE-OPEN items/entities (eg. parameters, ports, …)."
This collection is of "CapeParameter" type. It can be converted
to a template that takes other types of entities.
*/

#pragma once
#include <COBIA.h>

using namespace COBIA;

class ParameterCollection :
	public CapeOpenObject<ParameterCollection>,
	public CAPEOPEN_1_2::CapeCollectionAdapter<CAPEOPEN_1_2::CapeParameter,ParameterCollection>,
	public CAPEOPEN_1_2::CapeIdentificationAdapter<ParameterCollection> {

	// Members
	CapeStringImpl& unitName;
	std::vector<CAPEOPEN_1_2::CapeParameter> parameters;

public:

	const CapeStringImpl getDescriptionForErrorSource() {
		return COBIATEXT("ParameterCollection of ") + unitName;
	}

	ParameterCollection(CapeStringImpl& _unitName) : unitName(_unitName) {
	}

	~ParameterCollection() {
	}

	// Method to adding parameters to the collection
	void addParameter(CAPEOPEN_1_2::CapeParameter param) {
		parameters.emplace_back(param);
	}

	//CAPEOPEN_1_2::ICapeIdentification

	void getComponentName(/*out*/ CapeString name) {
		name = COBIATEXT("Parameter collection");
	}

	void putComponentName(/*in*/ CapeString name) {
		throw cape_open_error(COBIAERR_Denied);
	}

	void getComponentDescription(/*out*/ CapeString desc) {
		desc = COBIATEXT("Parameter collection");
	}

	void putComponentDescription(/*in*/ CapeString desc) {
		throw cape_open_error(COBIAERR_Denied);
	}

	//CAPEOPEN_1_2::ICapeCollection<CAPEOPEN_1_2::ICapeParameter>
	
	// Lookup by index
	CAPEOPEN_1_2::CapeParameter Item(/*in*/ CapeInteger index) {
		if ((index < 0) || (index >= parameters.size()))
		{
			throw cape_open_error(COBIAERR_NoSuchItem);
		}
		return parameters[index];
	}
	
	// Lookup by name
	CAPEOPEN_1_2::CapeParameter Item(/*in*/ CapeString name) {
		CapeString paramName;
		for (CAPEOPEN_1_2::CapeParameter& p : parameters)
		{
			CAPEOPEN_1_2::CapeIdentification paramIdentification(p);
			paramIdentification.getComponentName(paramName);
			if (name == paramName)
			{
				return p;
			}
		}
		// If not found (no return)
		throw cape_open_error(COBIAERR_NoSuchItem);
	}
	
	CapeInteger getCount() {
		return (CapeInteger)parameters.size();
	}
};

using ParameterCollectionPtr = CapeOpenObjectSmartPointer<ParameterCollection>;
