/*
Implementation of Collection Common Interface for the developed unit:
"The ICapeCollection interface provides a means of collecting together
lists of CAPE-OPEN items/entities (eg. parameters, ports, …)."
This collection is of "CapeUnitPort" type. It can be converted
to a template that takes other types of entities.
*/

#pragma once
#include <COBIA.h>
#include <vector>

using namespace COBIA;

class PortCollection :
	public CapeOpenObject<PortCollection>,
	public CAPEOPEN_1_2::CapeIdentificationAdapter<PortCollection>,
	public CAPEOPEN_1_2::CapeCollectionAdapter<CAPEOPEN_1_2::CapeUnitPort,PortCollection> {

	// Members
	CapeStringImpl& unitName;
	std::vector<CAPEOPEN_1_2::CapeUnitPort> ports;

public:

	const CapeStringImpl getDescriptionForErrorSource() {
		return COBIATEXT("PortCollection of ") + unitName;
	}

	PortCollection(CapeStringImpl& _unitName) : unitName(_unitName) {
	}

	~PortCollection() {
	}

	// Method to adding ports to the collection
	void addPort(CAPEOPEN_1_2::CapeUnitPort port) {
		ports.emplace_back(port);
	}

	// Method to removing ports from the collection is missing

	// CAPEOPEN_1_2::ICapeIdentification

	void getComponentName(/*out*/ CapeString name) {
		name = COBIATEXT("Port collection");
	}

	void putComponentName(/*in*/ CapeString name) {
		throw cape_open_error(COBIAERR_Denied);
	}

	void getComponentDescription(/*out*/ CapeString desc) {
		desc = COBIATEXT("Port collection");
	}

	void putComponentDescription(/*in*/ CapeString desc) {
		throw cape_open_error(COBIAERR_Denied);
	}
	
	// CAPEOPEN_1_2::ICapeCollection<CAPEOPEN_1_2::ICapeUnitPort>
	
	// Lookup by index
	CAPEOPEN_1_2::CapeUnitPort Item(/*in*/ CapeInteger index) {
		if ((index < 0) || (index >= (CapeInteger)ports.size())) {
			throw cape_open_error(COBIAERR_NoSuchItem);
		}
		return ports[index];
	}

	// Lookup by name
	CAPEOPEN_1_2::CapeUnitPort Item(/*in*/ CapeString name) {
		CapeString portName;
		for (CAPEOPEN_1_2::CapeUnitPort& portPtr : ports) {
			CAPEOPEN_1_2::CapeIdentification portIdentification(portPtr);
			portIdentification.getComponentName(portName);
			if (portName == name)
			{
				return portPtr;
			}
		}
		// If not found (no return)
		throw cape_open_error(COBIAERR_NoSuchItem);
	}

	CapeInteger getCount() {
		return (CapeInteger)ports.size();
	}
};

using PortCollectionPtr = CapeOpenObjectSmartPointer<PortCollection>;
