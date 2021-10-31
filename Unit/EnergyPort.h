#pragma once
#include <COBIA.h>
#include <vector>

using namespace COBIA;

class EnergyPort :
	public CapeOpenObject<EnergyPort>,
	public CAPEOPEN_1_2::CapeIdentificationAdapter<EnergyPort>,
	public CAPEOPEN_1_2::CapeUnitPortAdapter<EnergyPort>,
	public CAPEOPEN_1_2::CapeCollectionAdapter<CAPEOPEN_1_2::CapeParameter,EnergyPort> {

	// Members
	CapeStringImpl& unitName;
	CapeStringImpl portName;

	CAPEOPEN_1_2::CapePortDirection direction;
	CAPEOPEN_1_2::CapeCollection<CAPEOPEN_1_2::CapeParameter> connectedEnergyObject;
	
	CAPEOPEN_1_2::CapeValidationStatus unitValidationStatus;

	std::vector<CAPEOPEN_1_2::CapeParameter> parameters;

public:

	const CapeStringImpl getDescriptionForErrorSource() {
		return portName + COBIATEXT(" port of ") + unitName;
	}

	EnergyPort(const COBIACHAR* _portName, CAPEOPEN_1_2::CapePortDirection _direction,
		CapeStringImpl& _unitName, CAPEOPEN_1_2::CapeValidationStatus& _unitValidationStatus) :
		portName(_portName), direction(_direction),
		unitName(_unitName), unitValidationStatus(_unitValidationStatus) {
	}

	~EnergyPort() {
	}

	//CAPEOPEN_1_2::ICapeIdentification
	
	void getComponentName(/*out*/ CapeString name) {
		name = this->portName;
	}

	void putComponentName(/*in*/ CapeString name) {
		throw cape_open_error(COBIAERR_Denied);
	}

	void getComponentDescription(/*out*/ CapeString desc) {
		desc = L"EnergyPort";
	}

	void putComponentDescription(/*in*/ CapeString desc) {
		throw cape_open_error(COBIAERR_Denied);
	}

	// Method to adding parameters to the collection
	void addParameter(CAPEOPEN_1_2::CapeParameter param) {
		parameters.emplace_back(param);
	}
	
	// Method to expose the energy stream as a real parameter collection
	CAPEOPEN_1_2::CapeCollection<CAPEOPEN_1_2::CapeRealParameter> getCollection() {
		return (CAPEOPEN_1_2::CapeCollection<CAPEOPEN_1_2::CapeRealParameter>) connectedEnergyObject;
	}

	//CAPEOPEN_1_2::ICapeUnitPort
	
	CAPEOPEN_1_2::CapePortType getPortType() {
		return CAPEOPEN_1_2::CAPE_ENERGY;
	}
	
	CAPEOPEN_1_2::CapePortDirection getDirection() {
		return direction;
	}
	
	CapeInterface getConnectedObject() {
		return connectedEnergyObject;
	}
	
	void Connect(/*in*/ CapeInterface objectToConnect) {
		CAPEOPEN_1_2::CapeCollection<CAPEOPEN_1_2::CapeParameter> newEnergyObject = (ICapeInterface*)objectToConnect;
		// CAPEOPEN_1_2::CapeParameter zero = newEnergyObject[0];
		if (!newEnergyObject) {
			//expected an energy object
			throw cape_open_error(COBIAERR_NoSuchInterface);
		}
		unitValidationStatus = CAPEOPEN_1_2::CAPE_NOT_VALIDATED;
		connectedEnergyObject = newEnergyObject;
	}
	
	void Disconnect() {
		unitValidationStatus = CAPEOPEN_1_2::CAPE_NOT_VALIDATED;
		connectedEnergyObject.clear();
	}
	
	//CAPEOPEN_1_2::ICapeCollection<CAPEOPEN_1_2::ICapeParameter>
	
	CAPEOPEN_1_2::CapeParameter Item(/*in*/ CapeInteger index) {
		if ((index < 0) || (index >= parameters.size()))
		{
			throw cape_open_error(COBIAERR_NoSuchItem);
		}
		return parameters[index];
	}
	
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

using EnergyPortPtr = CapeOpenObjectSmartPointer<EnergyPort>;
