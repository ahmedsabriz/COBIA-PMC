#pragma once
#include <COBIA.h>

using namespace COBIA;

class MaterialPort :
	public CapeOpenObject<MaterialPort>,
	public CAPEOPEN_1_2::CapeIdentificationAdapter<MaterialPort>,
	public CAPEOPEN_1_2::CapeUnitPortAdapter<MaterialPort> {

	// Members
	CapeStringImpl& unitName;
	CapeStringImpl portName;

	CAPEOPEN_1_2::CapePortDirection direction;
	CAPEOPEN_1_2::CapeThermoMaterial connectedMaterial;

	CAPEOPEN_1_2::CapeValidationStatus unitValidationStatus;

public:

	const CapeStringImpl getDescriptionForErrorSource() {
		return portName + COBIATEXT(" port of ") + unitName;
	}

	MaterialPort(const COBIACHAR* _portName, CAPEOPEN_1_2::CapePortDirection _direction,
		CapeStringImpl& _unitName, CAPEOPEN_1_2::CapeValidationStatus& _unitValidationStatus) :
		portName(_portName), direction(_direction),
		unitName(_unitName), unitValidationStatus(_unitValidationStatus) {
	}

	~MaterialPort() {
	}

	// CAPEOPEN_1_2::ICapeIdentification
	
	void getComponentName(/*out*/ CapeString name) {
		name = this->portName;
	}

	void putComponentName(/*in*/ CapeString name) {
		throw cape_open_error(COBIAERR_Denied);
	}

	void getComponentDescription(/*out*/ CapeString desc) {
		desc = L"Materialport";
	}

	void putComponentDescription(/*in*/ CapeString desc) {
		throw cape_open_error(COBIAERR_Denied);
	}
	
	// CAPEOPEN_1_2::ICapeUnitPort
	
	void Connect(/*in*/ CapeInterface objectToConnect) {
		CAPEOPEN_1_2::CapeThermoMaterial newConnectedMaterial = objectToConnect;
		if (!newConnectedMaterial) {
			//expected a material object
			throw cape_open_error(COBIAERR_NoSuchInterface);
		}
		unitValidationStatus = CAPEOPEN_1_2::CAPE_NOT_VALIDATED;
		connectedMaterial = newConnectedMaterial;
	}
	void Disconnect() {
		unitValidationStatus = CAPEOPEN_1_2::CAPE_NOT_VALIDATED;
		connectedMaterial.clear();
	}

	CAPEOPEN_1_2::CapePortType getPortType() {
		return CAPEOPEN_1_2::CAPE_MATERIAL;
	}

	CAPEOPEN_1_2::CapePortDirection getDirection() {
		return direction;
	}

	CapeInterface getConnectedObject() {
		return connectedMaterial;
	}

	// ICapeThermoMaterial
	CAPEOPEN_1_2::CapeThermoMaterial getMaterial() {
		return connectedMaterial;
	}
};

using MaterialPortPtr = CapeOpenObjectSmartPointer<MaterialPort>;
