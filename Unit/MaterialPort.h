#pragma once
#include <COBIA.h>
#include "Helpers.h"

using namespace COBIA;

class MaterialPort :
	public CapeOpenObject<MaterialPort>,
	public CAPEOPEN_1_2::CapeIdentificationAdapter<MaterialPort>,
	public CAPEOPEN_1_2::CapeUnitPortAdapter<MaterialPort> {

	// Members
	CapeStringImpl &unitName;
	CAPEOPEN_1_2::CapeValidationStatus &unitValidationStatus;

	CapeStringImpl portName;
	CapeBoolean primary;
	CAPEOPEN_1_2::CapePortDirection direction;	

	CAPEOPEN_1_2::CapeThermoMaterial connectedMaterial;
	CapeBoolean connected;
	CapeStreamSide side;

public:

	const CapeStringImpl getDescriptionForErrorSource() {
		return portName + COBIATEXT(" port of ") + unitName;
	}

	MaterialPort(CapeStringImpl& _unitName, CAPEOPEN_1_2::CapeValidationStatus& _unitValidationStatus,
		const COBIACHAR* _portName, CapeBoolean _primary,
		CAPEOPEN_1_2::CapePortDirection _direction) :
		unitName(_unitName), unitValidationStatus(_unitValidationStatus),
		portName(_portName), primary(_primary), direction(_direction) {
		connected = false;
		side = CapeStreamSide::CAPE_UNSELECTED;
	}

	~MaterialPort() {
	}

	CapeBoolean isPrimary() {
		return primary;
	}

	CapeBoolean isConnected() {
		return connected;
	}

	void ignoreOptionalStream() {
		this->connected = false;
	}

	CAPEOPEN_1_2::CapeThermoMaterial getMaterial() {
		return connectedMaterial;
	}

	// CAPEOPEN_1_2::ICapeIdentification
	void getComponentName(/*out*/ CapeString name) {
		name = this->portName;
	}
	void putComponentName(/*in*/ CapeString name) {
		throw cape_open_error(COBIAERR_Denied);
	}
	void getComponentDescription(/*out*/ CapeString desc) {
		desc = COBIATEXT("Material Port");
	}
	void putComponentDescription(/*in*/ CapeString desc) {
		throw cape_open_error(COBIAERR_Denied);
	}
	
	// CAPEOPEN_1_2::ICapeUnitPort
	CAPEOPEN_1_2::CapePortType getPortType() {
		return CAPEOPEN_1_2::CAPE_MATERIAL;
	}
	CAPEOPEN_1_2::CapePortDirection getDirection() {
		return direction;
	}
	CapeInterface getConnectedObject() {
		return connectedMaterial;
	}
	void Connect(/*in*/ CapeInterface objectToConnect) {
		CAPEOPEN_1_2::CapeThermoMaterial newConnectedMaterial = objectToConnect;
		if (!newConnectedMaterial) {
			//expected a material object
			throw cape_open_error(COBIAERR_NoSuchInterface);
		}
		unitValidationStatus = CAPEOPEN_1_2::CAPE_NOT_VALIDATED;
		connectedMaterial = newConnectedMaterial;
		this->connected = true;
	}
	void Disconnect() {
		unitValidationStatus = CAPEOPEN_1_2::CAPE_NOT_VALIDATED;
		connectedMaterial.clear();
		this->connected = false;
	}
};

using MaterialPortPtr = CapeOpenObjectSmartPointer<MaterialPort>;