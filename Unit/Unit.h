#pragma once
#include <COBIA.h>
#include "MaterialPort.h"
#include "OptionParameter.h"
#include "PortCollection.h"
#include "ParameterCollection.h"
#include "Validator.h"
#include "Solver.h"


#ifdef _DEBUG
	#ifdef _WIN64
		#define unitName COBIATEXT("CO MSHEX x64 Debug")
		// Class UUID = AAF02E89-291C-4D7C-836F-10EC28A705A9
		#define unitUUID 0xaa,0xf0,0x2e,0x89,0x29,0x1c,0x4d,0x7c,0x83,0x6f,0x10,0xec,0x28,0xa7,0x05,0xa9
	#else
		#define unitName COBIATEXT("CO MSHEX x86 Debug")
		// Class UUID = AAF02E89-291C-4D7C-836F-10EC28A705A8
		#define unitUUID 0xaa,0xf0,0x2e,0x89,0x29,0x1c,0x4d,0x7c,0x83,0x6f,0x10,0xec,0x28,0xa7,0x05,0xa8
	#endif
#endif

#ifndef _DEBUG
	#ifdef _WIN64
		#define unitName COBIATEXT("CO MSHEX x64")
		// Class UUID = AAF02E89-291C-4D7C-836F-10EC28A705FF
		#define unitUUID 0xaa,0xf0,0x2e,0x89,0x29,0x1c,0x4d,0x7c,0x83,0x6f,0x10,0xec,0x28,0xa7,0x05,0xff
	#else
		#define unitName COBIATEXT("CO MSHEX x86")
		// Class UUID = AAF02E89-291C-4D7C-836F-10EC28A705AA
		#define unitUUID 0xaa,0xf0,0x2e,0x89,0x29,0x1c,0x4d,0x7c,0x83,0x6f,0x10,0xec,0x28,0xa7,0x05,0xaa
	#endif
#endif

#define unitDescription COBIATEXT("MultiStream Heat Exchanger")

using namespace COBIA;

class Unit :
	public CapeOpenObject<Unit>,
	public CAPEOPEN_1_2::CapeIdentificationAdapter<Unit>,
	public CAPEOPEN_1_2::CapeUnitAdapter<Unit>,
	public CAPEOPEN_1_2::CapeUtilitiesAdapter<Unit>,
	public CAPEOPEN_1_2::CapePersistAdapter<Unit> {

	// Members: Identification
	CapeStringImpl name, description;

	// Members: Persistence and Validation
	CapeBoolean dirty;
	CAPEOPEN_1_2::CapeValidationStatus validationStatus;

	// Members: Ports and Parameters
	MaterialPortPtr in1, in2, in3, in4, in5, out1, out2, out3, out4, out5;

	CapeArrayStringImpl sideOptions;
	OptionParameterPtr in1side, in2side, in3side, in4side, in5side;

	// Members: Collections
	PortCollectionPtr portCollection;
	ParameterCollectionPtr paramCollection;

	// Flash arguments
	std::vector<CapeArrayStringImpl> productsPhaseIDs;
	std::vector<CapeArrayEnumerationImpl<CAPEOPEN_1_2::CapePhaseStatus>> productsPhaseStatus;
	CapeArrayStringImpl flashCond1, flashCond2;

public:

	// Returns a description of the current object for error handling
	const CapeStringImpl getDescriptionForErrorSource() {
		return COBIATEXT("Unit: ") + name;
	}

	// Registration info
	static const CapeUUID getObjectUUID() {
		return CapeUUID{unitUUID};
	}

	static void Register(CapePMCRegistrar registrar) {
		registrar.putName(unitName);
		registrar.putDescription(unitDescription);
		registrar.putCapeVersion(COBIATEXT("1.1"));
		registrar.putComponentVersion(COBIATEXT("0.5.0"));
		registrar.putAbout(COBIATEXT("Sample Unit Operation using COBIA."));
		registrar.putVendorURL(COBIATEXT("www.polimi.it"));
		// registrar.putProgId(COBIATEXT("Polimi.Unit"));
		registrar.addCatID(CAPEOPEN::categoryId_UnitOperation);
		registrar.addCatID(CAPEOPEN_1_2::categoryId_Component_1_2);
		//registrar.putCreationFlags(CapePMCRegistationFlag_None);
	}

	// Constructor and members initialisation
	Unit()  :
		name(unitName),
		description(unitDescription),
		validationStatus(CAPEOPEN_1_2::CAPE_NOT_VALIDATED),
		dirty(false),
		in1(new MaterialPort(name, validationStatus, COBIATEXT("Inlet 1"), true, CAPEOPEN_1_2::CAPE_INLET)),
		in2(new MaterialPort(name, validationStatus, COBIATEXT("Inlet 2"), true, CAPEOPEN_1_2::CAPE_INLET)),
		in3(new MaterialPort(name, validationStatus, COBIATEXT("Inlet 3"), false, CAPEOPEN_1_2::CAPE_INLET)),
		in4(new MaterialPort(name, validationStatus, COBIATEXT("Inlet 4"), false, CAPEOPEN_1_2::CAPE_INLET)),
		in5(new MaterialPort(name, validationStatus, COBIATEXT("Inlet 5"), false, CAPEOPEN_1_2::CAPE_INLET)),
		out1(new MaterialPort(name, validationStatus, COBIATEXT("Outlet 1"), true, CAPEOPEN_1_2::CAPE_OUTLET)),
		out2(new MaterialPort(name, validationStatus, COBIATEXT("Outlet 2"), true, CAPEOPEN_1_2::CAPE_OUTLET)),
		out3(new MaterialPort(name, validationStatus, COBIATEXT("Outlet 3"), false, CAPEOPEN_1_2::CAPE_OUTLET)),
		out4(new MaterialPort(name, validationStatus, COBIATEXT("Outlet 4"), false, CAPEOPEN_1_2::CAPE_OUTLET)),
		out5(new MaterialPort(name, validationStatus, COBIATEXT("Outlet 5"), false, CAPEOPEN_1_2::CAPE_OUTLET)),
		sideOptions(3),
		in1side(new OptionParameter(name, validationStatus, dirty, sideOptions, COBIATEXT("Inlet 1 Side"))),
		in2side(new OptionParameter(name, validationStatus, dirty, sideOptions, COBIATEXT("Inlet 2 Side"))),
		in3side(new OptionParameter(name, validationStatus, dirty, sideOptions, COBIATEXT("Inlet 3 Side"))),
		in4side(new OptionParameter(name, validationStatus, dirty, sideOptions, COBIATEXT("Inlet 4 Side"))),
		in5side(new OptionParameter(name, validationStatus, dirty, sideOptions, COBIATEXT("Inlet 5 Side"))),
		portCollection(new PortCollection(name)),
		paramCollection(new ParameterCollection(name)) {

		// Stream Side Options
		sideOptions[0] = COBIATEXT("Ignored");
		sideOptions[1] = COBIATEXT("Hot");
		sideOptions[2] = COBIATEXT("Cold");

		// Add ports to port collection
		// Inlet Stream should be followed by its outlet to reserve indexing.
		// Energy streams should be placed at the end
		portCollection->addPort(in1);
		portCollection->addPort(out1);
		portCollection->addPort(in2);
		portCollection->addPort(out2);
		portCollection->addPort(in3);
		portCollection->addPort(out3);
		portCollection->addPort(in4);
		portCollection->addPort(out4);
		portCollection->addPort(in5);
		portCollection->addPort(out5);

		// Add parameters to port collection
		paramCollection->addParameter(in1side);
		paramCollection->addParameter(in2side);
		paramCollection->addParameter(in3side);
		paramCollection->addParameter(in4side);
		paramCollection->addParameter(in5side);

		// Prepare T & P flash specifications for products flash
		// specification format:
		// CapeArrayRealImpl = { propertyIdentifier, basis, phaseLabel [, compoundIdentifier] }
		// basis is undefined when it is not a dependency of the property (e.g. T, P)
		flashCond1.resize(3);
		flashCond1[0] = COBIATEXT("temperature");
		flashCond1[2] = COBIATEXT("overall");
		flashCond2.resize(3);
		flashCond2[0] = COBIATEXT("pressure");
		flashCond2[2] = COBIATEXT("overall");
	}

	~Unit() {
	}

	// CAPEOPEN_1_2::ICapeIdentification
	void getComponentName(/*out*/ CapeString name) {
		name = this->name;
	}
	void putComponentName(/*in*/ CapeString name) {
		this->name = name;
		dirty = true;
	}
	void getComponentDescription(/*out*/ CapeString desc) {
		desc = description;
	}
	void putComponentDescription(/*in*/ CapeString desc) {
		description = desc;
		dirty = true;
	}

	// CAPEOPEN_1_2::ICapeUnit
	CAPEOPEN_1_2::CapeCollection<CAPEOPEN_1_2::CapeUnitPort> ports() {
		return portCollection;
	}
	CAPEOPEN_1_2::CapeValidationStatus getValStatus() {
		return validationStatus;
	}
	CapeBoolean Validate(/*out*/ CapeString message) {
		ValidatorPtr validator = new Validator(portCollection);
		if (validator->validateMaterialPorts(message) && validator->validateHeatExchangerInletOutlet(message)) {
			validator->preparePhaseIDs(productsPhaseIDs, productsPhaseStatus);
			validationStatus = CAPEOPEN_1_2::CAPE_VALID;
			return true;
		}
		else {
			validationStatus = CAPEOPEN_1_2::CAPE_INVALID;
			return false;
		}
	}
	void Calculate() {	
		// Check validation status before calculation
		if (validationStatus != CAPEOPEN_1_2::CAPE_VALID) {
			throw cape_open_error(COBIATEXT("Unit is not in a valid state"));
		}
		// Initiate Solver
		SolverPtr solver = new Solver(portCollection);
		solver->flashProduct(productsPhaseIDs, productsPhaseStatus, flashCond1, flashCond2);
	}

	// CAPEOPEN_1_2::ICapeUtilities
	CAPEOPEN_1_2::CapeCollection<CAPEOPEN_1_2::CapeParameter> getParameters() {
		if (paramCollection != NULL)
		{
			return paramCollection;
		}
		throw cape_open_error(COBIAERR_NotImplemented);
	}
	void putSimulationContext(/*in*/ CAPEOPEN_1_2::CapeSimulationContext context) {
		throw cape_open_error(COBIAERR_NotImplemented);
	}

	void Initialize() {
		// 1. The PME will order the PMC to get initialized through this method.
		// 2. Any initialisation that could fail must be placed here.
		// 3. Initialize is guaranteed to be the first method called by the client
		// (except low level methods such as class constructors or initialization persistence methods).
		// 4. Initialize has to be called once when the PMC is instantiated in a particular flowsheet.
		// 5. When the initialization fails, before signalling an error,
		// the PMC must free all the resources that were allocated before the failure occurred.
		// When the PME receives this error, it may not use the PMC anymore.
		// The method terminate of the current interface must not either be called.
		// Hence, the PME may only release the PMC through the middleware native mechanisms.
	}
	void Terminate() {
		// Disconnect ports
		CAPEOPEN_1_2::CapeCollection<CAPEOPEN_1_2::CapeUnitPort> collection(portCollection);
		for (CAPEOPEN_1_2::CapeUnitPort p : collection) {
			p.Disconnect();
		}
		// In case a reference to the simulation context is stored, it too must be released at Terminate
	}
	CAPEOPEN_1_2::CapeEditResult Edit(CapeWindowId parent) {
		throw cape_open_error(COBIAERR_NotImplemented);
	}

	//CAPEOPEN_1_2::ICapePersist
	void Save(/*in*/ CAPEOPEN_1_2::CapePersistWriter writer,/*in*/ CapeBoolean clearDirty) {
		writer.Add(ConstCapeString(COBIATEXT("name")), name);
		writer.Add(ConstCapeString(COBIATEXT("description")), description);
		// TODO: Iterate over parameter collection
		// and import parameter name using identification interface
		// writer.Add(ConstCapeString(COBIATEXT("param")), param->getValue());

		if (clearDirty) {
			dirty = false;
		}
	}
	void Load(/*in*/ CAPEOPEN_1_2::CapePersistReader reader) {
		reader.GetString(ConstCapeString(COBIATEXT("name")), name);
		reader.GetString(ConstCapeString(COBIATEXT("description")), description);
		// TODO: Iterate over parameter collection
		// and import parameter name using identification interface
		// param->putValue(reader.GetReal(ConstCapeString(COBIATEXT("param"))));
	}
	CapeBoolean getIsDirty() {
		return dirty;
	}
};
