#pragma once
#include <COBIA.h>

using namespace COBIA;


// Reaction package enums
enum class CapeReactionType : CapeEnumeration {
	// Focus on single option for Proof Of Concept
	CAPE_KINETIC = 0,
};

enum class CapeReactionBasis : CapeEnumeration {
	// Focus on single option for Proof Of Concept
	CAPE_MOLE_FRACTION = 0,
};

enum class CapeReactionPhase : CapeEnumeration {
	// Focus on single option for Proof Of Concept
	CAPE_VAPOR = 0,
};

// MiltiStream Heat Exchanger side
enum class CapeStreamSide : CapeEnumeration {
	CAPE_UNSELECTED = 0,
	CAPE_HOT = 1,
	CAPE_COLD = 2,
};

CapeReal calcOverallFromPhaseProps(/*in*/ CAPEOPEN_1_2::CapeThermoMaterial thermoMaterial,
	/*in*/ CapeStringImpl prop) {

	/* This function will calculate phase properties and phase fraction for a given stream
	and return the overall value of the property */

	CapeArrayStringImpl props(1);
	props[0] = prop;
	CapeReal overallValue = 0;

	// Get present phaseLabels to calculate their properties
	CapeArrayStringImpl phaseLabels;
	CapeArrayEnumerationImpl<CAPEOPEN_1_2::CapePhaseStatus> phaseStatus;
	thermoMaterial.GetPresentPhases(phaseLabels, phaseStatus);

	// Iterate over present phases and calculate selected phase properties
	for (ConstCapeString phaseLabel : phaseLabels) {

		// Calculate phase fraction
		CapeArrayRealImpl phaseFraction;
		thermoMaterial.GetSinglePhaseProp(ConstCapeString(COBIATEXT("phaseFraction")),
			phaseLabel, ConstCapeString(COBIATEXT("mole")), phaseFraction);

		/* Before getting a property from a material object,
		property must be calculated using CalcSinglePhaseProp.
		Exceptions are for properties that are available at phase equilibrium:
		flow rate, composition, pressure, temperature and phase fraction */
		CAPEOPEN_1_2::CapeThermoPropertyRoutine routine(thermoMaterial);
		routine.CalcSinglePhaseProp(props, phaseLabel);
		
		CapeArrayRealImpl value(1);
		thermoMaterial.GetSinglePhaseProp(prop, phaseLabel,
			ConstCapeString(COBIATEXT("mole")), value);
		overallValue += phaseFraction[0] * value[0];
		return overallValue;
	}
}