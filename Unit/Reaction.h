#pragma once
#include <COBIA.h>
#include "Helpers.h"

using namespace COBIA;

class Reaction :
	public CapeOpenObject<Reaction>,
	public CAPEOPEN_1_2::CapeIdentificationAdapter<Reaction> {

public:

	// Memebers: Identification
	CapeStringImpl reactionName, reactionDesc;

	CapeArrayIntegerImpl stoichiometry;
	const CapeReal forwarRateConstant, forward_ArrheniusEnergy, heatOfReaction;
	CapeInteger basisComponentIndex;

	CapeReactionType type;
	CapeReactionBasis basis;
	CapeReactionPhase phase;

	const CapeStringImpl getDescriptionForErrorSource() {
		return COBIATEXT("Reaction") + reactionName;
	}

	Reaction(CapeStringImpl _reactionName, CapeStringImpl _reactionDesc,
		CapeReactionType _type, CapeReactionPhase _phase,
		CapeArrayIntegerImpl _stoichiometry, CapeInteger _basisComponentIndex, CapeReactionBasis _basis,
		CapeReal _forwarRateConstant, CapeReal _forward_ArrheniusEnergy, CapeReal _heatOfReaction) :
		reactionName(_reactionName), reactionDesc(_reactionDesc), type(_type), phase(_phase),
		stoichiometry(_stoichiometry), basisComponentIndex(_basisComponentIndex), basis(_basis),
		forwarRateConstant(_forwarRateConstant), forward_ArrheniusEnergy(_forward_ArrheniusEnergy),
		heatOfReaction(_heatOfReaction) {

	}

	~Reaction() {
	}

	//CAPEOPEN_1_2::ICapeIdentification
	
	void getComponentName(/*out*/ CapeString name) {
		name = this->reactionName;
	}
	
	void putComponentName(/*in*/ CapeString name) {
		throw cape_open_error(COBIAERR_Denied);
	}
	
	void getComponentDescription(/*out*/ CapeString desc) {
		desc = this->reactionDesc;
	}
	
	void putComponentDescription(/*in*/ CapeString desc) {
		throw cape_open_error(COBIAERR_Denied);
	}
	
};

using ReactionPtr = CapeOpenObjectSmartPointer<Reaction>;
