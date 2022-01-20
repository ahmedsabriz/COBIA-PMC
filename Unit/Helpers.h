#pragma once
#include <COBIA.h>

using namespace COBIA;

// Get component name
CapeStringImpl getName(CapeInterface param) {
	CapeStringImpl paramName;
	CAPEOPEN_1_2::CapeIdentification identification(param);
	identification.getComponentName(paramName);
	return paramName;
}