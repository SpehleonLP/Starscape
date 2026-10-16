#pragma once
#include <glm/mat3x3.hpp>
#include <glm/vec3.hpp>
#include "field/galaxy.glsl.h"

namespace starscape
{

struct PreparedGalaxy
{
	Galaxy     galaxy;
	glm::dmat3 toLocal{1.0};
	glm::dvec3 position{0.0};
	double     boundRadius = 0.0;
};

double BoundRadius(Galaxy const& galaxy);
PreparedGalaxy PrepareGalaxy(Galaxy const& galaxy);

inline void Evaluate(PreparedGalaxy const& prepared, glm::dvec3 scenePosition, GalaxySample& out)
{
	EvaluateLocal(prepared.galaxy, prepared.toLocal * (scenePosition - prepared.position), out);
}

}
