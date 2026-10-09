#pragma once
#include "presets.h"
#include "scene.h"
#include <format>

namespace starscape::test
{

inline constexpr double kPi = 3.14159265358979323846;

inline Scene BaseScene()
{
	Scene scene;
	scene.populations = PresetPopulations();
	scene.output.faceSize = 16;
	scene.output.quality = Quality::Quick;
	scene.output.limitingMagnitude = 6.0;
	return scene;
}

// Every star of this population sits in one luminosity bin, so counts and fluxes have closed forms.
inline int AddSingleBinPopulation(Scene& scene, int bin, double kelvin = 5800.0)
{
	Population population;
	population.name = std::format("single-{}", bin);
	population.fraction[size_t(bin)] = 1.0;
	population.temperature.fill(kelvin);
	scene.populations.push_back(std::move(population));
	return int(scene.populations.size()) - 1;
}

// (1 + r / 10¹²)^-3.5 is 1 to nine digits, so this is a uniform ball.
inline Galaxy UniformBall(double density, double cutoff, int population, glm::dvec3 position)
{
	Galaxy galaxy;
	galaxy.position = position;
	galaxy.halo = {.enabled = true, .density = density, .scale = 1e12, .cutoff = cutoff, .population = population};
	return galaxy;
}

inline Galaxy UniformDust(double magnitudesPerParsec, double cutoff)
{
	Galaxy galaxy;
	galaxy.dust = {.enabled = true, .extinction = magnitudesPerParsec, .scaleLength = 1e12, .scaleHeight = 1e12, .cutoff = cutoff,
		.laneOffset = 0.0, .contrast = 1.0, .filament = 0.0};
	return galaxy;
}

inline Scene MilkyWayScene(glm::dvec3 observer)
{
	Scene scene = BaseScene();
	scene.galaxies.push_back(*PresetGalaxy("milky-way", scene.populations));
	scene.observer.position = observer;
	return scene;
}

}
