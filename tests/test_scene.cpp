#include "scene.h"
#include <doctest.h>
#include <glm/geometric.hpp>
#include <string>

using namespace starscape;

namespace
{

constexpr char const* kMinimal = R"(
[output]
face_size = 64
target_stars = 1000

[observer]
position = [8200.0, 0.0, 20.0]

[[galaxy]]
preset = "milky-way"

  [galaxy.disk]
  scale_height = 250
)";

std::string ErrorOf(std::string_view text)
{
	Expected<Scene> const scene = ParseScene(text, "test.toml");
	REQUIRE_FALSE(scene.has_value());
	return scene.error().message;
}

}

TEST_CASE("A minimal scene takes its preset and applies overrides")
{
	Expected<Scene> const scene = ParseScene(kMinimal, "test.toml");
	REQUIRE(scene.has_value());
	CHECK(scene->output.faceSize == 64);
	CHECK(scene->output.targetStars == 1000u);
	CHECK_FALSE(scene->output.limitingMagnitude.has_value());
	CHECK(scene->output.quality == Quality::Final);
	REQUIRE(scene->galaxies.size() == 1);
	CHECK(scene->galaxies[0].disk.scaleLength == 2600.0);   // from the preset
	CHECK(scene->galaxies[0].disk.scaleHeight == 250.0);    // overridden
	CHECK(scene->observer.position.x == 8200.0);
	CHECK(scene->populations.size() == 4);
	CHECK(scene->hash != 0);
	CHECK(scene->hash == ParseScene(kMinimal, "other.toml")->hash);
}

TEST_CASE("Unknown keys are errors with a position")
{
	// Mutation: disabling unknown-key rejection reddens this test.
	std::string const message = ErrorOf("[output]\nface_sise = 256\ntarget_stars = 10\n[observer]\nposition = [0.0, 0.0, 0.0]\n[[galaxy]]\npreset = \"dwarf\"\n");
	CHECK(message.find("test.toml:2:") != std::string::npos);
	CHECK(message.find("output.face_sise") != std::string::npos);
	CHECK(message.find("unknown key") != std::string::npos);
}

TEST_CASE("Exactly one of target_stars and limiting_magnitude")
{
	CHECK(ErrorOf("[output]\ntarget_stars = 10\nlimiting_magnitude = 6.0\n[observer]\nposition = [0.0, 0.0, 0.0]\n[[galaxy]]\npreset = \"dwarf\"\n").find("exactly one") != std::string::npos);
	CHECK(ErrorOf("[output]\n[observer]\nposition = [0.0, 0.0, 0.0]\n[[galaxy]]\npreset = \"dwarf\"\n").find("exactly one") != std::string::npos);
}

TEST_CASE("Values are range- and type-checked")
{
	char const* tail = "[observer]\nposition = [0.0, 0.0, 0.0]\n[[galaxy]]\npreset = \"dwarf\"\n";
	CHECK(ErrorOf(std::string("[output]\ntarget_stars = 10\nface_size = 300\n") + tail).find("power of two") != std::string::npos);
	CHECK(ErrorOf(std::string("[output]\ntarget_stars = 10\nface_size = \"big\"\n") + tail).find("expected an integer") != std::string::npos);
	CHECK(ErrorOf(std::string("[output]\ntarget_stars = 10\nquality = \"best\"\n") + tail).find("quality") != std::string::npos);
	CHECK(ErrorOf("[output]\ntarget_stars = 10\n[observer]\nposition = [0.0, 0.0, 0.0]\n[[galaxy]]\npreset = \"andromeda\"\n").find("unknown preset") != std::string::npos);
	CHECK(ErrorOf("[output]\ntarget_stars = 10\n[observer]\nposition = [0.0, 0.0, 0.0]\n[[galaxy]]\npreset = \"dwarf\"\n[galaxy.bulge]\npopulation = \"nebula\"\n").find("unknown population") != std::string::npos);
	CHECK(ErrorOf("[output]\ntarget_stars = 10\n[observer]\nposition = [0.0, 0.0, 0.0]\nspin_axis = [1.0, 0.0, 0.0]\nzero_meridian = [2.0, 0.0, 0.0]\n[[galaxy]]\npreset = \"dwarf\"\n").find("parallel") != std::string::npos);
	CHECK(ErrorOf("[output]\ntarget_stars = 10\n[observer]\nposition = [0.0, 0.0, 0.0]\n").find("galaxy") != std::string::npos);
	CHECK(ErrorOf("[output]\ntarget_stars = 10\n[[galaxy]]\npreset = \"dwarf\"\n").find("observer") != std::string::npos);
}

TEST_CASE("A component table enables its component")
{
	Expected<Scene> const scene = ParseScene(
		"[output]\nlimiting_magnitude = 6.0\n[observer]\nposition = [0.0, 0.0, 0.0]\n"
		"[[galaxy]]\npreset = \"dwarf\"\n[galaxy.thick_disk]\ndensity = 0.1\nscale_length = 1000\nscale_height = 300\ncutoff = 5000\npopulation = \"old-disk\"\n"
		"[galaxy.halo]\nenabled = false\n",
		"test.toml");
	REQUIRE(scene.has_value());
	CHECK(scene->galaxies[0].thickDisk.enabled);
	CHECK(scene->galaxies[0].thickDisk.population == 1);
	CHECK_FALSE(scene->galaxies[0].halo.enabled);
}

TEST_CASE("The sky frame is orthonormal with +z on the spin axis")
{
	Observer observer;
	observer.spinAxis = {0.0, -0.46, 0.89};
	observer.zeroMeridian = {1.0, 0.3, 0.0};
	glm::dmat3 const sky = SkyToScene(observer);
	CHECK(glm::length(sky[0]) == doctest::Approx(1.0));
	CHECK(glm::length(sky[1]) == doctest::Approx(1.0));
	CHECK(glm::dot(sky[0], sky[1]) == doctest::Approx(0.0));
	CHECK(glm::dot(sky[0], sky[2]) == doctest::Approx(0.0));
	CHECK(glm::dot(sky[2], glm::normalize(observer.spinAxis)) == doctest::Approx(1.0));
	CHECK(glm::dot(glm::cross(sky[0], sky[1]), sky[2]) == doctest::Approx(1.0));   // right-handed
	CHECK(glm::dot(sky[1], sky[2]) == doctest::Approx(0.0));
	CHECK(glm::dot(sky[0], glm::normalize(observer.zeroMeridian)) > 0.0);   // a flipped meridian keeps every other property
}

namespace
{

constexpr char const* kHead = "[output]\ntarget_stars = 10\n";
constexpr char const* kObs = "[observer]\nposition = [0.0, 0.0, 0.0]\n";

}

TEST_CASE("Scale lengths and arm pitch are range-checked")
{
	CHECK(ErrorOf(std::string(kHead) + kObs + "[[galaxy]]\npreset = \"milky-way\"\n[galaxy.disk]\nscale_length = -5.0\n").find("scale_length") != std::string::npos);
	CHECK(ErrorOf(std::string(kHead) + kObs + "[[galaxy]]\npreset = \"milky-way\"\n[galaxy.disk]\narms = {pitch_deg = 0.0}\n").find("pitch_deg") != std::string::npos);
	CHECK(ErrorOf(std::string(kHead) + kObs + "[[galaxy]]\npreset = \"milky-way\"\n[galaxy.disk]\narms = {pitch_deg = 90.0}\n").find("pitch_deg") != std::string::npos);
}

TEST_CASE("Observer axes must be usable")
{
	CHECK(ErrorOf(std::string(kHead) + "[observer]\nposition = [0.0, 0.0, 0.0]\nspin_axis = [0.0, 0.0, 0.0]\n[[galaxy]]\npreset = \"dwarf\"\n").find("spin_axis") != std::string::npos);
	CHECK(ErrorOf(std::string(kHead) + "[observer]\nposition = [0.0, 0.0, 0.0]\nspin_axis = [0.0, 0.0, 1.0]\nzero_meridian = [0.0, 0.0, -3.0]\n[[galaxy]]\npreset = \"dwarf\"\n").find("parallel") != std::string::npos);
	CHECK(ErrorOf(std::string(kHead) + "[observer]\nposition = [0.0, 0.0, 0.0]\nzero_meridian = [0.0, 0.0, 0.0]\n[[galaxy]]\npreset = \"dwarf\"\n").find("zero_meridian") != std::string::npos);
}

TEST_CASE("Vectors are bounded")
{
	CHECK(ErrorOf(std::string(kHead) + "[observer]\nposition = [1e300, 0.0, 0.0]\n[[galaxy]]\npreset = \"dwarf\"\n").find("within") != std::string::npos);
	CHECK(ErrorOf(std::string(kHead) + kObs + "[[galaxy]]\npreset = \"dwarf\"\nposition = [0.0, 1e300, 0.0]\n").find("within") != std::string::npos);
	CHECK(ErrorOf(std::string(kHead) + kObs + "[[galaxy]]\npreset = \"dwarf\"\nrotation = [0.0, 0.0, 1e11]\n").find("within") != std::string::npos);
}

TEST_CASE("Unknown keys inside nested tables are errors")
{
	CHECK(ErrorOf(std::string(kHead) + kObs + "[[galaxy]]\npreset = \"milky-way\"\n[galaxy.disk]\nscale_hieght = 3.0\n").find("unknown key") != std::string::npos);
	CHECK(ErrorOf(std::string(kHead) + kObs + "[[galaxy]]\npreset = \"milky-way\"\n[galaxy.disk]\narms = {cuont = 2}\n").find("arms.cuont") != std::string::npos);
}

TEST_CASE("An enabled component must be usable")
{
	std::string const tail = std::string(kHead) + kObs + "[[galaxy]]\npreset = \"dwarf\"\n";
	CHECK(ErrorOf(tail + "[galaxy.thick_disk]\ndensity = 0.1\ncutoff = 5000\n").find("no population") != std::string::npos);
	CHECK(ErrorOf(tail + "[galaxy.thick_disk]\ndensity = 0.1\npopulation = \"old-disk\"\n").find("no cutoff") != std::string::npos);
	CHECK(ErrorOf(tail + "[galaxy.dust]\nextinction = 0.1\n").find("no cutoff") != std::string::npos);
	CHECK(ErrorOf(std::string(kHead) + kObs + "[[galaxy]]\npreset = \"dwarf\"\n[galaxy.disk]\narm_young_fraction = 0.5\n").find("arm_population") != std::string::npos);
}

TEST_CASE("A dust bubble is read, and needs both semi-axes")
{
	std::string const head = std::string(kHead) + kObs + "[[galaxy]]\npreset = \"milky-way\"\n[galaxy.dust]\n";
	Expected<Scene> const scene = ParseScene(head + "bubble = { center = [8200.0, 0.0, 20.0], radius = 300.0, height = 800.0 }\n", "test.toml");
	REQUIRE(scene.has_value());
	Dust const& dust = scene->galaxies[0].dust;
	CHECK(dust.bubbleCenter == glm::dvec3(8200.0, 0.0, 20.0));
	CHECK(dust.bubbleRadius == 300.0);
	CHECK(dust.bubbleHeight == 800.0);
	CHECK(ErrorOf(head + "bubble = { center = [8200.0, 0.0, 20.0], radius = 300.0 }\n").find("dust.bubble: needs a positive radius and height") != std::string::npos);
	CHECK(ErrorOf(head + "bubble = { radius = 300.0, height = 800.0, size = 2.0 }\n").find("bubble.size") != std::string::npos);
}

TEST_CASE("Style is read and range-checked")
{
	// Mutation: a default star_saturation other than 1 reddens the first check.
	std::string const scene = std::string(kHead) + kObs + "[[galaxy]]\npreset = \"dwarf\"\n";
	CHECK(ParseScene(scene, "test.toml")->style.starSaturation == 1.0);
	CHECK(ParseScene(scene + "[style]\nstar_saturation = 1.5\n", "test.toml")->style.starSaturation == 1.5);
	CHECK(ErrorOf(scene + "[style]\nstar_saturation = 5.0\n").find("style.star_saturation") != std::string::npos);
	CHECK(ErrorOf(scene + "[style]\nstar_saturation = -1.0\n").find("must be in") != std::string::npos);
	CHECK(ErrorOf(scene + "[style]\nsaturation = 1.5\n").find("unknown key") != std::string::npos);
}

TEST_CASE("Gas keys are read and checked")
{
	// Mutation: deleting the arm_population check in ReadDisk reddens the last check.
	std::string const head = std::string(kHead) + kObs + "[[galaxy]]\npreset = \"milky-way\"\n[galaxy.disk.gas]\n";
	Expected<Scene> const scene = ParseScene(head + "glow = 1.5\ncolour = [1.0, 0.2, 0.3]\ncore_colour = [0.2, 1.0, 0.8]\n", "test.toml");
	REQUIRE(scene.has_value());
	Gas const& gas = scene->galaxies[0].disk.gas;
	CHECK(gas.glow == 1.5);
	CHECK(gas.colour == glm::dvec3(1.0, 0.2, 0.3));
	CHECK(gas.coreColour == glm::dvec3(0.2, 1.0, 0.8));
	CHECK(ParseScene(std::string(kHead) + kObs + "[[galaxy]]\npreset = \"milky-way\"\n", "test.toml")->galaxies[0].disk.gas.glow == 0.0);
	CHECK(ErrorOf(head + "glow = -1.0\n").find("disk.gas.glow") != std::string::npos);
	CHECK(ErrorOf(head + "glow = 1001.0\n").find("must be in") != std::string::npos);
	CHECK(ErrorOf(head + "colour = [-0.1, 1.0, 1.0]\n").find("disk.gas.colour") != std::string::npos);
	CHECK(ErrorOf(head + "core_colour = [0.0, 0.0, 0.0]\n").find("positive luminance") != std::string::npos);
	CHECK(ErrorOf(head + "colour = \"red\"\n").find("three numbers") != std::string::npos);
	CHECK(ErrorOf(head + "glw = 1.0\n").find("gas.glw") != std::string::npos);
	CHECK(ErrorOf(std::string(kHead) + kObs + "[[galaxy]]\npreset = \"milky-way\"\n[galaxy.thick_disk.gas]\nglow = 1.0\n").find("unknown key") != std::string::npos);
	CHECK(ErrorOf(std::string(kHead) + kObs + "[[galaxy]]\npreset = \"dwarf\"\n[galaxy.disk]\ndensity = 0.1\ncutoff = 5000\npopulation = \"old-disk\"\n[galaxy.disk.gas]\nglow = 1.0\n")
		.find("arm_population") != std::string::npos);
}

TEST_CASE("A gas palette is read and checked")
{
	// Mutation: deleting the luminance check in ReadPalette lets a black entry through to a division by zero.
	std::string const head = std::string(kHead) + kObs + "[[galaxy]]\npreset = \"milky-way\"\n[galaxy.disk.gas]\nglow = 1.0\n";
	Expected<Scene> const scene = ParseScene(head + "palette = [[1.0, 0.2, 0.3], [0.2, 1.0, 0.8], [1.0, 0.7, 0.2]]\npalette_wavelength = 900.0\n", "test.toml");
	REQUIRE(scene.has_value());
	Gas const& gas = scene->galaxies[0].disk.gas;
	CHECK(gas.paletteCount == 3);
	CHECK(gas.palette[2] == glm::dvec3(1.0, 0.7, 0.2));
	CHECK(gas.paletteWavelength == 900.0);
	CHECK(ParseScene(std::string(kHead) + kObs + "[[galaxy]]\npreset = \"milky-way\"\n", "test.toml")->galaxies[0].disk.gas.paletteCount == 0);
	CHECK(ParseScene(head + "palette = []\n", "test.toml")->galaxies[0].disk.gas.paletteCount == 0);
	CHECK(ErrorOf(head + "palette = [[1.0, 0.2, 0.3]]\n").find("0 or 2 to 8") != std::string::npos);
	CHECK(ErrorOf(head + "palette = [[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0]]\n").find("0 or 2 to 8") != std::string::npos);
	CHECK(ErrorOf(head + "palette = [[1.0, 0.2, 0.3], [0.0, 0.0, 0.0]]\n").find("palette[1]") != std::string::npos);
	CHECK(ErrorOf(head + "palette = [[1.0, 0.2, 0.3], [-1.0, 1.0, 1.0]]\n").find("positive luminance") != std::string::npos);
	CHECK(ErrorOf(head + "palette = [[1.0, 0.2, 0.3], [1.0, 1.0]]\n").find("three numbers") != std::string::npos);
	CHECK(ErrorOf(head + "palette = [1.0, 0.2, 0.3]\n").find("three numbers") != std::string::npos);
	CHECK(ErrorOf(head + "palette = \"warm\"\n").find("expected an array") != std::string::npos);
	CHECK(ErrorOf(head + "palette_wavelength = 5.0\n").find("disk.gas.palette_wavelength") != std::string::npos);
	CHECK(ErrorOf(std::string(kHead) + kObs + "[[galaxy]]\npreset = \"milky-way\"\n[galaxy.thick_disk.gas]\npalette = []\n").find("unknown key") != std::string::npos);
}

TEST_CASE("Large-scale gas keys are read and checked")
{
	// Mutation: deleting the arm_population check for large_glow in ReadDisk reddens the last check.
	std::string const head = std::string(kHead) + kObs + "[[galaxy]]\npreset = \"milky-way\"\n[galaxy.disk.gas]\n";
	Expected<Scene> const scene = ParseScene(head + "large_glow = 1.5\nlarge_wavelength = 120.0\n", "test.toml");
	REQUIRE(scene.has_value());
	CHECK(scene->galaxies[0].disk.gas.largeGlow == 1.5);
	CHECK(scene->galaxies[0].disk.gas.largeWavelength == 120.0);
	CHECK(ParseScene(std::string(kHead) + kObs + "[[galaxy]]\npreset = \"milky-way\"\n", "test.toml")->galaxies[0].disk.gas.largeGlow == 0.0);
	CHECK(ErrorOf(head + "large_glow = -1.0\n").find("disk.gas.large_glow") != std::string::npos);
	CHECK(ErrorOf(head + "large_wavelength = 200000.0\n").find("must be in") != std::string::npos);
	CHECK(ErrorOf(std::string(kHead) + kObs + "[[galaxy]]\npreset = \"dwarf\"\n[galaxy.disk]\ndensity = 0.1\ncutoff = 5000\npopulation = \"old-disk\"\n[galaxy.disk.gas]\nlarge_glow = 1.0\n")
		.find("large_glow needs arm_population") != std::string::npos);
}

TEST_CASE("The bubble's rim is read and needs young stars")
{
	// Mutation: deleting the arm_population check in ReadBubble reddens the last check.
	std::string const head = std::string(kHead) + kObs + "[[galaxy]]\npreset = \"milky-way\"\n[galaxy.dust]\n";
	Expected<Scene> const scene = ParseScene(head + "bubble = { center = [8200.0, 0.0, 20.0], radius = 300.0, height = 800.0, rim = 4.0 }\n", "test.toml");
	REQUIRE(scene.has_value());
	CHECK(scene->galaxies[0].dust.bubbleRim == 4.0);
	CHECK(ParseScene(head + "bubble = { center = [8200.0, 0.0, 20.0], radius = 300.0, height = 800.0 }\n", "test.toml")->galaxies[0].dust.bubbleRim == 0.0);
	CHECK(ErrorOf(head + "bubble = { center = [0.0, 0.0, 0.0], radius = 300.0, height = 800.0, rim = 1001.0 }\n").find("dust.bubble.rim") != std::string::npos);
	CHECK(ErrorOf(std::string(kHead) + kObs + "[[galaxy]]\npreset = \"dwarf\"\n[galaxy.disk]\ndensity = 0.1\ncutoff = 5000\npopulation = \"old-disk\"\n[galaxy.dust]\nbubble = { center = [0.0, 0.0, 0.0], radius = 300.0, height = 800.0, rim = 1.0 }\n")
		.find("rim needs the thin disk's arm_population") != std::string::npos);
}
