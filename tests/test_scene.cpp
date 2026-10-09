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
