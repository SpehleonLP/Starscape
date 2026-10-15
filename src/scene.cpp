#include "scene.h"
#include "presets.h"
#include <toml.hpp>
#include <algorithm>
#include <bit>
#include <cmath>
#include <format>
#include <fstream>
#include <limits>
#include <sstream>
#include <glm/geometric.hpp>

namespace starscape
{
namespace
{

uint64_t Fnv1a64(std::string_view bytes)
{
	uint64_t hash = 0xCBF29CE484222325ull;
	for (char const c : bytes)
	{
		hash ^= uint8_t(c);
		hash *= 0x100000001B3ull;
	}
	return hash;
}

constexpr double kHuge = 1e12;
// Positions and angles beyond this overflow squared distances in the galaxy math.
constexpr double kVectorLimit = 1e10;

// Reads one TOML table, remembering every key it was asked for so Finish() can reject the rest:
// a misspelt key must not silently fall back to a preset value.
class TableReader
{
public:
	TableReader(toml::table const& table, std::string_view source, std::string&& path)
	: m_table(table)
	, m_source(source)
	, m_path(std::move(path))
	{
	}

	std::string_view Source() const { return m_source; }
	std::string const& Path() const { return m_path; }

	std::unexpected<Error> Fail(toml::source_region const& where, std::string_view key, std::string_view message) const
	{
		std::string name = m_path;
		if (!key.empty())
		{
			name = name.empty() ? std::string(key) : std::format("{}.{}", name, key);
		}
		if (name.empty())
		{
			name = "scene";
		}
		return InputError(std::format("{}:{}:{}: {}: {}", m_source, where.begin.line, where.begin.column, name, message));
	}

	std::unexpected<Error> Fail(std::string_view key, std::string_view message) const
	{
		return Fail(m_table.source(), key, message);
	}

	bool Has(std::string_view key) const { return m_table.contains(key); }

	Expected<void> Number(std::string_view key, double& out, double min, double max)
	{
		toml::node const* node = Take(key);
		if (node == nullptr)
		{
			return {};
		}
		std::optional<double> const value = node->is_boolean() ? std::nullopt : node->value<double>();
		if (!value)
		{
			return Fail(node->source(), key, "expected a number");
		}
		if (!std::isfinite(*value) || *value < min || *value > max)
		{
			return Fail(node->source(), key, std::format("must be in [{}, {}]", min, max));
		}
		out = *value;
		return {};
	}

	Expected<void> Integer(std::string_view key, int64_t& out, int64_t min, int64_t max)
	{
		toml::node const* node = Take(key);
		if (node == nullptr)
		{
			return {};
		}
		std::optional<int64_t> const value = node->value_exact<int64_t>();
		if (!value)
		{
			return Fail(node->source(), key, "expected an integer");
		}
		if (*value < min || *value > max)
		{
			return Fail(node->source(), key, std::format("must be in [{}, {}]", min, max));
		}
		out = *value;
		return {};
	}

	Expected<void> Boolean(std::string_view key, bool& out)
	{
		toml::node const* node = Take(key);
		if (node == nullptr)
		{
			return {};
		}
		std::optional<bool> const value = node->value_exact<bool>();
		if (!value)
		{
			return Fail(node->source(), key, "expected true or false");
		}
		out = *value;
		return {};
	}

	Expected<void> Text(std::string_view key, std::string& out)
	{
		toml::node const* node = Take(key);
		if (node == nullptr)
		{
			return {};
		}
		std::optional<std::string> value = node->value_exact<std::string>();
		if (!value)
		{
			return Fail(node->source(), key, "expected a string");
		}
		out = std::move(*value);
		return {};
	}

	Expected<void> Vector(std::string_view key, glm::dvec3& out)
	{
		toml::node const* node = Take(key);
		if (node == nullptr)
		{
			return {};
		}
		toml::array const* array = node->as_array();
		if (array == nullptr || array->size() != 3)
		{
			return Fail(node->source(), key, "expected an array of three numbers");
		}
		for (size_t i = 0; i < 3; ++i)
		{
			toml::node const& element = (*array)[i];
			std::optional<double> const value = element.is_boolean() ? std::nullopt : element.value<double>();
			if (!value || !std::isfinite(*value))
			{
				return Fail(element.source(), key, "expected an array of three numbers");
			}
			if (std::abs(*value) > kVectorLimit)
			{
				return Fail(element.source(), key, std::format("components must be within +-{}", kVectorLimit));
			}
			out[int(i)] = *value;
		}
		return {};
	}

	Expected<toml::table const*> Table(std::string_view key)
	{
		toml::node const* node = Take(key);
		if (node == nullptr)
		{
			return nullptr;
		}
		if (!node->is_table())
		{
			return Fail(node->source(), key, "expected a table");
		}
		return node->as_table();
	}

	Expected<toml::array const*> Array(std::string_view key)
	{
		toml::node const* node = Take(key);
		if (node == nullptr)
		{
			return nullptr;
		}
		if (!node->is_array())
		{
			return Fail(node->source(), key, "expected an array");
		}
		return node->as_array();
	}

	Expected<void> Finish() const
	{
		for (auto const& [key, node] : m_table)
		{
			if (std::find(m_taken.begin(), m_taken.end(), key.str()) == m_taken.end())
			{
				return Fail(node.source(), key.str(), "unknown key");
			}
		}
		return {};
	}

private:
	toml::node const* Take(std::string_view key)
	{
		m_taken.emplace_back(key);
		return m_table.get(key);
	}

	toml::table const&       m_table;
	std::string_view         m_source;
	std::string              m_path;
	std::vector<std::string> m_taken;
};

Expected<TableReader> ChildReader(TableReader& parent, std::string_view key, toml::table const& table)
{
	std::string path = parent.Path().empty() ? std::string(key) : std::format("{}.{}", parent.Path(), key);
	return TableReader(table, parent.Source(), std::move(path));
}

Expected<void> ReadPopulation(TableReader& reader, std::string_view key, std::span<Population const> populations, int& out)
{
	std::string name;
	STARSCAPE_TRY(reader.Text(key, name));
	if (name.empty())
	{
		return {};
	}
	std::optional<int> const found = FindPopulation(populations, name);
	if (!found)
	{
		return reader.Fail(key, std::format("unknown population '{}'", name));
	}
	out = *found;
	return {};
}

Expected<void> ReadArms(TableReader& parent, Arms& arms)
{
	STARSCAPE_TRY_ASSIGN(table, parent.Table("arms"));
	if (table == nullptr)
	{
		return {};
	}
	STARSCAPE_TRY_ASSIGN(reader, ChildReader(parent, "arms", *table));
	int64_t count = arms.count;
	STARSCAPE_TRY(reader.Integer("count", count, 0, 16));
	arms.count = int(count);
	STARSCAPE_TRY(reader.Number("pitch_deg", arms.pitchDeg, 1.0, 89.0));
	STARSCAPE_TRY(reader.Number("width", arms.width, 1e-3, kHuge));
	STARSCAPE_TRY(reader.Number("contrast", arms.contrast, 0.0, 1000.0));
	STARSCAPE_TRY(reader.Number("clumping", arms.clumping, 0.0, 10.0));
	STARSCAPE_TRY(reader.Number("phase_deg", arms.phaseDeg, -360.0, 360.0));
	STARSCAPE_TRY(reader.Number("r0", arms.r0, 1e-3, kHuge));
	return reader.Finish();
}

Expected<void> ReadColour(TableReader& reader, std::string_view key, glm::dvec3& colour)
{
	STARSCAPE_TRY(reader.Vector(key, colour));
	// Normalisation divides by the luminance, so a black or negative colour must stop here.
	if (colour.x < 0.0 || colour.y < 0.0 || colour.z < 0.0 || !(Luminance(colour) > 0.0))
	{
		return reader.Fail(key, "needs non-negative components and a positive luminance");
	}
	return {};
}

// An array of colours, each checked like `colour`: normalisation divides by the luminance.
Expected<void> ReadPalette(TableReader& reader, Gas& gas)
{
	STARSCAPE_TRY_ASSIGN(array, reader.Array("palette"));
	if (array == nullptr)
	{
		return {};
	}
	if (array->size() == 1 || array->size() > size_t(kMaxPalette))
	{
		return reader.Fail(array->source(), "palette", std::format("needs 0 or 2 to {} colours", kMaxPalette));
	}
	for (size_t i = 0; i < array->size(); ++i)
	{
		toml::node const& element = (*array)[i];
		std::string const key = std::format("palette[{}]", i);
		toml::array const* rgb = element.as_array();
		if (rgb == nullptr || rgb->size() != 3)
		{
			return reader.Fail(element.source(), key, "expected an array of three numbers");
		}
		glm::dvec3 colour(0.0);
		for (size_t c = 0; c < 3; ++c)
		{
			std::optional<double> const value = (*rgb)[c].is_boolean() ? std::nullopt : (*rgb)[c].value<double>();
			if (!value || !std::isfinite(*value))
			{
				return reader.Fail((*rgb)[c].source(), key, "expected an array of three numbers");
			}
			colour[int(c)] = *value;
		}
		if (colour.x < 0.0 || colour.y < 0.0 || colour.z < 0.0 || !(Luminance(colour) > 0.0))
		{
			return reader.Fail(element.source(), key, "needs non-negative components and a positive luminance");
		}
		gas.palette[i] = colour;
	}
	gas.paletteCount = int(array->size());
	return {};
}

Expected<void> ReadGas(TableReader& parent, Gas& gas)
{
	STARSCAPE_TRY_ASSIGN(table, parent.Table("gas"));
	if (table == nullptr)
	{
		return {};
	}
	STARSCAPE_TRY_ASSIGN(reader, ChildReader(parent, "gas", *table));
	STARSCAPE_TRY(reader.Number("glow", gas.glow, 0.0, 1000.0));
	STARSCAPE_TRY(ReadColour(reader, "colour", gas.colour));
	STARSCAPE_TRY(ReadColour(reader, "core_colour", gas.coreColour));
	STARSCAPE_TRY(ReadPalette(reader, gas));
	STARSCAPE_TRY(reader.Number("palette_wavelength", gas.paletteWavelength, 10.0, 100000.0));
	STARSCAPE_TRY(reader.Number("large_glow", gas.largeGlow, 0.0, 1000.0));
	STARSCAPE_TRY(reader.Number("large_wavelength", gas.largeWavelength, 10.0, 100000.0));
	return reader.Finish();
}

// A preset leaves unused components with no population and no cutoff; enabling one from a scene
// must not let that reach the sampler as an out-of-range population index.
Expected<void> CheckUsable(TableReader const& reader, std::string_view name, bool enabled, int population, double cutoff)
{
	if (!enabled)
	{
		return {};
	}
	if (population < 0)
	{
		return reader.Fail("", std::format("{} is enabled but has no population", name));
	}
	if (cutoff <= 0.0)
	{
		return reader.Fail("", std::format("{} is enabled but has no cutoff", name));
	}
	return {};
}

Expected<void> ReadDisk(TableReader& parent, std::string_view key, Disk& disk, std::span<Population const> populations, bool thin)
{
	STARSCAPE_TRY_ASSIGN(table, parent.Table(key));
	if (table == nullptr)
	{
		return {};
	}
	STARSCAPE_TRY_ASSIGN(reader, ChildReader(parent, key, *table));
	disk.enabled = true;
	STARSCAPE_TRY(reader.Boolean("enabled", disk.enabled));
	STARSCAPE_TRY(reader.Number("density", disk.density, 0.0, kHuge));
	STARSCAPE_TRY(reader.Number("scale_length", disk.scaleLength, 1e-3, kHuge));
	STARSCAPE_TRY(reader.Number("scale_height", disk.scaleHeight, 1e-3, kHuge));
	STARSCAPE_TRY(reader.Number("cutoff", disk.cutoff, 0.0, kHuge));
	STARSCAPE_TRY(ReadPopulation(reader, "population", populations, disk.population));
	if (thin)
	{
		STARSCAPE_TRY(ReadPopulation(reader, "arm_population", populations, disk.armPopulation));
		STARSCAPE_TRY(reader.Number("arm_young_fraction", disk.armYoungFraction, 0.0, 1.0));
		STARSCAPE_TRY(reader.Number("arm_scale_height", disk.armScaleHeight, 0.0, kHuge));
		STARSCAPE_TRY(reader.Number("young_fraction", disk.youngFraction, 0.0, 1.0));
		STARSCAPE_TRY(reader.Number("association_size", disk.associationSize, 0.0, kHuge));
		STARSCAPE_TRY(reader.Number("association_clumping", disk.associationClumping, 0.0, 10.0));
		STARSCAPE_TRY(ReadArms(reader, disk.arms));
		if (disk.enabled && (disk.armYoungFraction > 0.0 || disk.youngFraction > 0.0) && disk.armPopulation < 0)
		{
			return reader.Fail("arm_young_fraction", "needs arm_population to be set");
		}
		STARSCAPE_TRY(ReadGas(reader, disk.gas));
		if (disk.enabled && disk.gas.glow > 0.0 && disk.armPopulation < 0)
		{
			return reader.Fail("gas", "glow needs arm_population to be set: young stars light the gas");
		}
		if (disk.enabled && disk.gas.largeGlow > 0.0 && disk.armPopulation < 0)
		{
			return reader.Fail("gas", "large_glow needs arm_population to be set: young stars light the gas");
		}
	}
	STARSCAPE_TRY(CheckUsable(reader, key, disk.enabled, disk.population, disk.cutoff));
	return reader.Finish();
}

Expected<void> ReadBulge(TableReader& parent, Bulge& bulge, std::span<Population const> populations)
{
	STARSCAPE_TRY_ASSIGN(table, parent.Table("bulge"));
	if (table == nullptr)
	{
		return {};
	}
	STARSCAPE_TRY_ASSIGN(reader, ChildReader(parent, "bulge", *table));
	bulge.enabled = true;
	STARSCAPE_TRY(reader.Boolean("enabled", bulge.enabled));
	STARSCAPE_TRY(reader.Number("density", bulge.density, 0.0, kHuge * 1e6));
	STARSCAPE_TRY(reader.Number("scale", bulge.scale, 1e-3, kHuge));
	STARSCAPE_TRY(reader.Number("axis_y", bulge.axisY, 1e-3, 1.0));
	STARSCAPE_TRY(reader.Number("axis_z", bulge.axisZ, 1e-3, 1.0));
	STARSCAPE_TRY(reader.Number("bar_angle_deg", bulge.barAngleDeg, -180.0, 180.0));
	STARSCAPE_TRY(reader.Number("cutoff", bulge.cutoff, 0.0, kHuge));
	STARSCAPE_TRY(ReadPopulation(reader, "population", populations, bulge.population));
	STARSCAPE_TRY(CheckUsable(reader, "bulge", bulge.enabled, bulge.population, bulge.cutoff));
	return reader.Finish();
}

Expected<void> ReadHalo(TableReader& parent, Halo& halo, std::span<Population const> populations)
{
	STARSCAPE_TRY_ASSIGN(table, parent.Table("halo"));
	if (table == nullptr)
	{
		return {};
	}
	STARSCAPE_TRY_ASSIGN(reader, ChildReader(parent, "halo", *table));
	halo.enabled = true;
	STARSCAPE_TRY(reader.Boolean("enabled", halo.enabled));
	STARSCAPE_TRY(reader.Number("density", halo.density, 0.0, kHuge));
	STARSCAPE_TRY(reader.Number("scale", halo.scale, 1e-3, kHuge));
	STARSCAPE_TRY(reader.Number("cutoff", halo.cutoff, 0.0, kHuge));
	STARSCAPE_TRY(ReadPopulation(reader, "population", populations, halo.population));
	STARSCAPE_TRY(CheckUsable(reader, "halo", halo.enabled, halo.population, halo.cutoff));
	return reader.Finish();
}

Expected<void> ReadBubble(TableReader& parent, Dust& dust)
{
	STARSCAPE_TRY_ASSIGN(table, parent.Table("bubble"));
	if (table == nullptr)
	{
		return {};
	}
	STARSCAPE_TRY_ASSIGN(reader, ChildReader(parent, "bubble", *table));
	STARSCAPE_TRY(reader.Vector("center", dust.bubbleCenter));
	STARSCAPE_TRY(reader.Number("radius", dust.bubbleRadius, 0.0, kHuge));
	STARSCAPE_TRY(reader.Number("height", dust.bubbleHeight, 0.0, kHuge));
	// Either semi-axis left at zero would make a bubble that silently is not there, or divides by zero.
	if (dust.bubbleRadius <= 0.0 || dust.bubbleHeight <= 0.0)
	{
		return reader.Fail("", "needs a positive radius and height");
	}
	return reader.Finish();
}

Expected<void> ReadDust(TableReader& parent, Dust& dust)
{
	STARSCAPE_TRY_ASSIGN(table, parent.Table("dust"));
	if (table == nullptr)
	{
		return {};
	}
	STARSCAPE_TRY_ASSIGN(reader, ChildReader(parent, "dust", *table));
	dust.enabled = true;
	STARSCAPE_TRY(reader.Boolean("enabled", dust.enabled));
	STARSCAPE_TRY(reader.Number("extinction", dust.extinction, 0.0, 10.0));
	STARSCAPE_TRY(reader.Number("scale_length", dust.scaleLength, 1e-3, kHuge));
	STARSCAPE_TRY(reader.Number("scale_height", dust.scaleHeight, 1e-3, kHuge));
	STARSCAPE_TRY(reader.Number("cutoff", dust.cutoff, 0.0, kHuge));
	STARSCAPE_TRY(reader.Number("lane_offset", dust.laneOffset, -5.0, 5.0));
	STARSCAPE_TRY(reader.Number("contrast", dust.contrast, 0.0, 1000.0));
	STARSCAPE_TRY(reader.Number("filament", dust.filament, 0.0, 10.0));
	STARSCAPE_TRY(ReadBubble(reader, dust));
	STARSCAPE_TRY(CheckUsable(reader, "dust", dust.enabled, 0, dust.cutoff));
	return reader.Finish();
}

Expected<Galaxy> ReadGalaxy(toml::table const& table, std::string_view source, size_t index, std::span<Population const> populations)
{
	TableReader reader(table, source, std::format("galaxy[{}]", index));
	std::string preset;
	STARSCAPE_TRY(reader.Text("preset", preset));
	if (preset.empty())
	{
		return reader.Fail("preset", "every galaxy needs a preset");
	}
	std::optional<Galaxy> base = PresetGalaxy(preset, populations);
	if (!base)
	{
		return reader.Fail("preset", std::format("unknown preset '{}'", preset));
	}
	Galaxy galaxy = std::move(*base);
	STARSCAPE_TRY(reader.Vector("position", galaxy.position));
	STARSCAPE_TRY(reader.Vector("rotation", galaxy.rotationDeg));
	int64_t seed = int64_t(galaxy.seed);
	STARSCAPE_TRY(reader.Integer("seed", seed, 0, std::numeric_limits<int64_t>::max()));
	galaxy.seed = uint64_t(seed);
	STARSCAPE_TRY(ReadDisk(reader, "disk", galaxy.disk, populations, true));
	STARSCAPE_TRY(ReadDisk(reader, "thick_disk", galaxy.thickDisk, populations, false));
	STARSCAPE_TRY(ReadBulge(reader, galaxy.bulge, populations));
	STARSCAPE_TRY(ReadHalo(reader, galaxy.halo, populations));
	STARSCAPE_TRY(ReadDust(reader, galaxy.dust));
	STARSCAPE_TRY(reader.Finish());
	return galaxy;
}

Expected<void> ReadOutput(TableReader& root, Output& output)
{
	STARSCAPE_TRY_ASSIGN(table, root.Table("output"));
	if (table == nullptr)
	{
		return root.Fail("output", "missing [output] table");
	}
	STARSCAPE_TRY_ASSIGN(reader, ChildReader(root, "output", *table));
	int64_t faceSize = output.faceSize;
	STARSCAPE_TRY(reader.Integer("face_size", faceSize, 1, 4096));
	if (!std::has_single_bit(uint64_t(faceSize)))
	{
		return reader.Fail("face_size", "must be a power of two");
	}
	output.faceSize = uint32_t(faceSize);
	if (reader.Has("target_stars") == reader.Has("limiting_magnitude"))
	{
		return reader.Fail("", "set exactly one of target_stars and limiting_magnitude");
	}
	if (reader.Has("target_stars"))
	{
		int64_t target = 0;
		STARSCAPE_TRY(reader.Integer("target_stars", target, 1, 10000000));
		output.targetStars = uint32_t(target);
	}
	else
	{
		double magnitude = 0.0;
		STARSCAPE_TRY(reader.Number("limiting_magnitude", magnitude, -30.0, 30.0));
		output.limitingMagnitude = magnitude;
	}
	std::string quality;
	STARSCAPE_TRY(reader.Text("quality", quality));
	if (quality == "quick")
	{
		output.quality = Quality::Quick;
	}
	else if (quality.empty() || quality == "final")
	{
		output.quality = Quality::Final;
	}
	else
	{
		return reader.Fail("quality", "must be \"final\" or \"quick\"");
	}
	int64_t seed = int64_t(output.seed);
	STARSCAPE_TRY(reader.Integer("seed", seed, 0, std::numeric_limits<int64_t>::max()));
	output.seed = uint64_t(seed);
	STARSCAPE_TRY_ASSIGN(preview, reader.Table("preview"));
	if (preview != nullptr)
	{
		STARSCAPE_TRY_ASSIGN(previewReader, ChildReader(reader, "preview", *preview));
		int64_t width = output.previewWidth;
		STARSCAPE_TRY(previewReader.Integer("width", width, 16, 65536));
		if (width % 2 != 0)
		{
			return previewReader.Fail("width", "must be even (the preview is width × width/2)");
		}
		output.previewWidth = uint32_t(width);
		STARSCAPE_TRY(previewReader.Boolean("march", output.previewMarch));
		STARSCAPE_TRY(previewReader.Finish());
	}
	return reader.Finish();
}

Expected<void> ReadObserver(TableReader& root, Observer& observer)
{
	STARSCAPE_TRY_ASSIGN(table, root.Table("observer"));
	if (table == nullptr)
	{
		return root.Fail("observer", "missing [observer] table");
	}
	STARSCAPE_TRY_ASSIGN(reader, ChildReader(root, "observer", *table));
	if (!reader.Has("position"))
	{
		return reader.Fail("position", "the observer needs a position");
	}
	STARSCAPE_TRY(reader.Vector("position", observer.position));
	STARSCAPE_TRY(reader.Vector("spin_axis", observer.spinAxis));
	STARSCAPE_TRY(reader.Vector("zero_meridian", observer.zeroMeridian));
	STARSCAPE_TRY(reader.Number("s_min", observer.sMin, 1e-6, 1e6));
	if (glm::length(observer.spinAxis) < 1e-6)
	{
		return reader.Fail("spin_axis", "must have length of at least 1e-6");
	}
	double const meridianLength = glm::length(observer.zeroMeridian);
	if (meridianLength < 1e-6)
	{
		return reader.Fail("zero_meridian", "must have length of at least 1e-6");
	}
	glm::dvec3 const z = glm::normalize(observer.spinAxis);
	glm::dvec3 const across = observer.zeroMeridian - glm::dot(observer.zeroMeridian, z) * z;
	if (glm::length(across) < 1e-6 * meridianLength)
	{
		return reader.Fail("zero_meridian", "must not be parallel to spin_axis");
	}
	return reader.Finish();
}

Expected<void> ReadStyle(TableReader& root, Style& style)
{
	STARSCAPE_TRY_ASSIGN(table, root.Table("style"));
	if (table == nullptr)
	{
		return {};
	}
	STARSCAPE_TRY_ASSIGN(reader, ChildReader(root, "style", *table));
	STARSCAPE_TRY(reader.Number("star_saturation", style.starSaturation, 0.0, 4.0));
	return reader.Finish();
}

}

glm::dmat3 SkyToScene(Observer const& observer)
{
	glm::dvec3 const z = glm::normalize(observer.spinAxis);
	glm::dvec3 const x = glm::normalize(observer.zeroMeridian - glm::dot(observer.zeroMeridian, z) * z);
	glm::dvec3 const y = glm::cross(z, x);
	return glm::dmat3(x, y, z);
}

Expected<Scene> ParseScene(std::string_view text, std::string_view sourceName)
{
	toml::parse_result parsed = toml::parse(text, sourceName);
	if (!parsed)
	{
		toml::parse_error const& error = parsed.error();
		return InputError(std::format("{}:{}:{}: {}", sourceName, error.source().begin.line, error.source().begin.column, error.description()));
	}
	toml::table const& root = parsed.table();

	Scene scene;
	scene.populations = PresetPopulations();
	scene.hash = Fnv1a64(text);

	TableReader reader(root, sourceName, "");
	STARSCAPE_TRY(ReadOutput(reader, scene.output));
	STARSCAPE_TRY(ReadObserver(reader, scene.observer));
	STARSCAPE_TRY(ReadStyle(reader, scene.style));
	STARSCAPE_TRY_ASSIGN(galaxies, reader.Array("galaxy"));
	if (galaxies == nullptr || galaxies->empty())
	{
		return reader.Fail("galaxy", "the scene needs at least one [[galaxy]]");
	}
	for (size_t i = 0; i < galaxies->size(); ++i)
	{
		toml::table const* table = (*galaxies)[i].as_table();
		if (table == nullptr)
		{
			return reader.Fail((*galaxies)[i].source(), "galaxy", "expected [[galaxy]] tables");
		}
		STARSCAPE_TRY_ASSIGN(galaxy, ReadGalaxy(*table, sourceName, i, scene.populations));
		scene.galaxies.push_back(std::move(galaxy));
	}
	STARSCAPE_TRY(reader.Finish());
	return scene;
}

Expected<Scene> LoadScene(std::filesystem::path const& path)
{
	std::ifstream file(path, std::ios::binary);
	if (!file)
	{
		return IoError(std::format("{}: cannot open scene file", path.string()));
	}
	std::ostringstream text;
	text << file.rdbuf();
	if (file.bad())
	{
		return IoError(std::format("{}: read failed", path.string()));
	}
	return ParseScene(text.str(), path.string());
}

}
