#include "cook.h"
#include "march.h"
#include "preview.h"
#include "scene.h"
#include "writer.h"
#include <kre_starscape.h>
#include <algorithm>
#include <array>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <format>
#include <optional>
#include <string>
#include <string_view>
#include <system_error>
#include <thread>
#include <vector>

using namespace starscape;

namespace
{

constexpr char const* kUsage =
	"usage:\n"
	"  starscape cook    <scene.toml> -o <out.kre_starscape> [--preview <out.hdr>] [--threads N] [--quick]\n"
	"  starscape preview <in.kre_starscape> -o <out.hdr> [--width N] [--exposure EV] [--threads N]\n"
	"  starscape inspect <in.kre_starscape>\n";

auto const gStart = std::chrono::steady_clock::now();

void Log(std::string_view message)
{
	auto const elapsed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - gStart).count();
	std::printf("[+%02lld:%02lld:%02lld] %.*s\n", (long long)(elapsed / 3600), (long long)(elapsed / 60 % 60), (long long)(elapsed % 60),
		int(message.size()), message.data());
	std::fflush(stdout);
}

int Report(Error const& error)
{
	std::fprintf(stderr, "starscape: %s\n", error.message.c_str());
	return error.kind == Error::Kind::Io ? 2 : 1;
}

int Usage(std::string_view problem)
{
	std::fprintf(stderr, "starscape: %.*s\n%s", int(problem.size()), problem.data(), kUsage);
	return 1;
}

struct Arguments
{
	std::vector<std::string_view> positional;
	std::optional<std::string_view> output;
	std::optional<std::string_view> preview;
	unsigned threads = std::max(1u, std::thread::hardware_concurrency());
	bool quick = false;
	std::optional<uint32_t> width;
	double exposure = 0.0;
};

template<class T>
bool ParseNumber(std::string_view text, T& out)
{
	char const* end = text.data() + text.size();
	auto const [pointer, error] = std::from_chars(text.data(), end, out);
	return error == std::errc() && pointer == end;
}

std::optional<std::string> ParseArguments(int argc, char** argv, Arguments& out)
{
	for (int i = 2; i < argc; ++i)
	{
		std::string_view const arg = argv[i];
		auto value = [&]() -> std::optional<std::string_view>
		{
			if (i + 1 >= argc)
			{
				return std::nullopt;
			}
			return std::string_view(argv[++i]);
		};
		if (arg == "--quick")
		{
			out.quick = true;
		}
		else if (arg == "-o" || arg == "--preview" || arg == "--threads" || arg == "--width" || arg == "--exposure")
		{
			std::optional<std::string_view> const v = value();
			if (!v)
			{
				return std::format("{} needs a value", arg);
			}
			if (arg == "-o")
			{
				out.output = *v;
			}
			else if (arg == "--preview")
			{
				out.preview = *v;
			}
			else if (arg == "--threads")
			{
				if (!ParseNumber(*v, out.threads) || out.threads == 0)
				{
					return "--threads needs a positive integer";
				}
			}
			else if (arg == "--width")
			{
				uint32_t width = 0;
				if (!ParseNumber(*v, width) || width < 16 || width % 2 != 0)
				{
					return "--width needs an even integer >= 16";
				}
				out.width = width;
			}
			else
			{
				if (!ParseNumber(*v, out.exposure) || !std::isfinite(out.exposure))
				{
					return "--exposure needs a number";
				}
			}
		}
		else if (arg.starts_with("-"))
		{
			return std::format("unknown option {}", arg);
		}
		else
		{
			out.positional.push_back(arg);
		}
	}
	return std::nullopt;
}

Expected<std::vector<uint8_t>> LoadSky(std::string_view path, kre_starscape::View& view)
{
	STARSCAPE_TRY_ASSIGN(bytes, ReadFileBytes(std::string(path)));
	// `view` points into `bytes`, so the caller keeps the returned buffer alive while using it.
	kre_starscape::Result const result = kre_starscape::Parse(bytes, view);
	if (result != kre_starscape::Result::Ok)
	{
		return InputError(std::format("{}: {}", path, kre_starscape::Describe(result)));
	}
	return bytes;
}

int Cook(Arguments const& args)
{
	if (args.positional.size() != 1 || !args.output)
	{
		return Usage("cook needs one scene file and -o");
	}
	Expected<Scene> const scene = LoadScene(std::string(args.positional[0]));
	if (!scene)
	{
		return Report(scene.error());
	}
	Log(std::format("cooking {} on {} threads", args.positional[0], args.threads));
	int lastPercent = -1;
	std::string lastPhase;
	CookOptions options;
	options.threads = args.threads;
	if (args.quick)
	{
		options.quality = Quality::Quick;
	}
	auto const marchStart = std::chrono::steady_clock::now();
	options.progress = [&](std::string_view phase, double done)
	{
		int const percent = int(done * 100.0);
		if (phase == lastPhase && percent == lastPercent)
		{
			return;
		}
		lastPhase = std::string(phase);
		lastPercent = percent;
		double const seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - marchStart).count();
		std::string eta;
		if (phase == "march" && done > 0.01)
		{
			eta = std::format(", ETA {:.0f} s", seconds / done * (1.0 - done));
		}
		Log(std::format("{} {}%{}", phase, percent, eta));
	};
	Cooked const cooked = starscape::Cook(*scene, options);
	std::vector<uint8_t> const bytes = BuildFile(cooked, scene->hash);
	if (Expected<void> const written = WriteFileAtomic(std::string(*args.output), bytes); !written)
	{
		return Report(written.error());
	}

	double const total = cooked.stats.pointFlux + cooked.stats.diffuseFlux;
	Log(std::format("wrote {} ({} bytes)", *args.output, bytes.size()));
	Log(std::format("threshold F_res = {:.4g} (magnitude {:.2f}); {} stars; points carry {:.2f}% of the light",
		cooked.fluxThreshold, -2.5 * std::log10(cooked.fluxThreshold), cooked.stats.starCount,
		total > 0.0 ? 100.0 * cooked.stats.pointFlux / total : 0.0));
	if (cooked.stats.budgetUnreachable)
	{
		Log("warning: the scene holds fewer stars than target_stars; every star is a point");
	}
	if (cooked.stats.droppedStars > 0)
	{
		Log(std::format("warning: {} drawn stars found nothing resolved on re-march and were dropped", cooked.stats.droppedStars));
	}

	if (args.preview)
	{
		kre_starscape::View view;
		kre_starscape::Parse(bytes, view);
		Image const image = scene->output.previewMarch
			? RenderMarchedPreview(PrepareScene(*scene), ParamsFor(options.quality.value_or(scene->output.quality)), view,
				scene->output.previewWidth, 0.0, args.threads)
			: RenderPreview(view, scene->output.previewWidth, 0.0, args.threads);
		if (Expected<void> const written = WriteFileAtomic(std::string(*args.preview), EncodeRadiance(image)); !written)
		{
			return Report(written.error());
		}
		Log(std::format("wrote preview {} ({}x{})", *args.preview, image.width, image.height));
	}
	return 0;
}

int Preview(Arguments const& args)
{
	if (args.positional.size() != 1 || !args.output)
	{
		return Usage("preview needs one .kre_starscape file and -o");
	}
	kre_starscape::View view;
	Expected<std::vector<uint8_t>> const bytes = LoadSky(args.positional[0], view);
	if (!bytes)
	{
		return Report(bytes.error());
	}
	Image const image = RenderPreview(view, args.width.value_or(8192), args.exposure, args.threads);
	if (Expected<void> const written = WriteFileAtomic(std::string(*args.output), EncodeRadiance(image)); !written)
	{
		return Report(written.error());
	}
	Log(std::format("wrote preview {} ({}x{})", *args.output, image.width, image.height));
	return 0;
}

int Inspect(Arguments const& args)
{
	if (args.positional.size() != 1)
	{
		return Usage("inspect needs one .kre_starscape file");
	}
	kre_starscape::View view;
	Expected<std::vector<uint8_t>> const bytes = LoadSky(args.positional[0], view);
	if (!bytes)
	{
		return Report(bytes.error());
	}
	kre_starscape::Header const& h = view.header;
	double const scale = h.valueScale;
	std::printf("scene hash      %016llx\n", (unsigned long long)h.sceneHash);
	std::printf("diffuse         RGB9E5 cubemap, %u^2 per face, %u mips, %llu bytes\n", h.diffuseFaceSize, h.diffuseMipCount, (unsigned long long)h.diffuseBytes);
	std::printf("stars           %u\n", h.starCount);
	std::printf("flux threshold  %.4g (magnitude %.2f)\n", h.fluxThreshold, -2.5 * std::log10(h.fluxThreshold));
	std::printf("value scale     %g\n", scale);

	double diffuse = 0.0;
	for (uint32_t face = 0; face < 6; ++face)
	{
		for (uint32_t y = 0; y < h.diffuseFaceSize; ++y)
		{
			for (uint32_t x = 0; x < h.diffuseFaceSize; ++x)
			{
				auto const v = kre_starscape::DecodeRGB9E5(kre_starscape::DiffuseTexel(view, 0, face, x, y));
				diffuse += Luminance(glm::dvec3(v[0], v[1], v[2])) / scale * kre_starscape::TexelSolidAngle(h.diffuseFaceSize, x, y);
			}
		}
	}
	double points = 0.0;
	std::array<uint32_t, 40> histogram{};   // magnitude bins from -10 to +30
	for (uint32_t i = 0; i < h.starCount; ++i)
	{
		auto const v = kre_starscape::DecodeRGB9E5(kre_starscape::StarAt(view, i).fluxRGB9E5);
		double const flux = Luminance(glm::dvec3(v[0], v[1], v[2])) / scale;
		points += flux;
		if (flux > 0.0)
		{
			int const bin = std::clamp(int(std::floor(-2.5 * std::log10(flux))) + 10, 0, 39);
			++histogram[size_t(bin)];
		}
	}
	std::printf("total light     %.6g (points %.2f%%, diffuse %.2f%%)\n", points + diffuse,
		100.0 * points / std::max(points + diffuse, 1e-300), 100.0 * diffuse / std::max(points + diffuse, 1e-300));
	std::printf("stars by magnitude:\n");
	for (size_t bin = 0; bin < histogram.size(); ++bin)
	{
		if (histogram[bin] > 0)
		{
			std::printf("  %+3d .. %+3d  %u\n", int(bin) - 10, int(bin) - 9, histogram[bin]);
		}
	}
	return 0;
}

}

int main(int argc, char** argv)
{
	if (argc < 2)
	{
		return Usage("no command");
	}
	Arguments args;
	if (std::optional<std::string> const problem = ParseArguments(argc, argv, args))
	{
		return Usage(*problem);
	}
	std::string_view const command = argv[1];
	if (command == "cook")
	{
		return Cook(args);
	}
	if (command == "preview")
	{
		return Preview(args);
	}
	if (command == "inspect")
	{
		return Inspect(args);
	}
	return Usage(std::format("unknown command {}", command));
}
