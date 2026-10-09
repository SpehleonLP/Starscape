#pragma once
#include <expected>
#include <string>
#include <utility>

namespace starscape
{

struct Error
{
	// Input = bad scene or file contents (exit 1); Io = the filesystem refused (exit 2).
	enum class Kind { Input, Io };
	Kind        kind = Kind::Input;
	std::string message;
};

template<class T>
using Expected = std::expected<T, Error>;

inline std::unexpected<Error> InputError(std::string&& message)
{
	return std::unexpected(Error{Error::Kind::Input, std::move(message)});
}

inline std::unexpected<Error> IoError(std::string&& message)
{
	return std::unexpected(Error{Error::Kind::Io, std::move(message)});
}

}

#define STARSCAPE_TRY(expr) \
	do \
	{ \
		if (auto starscapeTry_ = (expr); !starscapeTry_) \
		{ \
			return std::unexpected(std::move(starscapeTry_.error())); \
		} \
	} while (0)

#define STARSCAPE_TRY_ASSIGN(name, expr) \
	auto name##Result_ = (expr); \
	if (!name##Result_) \
	{ \
		return std::unexpected(std::move(name##Result_.error())); \
	} \
	auto name = std::move(*name##Result_)
