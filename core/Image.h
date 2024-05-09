#pragma once

#include <filesystem>
#include <span>

#include <cstddef>

namespace Pivot {
	class Image {
	public:
		static void WriteBytes(std::filesystem::path const &filename, std::span<std::byte const> bytes, int width, int height, int numChannels, bool flipped);
	};
}
