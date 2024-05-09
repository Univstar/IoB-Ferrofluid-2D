#include "Image.h"

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

namespace Pivot {
	void Image::WriteBytes(std::filesystem::path const &filename, std::span<std::byte const> bytes, int width, int height, int numChannels, bool flipped) {
		stbi_flip_vertically_on_write(flipped);
		if (filename.extension() == ".bmp") {
			stbi_write_bmp(filename.string().c_str(), width, height, numChannels, bytes.data());
		} else if (filename.extension() == ".jpg") {
			stbi_write_jpg(filename.string().c_str(), width, height, numChannels, bytes.data(), 100);
		} else {
			stbi_write_png(filename.string().c_str(), width, height, numChannels, bytes.data(), width * numChannels);
		}
	}
}
