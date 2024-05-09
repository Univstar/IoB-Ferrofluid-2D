#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <fmt/chrono.h>
#include <fmt/color.h>
#include <fmt/core.h>

#include <spdlog/spdlog.h>

#include <tbb/tbb.h>

#include <yaml-cpp/yaml.h>

#include <array>
#include <atomic>
#include <chrono>
#include <concepts>
#include <filesystem>
#include <fstream>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <memory>
#include <numbers>
#include <span>
#include <string>
#include <string_view>
#include <thread>
#include <type_traits>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>

namespace YAML {
	template <typename Derived, int Rows>
	struct convert<Eigen::Matrix<Derived, Rows, 1>> {
		static Node encode(Eigen::Matrix<Derived, Rows, 1> const &rhs) {
			Node node;
			for (std::size_t i = 0; i < rhs.size(); i++) {
				node.push_back(rhs(i));
			}
			node.SetStyle(EmitterStyle::Flow);
			return node;
		}
	};
}

namespace Pivot {
	using namespace Eigen;

	template <typename Type>
	inline Type Zero() {
		if constexpr (requires { { Type::Zero() } -> std::convertible_to<Type>; }) {
			return Type::Zero();
		} else {
			return 0;
		}
	}

	class IO {
	public:
		template <typename T>
		static void Read(std::istream &in, T &val) {
			if constexpr (requires { std::span(val); }) {
				auto buf = std::span(val);
				in.read(reinterpret_cast<char *>(buf.data()), buf.size_bytes());
			} else {
				in.read(reinterpret_cast<char *>(&val), sizeof(val));
			}
		}

		template <typename T>
		static void Write(std::ostream &out, T const &val) {
			if constexpr (requires { std::span(val); }) {
				auto buf = std::span(val);
				out.write(reinterpret_cast<char const *>(buf.data()), buf.size_bytes());
			} else {
				out.write(reinterpret_cast<char const *>(&val), sizeof(val));
			}
		}
	};

	struct CellTypeBits {
		enum : std::uint8_t {
			Air    = 1 << 0,
			Liquid = 1 << 1,
			Solid  = 1 << 2,
		};
	};
}
