#pragma once

#include "Common.h"

namespace Pivot {
	class Grid {
	public:
		Grid(double spacing, Vector2i const &size, Vector2d const &origin) :
			m_Spacing { spacing },
			m_InvSpacing { 1 / m_Spacing },
			m_Size { size },
			m_Origin { origin } {
		}

		bool operator==(Grid const &rhs) const { return m_Spacing == rhs.m_Spacing && m_Size == rhs.m_Size && m_Origin == rhs.m_Origin; }

		double   GetSpacing   () const { return m_Spacing; }
		double   GetInvSpacing() const { return m_InvSpacing; }
		Vector2i GetSize      () const { return m_Size; }
		Vector2d GetOrigin    () const { return m_Origin; }

		int GetNumVertices() const { return m_Size.prod(); }

		bool IsInside(Vector2i const &coord, int offset = 0) const { return (coord.array() >= offset).all() && (coord.array() < m_Size.array() - offset).all(); }
		bool IsValid (Vector2i const &coord)                 const { return IsInside(coord, 0); }

		Vector2i Clamp     (Vector2i const &coord) const { return coord.cwiseMax(0).cwiseMin(m_Size - Vector2i::Ones()); }
		int      IndexOf   (Vector2i const &coord) const { return coord.y() + m_Size.y() * coord.x(); }
		Vector2d PositionOf(Vector2i const &coord) const { return m_Origin + coord.cast<double>() * m_Spacing; }
		
		Vector2i CoordOf(int index) const { return { index / m_Size.y(), index % m_Size.y() }; }

		template <int Order> Vector2i CalcLower(Vector2d const &pos) const { return ((pos - m_Origin) * m_InvSpacing - Vector2d::Ones() * (Order - 1) / 2).array().floor().cast<int>().matrix(); }

		Array2d CalcLowerFrac(Vector2d const &pos, Vector2i const &lower) const { return (pos - m_Origin - lower.cast<double>() * m_Spacing).array() * m_InvSpacing; }

		static constexpr int GetNumNeighbors() { return 4; }
		static constexpr int NeighborAxisOf(int ord) { return ord >> 1; }
		static constexpr int NeighborSideOf(int ord) { return ord & 1 ? 1 : -1; }
		static Vector2i NeighborOf(Vector2i const &coord, int ord) { return coord + Vector2i::Unit(NeighborAxisOf(ord)) * NeighborSideOf(ord); }

	private:
		double   m_Spacing;
		double   m_InvSpacing;
		Vector2i m_Size;
		Vector2d m_Origin;
	};

	template <typename Func>
		requires std::is_convertible_v<Func, std::function<void(Vector2i const &)>>
	inline void ForEach(Grid const &grid, Func &&func) {
		for (int i = 0; i < grid.GetSize().x(); i++) {
			for (int j = 0; j < grid.GetSize().y(); j++) {
				func(Vector2i(i, j));
			}
		}
	}

	template <typename Func>
		requires std::is_convertible_v<Func, std::function<void(Vector2i const &)>>
	inline void ParallelForEach(Grid const &grid, Func &&func) {
		tbb::parallel_for(tbb::blocked_range2d<int>(0, grid.GetSize().x(), 0, grid.GetSize().y()), [&](tbb::blocked_range2d<int> const &r) {
			for (int i = r.rows().begin(); i != r.rows().end(); i++) {
				for (int j = r.cols().begin(); j != r.cols().end(); j++) {
					func(Vector2i(i, j));
				}
			}
		});
	}

	template <typename Func>
		requires std::is_convertible_v<Func, std::function<void(int, Vector2i const &)>>
	inline void ForEach(std::array<Grid, 2> const &grids, Func &&func) {
		for (int axis = 0; axis < 2; axis++) {
			ForEach(grids[axis], [&](Vector2i const &coord) {
				func(axis, coord);
			});
		}
	}

	template <typename Func>
		requires std::is_convertible_v<Func, std::function<void(int, Vector2i const &)>>
	inline void ParallelForEach(std::array<Grid, 2> const &grids, Func &&func) {
		tbb::parallel_for(0, 2, [&](int axis) {
			ParallelForEach(grids[axis], [&](Vector2i const &coord) {
				func(axis, coord);
			});
		});
	}
}
