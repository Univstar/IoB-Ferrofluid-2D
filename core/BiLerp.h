#pragma once

#include "GridData.h"

namespace Pivot {
	class BiLerp {
	private:
		using WtPoint     = std::pair<Vector2i, double>;
		using GradWtPoint = std::pair<Vector2i, Vector2d>;

	public:
		static std::array<Vector2i, 4> GetPoints(Grid const &grid, Vector2d const &pos) {
			Vector2i const lower = grid.CalcLower<1>(pos);
			return {
				lower + Vector2i(0, 0), lower + Vector2i(1, 0),
				lower + Vector2i(0, 1), lower + Vector2i(1, 1)
			};
		}

		static std::array<WtPoint, 4> GetWtPoints(Grid const &grid, Vector2d const &pos) {
			Vector2i const lower = grid.CalcLower<1>(pos);
			Array2d  const frac = grid.CalcLowerFrac(pos, lower);
			std::array<Array2d, 2> const w = { 1. - frac, frac };
			return {
				WtPoint(lower + Vector2i(0, 0), w[0][0] * w[0][1]),
				WtPoint(lower + Vector2i(1, 0), w[1][0] * w[0][1]),
				WtPoint(lower + Vector2i(0, 1), w[0][0] * w[1][1]),
				WtPoint(lower + Vector2i(1, 1), w[1][0] * w[1][1]),
			};
		}

		static std::array<GradWtPoint, 4> GetGradWtPoints(Grid const &grid, Vector2d const &pos) {
			Vector2i const lower = grid.CalcLower<1>(pos);
			Array2d  const frac = grid.CalcLowerFrac(pos, lower);
			std::array<Array2d, 2> const w = { 1. - frac, frac };
			return {
				GradWtPoint(lower + Vector2i(0, 0), Vector2d(-w[0][1], -w[0][0]) * grid.GetInvSpacing()),
				GradWtPoint(lower + Vector2i(1, 0), Vector2d( w[0][1], -w[1][0]) * grid.GetInvSpacing()),
				GradWtPoint(lower + Vector2i(0, 1), Vector2d(-w[1][1],  w[0][0]) * grid.GetInvSpacing()),
				GradWtPoint(lower + Vector2i(1, 1), Vector2d( w[1][1],  w[1][0]) * grid.GetInvSpacing()),
			};
		}

		template <typename Type>
		static Type Interpolate(GridData<Type> const &grData, Vector2d const &pos) {
			Type val = Zero<Type>();
			for (auto const [coord, weight] : GetWtPoints(grData.GetGrid(), pos)) {
				val += grData.At(coord) * weight;
			}
			return val;
		}

		template <typename Type> 
			requires (std::is_arithmetic_v<Type>)
		static Vector<Type, 2> Interpolate(SGridData<Type> const &sgrData, Vector2d const &pos) {
			return { Interpolate(sgrData[0], pos), Interpolate(sgrData[1], pos) };
		}
	};
}
