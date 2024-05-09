#pragma once

#include "GridData.h"

namespace Pivot {
class BiCuInterp {
  private:
    using WtPoint = std::pair<Vector2i, double>;

  public:
    static std::array<WtPoint, 16> GetWtPoints(Grid const &grid,
                                               Vector2d const &pos) {
        Vector2i const lower = grid.CalcLower<1>(pos);
        Array2d const frac = grid.CalcLowerFrac(pos, lower);
        std::array<Array2d, 3> const s = {frac, frac * frac,
                                          frac * frac * frac};
        std::array<Array2d, 4> const w = {
            -s[0] / 3 + s[1] / 2 - s[2] / 6, 1 - s[1] + (s[2] - s[0]) / 2,
            s[0] + (s[1] - s[2]) / 2, (s[2] - s[0]) / 6};
        return {
            WtPoint(lower + Vector2i(-1, -1), w[0][0] * w[0][1]),
            WtPoint(lower + Vector2i(+0, -1), w[1][0] * w[0][1]),
            WtPoint(lower + Vector2i(+1, -1), w[2][0] * w[0][1]),
            WtPoint(lower + Vector2i(+2, -1), w[3][0] * w[0][1]),
            WtPoint(lower + Vector2i(-1, +0), w[0][0] * w[1][1]),
            WtPoint(lower + Vector2i(+0, +0), w[1][0] * w[1][1]),
            WtPoint(lower + Vector2i(+1, +0), w[2][0] * w[1][1]),
            WtPoint(lower + Vector2i(+2, +0), w[3][0] * w[1][1]),
            WtPoint(lower + Vector2i(-1, +1), w[0][0] * w[2][1]),
            WtPoint(lower + Vector2i(+0, +1), w[1][0] * w[2][1]),
            WtPoint(lower + Vector2i(+1, +1), w[2][0] * w[2][1]),
            WtPoint(lower + Vector2i(+2, +1), w[3][0] * w[2][1]),
            WtPoint(lower + Vector2i(-1, +2), w[0][0] * w[3][1]),
            WtPoint(lower + Vector2i(+0, +2), w[1][0] * w[3][1]),
            WtPoint(lower + Vector2i(+1, +2), w[2][0] * w[3][1]),
            WtPoint(lower + Vector2i(+2, +2), w[3][0] * w[3][1]),
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
        requires(std::is_arithmetic_v<Type>)
    static Vector<Type, 2> Interpolate(SGridData<Type> const &sgrData,
                                       Vector2d const &pos) {
        return {Interpolate(sgrData[0], pos), Interpolate(sgrData[1], pos)};
    }
};
} // namespace Pivot
