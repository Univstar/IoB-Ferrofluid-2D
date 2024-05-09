#pragma once

#include "Surface.h"

namespace Pivot {
class SurfaceMesh : public Surface {
  public:
    SurfaceMesh() = default;

    virtual Vector2d ClosestNormalOf(Vector2d const &pos) const override {
        return Vector2d::Zero();
    } // FIXME
    virtual double SignedDistanceTo(Vector2d const &pos) const override {
        return 0;
    } // FIXME

    void Clear();
    void Export(std::ostream &out) const;

    void ComputeAreas();
    void ComputeVolume();
    void ComputeMeanCurvatures();

    int size() { return Positions.size(); }

  public:
    std::vector<Vector2d> Positions;
    std::vector<Vector2d> Normals;
    std::vector<std::uint32_t> Indices;

    std::vector<double> Areas;
    double TotalArea;
    double TotalVolume;
    std::vector<double> MeanCurvatures;
};
} // namespace Pivot
