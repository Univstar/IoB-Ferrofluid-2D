#include "SurfaceMesh.h"

namespace Pivot {
void SurfaceMesh::Clear() {
    Positions.clear();
    Normals.clear();
    Indices.clear();
}

void SurfaceMesh::Export(std::ostream &out) const {
    IO::Write(out, static_cast<std::uint32_t>(Positions.size()));
    for (auto const &pos : Positions) {
        IO::Write(out, pos.cast<float>().eval());
    }
    // for (auto const curv : MeanCurvatures) {
    // 	IO::Write(out, static_cast<float>(std::abs(curv)));
    // }
    IO::Write(out, static_cast<std::uint32_t>(Indices.size()));
    IO::Write(out, Indices);
}

void SurfaceMesh::ComputeAreas() {
    TotalArea = 0;
    Areas.resize(Positions.size());
    std::fill(Areas.begin(), Areas.end(), 0.);
    for (std::size_t i = 0; i < Indices.size(); i += 2) {
        auto const i0 = Indices[i + 0];
        auto const i1 = Indices[i + 1];
        double const faceArea = (Positions[i1] - Positions[i0]).norm();
        TotalArea += faceArea;
        Areas[i0] += faceArea / 2;
        Areas[i1] += faceArea / 2;
    }
}

void SurfaceMesh::ComputeVolume() {
    TotalVolume = 0;
    for (std::size_t i = 0; i < Indices.size(); i += 2) {
        auto const i0 = Indices[i + 0];
        auto const i1 = Indices[i + 1];
        Vector2d v0 = Positions[i0];
        Vector2d v1 = Positions[i1];
        TotalVolume += (v0.x() * v1.y() - v0.y() * v1.x()) / 2;
    }
}

void SurfaceMesh::ComputeMeanCurvatures() {
    ComputeAreas();

    std::vector<Vector2d> inEdge(Positions.size());
    std::vector<Vector2d> outEdge(Positions.size());
    for (std::size_t i = 0; i < Indices.size(); i += 2) {
        auto const i0 = Indices[i + 0];
        auto const i1 = Indices[i + 1];
        Vector2d const vec = Positions[i1] - Positions[i0];
        inEdge[i1] = vec;
        outEdge[i0] = vec;
    }
    MeanCurvatures.resize(Positions.size());
    tbb::parallel_for(
        static_cast<std::size_t>(0), Positions.size(), [&](std::size_t i) {
            double const cosTheta = inEdge[i].dot(outEdge[i]) /
                                    (inEdge[i].norm() * outEdge[i].norm());
            int const sign = inEdge[i].x() * outEdge[i].y() -
                                         inEdge[i].y() * outEdge[i].x() >
                                     0
                                 ? +1
                                 : -1;
            MeanCurvatures[i] = std::acos(cosTheta) * sign / Areas[i];
        });
}
} // namespace Pivot
