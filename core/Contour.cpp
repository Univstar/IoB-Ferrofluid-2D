#include "Contour.h"

#include "StaggeredGrid.h"

#include "FiniteDiff.h"

#include "fraction.hpp"

namespace Pivot {
static inline std::uint8_t GetCellType(GridData<double> const &grData,
                                       Vector2i const &cell, double value) {
    std::uint8_t type = 0;
    for (int i = 0; i < StaggeredGrid::GetNumNodesPerCell(); i++) {
        Vector2i const node = StaggeredGrid::NodeOfCell(cell, i);
        if (grData[node] <= value)
            type |= 1 << i;
    }
    // Handle diagonal cases for 2D
    if (type == 6 || type == 9) {
        double sumPhi = 0;
        for (int i = 0; i < StaggeredGrid::GetNumNodesPerCell(); i++) {
            sumPhi += grData[StaggeredGrid::NodeOfCell(cell, i)];
        }
        if (sumPhi <= 4 * value)
            type = 16 + (type > 8);
    }
    return type;
}

Contour::Contour(Grid const &nodeGrid)
    : m_NodeGrid{nodeGrid},
      m_CellGrid(m_NodeGrid.GetSpacing(),
                 m_NodeGrid.GetSize() - Vector2i::Ones(),
                 m_NodeGrid.GetOrigin() +
                     Vector2d::Constant(m_NodeGrid.GetSpacing() / 2)),
      m_EdgeGrids{
          Grid(m_NodeGrid.GetSpacing(),
               m_NodeGrid.GetSize() - Vector2i::Unit(0),
               m_NodeGrid.GetOrigin() +
                   Vector2d::Unit(0) * m_NodeGrid.GetSpacing() / 2),
          Grid(m_NodeGrid.GetSpacing(),
               m_NodeGrid.GetSize() - Vector2i::Unit(1),
               m_NodeGrid.GetOrigin() +
                   Vector2d::Unit(1) * m_NodeGrid.GetSpacing() / 2),
      },
      m_EdgeMark{
          GridData<int>(m_EdgeGrids[0]),
          GridData<int>(m_EdgeGrids[1]),
      } {}

void Contour::Generate(GridData<double> const &grData, double value) {
    if (m_NodeGrid != grData.GetGrid()) {
        spdlog::critical(
            "Failed to generate contour because of incompatible grids");
        std::exit(EXIT_FAILURE);
    }

    m_EdgeMark[0].SetConstant(-1);
    m_EdgeMark[1].SetConstant(-1);
    m_Mesh.Clear();

    ForEach(m_CellGrid, [&](Vector2i const &cell) {
        auto const cellType = GetCellType(grData, cell, value);
        auto edgeState = c_EdgeStateTable2[cellType];
        for (int i = 0; i < StaggeredGrid::GetNumEdgesPerCell(); i++) {
            if (edgeState >> i & 1) {
                auto const [axis, edge] = StaggeredGrid::EdgeOfCell(cell, i);
                if (m_EdgeMark[axis][edge] < 0) {
                    Vector2i const node0 =
                        StaggeredGrid::NodeOfEdge(axis, edge, 0);
                    Vector2i const node1 =
                        StaggeredGrid::NodeOfEdge(axis, edge, 1);
                    double const theta = (grData[node0] - value) /
                                         (grData[node0] - grData[node1]);
                    Vector2d const pos =
                        (1 - theta) * m_NodeGrid.PositionOf(node0) +
                        theta * m_NodeGrid.PositionOf(node1);
                    m_EdgeMark[axis][edge] =
                        static_cast<int>(m_Mesh.Positions.size());
                    m_Mesh.Positions.push_back(pos);
                }
            }
        }
        for (auto *it = c_EdgeOrdsTable2[cellType]; *it != -1; it++) {
            auto const [axis, edge] = StaggeredGrid::EdgeOfCell(cell, *it);
            m_Mesh.Indices.push_back(
                static_cast<std::uint32_t>(m_EdgeMark[axis][edge]));
        }
    });
}

void Contour::ComputeVertexInfos() { m_Mesh.ComputeMeanCurvatures(); }

void Contour::ComputeVertexInfosFromLS(GridData<double> const &levelSet) {
    m_Mesh.ComputeAreas();

    m_Mesh.MeanCurvatures.resize(m_Mesh.Positions.size());
    m_Mesh.Normals.resize(m_Mesh.Positions.size());
    ForEach(m_CellGrid, [&](Vector2i const &cell) {
        for (int i = 0; i < StaggeredGrid::GetNumEdgesPerCell(); i++) {
            auto const [axis, edge] = StaggeredGrid::EdgeOfCell(cell, i);
            int index = VertexIndexOf(axis, edge);
            if (index >= 0) {
                Vector2i const cell0 = StaggeredGrid::NodeOfEdge(axis, edge, 0);
                Vector2i const cell1 = StaggeredGrid::NodeOfEdge(axis, edge, 1);
                double const theta =
                    levelSet[cell0] / (levelSet[cell0] - levelSet[cell1]);
                Vector2d grad0 =
                    Vector2d(FiniteDiff::CalcFirstDrv(levelSet, cell0, 0),
                             FiniteDiff::CalcFirstDrv(levelSet, cell0, 1))
                        .normalized();
                Vector2d grad1 =
                    Vector2d(FiniteDiff::CalcFirstDrv(levelSet, cell1, 0),
                             FiniteDiff::CalcFirstDrv(levelSet, cell1, 1))
                        .normalized();
                Vector2d grad =
                    ((1 - theta) * grad0 + theta * grad1).normalized();
                m_Mesh.Normals[index] = grad;
                double const kappa0 =
                    FiniteDiff::CalcCurvature(levelSet, cell0);
                double const kappa1 =
                    FiniteDiff::CalcCurvature(levelSet, cell1);
                double const kappa = (1 - theta) * kappa0 + theta * kappa1;
                m_Mesh.MeanCurvatures[index] = kappa;
            }
        }
    });
}
void Contour::ComputeVolumeFromLS(GridData<double> const &levelSet) {
    double vol = 0;
    ForEach(m_CellGrid, [&](Vector2i const &cell) {
        std::array<double, 4> phi2d{
            levelSet.At(cell + Vector2i(0, 0)),
            levelSet.At(cell + Vector2i(1, 0)),
            levelSet.At(cell + Vector2i(1, 1)),
            levelSet.At(cell + Vector2i(0, 1)),
        };
        vol += Fraction::get_ms_area(phi2d);
    });
    m_Mesh.TotalVolume =
        vol * levelSet.GetGrid().GetSpacing() * levelSet.GetGrid().GetSpacing();
}
} // namespace Pivot
