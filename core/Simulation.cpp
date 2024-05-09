#include "Simulation.h"

#include "Advection.h"
#include "BiLerp.h"
#include "CSG.h"
#include "Extrapolation.h"
#include "FiniteDiff.h"
#include "Image.h"
#include "Reinitialization.h"

namespace Pivot {
Simulation::Simulation(StaggeredGrid const &sgrid)
    : m_SGrid{sgrid}, m_Collider(m_SGrid), m_Pressure(m_SGrid),
      m_Velocity(m_SGrid.GetFaceGrids()),
      m_LevelSet(m_SGrid.GetCellGrid(),
                 std::numeric_limits<double>::infinity()),
      m_Contour(m_SGrid.GetCellGrid()) {}

void Simulation::Export(std::filesystem::path const &filename) const {
    { // Export the image
        Vector2i const size = m_SGrid.GetCellGrid().GetSize();
        std::vector<std::uint8_t> pixels(size.x() * size.y() * 16);
        for (int i = 0; i < size.x() * 4; i++) {
            for (int j = 0; j < size.y() * 4; j++) {
                Vector2d const pos = m_SGrid.GetDomainOrigin() + Vector2d(i + .5, j + .5).cwiseQuotient(size.cast<double>() * 4).cwiseProduct(m_SGrid.GetDomainLengths());
                double const phi = BiLerp::Interpolate(m_LevelSet, pos);
                pixels[j * size.x() * 4 + i] = phi <= 0 ? 0 : 255;
            }
        }
        std::span<std::uint8_t const> pixelsSpan(pixels);
        Image::WriteBytes(filename, std::as_bytes(pixelsSpan), size.x() * 4, size.y() * 4, 1, true);
    }
}

void Simulation::Initialize() {
    m_Collider.Finish(m_SGrid);
    CSG::Intersect(m_LevelSet, m_Collider.GetDomainBox());

    ReinitializeLevelSet(true);
    if (m_MagneticEnabled) {
        m_Magnetic.Solve(m_Contour.GetMesh());
    }
}

void Simulation::Advance(double deltaTime) {
    AdvectFields(deltaTime);
    ApplyBodyForces(deltaTime);
    ApplySurfacePressure(deltaTime);
    ProjectVelocity(deltaTime);
}

void Simulation::AdvectFields(double dt) {
    Advection::Solve<2>(m_LevelSet, m_Velocity, dt);
    Advection::Solve<2>(m_Velocity, m_Velocity, dt);

    ReinitializeLevelSet();
}

void Simulation::ApplyBodyForces(double dt) {
    if (m_GravityEnabled) {
        ParallelForEach(m_Velocity[1].GetGrid(), [&](Vector2i const &face) {
            m_Velocity[1][face] -= 9.8 * dt;
        });
    }
}

void Simulation::ApplySurfacePressure(double dt) {
    if (m_MagneticEnabled) {
        m_Magnetic.Solve(m_Contour.GetMesh());
    }
    m_Pressure.SetPressureJump([&](int axis, Vector2i const &face,
                                   double theta) -> double {
        int index = m_Contour.VertexIndexOf(axis, face - Vector2i::Unit(axis));
        if (index < 0) {
            return 0;
        } else {
            double pressure = 0;
            if (m_SurfaceTensionEnabled) {
                double kappa = m_Contour.GetMesh().MeanCurvatures[index];
                pressure += kappa * m_SurfaceTensionCoeff / m_LiquidDensity *
                            m_SGrid.GetInvSpacing() * dt;
            }
            if (m_MagneticEnabled) {
                pressure += m_Magnetic.m_MagneticPressure[index] /
                            m_LiquidDensity * m_SGrid.GetInvSpacing() * dt;
            }
            return pressure;
        }
    });
}

void Simulation::ProjectVelocity(double dt) {
    double x = (m_CurrentVolume - m_InitVolume) / (m_InitVolume);
    m_CumulVolError += x * dt;
    double kp = 0.1 / dt;
    double ki = kp * kp / 16;
    double c = 1 / (x + 1) * (-kp * x - ki * m_CumulVolError);
    m_Pressure.Project(m_Velocity, m_LevelSet, m_Collider,
                       c * m_SGrid.GetSpacing());
    // m_Pressure.Project(m_Velocity, m_LevelSet, m_Collider);
    Extrapolation::Solve(
        m_Velocity, 0., 6, [&](int axis, Vector2i const &face) {
            Vector2i const cell0 = StaggeredGrid::AdjCellOfFace(axis, face, 0);
            Vector2i const cell1 = StaggeredGrid::AdjCellOfFace(axis, face, 1);
            return m_Collider.GetFraction()[axis][face] < 1 &&
                   (m_LevelSet[cell0] <= 0 || m_LevelSet[cell1] <= 0);
        });
    m_Collider.Enforce(m_Velocity);
}

void Simulation::ReinitializeLevelSet(bool initial) {
    Extrapolation::Solve(
        m_LevelSet, 1.5 * m_SGrid.GetSpacing(), 1,
        [&](Vector2i const &cell) { return !m_Collider.IsInside(cell); });
    Reinitialization::Solve(m_LevelSet, 5);

    auto opLevelSet = m_LevelSet;
    CSG::Except(opLevelSet, m_Collider.GetAuxLevelSet());
    m_Contour.Generate(opLevelSet);
    // m_Contour.ComputeVertexInfos();
    m_Contour.ComputeVertexInfosFromLS(opLevelSet);
    m_Contour.ComputeVolumeFromLS(opLevelSet);
    m_CurrentVolume = m_Contour.GetMesh().TotalVolume;
    if (initial) {
        m_InitVolume = m_CurrentVolume;
    }
    fmt::print("volume {:.3e}/{:.3e} ", m_CurrentVolume, m_InitVolume);
}
} // namespace Pivot
