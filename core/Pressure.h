#pragma once

#include "Collider.h"

namespace Pivot {
class Pressure {
  public:
    explicit Pressure(StaggeredGrid const &sgrid);

    void Project(SGridData<double> &velocity, GridData<double> const &levelSet,
                 Collider const &collider, double volError = 0);

    template <typename Func>
        requires(std::is_convertible_v<
                 Func, std::function<double(int, Vector2i const &, double)>>)
    void SetPressureJump(Func &&pressureJump) {
        m_PressureJump = pressureJump;
    }

  private:
    void BuildProjectionMatrix(SGridData<double> const &velocity,
                               GridData<double> const &levelSet,
                               Collider const &collider, double volError = 0);

    void SetUnKnowns(GridData<double> const &levelSet);

    void SolveLinearSystem();

    void ApplyProjection(SGridData<double> &velocity,
                         GridData<double> const &levelSet,
                         Collider const &collider);

  private:
    GridData<int> m_Grid2Mat;
    std::vector<int> m_Mat2Grid;

    SparseMatrix<double, RowMajor> m_MatL; // A matrix of the Laplacian operator

    VectorXd m_RdP; // reduced pressure
    VectorXd m_Rhs;

    // Pressure jump: p_liquid - p_air
    std::function<double(int, Vector2i const &, double)> m_PressureJump =
        nullptr;
};
} // namespace Pivot
