#include "Collider.h"

#include "BiLerp.h"
#include "FiniteDiff.h"
#include "Reinitialization.h"

namespace Pivot {
	Collider::Collider(StaggeredGrid const &sgrid) :
		LevelSet(sgrid.GetNodeGrid()),
		Velocity(sgrid.GetFaceGrids()),
		m_DomainBox(sgrid.GetDomainOrigin(), sgrid.GetDomainLengths()),
		m_Fraction(sgrid.GetFaceGrids()),
		m_Normal(sgrid.GetFaceGrids()),
		m_AuxLevelSet(sgrid.GetCellGrid()) {
		ParallelForEach(LevelSet.GetGrid(), [&](Vector2i const &node) {
			Vector2d const pos = LevelSet.GetGrid().PositionOf(node);
			LevelSet[node] = -.01 * sgrid.GetSpacing() - m_DomainBox.SignedDistanceTo(pos);
		});
	}

	void Collider::Finish(StaggeredGrid const &sgrid) {
		Reinitialization::Solve(LevelSet, -1);
		ParallelForEach(m_AuxLevelSet.GetGrid(), [&](Vector2i const &cell) {
			for (int i = 0; i < StaggeredGrid::GetNumNodesPerCell(); i++) {
				m_AuxLevelSet[cell] += LevelSet[StaggeredGrid::NodeOfCell(cell, i)];
			}
			m_AuxLevelSet[cell] /= StaggeredGrid::GetNumNodesPerCell();
		});
		auto levelSetDrvs = std::to_array({ GridData<double>(sgrid.GetNodeGrid()), GridData<double>(sgrid.GetNodeGrid()) });
		tbb::parallel_for(0, 2, [&](int axis) {
			ParallelForEach(sgrid.GetNodeGrid(), [&](Vector2i const &node) {
				levelSetDrvs[axis][node] = FiniteDiff::CalcFirstDrv(LevelSet, node, axis);
			});
		});

		ParallelForEach(m_Fraction.GetGrids(), [&](int axis, Vector2i const &face) {
			m_Fraction[axis][face] = CalcFaceFraction(axis, face);
			double sum = 0;
			for (int i = 0; i < StaggeredGrid::GetNumNodesPerFace(); i++) {
				sum += levelSetDrvs[axis][StaggeredGrid::NodeOfFace(axis, face, i)];
			}
			m_Normal[axis][face] = sum / StaggeredGrid::GetNumNodesPerFace();
		});
	}

	double Collider::CalcFaceFraction(int axis, Vector2i const &face) const {
		constexpr auto theta = [](double phi0, double phi1) { return phi0 / (phi0 - phi1); };

		constexpr auto fraction = [](const double phi0, const double phi1) {
			if (phi0 <= 0 && phi1 <= 0) {
				return 1.;
			} else if (phi0 <= 0 && phi1 > 0) {
				return theta(phi0, phi1);
			} else if (phi0 > 0 && phi1 <= 0) {
				return theta(phi1, phi0);
			} else {
				return 0.;
			}
		};

		double const phi0 = LevelSet[StaggeredGrid::NodeOfFace(axis, face, 0)];
		double const phi1 = LevelSet[StaggeredGrid::NodeOfFace(axis, face, 1)];
		double faceFraction = fraction(phi0, phi1);

		return faceFraction > .9 ? 1. : faceFraction;
	}

	void Collider::Enforce(SGridData<double> &fluidVelocity) const {
		auto const oldFluidVelocity = fluidVelocity;
		ParallelForEach(fluidVelocity.GetGrids(), [&](int axis, Vector2i const &face) {
			if (m_Fraction[axis][face] == 1.) {
				Vector2d const pos = fluidVelocity[axis].GetGrid().PositionOf(face);
				Vector2d const vel = BiLerp::Interpolate(oldFluidVelocity, pos);
				Vector2d const vel0 = BiLerp::Interpolate(Velocity, pos);
				Vector2d const n = BiLerp::Interpolate(m_Normal, pos).normalized();
				fluidVelocity[axis][face] = (vel - (vel - vel0).dot(n) * n)[axis];
			}
		});
	}
}
