#include "CSG.h"

namespace Pivot {
	void CSG::Union(GridData<double> &levelSet, Surface const &surface) {
		ParallelForEach(levelSet.GetGrid(), [&](Vector2i const &coord) {
			Vector2d const pos = levelSet.GetGrid().PositionOf(coord);
			levelSet[coord] = std::min(levelSet[coord], surface.SignedDistanceTo(pos));
		});
	}

	void CSG::Intersect(GridData<double> &levelSet, Surface const &surface) {
		ParallelForEach(levelSet.GetGrid(), [&](Vector2i const &coord) {
			Vector2d const pos = levelSet.GetGrid().PositionOf(coord);
			levelSet[coord] = std::max(levelSet[coord], surface.SignedDistanceTo(pos));
		});
	}

	void CSG::Except(GridData<double> &levelSet, Surface const &surface) {
		ParallelForEach(levelSet.GetGrid(), [&](Vector2i const &coord) {
			Vector2d const pos = levelSet.GetGrid().PositionOf(coord);
			levelSet[coord] = std::max(levelSet[coord], -surface.SignedDistanceTo(pos));
		});
	}

	void CSG::Union(GridData<double> &lhs, GridData<double> const &rhs) {
		ParallelForEach(lhs.GetGrid(), [&](Vector2i const &coord) {
			lhs[coord] = std::min(lhs[coord], rhs[coord]);
		});
	}

	void CSG::Intersect(GridData<double> &lhs, GridData<double> const &rhs) {
		ParallelForEach(lhs.GetGrid(), [&](Vector2i const &coord) {
			lhs[coord] = std::max(lhs[coord], rhs[coord]);
		});
	}

	void CSG::Except(GridData<double> &lhs, GridData<double> const &rhs) {
		ParallelForEach(lhs.GetGrid(), [&](Vector2i const &coord) {
			lhs[coord] = std::max(lhs[coord], -rhs[coord]);
		});
	}
}
