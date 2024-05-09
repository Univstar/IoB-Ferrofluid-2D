#include "Reinitialization.h"

namespace Pivot {
	static double SolveQuadratic(double p0, double dx) { return p0 + dx; }

	static double SolveQuadratic(double p0, double p1, double dx) {
		if (p0 > p1) std::swap(p0, p1);
		if (std::isinf(p1) || p1 - p0 > dx) {
			return SolveQuadratic(p0, dx);
		} else {
			return ((p0 + p1) + std::sqrt(2 * dx * dx - (p0 - p1) * (p0 - p1))) * .5;
		}
	}

	void Reinitialization::Solve(GridData<double> &phi, int maxSteps) {
		double const bandWidth = maxSteps * phi.GetGrid().GetSpacing();
		GridData<std::int8_t> visited(phi.GetGrid());
		GridData<double> tent(phi.GetGrid(), bandWidth > 0 ? bandWidth : std::numeric_limits<double>::infinity());
		std::vector<int> intfIndices;
		Heap heap;
		// Initialize interface cells
		ForEach(phi.GetGrid(), [&](Vector2i const &coord) {
			Vector2d tempPhi = Vector2d::Ones() * std::numeric_limits<double>::infinity();
			for (int i = 0; i < Grid::GetNumNeighbors(); i++) {
				Vector2i nbCoord = Grid::NeighborOf(coord, i);
				if (phi.GetGrid().IsValid(nbCoord) && phi[coord] * phi[nbCoord] <= 0) {
					int const axis = Grid::NeighborAxisOf(i);
					tempPhi[axis] = std::min(tempPhi[axis], phi[coord] / (phi[coord] - phi[nbCoord]) * phi.GetGrid().GetSpacing());
				}
			}
			if (tempPhi.array().isFinite().any()) {
				tent[coord] = 1. / tempPhi.cwiseInverse().norm();
				visited[coord] = true;
				intfIndices.push_back(phi.GetGrid().IndexOf(coord));
			}
		});
		// Perform the algorithm
		for (auto const index : intfIndices) {
			UpdateNeighbors(tent.GetGrid().CoordOf(index), visited, tent, heap);
		}
		while (!heap.empty()) {
			double const val = heap.top().first;
			Vector2i const coord = tent.GetGrid().CoordOf(heap.top().second);
			heap.pop();
			if (tent[coord] != val) continue;
			visited[coord] = true;
			UpdateNeighbors(coord, visited, tent, heap);
		}
		ParallelForEach(phi.GetGrid(), [&](Vector2i const &coord) {
			phi[coord] = (phi[coord] <= 0 ? -1 : 1) * tent[coord];
		});
	}

	void Reinitialization::UpdateNeighbors(Vector2i const &coord, GridData<std::int8_t> const &visited, GridData<double> &tent, Heap &heap) {
		for (int i = 0; i < Grid::GetNumNeighbors(); i++) {
			Vector2i const nbCoord = Grid::NeighborOf(coord, i);
			if (!tent.GetGrid().IsValid(nbCoord) || visited[nbCoord]) continue;
			if (auto const temp = SolveEikonalEquation(nbCoord, visited, tent); temp < tent[nbCoord]) {
				tent[nbCoord] = temp;
				heap.push(HeapElement(temp, tent.GetGrid().IndexOf(nbCoord)));
			}
		}
	}

	double Reinitialization::SolveEikonalEquation(Vector2i const &coord, GridData<std::int8_t> const &visited, GridData<double> const &tent) {
		Vector2d tempPhi = Vector2d::Ones() * std::numeric_limits<double>::infinity();
		Vector2i mark = Vector2i::Zero();
		for (int i = 0; i < Grid::GetNumNeighbors(); i++) {
			Vector2i const nbCoord = Grid::NeighborOf(coord, i);
			if (tent.GetGrid().IsValid(nbCoord) && visited[nbCoord]) {
				int const axis = Grid::NeighborAxisOf(i);
				tempPhi[axis] = std::min(tempPhi[axis], tent[nbCoord]);
			}
		}
		double newPhi;
		newPhi = SolveQuadratic(tempPhi.x(), tempPhi.y(),tent.GetGrid().GetSpacing());
		if (!std::isfinite(newPhi)) {
			spdlog::critical("Failed to solve Eikonal equation");
			std::exit(EXIT_FAILURE);
		}
		return newPhi;
	}
}
