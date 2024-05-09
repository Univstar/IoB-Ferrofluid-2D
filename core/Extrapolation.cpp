#include "Extrapolation.h"

namespace Pivot {
	void Extrapolation::Solve(GridData<double> &grData, double clearVal, int maxSteps, GridData<std::uint8_t> &valid) {
		auto valid1 = valid;
		for (int iter = 0; iter < maxSteps; iter++) {
			ParallelForEach(grData.GetGrid(), [&](Vector2i const &coord) {
				if (!valid[coord]) {
					int    cnt = 0;
					double sum = 0;
					for (int i = 0; i < Grid::GetNumNeighbors(); i++) {
						Vector2i const nbCoord = Grid::NeighborOf(coord, i);
						if (grData.GetGrid().IsValid(nbCoord) && valid[nbCoord]) {
							sum += grData[nbCoord];
							cnt++;
						}
					}
					if (cnt > 0) {
						grData[coord] = sum / cnt;
						valid1[coord] = true;
					} else {
						grData[coord] = clearVal;
					}
				} else {
					valid1[coord] = true;
				}
			});
			valid.GetData().swap(valid1.GetData());
		}
	}
}
