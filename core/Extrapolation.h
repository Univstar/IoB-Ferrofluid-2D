#include "SGridData.h"

namespace Pivot {
	class Extrapolation {
	public:
		static void Solve(GridData<double> &grData, double clearVal, int maxSteps, GridData<std::uint8_t> &valid);

		template <typename Func>
			requires std::is_convertible_v<Func, std::function<bool(Vector2i const &)>>
		static void Solve(GridData<double> &grData, double clearVal, int maxSteps, Func &&isValid) {
			GridData<std::uint8_t> valid(grData.GetGrid());
			ParallelForEach(grData.GetGrid(), [&](Vector2i const &coord) {
				valid[coord] = isValid(coord);
			});
			Solve(grData, clearVal, maxSteps, valid);
		}

		template <typename Func>
			requires std::is_convertible_v<Func, std::function<bool(int, Vector2i const &)>>
		static void Solve(SGridData<double> &sgrData, double clearVal, int maxSteps, Func &&isValid) {
			SGridData<std::uint8_t> valid(sgrData.GetGrids());
			ParallelForEach(sgrData.GetGrids(), [&](int axis, Vector2i const &face) {
				valid[axis][face] = isValid(axis, face);
			});
			tbb::parallel_for(0, 2, [&](int axis) {
				Solve(sgrData[axis], clearVal, maxSteps, valid[axis]);
			});
		}
	};
}
