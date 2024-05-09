#pragma once

#include "BiLerp.h"
#include "BiCuInterp.h"
#include "SGridData.h"

namespace Pivot {
	class Advection {
	public:
		template <int RkOrder>
			requires (1 <= RkOrder && RkOrder <= 4)
		static Vector2d Trace(Vector2d const &startPos, SGridData<double> const &flow, double dt) {
			if constexpr (RkOrder == 1) {
				return startPos + BiLerp::Interpolate(flow, startPos) * dt;
			} else if constexpr (RkOrder == 2) { // the TVD Runge-Kutta scheme of second order
				Vector2d const vel0 = BiLerp::Interpolate(flow, startPos);
				Vector2d const vel1 = BiLerp::Interpolate(flow, startPos + vel0 * dt);
				return startPos + (vel0 + vel1) * dt / 2;
			} else if constexpr (RkOrder == 3) { // the TVD Runge-Kutta scheme of third order
				Vector2d const vel0 = BiLerp::Interpolate(flow, startPos);
				Vector2d const vel1 = BiLerp::Interpolate(flow, startPos + vel0 * dt);
				Vector2d const vel2 = BiLerp::Interpolate(flow, startPos + (vel0 + vel1) * dt / 4);
				return startPos + (vel0 + vel1 + 4 * vel2) * dt / 6;
			} else if constexpr (RkOrder == 4) { // the classical fourth order Runge-Kutta scheme
				Vector2d const vel0 = BiLerp::Interpolate(flow, startPos);
				Vector2d const vel1 = BiLerp::Interpolate(flow, startPos + vel0 * dt / 2);
				Vector2d const vel2 = BiLerp::Interpolate(flow, startPos + vel1 * dt / 2);
				Vector2d const vel3 = BiLerp::Interpolate(flow, startPos + vel2 * dt);
				return startPos + (vel0 + 2 * vel1 + 2 * vel2 + vel3) * dt / 6;
			}
		}

		template <int RkOrder, typename Type>
			requires (1 <= RkOrder && RkOrder <= 4)
		static void Solve(GridData<Type> &grData, SGridData<double> const &flow, double dt) {
			GridData<Type> newGrData(grData.GetGrid());
			ParallelForEach(grData.GetGrid(), [&](Vector2i const &coord) {
				Vector2d const pos = grData.GetGrid().PositionOf(coord);
				// newGrData[coord] = BiLerp::Interpolate(grData, Trace<RkOrder>(pos, flow, -dt));
				newGrData[coord] = BiCuInterp::Interpolate(grData, Trace<RkOrder>(pos, flow, -dt));
			});
			grData = newGrData;
		}

		template <int RkOrder, typename Type>
			requires (1 <= RkOrder && RkOrder <= 4)
		static void Solve(SGridData<Type> &sgrData, SGridData<double> const &flow, double dt) {
			SGridData<Type> newSgrData(sgrData.GetGrids());
			ParallelForEach(sgrData.GetGrids(), [&](int axis, Vector2i const &face) {
				Vector2d const pos = sgrData[axis].GetGrid().PositionOf(face);
				newSgrData[axis][face] = BiLerp::Interpolate(sgrData[axis], Trace<RkOrder>(pos, flow, -dt));
			});
			sgrData = newSgrData;
		}
	};
}
