#pragma once

#include "GridData.h"

namespace Pivot {
	class FiniteDiff {
	public:
		template <typename Type>
		static Type CalcFirstDrv(GridData<Type> const &grData, Vector2i const &coord, int axis) {
			double const invDx = grData.GetGrid().GetInvSpacing();
			auto const f = [&](int i)->Type { return grData[coord + Vector2i::Unit(axis) * i]; };
			if (coord[axis] == 0) {
				return CalcFirstForward(f(0), f(+1), f(+2)) * invDx * .5;
			} else if (coord[axis] + 1 == grData.GetGrid().GetSize()[axis]) {
				return CalcFirstBackward(f(-2), f(-1), f(0)) * invDx * .5;
			} else {
				return CalcFirstCentral(f(-1), f(+1)) * invDx * .5;
			}
		}

		template <typename Type>
		static Type CalcSecondDrv(GridData<Type> const &grData, Vector2i const &coord, int axis1, int axis2) {
			double const invDx = grData.GetGrid().GetInvSpacing();
			if (axis1 == axis2) {
				auto const f = [&](int i)->Type { return grData[coord + Vector2i::Unit(axis1) * i]; };
				if (coord[axis1] == 0) {
					return CalcSecondForward(f(0), f(+1), f(+2), f(+3)) * invDx * invDx;
				} else if (coord[axis1] + 1 == grData.GetGrid().GetSize()[axis1]) {
					return CalcSecondBackward(f(-3), f(-2), f(-1), f(0)) * invDx * invDx;
				} else {
					return CalcSecondCentral(f(-1), f(0), f(+1)) * invDx * invDx;
				}
			} 
			auto const f = [&](int i)->Type { return CalcFirstDrv(grData, coord + Vector2i::Unit(axis1) * i, axis2); };
			if (coord[axis1] == 0) {
				return CalcFirstForward(f(0), f(+1), f(+2)) * invDx * .5;
			} else if (coord[axis1] + 1 == grData.GetGrid().GetSize()[axis1]) {
				return CalcFirstBackward(f(-2), f(-1), f(0)) * invDx * .5;
			} else {
				return CalcFirstCentral(f(-1), f(+1)) * invDx * .5;
			}
		}

		static double CalcCurvature(GridData<double> const &grData, Vector2i const &coord);
	
	private:
		template <typename Type> static Type CalcFirstCentral (Type const &f0, Type const &f2)                 { return f2 - f0; }
		template <typename Type> static Type CalcFirstForward (Type const &f0, Type const &f1, Type const &f2) { return 4 * f1 - 3 * f0 - f2; }
		template <typename Type> static Type CalcFirstBackward(Type const &f0, Type const &f1, Type const &f2) { return 3 * f2 - 4 * f1 + f0; }

		template <typename Type> static Type CalcSecondCentral (Type const &f0, Type const &f1, Type const &f2)                 { return f0 - 2 * f1 + f2; }
		template <typename Type> static Type CalcSecondForward (Type const &f0, Type const &f1, Type const &f2, Type const &f3) { return 2 * f0 - 5 * f1 + 4 * f2 - f3; }
		template <typename Type> static Type CalcSecondBackward(Type const &f0, Type const &f1, Type const &f2, Type const &f3) { return 2 * f3 - 5 * f2 + 4 * f1 - f0; }
	};
}
