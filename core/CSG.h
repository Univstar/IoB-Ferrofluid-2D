#pragma once

#include "GridData.h"
#include "Surface.h"

namespace Pivot {
	class CSG {
	public:
		static void Union    (GridData<double> &levelSet, Surface const &surface);
		static void Intersect(GridData<double> &levelSet, Surface const &surface);
		static void Except   (GridData<double> &levelSet, Surface const &surface);

		static void Union    (GridData<double> &lhs, GridData<double> const &rhs);
		static void Intersect(GridData<double> &lhs, GridData<double> const &rhs);
		static void Except   (GridData<double> &lhs, GridData<double> const &rhs);
	};
}
