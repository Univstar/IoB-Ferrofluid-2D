#include "FiniteDiff.h"

namespace Pivot {
	double FiniteDiff::CalcCurvature(GridData<double> const &grData, Vector2i const &coord) {
		double const phi_x = CalcFirstDrv(grData, coord, 0);
		double const phi_y = CalcFirstDrv(grData, coord, 1);
		double const phi_xx = CalcSecondDrv(grData, coord, 0, 0);
		double const phi_yy = CalcSecondDrv(grData, coord, 1, 1);
		double const phi_xy = CalcSecondDrv(grData, coord, 0, 1);
		double const gradNorm2 = phi_x * phi_x + phi_y * phi_y;
		double const kappa = (phi_xx * phi_y * phi_y - 2 * phi_x * phi_y * phi_xy + phi_yy * phi_x * phi_x) / (gradNorm2 * std::sqrt(gradNorm2));
		double const invDx = grData.GetGrid().GetInvSpacing();

		// Vector2d nr = { CalcFirstDrv(grData, coord + Vector2i(1, 0), 0), CalcFirstDrv(grData, coord + Vector2i(1, 0), 1) };
		// Vector2d nl = { CalcFirstDrv(grData, coord - Vector2i(1, 0), 0), CalcFirstDrv(grData, coord - Vector2i(1, 0), 1) };
		// Vector2d nu = { CalcFirstDrv(grData, coord + Vector2i(0, 1), 0), CalcFirstDrv(grData, coord + Vector2i(0, 1), 1) };
		// Vector2d nd = { CalcFirstDrv(grData, coord - Vector2i(0, 1), 0), CalcFirstDrv(grData, coord - Vector2i(0, 1), 1) };
		// double kappa = (nr.normalized()[0] - nl.normalized()[0] + nu.normalized()[1] - nd.normalized()[1]) * invDx * .5;

		return kappa;
		// return std::abs(kappa) < invDx ? kappa: (kappa < 0 ? -1 : +1) * invDx;
	}
}
