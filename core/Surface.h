#pragma once

#include "Common.h"

namespace Pivot {
	class Surface {
	public:
		virtual ~Surface() = default;

		virtual Vector2d ClosestPositionOf(Vector2d const &pos) const { return pos - SignedDistanceTo(pos) * ClosestNormalOf(pos); }
		virtual Vector2d ClosestNormalOf  (Vector2d const &pos) const = 0;
		virtual double   DistanceTo       (Vector2d const &pos) const { return std::abs(SignedDistanceTo(pos)); }
		virtual double   SignedDistanceTo (Vector2d const &pos) const = 0;
		virtual bool     Surrounds        (Vector2d const &pos) const { return SignedDistanceTo(pos) <= 0; }
	};
}
