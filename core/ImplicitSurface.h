#pragma once

#include "Surface.h"

namespace Pivot {

	class ImplicitSphere : public Surface {
	public:

		ImplicitSphere(Vector2d const &center, double radius) : m_Center { center }, m_Radius { radius } { }

		virtual Vector2d ClosestPositionOf(Vector2d const &pos) const override { return m_Center + ClosestNormalOf(pos) * m_Radius; }
		virtual Vector2d ClosestNormalOf  (Vector2d const &pos) const override { return (pos - m_Center).normalized(); }
		virtual double   SignedDistanceTo (Vector2d const &pos) const override { return (pos - m_Center).norm() - m_Radius; }

	private:
		Vector2d m_Center;
		double   m_Radius;
	};

	class ImplicitBox : public Surface {
	public:

		ImplicitBox(Vector2d const &minCorner, Vector2d const &lengths) : m_Center { minCorner + lengths / 2 }, m_HalfLengths { lengths / 2 } { }

		virtual Vector2d ClosestNormalOf(Vector2d const &pos) const override {
			Vector2d const phi = (pos - m_Center).cwiseAbs() - m_HalfLengths;
			Vector2d normal;
			if ((phi.array() <= 0).all()) {
				int axis;
				phi.maxCoeff(&axis);
				normal = Vector2d::Unit(axis);
			} else {
				normal = phi.cwiseMax(0);
			}
			return normal.cwiseProduct((pos - m_Center).cwiseSign()).normalized();
		}

		virtual double SignedDistanceTo(Vector2d const &pos) const override {
			Vector2d const phi = (pos - m_Center).cwiseAbs() - m_HalfLengths;
			if ((phi.array() <= 0).all()) {
				return phi.maxCoeff();
			} else {
				return phi.cwiseMax(0).norm();
			}
		}

	private:
		Vector2d m_Center;
		Vector2d m_HalfLengths;
	};

	class ImplicitPlane : public Surface {
	public:
		ImplicitPlane(Vector2d const &position, Vector2d const &direction) : m_Position(position), m_Normal(direction.normalized()) { }

		virtual Vector2d ClosestNormalOf (Vector2d const &pos) const override { return m_Normal; }
		virtual double   SignedDistanceTo(Vector2d const &pos) const override { return (pos - m_Position).dot(m_Normal); }

	private:
		Vector2d const m_Position;
		Vector2d const m_Normal;
	};

	class ImplicitEllipsoid : public Surface {
	public:

		ImplicitEllipsoid(Vector2d const &center, Vector2d const &semiAxels) : m_Center(center), m_SemiAxels(semiAxels) { }

		virtual Vector2d ClosestPositionOf(Vector2d const &pos) const override { return pos / pos.cwiseQuotient(m_SemiAxels).norm(); } // not accurate solution
		virtual Vector2d ClosestNormalOf  (Vector2d const &pos) const override { return (pos - ClosestPositionOf(pos)).normalized() * (Surrounds(pos) ? -1 : 1); }
		virtual double   DistanceTo       (Vector2d const &pos) const override { return (pos - ClosestPositionOf(pos)).norm(); }
		virtual double   SignedDistanceTo (Vector2d const &pos) const override { return DistanceTo(pos) * (Surrounds(pos) ? -1 : 1); }
		virtual bool     Surrounds        (Vector2d const &pos) const override { return pos.cwiseQuotient(m_SemiAxels).squaredNorm() <= 1; }

	private:
		Vector2d m_Center;
		Vector2d m_SemiAxels;
	};
}
