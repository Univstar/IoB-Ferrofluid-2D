#pragma once

#include "ImplicitSurface.h"
#include "SGridData.h"
#include "StaggeredGrid.h"

namespace Pivot {
	class Collider {
	public:
		explicit Collider(StaggeredGrid const &sgrid);

		ImplicitBox       const &GetDomainBox  () const { return m_DomainBox; }
		SGridData<double> const &GetFraction   () const { return m_Fraction; }
		SGridData<double> const &GetNormal     () const { return m_Normal; }
		GridData<double>  const &GetAuxLevelSet() const { return m_AuxLevelSet; }

		bool IsInside(Vector2i const &cell) const { return m_AuxLevelSet[cell] <= 0; }

		void Finish(StaggeredGrid const &sgrid);

		void Enforce(SGridData<double> &fluidVelocity) const;

	private:
		double CalcFaceFraction(int axis, Vector2i const &face) const;

	public:
		GridData<double>  LevelSet;
		SGridData<double> Velocity;
	
	private:
		ImplicitBox       m_DomainBox;
		SGridData<double> m_Fraction;
		SGridData<double> m_Normal;
		GridData<double>  m_AuxLevelSet;
	};
}
