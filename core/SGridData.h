#pragma once

#include "GridData.h"

namespace Pivot {
	template <typename Type>
	class SGridData {
	public:
		explicit SGridData(std::array<Grid, 2> const &grids, Vector2<Type> const &value = Vector2<Type>::Zero()) :
			m_Grids { grids },
			m_Datas {
				GridData<Type>(m_Grids[0], value[0]),
				GridData<Type>(m_Grids[1], value[1]),
			} {
		}

		SGridData &operator=(SGridData const &rhs) {
			for (int axis = 0; axis < 2; axis++) {
				m_Datas[axis] = rhs.m_Datas[axis];
			}
			return *this;
		}

		std::array<Grid, 2> const &GetGrids() const { return m_Grids; }

		GridData<Type>       &operator[](int axis)       { return m_Datas[axis]; }
		GridData<Type> const &operator[](int axis) const { return m_Datas[axis]; }

		Type GetMaxAbsComponent() const requires (std::is_arithmetic_v<Type>) { return std::max(m_Datas[0].GetMaxAbsValue(), m_Datas[1].GetMaxAbsValue()); }

		void SetConstant(Vector2<Type> const &value) { for (int axis = 0; axis < 2; axis++) m_Datas[axis].SetConstant(value[axis]); }
		void SetZero    ()                           { SetConstant(Vector2<Type>::Zero()); }
		
	private:
		std::array<Grid, 2>           const &m_Grids;
		std::array<GridData<Type>, 2>        m_Datas;
	};
}
