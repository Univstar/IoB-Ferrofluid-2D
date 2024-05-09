#pragma once

#include "Grid.h"

namespace Pivot {
	class StaggeredGrid {
	public:
		StaggeredGrid(int boundaryWidth, double spacing, Vector2i resolution, Vector2d center = Vector2d::Zero()) :
			m_BoundaryWidth { boundaryWidth },
			m_Spacing { spacing },
			m_InvSpacing { 1 / m_Spacing },
			m_Resolution { resolution },
			m_Origin { center - m_Resolution.cast<double>() * m_Spacing / 2 },
			m_NodeGrid(m_Spacing, m_Resolution + Vector2i::Ones(), m_Origin),
			m_CellGrid(m_Spacing, m_Resolution, m_Origin + Vector2d::Constant(m_Spacing / 2)),
			m_FaceGrids {
				Grid(m_Spacing, m_Resolution + Vector2i(1, 0), m_Origin + Vector2d(0, 1) * m_Spacing / 2),
				Grid(m_Spacing, m_Resolution + Vector2i(0, 1), m_Origin + Vector2d(1, 0) * m_Spacing / 2),
			} {
		}

		int      GetBoundaryWidth() const { return m_BoundaryWidth; }
		double   GetSpacing      () const { return m_Spacing; }
		double   GetInvSpacing   () const { return m_InvSpacing; }
		Vector2i GetResolution   () const { return m_Resolution; }
		Vector2d GetOrigin       () const { return m_Origin; }

		int GetNumNodes() const { return m_NodeGrid.GetNumVertices(); }
		int GetNumCells() const { return m_CellGrid.GetNumVertices(); }
		int GetNumFaces() const { return m_FaceGrids[0].GetNumVertices() + m_FaceGrids[1].GetNumVertices(); }

		Grid                const &GetNodeGrid () const { return m_NodeGrid; }
		Grid                const &GetCellGrid () const { return m_CellGrid; }
		std::array<Grid, 2> const &GetFaceGrids() const { return m_FaceGrids; }

		Grid CreateRefinedCellGrid(int factor) const { return Grid(m_Spacing / factor, m_Resolution * factor, m_Origin + Vector2d::Constant(m_Spacing / factor / 2)); }

		Vector2d GetDomainOrigin () const { return m_Origin + Vector2d::Constant(m_BoundaryWidth * m_Spacing); }
		Vector2d GetDomainLengths() const { return (m_Resolution - Vector2i::Constant(m_BoundaryWidth * 2)).cast<double>() * m_Spacing; }
		double   GetDomainRadius () const { return m_Resolution.y() * m_Spacing * .5; }

		bool IsInsideNode  (Vector2i const &node)           const { return m_NodeGrid.IsInside(node, m_BoundaryWidth); }
		bool IsBoundaryNode(Vector2i const &node)           const { return !m_NodeGrid.IsInside(node, m_BoundaryWidth + 1); }
		bool IsInsideCell  (Vector2i const &cell)           const { return m_CellGrid.IsInside(cell, m_BoundaryWidth); }
		bool IsBoundaryCell(Vector2i const &cell)           const { return !IsInsideCell(cell); }
		bool IsInsideFace  (int axis, Vector2i const &face) const { return m_FaceGrids[axis].IsInside(face, m_BoundaryWidth); }
		bool IsBoundaryFace(int axis, Vector2i const &face) const { return face[axis] <= m_BoundaryWidth || face[axis] >= m_Resolution[axis] - m_BoundaryWidth || !IsInsideFace(axis, face); }

		static constexpr int GetNumFacesPerCell() { return 4; }
		static constexpr int GetNumEdgesPerCell() { return 4; }
		static constexpr int GetNumNodesPerCell() { return 4; }
		static constexpr int GetNumNodesPerFace() { return 2; }

		static constexpr int FaceAxisOfCell(int ord) { return ord >> 1; }
		static constexpr int FaceSideOfCell(int ord) { return ord & 1 ? 1: -1; }
		static std::pair<int, Vector2i> FaceOfCell(Vector2i const &cell, int ord) { return { FaceAxisOfCell(ord), cell + Vector2i::Unit(FaceAxisOfCell(ord)) * (ord & 1) }; }
		static Vector2i AdjCellOfFace(int axis, Vector2i const &face, int ord) { return face - Vector2i::Unit(axis) * (ord & 1 ^ 1); }

		static constexpr int EdgeAxisOfCell(int ord) { return ord >> 1; }
		static std::pair<int, Vector2i> EdgeOfCell(Vector2i const &cell, int ord) { return { EdgeAxisOfCell(ord), cell + Vector2i::Unit(EdgeAxisOfCell(ord) ^ 1) * (ord & 1) }; }

		static Vector2i NodeOfCell(Vector2i const &cell, int ord) { return cell + Vector2i(ord & 1, ord >> 1 & 1); }
		static Vector2i NodeOfFace(int axis, Vector2i const &face, int ord) { return face + Vector2i::Unit(axis ^ 1) * (ord & 1); }
		static Vector2i NodeOfEdge(int axis, Vector2i const &edge, int ord) { return edge + Vector2i::Unit(axis) * (ord & 1); }
		
	private:
		int      m_BoundaryWidth;
		double   m_Spacing;
		double   m_InvSpacing;
		Vector2i m_Resolution;
		Vector2d m_Origin;

		Grid                m_NodeGrid;
		Grid                m_CellGrid;
		std::array<Grid, 2> m_FaceGrids;
	};
}
