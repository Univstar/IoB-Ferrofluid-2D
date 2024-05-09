#pragma once

#include "Simulation.h"

namespace Pivot {
	struct DriverCreateOptions {
		std::filesystem::path Dirname;
		std::uint32_t         BeginFrame    = 0;
		std::uint32_t         EndFrame;
		double                FrameRate     = 25;
		double                CourantNumber = 1;
	};

	class Driver {
	public:
		explicit Driver(DriverCreateOptions const &options);

		void Run(Simulation *simulation) const;

	private:
		void ExportAndSaveFrame(Simulation *simulation, std::uint32_t frame) const;
		void AdvanceTimeBySteps(Simulation *simulation, std::uint32_t frame) const;

	private:
		std::filesystem::path m_Dirname;
		std::uint32_t         m_BeginFrame;
		std::uint32_t         m_EndFrame;
		double                m_SecondPerFrame;
		double                m_CourantNumber;
	};
}
