#pragma once

#include "Simulation.h"

namespace Pivot {
	struct SimBuildOptions {
		Simulation::Scene Scene;
		int               Scale;
	};

	class SimBuilder {
	public:
		static std::unique_ptr<Simulation> Build(SimBuildOptions const &options);
	
	private:
		static std::unique_ptr<Simulation> BuildBox(SimBuildOptions const &options);
	};
}
