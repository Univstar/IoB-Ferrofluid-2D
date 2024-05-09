#include "SimBuilder.h"

#include "CSG.h"

namespace Pivot {
std::unique_ptr<Simulation> SimBuilder::Build(SimBuildOptions const &options) {
    std::unique_ptr<Simulation> simulation;
    switch (options.Scene) {
    case Simulation::Scene::Box:
        simulation = BuildBox(options);
        break;
    }
    simulation->m_Scene = options.Scene;
    return simulation;
}

std::unique_ptr<Simulation>
SimBuilder::BuildBox(SimBuildOptions const &options) {
    constexpr double length = .15;
    constexpr int bw = 2;
    int const scale = options.Scale < 0 ? 128 : options.Scale;
    StaggeredGrid sgrid(2, length / (scale - bw * 2), Vector2i(1, 1) * scale);
    Vector2d center = Vector2d::Unit(1) * length * -.3;
    Vector2d halfLength = Vector2d(length * 1., length * .2);
    auto sim = std::make_unique<Simulation>(sgrid);
    sim->m_GravityEnabled = true;
    sim->m_SurfaceTensionEnabled = true;
    sim->m_MagneticEnabled = true;
    CSG::Union(sim->m_LevelSet,
               ImplicitBox(center - halfLength, halfLength * 2));
    return sim;
}
} // namespace Pivot
