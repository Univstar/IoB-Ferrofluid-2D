#include "Driver.h"

#include "StopWatch.h"

namespace Pivot {
	template <typename Duration>
	static std::string DurationFormat(Duration d) {
		auto ms = std::chrono::ceil<std::chrono::milliseconds>(d).count();
		auto const h = static_cast<int>(ms / 3600000);
		ms %= 3600000;
		auto const m = static_cast<int>(ms / 60000);
		ms %= 60000;
		auto const s = static_cast<int>(ms / 1000);
		ms %= 1000;

		if (h) {
			return fmt::format("{}h{}m{}.{:0>3}s", h, m, s, ms);
		} else if (m) {
			return fmt::format("{}m{}.{:0>3}s", m, s, ms);
		} else {
			return fmt::format("{}.{:0>3}s", s, ms);
		}
	}

	Driver::Driver(DriverCreateOptions const &options) :
		m_Dirname { std::filesystem::absolute(options.Dirname) },
		m_BeginFrame { options.BeginFrame },
		m_EndFrame { options.EndFrame },
		m_SecondPerFrame { 1. / options.FrameRate },
		m_CourantNumber { options.CourantNumber } {
	}

	void Driver::Run(Simulation *simulation) const {
		using Clock = std::chrono::steady_clock;
		auto const initTime = Clock::now();

		if (m_BeginFrame == 0) {
			fmt::print(fmt::fg(fmt::color::yellow_green), "[Initialize] ");
			{ // Initialize simulation
				auto sw = StopWatch("init.");
				simulation->Initialize();
				fmt::print("    ... {:>8.3f}s used\n", sw.Stop());
			}
			if (std::filesystem::is_directory(m_Dirname)) {
				spdlog::info("Output to the existing directory \"{}\"", m_Dirname.string());
			} else {
				std::filesystem::create_directories(m_Dirname);
				spdlog::info("Output to a new directory \"{}\"", m_Dirname.string());
			}
			std::filesystem::create_directory(m_Dirname);
			ExportAndSaveFrame(simulation, 0);
		} else {
			return;
		}

		// Initialize timing
		auto const beginTime = Clock::now();
		spdlog::info("Begin simulating (elapsed time: {})\n", DurationFormat(beginTime - initTime));
		auto lastTime = beginTime;
		auto beginFrame = std::max(m_BeginFrame, 1U);
		for (auto frame = beginFrame; frame < m_EndFrame; frame++) {
			// Simulate
			spdlog::info("Start to simulate Frame {}", frame);
			AdvanceTimeBySteps(simulation, frame);
			// Export and save files for the frame
			ExportAndSaveFrame(simulation, frame);
			// Output timing
			auto const currentTime = Clock::now();
			auto const frameTime = currentTime - lastTime;
			auto const totalTime = currentTime - initTime;
			spdlog::info("Finish simulating Frame {} (elapsed time: {}/{})", frame, DurationFormat(frameTime), DurationFormat(totalTime));
			auto const elapsedRatio = (m_EndFrame - beginFrame) / (frame - beginFrame + 1.);
			auto const prediTime = (currentTime - beginTime) * elapsedRatio + (beginTime - initTime);
			spdlog::info("Estimated total time: {}\n", DurationFormat(prediTime));
			lastTime = currentTime;
		}
		spdlog::info("Completed simulating! (elapsed time: {})", DurationFormat(lastTime - initTime));
		StopWatch::PrintStats();
	}

	void Driver::ExportAndSaveFrame(Simulation *simulation, std::uint32_t frame) const {
		spdlog::info("Export results of Frame {}", frame);
		{ // Export results
			auto const filename = m_Dirname / (std::to_string(frame) + ".png");
			simulation->Export(filename);
		}
	}

	void Driver::AdvanceTimeBySteps(Simulation *simulation, std::uint32_t frame) const {
		auto const startTime = (frame - 1) * m_SecondPerFrame;
		double time = 0;
		bool done = (time >= m_SecondPerFrame);
		while (!done) {
			simulation->SetTime(startTime + time);			
			// Calculate delta time
			auto deltaTime = std::min(m_SecondPerFrame, simulation->GetCourantTimeStep() * m_CourantNumber);
			if (time + deltaTime >= m_SecondPerFrame) {
				deltaTime = m_SecondPerFrame - time;
				done = true;
			} else if (time + 2 * deltaTime >= m_SecondPerFrame) {
				deltaTime = (m_SecondPerFrame - time) * .5; 
			}
			fmt::print(fmt::fg(fmt::color::yellow_green), "[{:>6.0f} SPF] ", m_SecondPerFrame / deltaTime);
			{ // Advance with timing
				auto sw = StopWatch("alg.");
				simulation->Advance(deltaTime);
				fmt::print("    ... {:>8.3f}s used\n", sw.Stop());
			}
			time += deltaTime;
		}
	}
}
