#ifndef FORESTFIRE_HPP_
#define FORESTFIRE_HPP_

#include <cmath>
#include <nlohmann/json.hpp>
#include <cadmium/modeling/celldevs/grid/cell.hpp>
#include <cadmium/modeling/celldevs/grid/config.hpp>
#include "cadmium/simulation/core/coordinator.hpp"
#include "cadmium/simulation/core/abs_simulator.hpp"
#include "forestfireState.hpp"

using namespace cadmium::celldevs;

//! Rules for the cells
class forestfire : public GridCell<forestfireState, double> {
	//protected:
	//	std::shared_ptr<cadmium::Coordinator> topCoordinator;  //!< Pointer to top coordinator.
	public:
	forestfire(const std::vector<int>& id, 
			const std::shared_ptr<const GridCellConfig<forestfireState, double>>& config
		  ): GridCell<forestfireState, double>(id, config) {}

	[[nodiscard]] forestfireState localComputation(forestfireState state,
		const std::unordered_map<std::vector<int>, NeighborData<forestfireState, double>>& neighborhood) const override {
		
			double new_temp = 0;
			double time;

			for (const auto& [neighborId, neighborData]: neighborhood) {
				auto nState = neighborData.state;
				double v = neighborData.vicinity;
				double nTemp = nState->temp;
				double t_change = nTemp*v;
				new_temp += t_change;
			}

		if (state.fire_status == 0) {
			state.temp = 300.0;
		}

		else if(state.fire_status == 1) {
			if(state.temp > 300.0) {
				double unburned_temp = new_temp+0.213;
				state.temp = unburned_temp;
				state.fire_status = 2;	
			}
			else if(state.temp == 300.0){
				double unburned_temp = new_temp+0.213;
				state.temp = unburned_temp;
			}
		}	
		else if(state.fire_status == 2) {
			if(state.temp < 573.0) {
				double unburned_temp = new_temp+0.213;		
				state.temp = unburned_temp;
				}	
			else if(state.temp >= 573.0) {
				state.fire_status = 3;
				state.t_ig = 1.0;
				double unburned_temp = new_temp+0.213;	
				state.temp = unburned_temp;
			}
		}
		else if(state.fire_status == 3) {
			if(state.temp > 333.0) {
				double burning_temp = new_temp+0.213+2.74*exp(-0.0019*(state.t_ig));		
				state.temp = burning_temp;
				state.t_ig += 1.0;
			}	
			else if(state.temp <= 333.0) {
				state.fire_status = 4;	
				state.temp = 300.0;
			}
		}
		else if(state.fire_status == 4) {
			double burned_temp = new_temp+0.213;		
			state.temp = burned_temp;
		}
		
		return state;
	}

	[[nodiscard]] double outputDelay(const forestfireState& state) const override {
		return 1.0;   
	}
};

#endif
