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
	double wind_dir;
	double wind_spd;
	public:
	
	forestfire(const std::vector<int>& id, 
			const std::shared_ptr<const GridCellConfig<forestfireState, double>>& config
		  ): GridCell<forestfireState, double>(id, config) {
			//wind direction, expressed in direction the wind comes from with in degrees with 0 as north, for example wind
			// from the east is 90, config in json file
			config->rawCellConfig.at("wind").get_to(wind_dir);
			//wind speed in mph from config in json file
			config->rawCellConfig.at("windspeed").get_to(wind_spd);
		  }
	
	[[nodiscard]] forestfireState localComputation(forestfireState state,
		const std::unordered_map<std::vector<int>, NeighborData<forestfireState, double>>& neighborhood) const override {
		
			double new_temp = 0;
			// Convert degrees to radians and adjust for difference between cardinal direction and cartesian
			double radians = (90-wind_dir) * (M_PI / 180.0);
			// Calculate x and y components using cosine and sine
			double x = cos(radians);
			double y = sin(radians);
			int myRow;
			int myCol;

			for (const auto& [neighborId, neighborData]: neighborhood) {
				double v = neighborData.vicinity;
					if (v == 0.98689){
						myRow = neighborId[0];
						myCol = neighborId[1];
					}
			}
			for (const auto& [neighborId, neighborData]: neighborhood) {
				auto nState = neighborData.state;
				double v = neighborData.vicinity;
				double nTemp = nState->temp;
				int x_rel = neighborId[1] - myCol;
				int y_rel = neighborId[0] - myRow;
				double t_change = (1+(wind_spd/30)*(x*x_rel-y*y_rel))*nTemp*v;
				new_temp += t_change;
			}

		if (state.fire_status == 0) {
			state.temp = 300.0;
			state.sigma = std::numeric_limits<double>::infinity();
		}

		else if(state.fire_status == 1) {
			if(state.temp == 300.0){
				double unburned_temp = new_temp+0.213;
				state.temp = unburned_temp;
				if(state.temp == 300){
					state.sigma = std::numeric_limits<double>::infinity();
				}
				else {
					double unburned_temp = new_temp+0.213;
					state.temp = unburned_temp;
					state.fire_status = 2;	
					state.sigma = 1.0;
				}
			}
		}	
		else if(state.fire_status == 2) {
			if(state.temp < 573.0) {
				double unburned_temp = new_temp+0.213;		
				state.temp = unburned_temp;
				if(state.temp == 300){
					state.sigma = std::numeric_limits<double>::infinity();
					state.fire_status = 1;	
				}	
				else {
					state.sigma = 1.0;
				}
			}	
			else if(state.temp >= 573.0) {
				double burning_temp = new_temp+0.213+2.74*exp(-0.0019*(state.t_ig));		
				state.fire_status = 3;
				state.t_ig = 1.0;
				state.temp = burning_temp;
				state.sigma = 1.0;
			}
		}
		else if(state.fire_status == 3) {
			if(state.temp > 333.0) {
				double burning_temp = new_temp+0.213+2.74*exp(-0.0019*(state.t_ig));		
				state.temp = burning_temp;
				state.t_ig += 1.0;
				state.sigma = 1.0;
			}	
			else if(state.temp <= 333.0) {
				double burned_temp = new_temp+0.213;
				state.fire_status = 4;	
				state.temp = burned_temp;
				state.sigma = 1.0;
			}
		}
		else if(state.fire_status == 4) {
			double burned_temp = new_temp+0.213;		
			state.temp = burned_temp;
			if(state.temp == 300){
				state.sigma = std::numeric_limits<double>::infinity();
			}
			else {
				state.sigma = 1.0;
			}
		}
		
		return state;
	}

	[[nodiscard]] double outputDelay(const forestfireState& state) const override {
		return state.sigma;   
	}
};

#endif
