#ifndef FLUID_OPTIONS_H
#define FLUID_OPTIONS_H

#include <eigenHeaders.h>
#include <string>
using std::string;

struct FluidOptions {
	// Scene physical quantity parameters
	double dt;
	double influx, influx_inc;
	double sigma, sigma_sa, sigma_sl;		// Surface tension coefficient (liquid - liquid, solid - air, solid - liquid)
	double rho;								// Density
	vec3 gravity;

	// Scene physical boundary parameters
	// bbox has precedence over floorOnly
	bool bbox;			// Fluid is contained within a box
	vec3 bbox_dim, bbox_center;

	bool floorOnly;		// Bounding box only consists of a floor
	double floorz;

	double mean_edge_length;
	double smoothing_coefficient;
	double view_range;

	bool ibr;			// Internal bubble removal (eventually should be turned off)
	bool tpcf;			// Triple junction positional constraint
	bool tdmc;			// Tiny droplet momentum conservation

	// Scene utility parameters
	int remesh_iter;
};

class OptionsLoader {
	public:
		OptionsLoader() {

		}
		~OptionsLoader() {

		}

		struct FluidOptions load(string fname) {

		}

	private:
		
};

#endif