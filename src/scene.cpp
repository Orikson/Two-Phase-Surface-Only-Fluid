#include "scene.h"

void createUVSphere(const vec3& center, double r, int N, int M, vector<vec3>& vs, vector<vec3i>& fs, vector<vec2i>& ls, const vec2i& label) {
	vs.push_back(center + r * vec3(0, 0, 1));
	for (int i = 1; i < N; i++)
		for (int j = 0; j < M; j++)
			vs.push_back(center + r * vec3(sin(i * M_PI / N) * cos(j * 2 * M_PI / M), sin(i * M_PI / N) * sin(j * 2 * M_PI / M), cos(i * M_PI / N)));
	vs.push_back(center + r * vec3(0, 0, -1));

	for (int j = 0; j < M; j++)
		fs.push_back(vec3i(0, j + 1, (j + 1) % M + 1)), ls.push_back(label);
	for (int i = 0; i < N - 2; i++)
		for (int j = 0; j < M; j++)
			fs.push_back(vec3i(i * M + j + 1, (i + 1) * M + j + 1, (i + 1) * M + (j + 1) % M + 1)), ls.push_back(label),
			fs.push_back(vec3i((i + 1) * M + (j + 1) % M + 1, i * M + (j + 1) % M + 1, i * M + j + 1)), ls.push_back(label);
	for (int j = 0; j < M; j++)
		fs.push_back(vec3i((N - 1) * M + 1, (N - 2) * M + (j + 1) % M + 1, (N - 2) * M + j + 1)), ls.push_back(label);
}

LosTopos::SurfTrack* s_crownsplash(vector<LosTopos::Vec3d>& vs, vector<LosTopos::Vec3st>& fs, vector<LosTopos::Vec2i>& ls, vector<vec3>& vels, vector<char>& solids, struct FluidOptions& options) {
	int M = 48; int N = 12;
	float dt = 0.01;
	double radius = 1.0; double initVelZ = -1.0;

	std::vector<vec3> v;
	std::vector<vec3i> f;
	std::vector<vec2i> l;

	// pool
	double r_pool = 40.0;
	double d_pool = 40.0;
	double r_droplet = radius;
	double z_droplet = radius - initVelZ * dt + 1e-6;

	// pool
	v.push_back(vec3(0, 0, 0));
	for (int i = 1; i < N; i++)
		for (int j = 0; j < M; j++)
			v.push_back(vec3(r_pool * i / N * cos(j * 2 * M_PI / M), r_pool * i / N * sin(j * 2 * M_PI / M), 0));
	for (int i = 0; i <= N; i++)
		for (int j = 0; j < M; j++)
			v.push_back(vec3(r_pool * cos(j * 2 * M_PI / M), r_pool * sin(j * 2 * M_PI / M), -d_pool * (double)i / N));
	v.push_back(vec3(0, 0, -d_pool));

	for (int j = 0; j < M; j++)
		f.push_back(vec3i(0, j + 1, (j + 1) % M + 1)), l.push_back(vec2i(1, 0));
	for (int i = 0; i < N * 2 - 1; i++)
		for (int j = 0; j < M; j++)
			f.push_back(vec3i(i * M + j + 1, (i + 1) * M + j + 1, (i + 1) * M + (j + 1) % M + 1)), l.push_back(vec2i(1, 0)),
			f.push_back(vec3i((i + 1) * M + (j + 1) % M + 1, i * M + (j + 1) % M + 1, i * M + j + 1)), l.push_back(vec2i(1, 0));
	for (int j = 0; j < M; j++)
		f.push_back(vec3i((2 * N) * M + 1, (2 * N - 1) * M + (j + 1) % M + 1, (2 * N - 1) * M + j + 1)), l.push_back(vec2i(1, 0));

	size_t v0 = vs.size();
	vs.reserve(vs.size() + v.size()); for (size_t i = 0; i < v.size(); i++) vs.push_back(LosTopos::Vec3d(v[i][0], v[i][1], v[i][2]));
	fs.reserve(fs.size() + f.size()); for (size_t i = 0; i < f.size(); i++) fs.push_back(LosTopos::Vec3st(v0 + f[i][0], v0 + f[i][1], v0 + f[i][2]));
	ls.reserve(ls.size() + l.size()); for (size_t i = 0; i < l.size(); i++) ls.push_back(LosTopos::Vec2i(l[i][0], l[i][1]));

	vels.reserve(vels.size() + v.size());
	for (size_t i = 0; i < v.size(); i++) vels.push_back(vec3(0, 0, 0));

	solids.reserve(solids.size() + v.size());
	for (size_t i = 0; i < v.size(); i++) solids.push_back(true);
	for (int i = 0; i < N - 1; i++) for (int j = 0; j < M; j++) solids[v0 + i * M + j + 1] = false;
	solids[v0] = false;

	// droplet
	v.clear();
	f.clear();
	l.clear();
	
	createUVSphere(vec3(0, 0, z_droplet), r_droplet, 16, 32, v, f, l);

	v0 = vs.size();
	vs.reserve(vs.size() + v.size()); for (size_t i = 0; i < v.size(); i++) vs.push_back(LosTopos::Vec3d(v[i][0], v[i][1], v[i][2]));
	fs.reserve(fs.size() + f.size()); for (size_t i = 0; i < f.size(); i++) fs.push_back(LosTopos::Vec3st(v0 + f[i][0], v0 + f[i][1], v0 + f[i][2]));
	ls.reserve(ls.size() + l.size()); for (size_t i = 0; i < l.size(); i++) ls.push_back(LosTopos::Vec2i(l[i][0], l[i][1]));

	vels.reserve(vels.size() + v.size());
	for (size_t i = 0; i < v.size(); i++) vels.push_back(vec3(0, 0, -1));

	solids.reserve(solids.size() + v.size());
	for (size_t i = 0; i < v.size(); i++) solids.push_back(false);

	// construct the surface tracker
	double target_edge_len = 0.3;
	//m_sim_options.mean_edge_length = target_edge_len;

	LosTopos::SurfTrackInitializationParameters params;
	params.m_proximity_epsilon = 1e-4 * target_edge_len;
	params.m_merge_proximity_epsilon = 0.1 * target_edge_len;
	params.m_merge_proximity_epsilon_for_liquid_sheet_puncture = params.m_merge_proximity_epsilon * 1.0;
	params.m_allow_vertex_movement_during_collapse = true;
	params.m_perform_smoothing = false;
	//    params.m_min_to_target_ratio = 0.5;   // default
	//    params.m_max_to_target_ratio = 1.5;   // default
	params.m_max_adjacent_target_edge_length_ratio = 1.5;   // default: 1.5
	params.m_target_edge_length_coef_curvature = 0.05;   // default: 0.05
	params.m_target_edge_length_coef_velocity = 0.01 / dt; // default: 0.01/dt
	params.m_min_target_edge_length = target_edge_len;
	params.m_max_target_edge_length = target_edge_len * 3;
	params.m_refine_solid = false;
	params.m_refine_triple_junction = false;
	params.m_max_volume_change = 1e-2 * pow(target_edge_len, 3);
	params.m_min_triangle_angle = 3;
	params.m_max_triangle_angle = 177;
	params.m_large_triangle_angle_to_split = 160;
	params.m_min_triangle_area = 0.02 * pow(target_edge_len, 2);
	params.m_verbose = false;
	params.m_allow_non_manifold = true;
	params.m_allow_topology_changes = true;
	params.m_collision_safety = true;
	params.m_remesh_boundaries = true;
	params.m_t1_transition_enabled = true;
	params.m_pull_apart_distance = 0.1 * target_edge_len;

	params.m_velocity_field_callback = NULL;

	// alternatively, use LosTopos::MidpointScheme()
	params.m_subdivision_scheme = new LosTopos::ButterflyScheme();

	params.m_use_curvature_when_collapsing = false;
	params.m_use_curvature_when_splitting = false;

	// finalize surface tracker
	vector<LosTopos::Vec3d> masses(vs.size(), LosTopos::Vec3d(1, 1, 1));
	for (size_t i = 0; i < vs.size(); i++)  // the solids labels supercede the constrained vertices; the latter has really no effect.
		masses[i] = LosTopos::Vec3d(1, 1, solids[i] ? std::numeric_limits<double>::infinity() : 1); //&&&& for now assume that the solid surface is always perpendicular to the z axis. For generality, LosTopos's impulse code will need a major overhaul anyway.
	
	// construct fluid options
	options.bbox = false;
	options.floorOnly = true;
	options.floorz = -1e+10;
	options.gravity = vec3(0, 0, 0);
	options.smoothing_coefficient = 20;
	options.sigma = 1;
	options.sigma_sa = 1;
	options.sigma_sl = 1;
	options.rho = 1.0;
	options.dt = dt;
	options.influx = 0;
	options.influx_inc = 0;
	options.remesh_iter = 5;
	options.view_range = 20;
	options.ibr = true;
	options.tpcf = true;
	options.tdmc = true;

	return new LosTopos::SurfTrack(vs, fs, ls, masses, params);
}

LosTopos::SurfTrack* s_ball(vector<LosTopos::Vec3d>& vs, vector<LosTopos::Vec3st>& fs, vector<LosTopos::Vec2i>& ls, vector<vec3>& vels, vector<char>& solids, struct FluidOptions& options) {
	int M = 24; int N = 12;
	double dt = 0.01;

	// Create sphere
	std::vector<vec3> v;
	std::vector<vec3i> f;
	std::vector<vec2i> l;

	double r = 1.0;
	double h = 1.1;
	//        createIcoSphere(Vec3d(0, 0, h), r, N, v, f, l);
	createUVSphere(vec3(0, 0, h), r, N, M, v, f, l);

	vs.resize(v.size()); for (size_t i = 0; i < v.size(); i++) vs[i] = LosTopos::Vec3d(v[i][0], v[i][1], v[i][2]);
	fs.resize(f.size()); for (size_t i = 0; i < f.size(); i++) fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
	ls.resize(l.size()); for (size_t i = 0; i < l.size(); i++) ls[i] = LosTopos::Vec2i(l[i][0], l[i][1]);

	vels.resize(v.size());
	for (size_t i = 0; i < v.size(); i++) vels[i] = vec3(0, 0, -2);
	solids.resize(v.size(), false);

	// construct the surface tracker
	double target_edge_len = 0;
	double mean_edge_len = 0;
	for (size_t i = 0; i < fs.size(); i++) {
		mean_edge_len += mag(vs[fs[i][0]] - vs[fs[i][1]]);
		mean_edge_len += mag(vs[fs[i][1]] - vs[fs[i][2]]);
		mean_edge_len += mag(vs[fs[i][2]] - vs[fs[i][0]]);
	}
	mean_edge_len /= (fs.size() * 3);
	target_edge_len = mean_edge_len;

	//m_sim_options.mean_edge_length = target_edge_len;

	LosTopos::SurfTrackInitializationParameters params;
	params.m_proximity_epsilon = 1e-4 * target_edge_len;
	params.m_merge_proximity_epsilon = 0.1 * target_edge_len;
	params.m_merge_proximity_epsilon_for_liquid_sheet_puncture = params.m_merge_proximity_epsilon * 1.0;
	params.m_allow_vertex_movement_during_collapse = true;
	params.m_perform_smoothing = false;
	//    params.m_min_to_target_ratio = 0.5;   // default
	//    params.m_max_to_target_ratio = 1.5;   // default
	params.m_max_adjacent_target_edge_length_ratio = 1.5;   // default: 1.5
	params.m_target_edge_length_coef_curvature = 0.05;   // default: 0.05
	params.m_target_edge_length_coef_velocity = 0.01 / dt; // default: 0.01/dt
	params.m_min_target_edge_length = target_edge_len;
	params.m_max_target_edge_length = target_edge_len * 100;
	params.m_refine_solid = true;
	params.m_refine_triple_junction = true;
	params.m_max_volume_change = 1e-2 * pow(target_edge_len, 3);
	params.m_min_triangle_angle = 3;
	params.m_max_triangle_angle = 177;
	params.m_large_triangle_angle_to_split = 160;
	params.m_min_triangle_area = 0.02 * pow(target_edge_len, 2);
	params.m_verbose = false;
	params.m_allow_non_manifold = false;
	params.m_allow_topology_changes = true;
	params.m_collision_safety = true;
	params.m_remesh_boundaries = true;
	params.m_t1_transition_enabled = true;
	params.m_pull_apart_distance = 0.1 * target_edge_len;

	params.m_velocity_field_callback = NULL;

	// alternatively, use LosTopos::MidpointScheme()
	params.m_subdivision_scheme = new LosTopos::ButterflyScheme();

	params.m_use_curvature_when_collapsing = false;
	params.m_use_curvature_when_splitting = false;

	// finalize surface tracker
	vector<LosTopos::Vec3d> masses(vs.size(), LosTopos::Vec3d(1, 1, 1));
	for (size_t i = 0; i < vs.size(); i++)  // the solids labels supercede the constrained vertices; the latter has really no effect.
		masses[i] = LosTopos::Vec3d(1, 1, solids[i] ? std::numeric_limits<double>::infinity() : 1); //&&&& for now assume that the solid surface is always perpendicular to the z axis. For generality, LosTopos's impulse code will need a major overhaul anyway.

	// construct fluid options
	options.bbox = false;
	options.floorOnly = true;
	options.floorz = 0;
	options.gravity = vec3(0, 0, 0);
	options.smoothing_coefficient = 20;
	options.sigma = 0.5;
	options.sigma_sa = 0.0;
	options.sigma_sl = 0.0;
	options.rho = 1.0;
	options.dt = dt;
	options.influx = 0;
	options.influx_inc = 0;
	options.remesh_iter = 5;
	options.view_range = 20;
	options.ibr = true;
	options.tpcf = true;
	options.tdmc = true;

	return new LosTopos::SurfTrack(vs, fs, ls, masses, params);
}

LosTopos::SurfTrack* s_collision(vector<LosTopos::Vec3d>& vs, vector<LosTopos::Vec3st>& fs, vector<LosTopos::Vec2i>& ls, vector<vec3>& vels, vector<char>& solids, struct FluidOptions& options) {
	int M = 20; int N = 10;
	double dt = 0.01;

	// Create spheres
	std::vector<vec3> v;
	std::vector<vec3i> f;
	std::vector<vec2i> l;

	double r = 1;
	double d = 2.2;
	double od = 0;

	//        createIcoSphere(Vec3d(-1.1, 0, 0), r, N, v, f, l);
	createUVSphere(vec3(-d / 2, 0, -od / 2), r, N, M, v, f, l);

	size_t v0 = vs.size();
	vs.reserve(vs.size() + v.size()); for (size_t i = 0; i < v.size(); i++) vs.push_back(LosTopos::Vec3d(v[i][0], v[i][1], v[i][2]));
	fs.reserve(fs.size() + f.size()); for (size_t i = 0; i < f.size(); i++) fs.push_back(LosTopos::Vec3st(v0 + f[i][0], v0 + f[i][1], v0 + f[i][2]));
	ls.reserve(ls.size() + l.size()); for (size_t i = 0; i < l.size(); i++) ls.push_back(LosTopos::Vec2i(l[i][0], l[i][1]));

	vels.reserve(vels.size() + v.size());
	for (size_t i = 0; i < v.size(); i++) vels.push_back(vec3(1, 0, 0));

	v.clear();
	f.clear();
	l.clear();

	//        createIcoSphere(Vec3d( 1.1, 0, 0), r, N, v, f, l);
	createUVSphere(vec3(d / 2, 0, od / 2), r, N, M, v, f, l);

	v0 = vs.size();
	vs.reserve(vs.size() + v.size()); for (size_t i = 0; i < v.size(); i++) vs.push_back(LosTopos::Vec3d(v[i][0], v[i][1], v[i][2]));
	fs.reserve(fs.size() + f.size()); for (size_t i = 0; i < f.size(); i++) fs.push_back(LosTopos::Vec3st(v0 + f[i][0], v0 + f[i][1], v0 + f[i][2]));
	ls.reserve(ls.size() + l.size()); for (size_t i = 0; i < l.size(); i++) ls.push_back(LosTopos::Vec2i(l[i][0], l[i][1]));

	vels.reserve(vels.size() + v.size());
	for (size_t i = 0; i < v.size(); i++) vels.push_back(vec3(-1, 0, 0));

	solids.resize(vs.size(), false);

	// construct the surface tracker
	double target_edge_len = 0;
	double mean_edge_len = 0;
	for (size_t i = 0; i < fs.size(); i++) {
		mean_edge_len += mag(vs[fs[i][0]] - vs[fs[i][1]]);
		mean_edge_len += mag(vs[fs[i][1]] - vs[fs[i][2]]);
		mean_edge_len += mag(vs[fs[i][2]] - vs[fs[i][0]]);
	}
	mean_edge_len /= (fs.size() * 3);
	target_edge_len = mean_edge_len;

	//m_sim_options.mean_edge_length = target_edge_len;

	LosTopos::SurfTrackInitializationParameters params;
	params.m_proximity_epsilon = 1e-4 * target_edge_len;
	params.m_merge_proximity_epsilon = 0.1 * target_edge_len;
	params.m_merge_proximity_epsilon_for_liquid_sheet_puncture = params.m_merge_proximity_epsilon * 1.0;
	params.m_allow_vertex_movement_during_collapse = true;
	params.m_perform_smoothing = false;
	//    params.m_min_to_target_ratio = 0.5;   // default
	//    params.m_max_to_target_ratio = 1.5;   // default
	params.m_max_adjacent_target_edge_length_ratio = 1.5;   // default: 1.5
	params.m_target_edge_length_coef_curvature = 0.05;   // default: 0.05
	params.m_target_edge_length_coef_velocity = 0.01 / dt; // default: 0.01/dt
	params.m_min_target_edge_length = target_edge_len;
	params.m_max_target_edge_length = target_edge_len * 100;
	params.m_refine_solid = true;
	params.m_refine_triple_junction = true;
	params.m_max_volume_change = 1e-2 * pow(target_edge_len, 3);
	params.m_min_triangle_angle = 3;
	params.m_max_triangle_angle = 177;
	params.m_large_triangle_angle_to_split = 160;
	params.m_min_triangle_area = 0.02 * pow(target_edge_len, 2);
	params.m_verbose = false;
	params.m_allow_non_manifold = true;
	params.m_allow_topology_changes = true;
	params.m_collision_safety = true;
	params.m_remesh_boundaries = true;
	params.m_t1_transition_enabled = true;
	params.m_pull_apart_distance = 0.1 * target_edge_len;

	params.m_velocity_field_callback = NULL;

	// alternatively, use LosTopos::MidpointScheme()
	params.m_subdivision_scheme = new LosTopos::ButterflyScheme();

	params.m_use_curvature_when_collapsing = false;
	params.m_use_curvature_when_splitting = false;

	// finalize surface tracker
	vector<LosTopos::Vec3d> masses(vs.size(), LosTopos::Vec3d(1, 1, 1));
	for (size_t i = 0; i < vs.size(); i++)  // the solids labels supercede the constrained vertices; the latter has really no effect.
		masses[i] = LosTopos::Vec3d(1, 1, solids[i] ? std::numeric_limits<double>::infinity() : 1); //&&&& for now assume that the solid surface is always perpendicular to the z axis. For generality, LosTopos's impulse code will need a major overhaul anyway.

	// construct fluid options
	options.bbox = false;
	options.floorOnly = true;
	options.floorz = -1e+10;
	options.gravity = vec3(0, 0, 0);
	options.smoothing_coefficient = 20;
	options.sigma = 0.25;
	options.sigma_sa = 0.25;
	options.sigma_sl = 1.0;
	options.rho = 1.0;
	options.dt = dt;
	options.influx = 0;
	options.influx_inc = 0;
	options.remesh_iter = 5;
	options.view_range = 20;
	options.ibr = true;
	options.tpcf = true;
	options.tdmc = true;

	return new LosTopos::SurfTrack(vs, fs, ls, masses, params);
}

LosTopos::SurfTrack* s_dripping(vector<LosTopos::Vec3d>& vs, vector<LosTopos::Vec3st>& fs, vector<LosTopos::Vec2i>& ls, vector<vec3>& vels, vector<char>& solids, struct FluidOptions& options) {
	int N = 12, M = 24;
	double dt = 0.01;

	// dripping
	std::vector<vec3> v;
	std::vector<vec3i> f;
	std::vector<vec2i> l;

	double alpha_end = M_PI / 3;
	double r = 2;
	v.push_back(vec3(0, 0, -r));
	for (int i = 1; i <= N; i++)
		for (int j = 0; j < M; j++)
			v.push_back(vec3(sin(i * alpha_end / N) * cos((j - 0.0 * i) * 2 * M_PI / M), -sin(i * alpha_end / N) * sin((j - 0.0 * i) * 2 * M_PI / M), -cos(i * alpha_end / N)) * r);
	v.push_back(vec3(0, 0, -cos(alpha_end) * r));

	for (int j = 0; j < M; j++)
		f.push_back(vec3i(0, j + 1, (j + 1) % M + 1)), l.push_back(vec2i(1, 0));
	for (int i = 0; i < N - 1; i++)
		for (int j = 0; j < M; j++)
			f.push_back(vec3i(i * M + j + 1, (i + 1) * M + j + 1, (i + 1) * M + (j + 1) % M + 1)), l.push_back(vec2i(1, 0)),
			f.push_back(vec3i((i + 1) * M + (j + 1) % M + 1, i * M + (j + 1) % M + 1, i * M + j + 1)), l.push_back(vec2i(1, 0));
	for (int j = 0; j < M; j++)
		f.push_back(vec3i(N * M + 1, (N - 1) * M + (j + 1) % M + 1, (N - 1) * M + j + 1)), l.push_back(vec2i(1, 0));

	vs.resize(v.size()); for (size_t i = 0; i < v.size(); i++) vs[i] = LosTopos::Vec3d(v[i][0], v[i][1], v[i][2]);
	fs.resize(f.size()); for (size_t i = 0; i < f.size(); i++) fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
	ls.resize(l.size()); for (size_t i = 0; i < l.size(); i++) ls[i] = LosTopos::Vec2i(l[i][0], l[i][1]);

	vels.resize(v.size());
	for (size_t i = 0; i < v.size(); i++) vels[i] = vec3(0, 0, 0);

	solids.resize(v.size(), false);
	for (size_t j = 0; j < M; j++)
		solids[(N - 1) * M + j + 1] = true;
	solids[N * M + 1] = true;

	// construct the surface tracker
	double target_edge_len = 0;
	double mean_edge_len = 0;
	for (size_t i = 0; i < fs.size(); i++) {
		mean_edge_len += mag(vs[fs[i][0]] - vs[fs[i][1]]);
		mean_edge_len += mag(vs[fs[i][1]] - vs[fs[i][2]]);
		mean_edge_len += mag(vs[fs[i][2]] - vs[fs[i][0]]);
	}
	mean_edge_len /= (fs.size() * 3);
	target_edge_len = mean_edge_len;

	//m_sim_options.mean_edge_length = target_edge_len;

	LosTopos::SurfTrackInitializationParameters params;
	params.m_proximity_epsilon = 1e-4 * target_edge_len;
	params.m_merge_proximity_epsilon = 0.1 * target_edge_len;
	params.m_merge_proximity_epsilon_for_liquid_sheet_puncture = params.m_merge_proximity_epsilon * 1.0;
	params.m_allow_vertex_movement_during_collapse = true;
	params.m_perform_smoothing = false;
	//    params.m_min_to_target_ratio = 0.5;   // default
	//    params.m_max_to_target_ratio = 1.5;   // default
	params.m_max_adjacent_target_edge_length_ratio = 1.5;   // default: 1.5
	params.m_target_edge_length_coef_curvature = 0.05;   // default: 0.05
	params.m_target_edge_length_coef_velocity = 0.01 / dt; // default: 0.01/dt
	params.m_min_target_edge_length = target_edge_len;
	params.m_max_target_edge_length = target_edge_len * 100;
	params.m_refine_solid = true;
	params.m_refine_triple_junction = true;
	params.m_max_volume_change = 1e-2 * pow(target_edge_len, 3);
	params.m_min_triangle_angle = 3;
	params.m_max_triangle_angle = 177;
	params.m_large_triangle_angle_to_split = 160;
	params.m_min_triangle_area = 0.02 * pow(target_edge_len, 2);
	params.m_verbose = false;
	params.m_allow_non_manifold = true;
	params.m_allow_topology_changes = true;
	params.m_collision_safety = true;
	params.m_remesh_boundaries = true;
	params.m_t1_transition_enabled = true;
	params.m_pull_apart_distance = 0.1 * target_edge_len;

	params.m_velocity_field_callback = NULL;

	// alternatively, use LosTopos::MidpointScheme()
	params.m_subdivision_scheme = new LosTopos::ButterflyScheme();

	params.m_use_curvature_when_collapsing = false;
	params.m_use_curvature_when_splitting = false;

	// finalize surface tracker
	vector<LosTopos::Vec3d> masses(vs.size(), LosTopos::Vec3d(1, 1, 1));
	for (size_t i = 0; i < vs.size(); i++)  // the solids labels supercede the constrained vertices; the latter has really no effect.
		masses[i] = LosTopos::Vec3d(1, 1, solids[i] ? std::numeric_limits<double>::infinity() : 1); //&&&& for now assume that the solid surface is always perpendicular to the z axis. For generality, LosTopos's impulse code will need a major overhaul anyway.

	// construct fluid options
	options.bbox = false;
	options.floorOnly = true;
	options.floorz = -1e+10;
	options.gravity = vec3(0, 0, -1);
	options.smoothing_coefficient = 20;
	options.sigma = 1.0;
	options.sigma_sa = 1.0;
	options.sigma_sl = 1.0;
	options.rho = 1.0;
	options.dt = dt;
	options.influx = 0;
	options.influx_inc = 0;
	options.remesh_iter = 3;
	options.view_range = 50;
	options.ibr = true;
	options.tpcf = true;
	options.tdmc = true;

	return new LosTopos::SurfTrack(vs, fs, ls, masses, params);
}