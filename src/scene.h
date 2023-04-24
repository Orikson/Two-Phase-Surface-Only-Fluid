#ifndef SCENE_H
#define SCENE_H

#include <subdivisionscheme.h>
#include <surftrack.h>
#include <vector>
using std::vector;
#include <eigenHeaders.h>
#include <fluidOptions.h>

void createUVSphere(const vec3& center, double r, int N, int M, vector<vec3>& vs, vector<vec3i>& fs, vector<vec2i>& ls, const vec2i& label = vec2i(1, 0));

LosTopos::SurfTrack* s_crownsplash(vector<LosTopos::Vec3d>& vs, vector<LosTopos::Vec3st>& fs, vector<LosTopos::Vec2i>& ls, vector<vec3>& vels, vector<char>& solids, struct FluidOptions& options);
LosTopos::SurfTrack* s_ball(vector<LosTopos::Vec3d>& vs, vector<LosTopos::Vec3st>& fs, vector<LosTopos::Vec2i>& ls, vector<vec3>& vels, vector<char>& solids, struct FluidOptions& options);
LosTopos::SurfTrack* s_collision(vector<LosTopos::Vec3d>& vs, vector<LosTopos::Vec3st>& fs, vector<LosTopos::Vec2i>& ls, vector<vec3>& vels, vector<char>& solids, struct FluidOptions& options);
LosTopos::SurfTrack* s_dripping(vector<LosTopos::Vec3d>& vs, vector<LosTopos::Vec3st>& fs, vector<LosTopos::Vec2i>& ls, vector<vec3>& vels, vector<char>& solids, struct FluidOptions& options);

#endif