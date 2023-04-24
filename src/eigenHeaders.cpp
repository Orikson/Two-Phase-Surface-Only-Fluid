#include "eigenHeaders.h"

vec3f df(const vec3& v) {
	return vec3f((float)v[0], (float)v[1], (float)v[2]);
}

vec3 fd(const vec3f& v) {
	return vec3((double)v[0], (double)v[1], (double)v[2]);
}

vec3 vc(const LosTopos::Vec3d& v) {
	return vec3(v[0], v[1], v[2]);
}

LosTopos::Vec3d vc(const vec3& v) {
	return LosTopos::Vec3d(v[0], v[1], v[2]);
}

vec3 clamp(vec3 v, vec3 c) {
	assert(c[0] > 0 && c[1] > 0 && c[2] > 0);
	return vec3(
		v[0] < 0 ? 0 : (v[0] > c[0] ? c[0] : v[0]),
		v[1] < 0 ? 0 : (v[1] > c[1] ? c[1] : v[1]),
		v[2] < 0 ? 0 : (v[2] > c[2] ? c[2] : v[2])
	);
}