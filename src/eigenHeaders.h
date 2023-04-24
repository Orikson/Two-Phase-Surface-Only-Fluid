#ifndef EIGEN_HEADERS_H
#define EIGEN_HEADERS_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <surftrack.h>

typedef Eigen::Matrix<double, 4, 4> mat4;
typedef Eigen::Matrix<double, 3, 3> mat3;
typedef Eigen::Matrix<double, 3, 2> mat32;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matX;

typedef Eigen::Matrix<double, 4, 1> vec4;
typedef Eigen::Matrix<double, 3, 1> vec3;
typedef Eigen::Matrix<double, 2, 1> vec2;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vecX;

typedef Eigen::Matrix<int, 3, 1> vec3i;
typedef Eigen::Matrix<int, 2, 1> vec2i;

typedef Eigen::Matrix<float, 3, 1> vec3f;
typedef Eigen::Matrix<unsigned int, 3, 1> vec3ui;

typedef Eigen::SparseMatrix<double> SparseMatrix;

vec3f df(const vec3& v);
vec3 fd(const vec3f& v);
vec3 vc(const LosTopos::Vec3d& v);
LosTopos::Vec3d vc(const vec3& v);

vec3 clamp(vec3 v, vec3 c);

class Vec2iComp {
public:
    bool operator () (const vec2i& v1, const vec2i& v2) const { return v1[0] < v2[0] || (v1[0] == v2[0] && v1[1] < v2[1]); }
};

#endif