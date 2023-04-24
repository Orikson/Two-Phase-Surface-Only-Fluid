//
//  BoundaryIntegral.cpp
//  Droplet3D
//
//  Created by Fang Da on 10/10/15.
//
//

#include "BoundaryIntegral.h"

const std::vector<vec2> & BoundaryIntegral::quadrature_line(int N)
{
    static std::map<int, std::vector<vec2> > s_quadratures;
    
    std::map<int, std::vector<vec2> >::iterator i = s_quadratures.find(N);
    if (i == s_quadratures.end())
    {
        // compute the quadrature weights now
        std::vector<vec2> q;
        switch (N)
        {
                // N = 1~2: Gaussian quadrature rules
            case 1:
                q.push_back(vec2(0, 2));
                break;
            case 2:
                q.push_back(vec2(-sqrt(1.0 / 3.0), 1));
                q.push_back(vec2( sqrt(1.0 / 3.0), 1));
                break;
            default:
                // N > 2: uniform grid of N points (midpoint rule)
                for (int j = 0; j < N; j++)
                    q.push_back(vec2((j + 0.5) / N * 2 - 1, 2.0 / N));
                break;
        }
        
        std::pair<std::map<int, std::vector<vec2> >::iterator, bool> res = s_quadratures.insert(std::pair<int, std::vector<vec2> >(N, q));
        assert(res.second);
        i = res.first;
    }
    
    return i->second;
}

const std::vector<vec3> & BoundaryIntegral::quadrature_square(int N)
{
    static std::map<int, std::vector<vec3> > s_quadratures;
    
    std::map<int, std::vector<vec3> >::iterator i = s_quadratures.find(N);
    if (i == s_quadratures.end())
    {
        // compute the quadrature weights now
        std::vector<vec3> q;
        
        // NxN quadrature, outer product of the corresponding N-point line quadrature (Gaussian or uniform depending on N)
        const std::vector<vec2> & lq = quadrature_line(N);
        assert(lq.size() == N);
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
                q.push_back(vec3(lq[k].x(), lq[j].x(), lq[j].y() * lq[k].y()));
        
        std::pair<std::map<int, std::vector<vec3> >::iterator, bool> res = s_quadratures.insert(std::pair<int, std::vector<vec3> >(N, q));
        assert(res.second);
        i = res.first;
    }
    
    return i->second;
}

const std::vector<vec3> & BoundaryIntegral::quadrature_triangle(int N)
{
    static std::map<int, std::vector<vec3> > s_quadratures;
    
    std::map<int, std::vector<vec3> >::iterator i = s_quadratures.find(N);
    if (i == s_quadratures.end())
    {
        // compute the quadrature weights now
        std::vector<vec3> q;

        switch (N)
        {
            case 1:
                q.push_back(vec3(.3333333333333330 + .3333333333333330, .3333333333333330, 1.0000000000000000 / 2));
                break;
            case 3:
                q.push_back(vec3(.1666666666666670 + .1666666666666660, .1666666666666660,  .3333333333333329 / 2));
                q.push_back(vec3(.6666666666666670 + .1666666666666660, .1666666666666660,  .3333333333333329 / 2));
                q.push_back(vec3(.1666666666666670 + .6666666666666661, .6666666666666661,  .3333333333333329 / 2));
                break;
            default:
                // for the other N numbers: uniform grid (N*(N+1)/2 + N*(N-1)/2 = N^2 congruent triangles)
                // note that for N values other than those above, the number of quadrature points is N^2 instead of N.
                for (int j = 0; j < N; j++)
                    for (int k = 0; k <= j; k++)
                        q.push_back(vec3(((double)j + 2.0 / 3.0) / N, ((double)k + 1.0 / 3.0) / N, 0.5 / (N * N)));
                for (int j = 0; j < N; j++)
                    for (int k = 0; k < j; k++)
                        q.push_back(vec3(((double)j + 1.0 / 3.0) / N, ((double)k + 2.0 / 3.0) / N, 0.5 / (N * N)));
                break;
        }
        
        std::pair<std::map<int, std::vector<vec3> >::iterator, bool> res = s_quadratures.insert(std::pair<int, std::vector<vec3> >(N, q));
        assert(res.second);
        i = res.first;
    }
    
    return i->second;
}

// Given the face geometry (x0, x1, x2) and the square domain quadrature point k, compute the world coordinates of
//  the quadrature point x (after applying Duffy transform), compute the jacobian, and compute the baricentric
//  coordinates of x in the face. Return value is the uv coordinates of x in the triangle ref domain.
vec2 BoundaryIntegral::duffyTransform(const vec3 & x0, const vec3 & x1, const vec3 & x2, const vec2 & k, vec3 & x, double & jacobian, vec3 & c)
{
    assert(k.x() >= -1 && k.x() <= 1 && k.y() >= -1 && k.y() <= 1);
    
    vec3 n = (x1 - x0).cross(x2 - x0);
    double a = n.norm();
    n /= a;     // now n is the triangle normal (oriented by the order (x0, x1, x2), not by outward/inward with respect to the liquid domain)
    a /= 2;     // now a is the triangle area
    assert(a != 0);
    assert(n == n);
    
    vec2 s = (k + vec2(1, 1)) / 2;        // convert the quadrature coordinates in the square ref domain (-1,1)^2 to the domain (0,1)^2
    double jacobian_s = 0.25;               // the jacobian ds/dk
    assert(s.x() >= 0 && s.x() <= 1 && s.y() >= 0 && s.y() <= 1);
    
    vec2 t = vec2(s.x(), s.x() * s.y());  // coordinates in the triangle ref domain ((0,0), (1,0), (1,1)), after Duffy transform
    double jacobian_t = s.x();              // the jacobian dt/ds
    
    Eigen::Matrix<double, 3, 2> defmap;     // deformation from the triangle ref domain to the world frame: mapping (0,0) to x0, (1,0) to x1, (1,1) to x2
    defmap.col(0) = x1 - x0;
    defmap.col(1) = x2 - x1;
    
    // the world coordinates of the quadrature point
    x = defmap * t + x0;
    double jacobian_x = a * 2;              // the jacobian dx/ds
    
    // the overall jacobian
    jacobian = jacobian_s * jacobian_t * jacobian_x;    // dx/dk
    
    // the baricentric coordinates of x in the face
    c[0] = 1 - t.x();
    c[1] = t.x() - t.y();
    c[2] = t.y();
    assert(c[0] >= 0 && c[0] <= 1);
    assert(c[1] >= 0 && c[1] <= 1);
    assert(c[2] >= 0 && c[2] <= 1);
    
    // the uv coordinates of x in the triangle ref domain
    return t;
}

vec2 BoundaryIntegral::duffyTransform(const mat3 & xs, int singularity, const vec2 & k, vec3 & x, double & jacobian, vec3 & c)
{
    mat3 xs_tmp = xs;
    if (singularity != 0)   // treat the specified vertex as the singularity point (move it to the first column in xs_tmp)
    {
        vec3 tmp = xs_tmp.col(0);
        xs_tmp.col(0) = xs_tmp.col(singularity);
        xs_tmp.col(singularity) = tmp;
    }
    vec2 ret = duffyTransform(xs_tmp.col(0), xs_tmp.col(1), xs_tmp.col(2), k, x, jacobian, c);
    if (singularity != 0)
        std::swap(c[0], c[singularity]);
    
    return ret;
}

// Given the edge geometry (x0, x1) and the quadrature point k in ref domain (0,1), compute the
//  world coordinates of the quadrature point x, compute the jacobian, and compute the baricentric coordinates
//  of x in the edge. No return value.
void BoundaryIntegral::lineTransform(const vec3 & x0, const vec3 & x1, double k, vec3 & x, double & jacobian, vec2 & c)
{
    assert(k >= -1 && k <= 1);
    
    vec3 t = x1 - x0;
    double l = t.norm();
    t /= l;
    assert(l != 0);
    assert(t == t);
    
    // the overall jacobian
    jacobian = l / 2;                  // dx/dk
    
    // the baricentric coordinates of x in the face
    c[0] = (1 - k) / 2;
    c[1] = (1 + k) / 2;
    assert(c[0] >= 0 && c[0] <= 1);
    assert(c[1] >= 0 && c[1] <= 1);
    
    // the world coordinates of the quadrature point
    x = x0 * c[0] + x1 * c[1];
}

void BoundaryIntegral::lineTransform(const Eigen::Matrix<double, 3, 2> & xs, double k, vec3 & x, double & jacobian, vec2 & c)
{
    lineTransform(xs.col(0), xs.col(1), k, x, jacobian, c);
}
