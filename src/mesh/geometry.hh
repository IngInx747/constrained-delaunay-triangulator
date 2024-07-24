#ifndef GEOMETRY_HH
#define GEOMETRY_HH

#include <cmath>

////////////////////////////////////////////////////////////////
/// Common
////////////////////////////////////////////////////////////////

template <typename T> constexpr T pi()
{
    return (T)(3.14159265358979323846264338327950288);
}

template <typename T> inline int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

template <typename T> inline T round(T val, T mtp)
{
    return mtp * std::round(val / mtp);
}

template <typename T> inline T radian(T deg)
{
    return (T)(deg * pi<T>() / (T)(180.));
}

template <typename T> inline T degree(T rad)
{
    return (T)(rad * (T)(180.) / pi<T>());
}

inline double cosine(double a, double b, double c)
{
    const double cs = (a*a + b*b - c*c) / (a*b*2);
    return cs < -1 ? -1 : cs > 1 ? 1 : cs;
}

inline double cos2cot(double cs)
{
    return cs / std::sqrt(1.0 - cs * cs);
}

////////////////////////////////////////////////////////////////
/// Vector
////////////////////////////////////////////////////////////////

#include "vector_n.hh"

inline double cross(const Vec2 &a, const Vec2 &b)
{
    return a[0]*b[1] - a[1]*b[0];
}

inline double dihedral_angle(const Vec3 &e, const Vec3 &n0, const Vec3 &n1)
{
    return atan2(dot(e.normalized(), cross(n0, n1)), dot(n0, n1)); // (-pi, pi]
}

inline Vec3 barycentric_coordinate(const Vec2 &a, const Vec2 &b, const Vec2 &c, const Vec2 &p)
{
    const auto ab = b - a;
    const auto bc = c - b;
    const auto ca = a - c;
    const auto ap = p - a;
    const auto bp = p - b;
    const auto cp = p - c;
    const double ma = cross(bc, cp);
    const double mb = cross(ca, ap);
    const double mc = cross(ab, bp);
    const Vec3 phi { ma, mb, mc };
    return phi / (ma + mb + mc);
}

//
//   2
//  / \
// 0---1
//  \ /
//   3
//
inline bool is_delaunay(const Vec2 &u0, const Vec2 &u1, const Vec2 &u2, const Vec2 &u3)
{
    const double l  = norm(u1 - u0);
    const double l0 = norm(u2 - u0);
    const double l1 = norm(u2 - u1);
    const double l2 = norm(u3 - u0);
    const double l3 = norm(u3 - u1);
    const double cs0 = cosine(l0, l1, l);
    const double cs1 = cosine(l2, l3, l);
    return cs0 + cs1 >= 0;
}

inline Vec2 circumcenter(const Vec2 &u0, const Vec2 &u1, const Vec2 &u2)
{
    const auto d0 = u1 - u0;
    const auto d1 = u2 - u0;
    const double a = cross(d0, d1);
    const double dd0 = dot(d0, d0);
    const double dd1 = dot(d1, d1);
    const Vec2 uc { d1[1]*dd0 - d0[1]*dd1, d0[0]*dd1 - d1[0]*dd0 };
    return uc / (a*2.) + u0;
}

inline double hausdorff_distance(const VecN<Vec2, 2> &s, const Vec2 &p)
{
    const auto &a = s[0];
    const auto &b = s[1];
    const auto d0 = p - a;
    const auto d1 = p - b;
    const auto e = (b - a).normalized();
    return dot(d0, e) < 0 ? norm(d0)
        :  dot(d1,-e) < 0 ? norm(d1)
        :  fabs(cross(d0, e));
}

inline double hausdorff_distance(const VecN<Vec3, 2> &s, const Vec3 &p)
{
    const auto &a = s[0];
    const auto &b = s[1];
    const auto d0 = p - a;
    const auto d1 = p - b;
    const auto e = (b - a).normalized();
    return dot(d0, e) < 0 ? norm(d0)
        :  dot(d1,-e) < 0 ? norm(d1)
        :  norm(cross(d0, e));
}

inline double hausdorff_distance(const VecN<Vec2, 3> &t, const Vec2 &p)
{
    const auto &a = t[0];
    const auto &b = t[1];
    const auto &c = t[2];
    const auto ab = b - a;
    const auto bc = c - b;
    const auto ca = a - c;
    const auto ap = p - a;
    const auto bp = p - b;
    const auto cp = p - c;
    const auto na = normalize(Vec2 { bc[1], -bc[0] });
    const auto nb = normalize(Vec2 { ca[1], -ca[0] });
    const auto nc = normalize(Vec2 { ab[1], -ab[0] });
    if (dot(ab, ap) <= 0 && dot(ca, ap) >= 0) return norm(ap);
    if (dot(bc, bp) <= 0 && dot(ab, bp) >= 0) return norm(bp);
    if (dot(ca, cp) <= 0 && dot(bc, cp) >= 0) return norm(cp);
    if (dot(nc, ap) >= 0 && dot(ab, ap) >= 0 && dot(ab, bp) <= 0) return dot(nc, ap);
    if (dot(na, bp) >= 0 && dot(bc, bp) >= 0 && dot(bc, cp) <= 0) return dot(na, bp);
    if (dot(nb, cp) >= 0 && dot(ca, cp) >= 0 && dot(ca, ap) <= 0) return dot(nb, cp);
    return 0; // inside triangle
}

inline double hausdorff_distance(const VecN<Vec3, 3> &t, const Vec3 &p)
{
    const auto &a = t[0];
    const auto &b = t[1];
    const auto &c = t[2];
    const auto ab = b - a;
    const auto bc = c - b;
    const auto ca = a - c;
    const auto ap = p - a;
    const auto bp = p - b;
    const auto cp = p - c;
    const auto n =  normalize(cross(ab,bc));
    const auto na = normalize(cross(bc, n));
    const auto nb = normalize(cross(ca, n));
    const auto nc = normalize(cross(ab, n));
    if (dot(ab, ap) <= 0 && dot(ca, ap) >= 0) return norm(ap);
    if (dot(bc, bp) <= 0 && dot(ab, bp) >= 0) return norm(bp);
    if (dot(ca, cp) <= 0 && dot(bc, cp) >= 0) return norm(cp);
    if (dot(nc, ap) >= 0 && dot(ab, ap) >= 0 && dot(ab, bp) <= 0) return norm(cross(ab, ap)) / norm(ab);
    if (dot(na, bp) >= 0 && dot(bc, bp) >= 0 && dot(bc, cp) <= 0) return norm(cross(bc, bp)) / norm(bc);
    if (dot(nb, cp) >= 0 && dot(ca, cp) >= 0 && dot(ca, ap) <= 0) return norm(cross(ca, cp)) / norm(ca);
    return fabs(dot(ap, n)); // inside triangle
}

#endif