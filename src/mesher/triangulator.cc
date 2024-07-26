#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <predicates.h>
#include "mesh.hh"
#include "topology.hh"
#include "geometry.hh"
#include "spatial_search.hh"
#include "mesh_io.hh"

using namespace OpenMesh;

////////////////////////////////////////////////////////////////
/// pred2D
////////////////////////////////////////////////////////////////

enum { CCW = 1, CW = -1, LINEAR = 0 };

inline int sign(const double r)
{
    return ((int)(r > 0) - (int)(r < 0)); // <0: -1, >0: +1, =0: 0
}

inline double determinant(const Vec2 &a, const Vec2 &b, const Vec2 &c)
{
    return orient2d(a.data(), b.data(), c.data());
}

inline int orientation(const Vec2 &a, const Vec2 &b, const Vec2 &c)
{
    return sign(orient2d(a.data(), b.data(), c.data())); // -1:CW, +1:CCW, 0:LINEAR
}

inline double determinant(const Vec2 &d0, const Vec2 &d1)
{
    return determinant({ 0,0 }, d0, d1);
}

inline Int4 intersection_info(const Vec2 &u0, const Vec2 &u1, const Vec2 &v0, const Vec2 &v1)
{
    const int ru0 = orientation(v0, v1, u0);
    const int ru1 = orientation(v0, v1, u1);
    const int rv0 = orientation(u0, u1, v0);
    const int rv1 = orientation(u0, u1, v1);
    return { ru0, ru1, rv0, rv1 };
}

inline Vec2 intersection_param(const Vec2 &u0, const Vec2 &u1, const Vec2 &v0, const Vec2 &v1)
{
    const double dt = determinant(u1 - u0, v1 - v0);
    const double t0 = determinant(v0 - u0, v1 - v0);
    const double t1 = determinant(v0 - u0, u1 - u0);
    return Vec2 { t0, t1 } / dt;
}

//inline bool is_intersecting(const Int4 &info)
//{
//    const int ru0 = info[0];
//    const int ru1 = info[1];
//    const int rv0 = info[2];
//    const int rv1 = info[3];
//
//    return
//        ru0 * ru1 < 0 && rv0 * rv1 < 0       || // intersecting
//        ru0 == 0 && ru1 != 0 && rv0*rv1 <= 0 || // u0 lies on [v0,v1]
//        ru1 == 0 && ru0 != 0 && rv0*rv1 <= 0 || // u1 lies on [v0,v1]
//        rv0 == 0 && rv1 != 0 && ru0*ru1 <= 0 || // v0 lies on [u0,u1]
//        rv1 == 0 && rv0 != 0 && ru0*ru1 <= 0 ;  // v1 lies on [u0,u1]
//}

//inline bool is_intersecting(const Vec2 &u0, const Vec2 &u1, const Vec2 &v0, const Vec2 &v1)
//{
//    return is_intersecting(intersection_info(u0, u1, v0, v1));
//}

//inline int projective_region(const Vec2 &u0, const Vec2 &u1, const Vec2 &u)
//{
//    return sign(dot(u1 - u0, u - u0));
//}

//inline Int8 overlapping_info(const Vec2 &u0, const Vec2 &u1, const Vec2 &v0, const Vec2 &v1)
//{
//    const auto du = u1 - u0;
//    const auto dv = v1 - v0;
//    const int pu0v0 = projective_region(v0, v1, u0);
//    const int pu1v0 = projective_region(v0, v1, u1);
//    const int pu0v1 = projective_region(v1, v0, u0);
//    const int pu1v1 = projective_region(v1, v0, u1);
//    const int pv0u0 = projective_region(u0, u1, v0);
//    const int pv1u0 = projective_region(u0, u1, v1);
//    const int pv0u1 = projective_region(u1, u0, v0);
//    const int pv1u1 = projective_region(u1, u0, v1);
//    return { pu0v0, pu1v0, pu0v1, pu1v1, pv0u0, pv1u0, pv0u1, pv1u1 };
//}

//inline bool is_overlapping(const Int8 &info)
//{
//    const int pu0v0 = info[0];
//    const int pu1v0 = info[1];
//    const int pu0v1 = info[2];
//    const int pu1v1 = info[3];
//    const int pv0u0 = info[4];
//    const int pv1u0 = info[5];
//    const int pv0u1 = info[6];
//    const int pv1u1 = info[7];
//    return
//        pu0v0 >= 0 && pu0v1 >= 0 || // u0 lies in [v0,v1]
//        pu1v0 >= 0 && pu1v1 >= 0 || // u1 lies in [v0,v1]
//        pv0u0 >= 0 && pv0u1 >= 0 || // v0 lies in [u0,u1]
//        pv1u0 >= 0 && pv1u1 >= 0;   // v1 lies in [u0,u1]
//}

//inline bool is_overlapping(const Vec2 &u0, const Vec2 &u1, const Vec2 &v0, const Vec2 &v1)
//{
//    return is_overlapping(overlapping_info(u0, u1, v0, v1));
//}

//inline bool is_parallel(const Vec2 &d0, const Vec2 &d1)
//{
//    return orientation({ 0,0 }, d0, d1) == 0;
//}

////////////////////////////////////////////////////////////////
/// Locations
////////////////////////////////////////////////////////////////

//
//     V2
//    /  \
//  E1    E0
//  /      \
// V0--E2--V1
//
enum
{
    OUTSIDE = 0,
    IN_TR,
    ON_E0,
    ON_E1,
    ON_E2,
    AT_V0,
    AT_V1,
    AT_V2
};

inline int locate(const Vec2 &u0, const Vec2 &u1, const Vec2 &u2, const Vec2 &u)
{
    const int r0 = orientation(u1, u2, u);
    const int r1 = orientation(u2, u0, u);
    const int r2 = orientation(u0, u1, u);

    return
        r0 == r1 && r0 == r2 ? IN_TR : // in triangle
        r0 == 0  && r1 == r2 ? ON_E0 : // on edge (u1, u2)
        r1 == 0  && r2 == r0 ? ON_E1 : // on edge (u2, u0)
        r2 == 0  && r0 == r1 ? ON_E2 : // on edge (u0, u1)
        r0 == 0  && r1 == 0  ? AT_V2 : // on vert u2
        r1 == 0  && r2 == 0  ? AT_V0 : // on vert u0
        r2 == 0  && r0 == 0  ? AT_V1 : // on vert u1
        OUTSIDE;
}

////////////////////////////////////////////////////////////////
/// utils
////////////////////////////////////////////////////////////////

inline double argument_angle(const Vec2 &du, const Vec2 &dw)
{ return atan2(cross(du, dw), dot(du, dw)); }

template <class MeshT>
inline Vec2 get_xy(const MeshT &mesh, Vh vh)
{ const auto p = mesh.point(vh); return { p[0], p[1] }; }

template <class MeshT>
inline Vec2 get_xy(const MeshT &mesh, Hh hh)
{ return get_xy(mesh, mesh.to_vertex_handle(hh)); }

template <class MeshT>
inline Vec2 get_dxy(const MeshT &mesh, Hh hh)
{ return get_xy(mesh, hh) - get_xy(mesh, mesh.opposite_halfedge_handle(hh)); }

template <class MeshT>
inline void set_xy(MeshT &mesh, const Vh &vh, const Vec2 &u)
{ const auto p = mesh.point(vh); mesh.set_point(vh, { u[0], u[1], p[2] }); }

template <class MeshT>
inline VecN<Vec2, 2> get_range(const MeshT &mesh)
{
    Vec2 bl { +1e20, +1e20 };
    Vec2 ur { -1e20, -1e20 };

    for (Vh vh : mesh.vertices())
    {
        const auto u = get_xy(mesh, vh);
        bl[0] = (u[0] < bl[0]) ? u[0] : bl[0];
        bl[1] = (u[1] < bl[1]) ? u[1] : bl[1];
        ur[0] = (u[0] > ur[0]) ? u[0] : ur[0];
        ur[1] = (u[1] > ur[1]) ? u[1] : ur[1];
    }

    return { bl, ur };
}

////////////////////////////////////////////////////////////////
/// Delaunay
////////////////////////////////////////////////////////////////

//   2
//  / \
// 0---1
//  \ /
//   3
inline bool is_delaunay(const TriMesh &mesh, const Eh &eh)
{
    Hh hh0 = mesh.halfedge_handle(eh, 0);
    Hh hh1 = mesh.halfedge_handle(eh, 1);
    const auto u0 = get_xy(mesh, hh1);
    const auto u1 = get_xy(mesh, hh0);
    const auto u2 = get_xy(mesh, mesh.next_halfedge_handle(hh0));
    const auto u3 = get_xy(mesh, mesh.next_halfedge_handle(hh1));
    return is_delaunay(u0, u1, u2, u3);
}

struct EuclideanDelaunay
{
    inline bool operator()(const TriMesh &mesh, const Eh &eh) const
    { return is_delaunay(mesh, eh); } // If true, do not flip
};

////////////////////////////////////////////////////////////////
/// Insertion
////////////////////////////////////////////////////////////////

inline int locate(TriMesh &mesh, const Fh &fh, const Vec2 &u, Hh &hh)
{
    constexpr double kEps = 1e-3;
    const auto hh0 = mesh.halfedge_handle(fh);
    const auto hh1 = mesh.next_halfedge_handle(hh0);
    const auto hh2 = mesh.next_halfedge_handle(hh1);
    const auto u0 = get_xy(mesh, mesh.to_vertex_handle(hh1));
    const auto u1 = get_xy(mesh, mesh.to_vertex_handle(hh2));
    const auto u2 = get_xy(mesh, mesh.to_vertex_handle(hh0));
    const int loc = locate(u0, u1, u2, u);
    if (loc == ON_E0) { hh = hh0; }
    if (loc == ON_E1) { hh = hh1; }
    if (loc == ON_E2) { hh = hh2; }
    return loc;
    //const auto bc = barycentric_coordinate(u0, u1, u2, u);
    //if (fabs(bc[1]) < kEps && fabs(bc[2]) < kEps) return AT_V0;
    //if (fabs(bc[2]) < kEps && fabs(bc[0]) < kEps) return AT_V1;
    //if (fabs(bc[0]) < kEps && fabs(bc[1]) < kEps) return AT_V2;
    //if (fabs(bc[0]) < kEps) { hh = hh0; return ON_E0; }
    //if (fabs(bc[1]) < kEps) { hh = hh1; return ON_E1; }
    //if (fabs(bc[2]) < kEps) { hh = hh2; return ON_E2; }
    //return IN_TR;
}

inline void split_edge(TriMesh &mesh, const Eh &eh, const Vh &vh)
{
    auto vh0 = mesh.from_vertex_handle(mesh.halfedge_handle(eh, 0));
    auto vh1 = mesh.to_vertex_handle  (mesh.halfedge_handle(eh, 0));

    mesh.split_edge_copy(eh, vh);

    for (auto hdge : mesh.voh_range(vh))
    if (hdge.to() != vh0)
    if (hdge.to() != vh1)
    { set_sharp(mesh, hdge.edge(), false); }
}

static Fh search_triangle(const TriMesh &mesh, const Vec2 &u)
{
    auto fh = search_triangle_guided_bfs(mesh, u, mesh.face_handle(0));
    if (!fh.is_valid()) fh = search_triangle_brute_force(mesh, u);
    return fh;
}

static void add_points(TriMesh &mesh, const LoopMesh &poly)
{
    assert(mesh.n_vertices() == 0);

    for (Vh vp : poly.vertices())
    {
        const auto u = get_xy(poly, vp);
        mesh.new_vertex({ u[0], u[1], 0 });
    }
}

static void set_domain(TriMesh &mesh, const LoopMesh &poly)
{
    // Get the size of the polygon
    const auto blur = get_range(poly);
    auto bl = blur[0], ur = blur[1];
    const auto dl = ur - bl;
    const double l = (dl[0]<dl[1]) ? dl[0] : dl[1];
    const Vec2 d { l, l };
    bl -= d;
    ur += d;

    Vh vh0 = mesh.new_vertex({ 0,0,0 });
    Vh vh1 = mesh.new_vertex({ 0,0,0 });
    Vh vh2 = mesh.new_vertex({ 0,0,0 });
    Vh vh3 = mesh.new_vertex({ 0,0,0 });
    set_xy(mesh, vh0, { bl[0], bl[1] });
    set_xy(mesh, vh1, { ur[0], bl[1] });
    set_xy(mesh, vh2, { ur[0], ur[1] });
    set_xy(mesh, vh3, { bl[0], ur[1] });
    mesh.add_face({ vh0, vh1, vh2 });
    mesh.add_face({ vh2, vh3, vh0 });
}

static int insert_vertices(TriMesh &mesh)
{
    const int max_num_edge_flip = 100;

    auto delaunifier = make_delaunifier(mesh, EuclideanDelaunay {});

    int nv {};

    for (Vh vh : mesh.vertices()) if (mesh.is_isolated(vh))
    {
        const auto u = get_xy(mesh, vh);

        // find the nearest triangle of new point
        auto fh = search_triangle(mesh, u);

        // skip any point out of domain
        if (!fh.is_valid()) continue;

        // on which primitive of the triangle the point is
        Hh hh {}; int loc = locate(mesh, fh, u, hh);
        if (loc == OUTSIDE || loc == AT_V0 || loc == AT_V1 || loc == AT_V2) continue;

        // insert the point in the triangle or on the edge
        if (loc == IN_TR) mesh.split(fh, vh);
        else split_edge(mesh, mesh.edge_handle(hh), vh);

        // edges to flip
        Eh ehs[4]; int ne {};

        for (auto hdge : mesh.voh_range(vh))
            if (!hdge.next().edge().is_boundary())
                ehs[ne++] = hdge.next().edge();

        // make delaunay
        delaunifier.reset(); delaunifier.to_flip(ehs, ne); int n_flip {};

        for (auto eh = delaunifier.flip(); eh.is_valid() && n_flip < max_num_edge_flip; eh = delaunifier.flip(), ++n_flip) {}

        ++nv;
    }

    return nv;
}

////////////////////////////////////////////////////////////////
/// Local search
////////////////////////////////////////////////////////////////

enum { NO_ITSC, IN_EDGE, ON_VERTEX };

inline int intersection_info(const TriMesh &mesh, const Vec2 &u0, const Vec2 &u1, Hh hh)
{
    const auto v0 = get_xy(mesh, mesh.from_vertex_handle(hh));
    const auto v1 = get_xy(mesh, mesh.to_vertex_handle  (hh));

    const auto ii = intersection_info(u0, u1, v0, v1);
    const int ru0 = ii[0];
    const int ru1 = ii[1];
    const int rv0 = ii[2];
    const int rv1 = ii[3];

    return
        (ru0 * ru1 < 0) && (rv0 * rv1 < 0)        ? IN_EDGE   : // intersecting
        (ru0 * ru1 < 0) && (rv1 == 0 && rv0 != 0) ? ON_VERTEX : // v1 lies on (u0,u1)
        (ru0 != 0 && ru1 == 0) && (rv1 == 0)      ? ON_VERTEX : // v1 overlaps u1
        NO_ITSC; // no intersecting, overlapping, v0 lies on (u0,u1), and other cases
}

inline double intersection_param(const TriMesh &mesh, const Vec2 &u0, const Vec2 &u1, Hh hh)
{
    const auto v0 = get_xy(mesh, mesh.from_vertex_handle(hh));
    const auto v1 = get_xy(mesh, mesh.to_vertex_handle  (hh));
    return intersection_param(u0, u1, v0, v1)[0]; // t in eq: u0 + (u1-u0)*t
}

inline int next_intersection(const TriMesh &mesh, const Vec2 &u0, const Vec2 &u1, Hh hho, Hh &hhc, Vh &vhc)
{
    int res { NO_ITSC };

    Hh hh0 = mesh.prev_halfedge_handle(hho);
    Hh hh1 = mesh.next_halfedge_handle(hho);

    for (Hh hh : { hh0, hh1 })
    {
        const int ii = intersection_info(mesh, u0, u1, hh);
        if (ii == NO_ITSC) continue;

        res = ii;

        if (ii == IN_EDGE)
        {
            hhc = hh;
        }
        else if (ii == ON_VERTEX)
        {
            vhc = mesh.to_vertex_handle(hh);
        }
    }

    return res;
}

inline int next_intersection(const TriMesh &mesh, const Vec2 &u0, const Vec2 &u1, Vh vho, Hh &hhc, Vh &vhc)
{
    int res { NO_ITSC };
    double tmax {}; // the parameter in line equation: u0 + (u1-u0)*t

    for (Hh hh : mesh.voh_range(vho)) if (!mesh.is_boundary(hh))
    {
        hh = mesh.next_halfedge_handle(hh); // apex edge of v0

        const int ii = intersection_info(mesh, u0, u1, hh);
        if (ii == NO_ITSC) continue;

        const double t = intersection_param(mesh, u0, u1, hh);
        if (tmax >= t) continue; // intersected but earlier

        tmax = t;
        res = ii;

        if (ii == IN_EDGE)
        {
            hhc = hh;
        }
        else if (ii == ON_VERTEX)
        {
            vhc = mesh.to_vertex_handle(hh);
        }
    }

    return res;
}

static int get_intersections(TriMesh &mesh, Vh vh0, Vh vh1, std::vector<Hh> &hhs, std::vector<Vh> &vhs)
{
    const auto u0 = get_xy(mesh, vh0);
    const auto u1 = get_xy(mesh, vh1);

    // search starts with v0
    int status = ON_VERTEX;
    Vh vhc = vh0;
    Hh hhc;

    const int max_num_iter = (int)mesh.n_edges();
    int iter {};

    for (int iter = 0; iter < max_num_iter; ++iter)
    {
        if (status == IN_EDGE)
        {
            status = next_intersection(mesh, u0, u1, hhc, hhc, vhc);
        }
        else if (status == ON_VERTEX)
        {
            status = next_intersection(mesh, u0, u1, vhc, hhc, vhc);
        }
        else // search in vain, for some reasons
        {
            printf("searching segment (%d, %d) lost on the way\n", vh0.idx()+1, vh1.idx()+1);
            break;
        }

        if (status == ON_VERTEX)
        {
            if (vhc == vh1) // search reaches v1 successfully
            {
                break;
            }
            else // (u0,u1) cannot go thru any vertex other than v0 and v1
            {
                vhs.push_back(vhc);
            }
        }
        else if (status == IN_EDGE)
        {
            hhs.push_back(hhc);
            hhc = mesh.opposite_halfedge_handle(hhc);
        }
    }

    if (iter >= max_num_iter)
    {
        printf("searching segment (%d, %d) reached maximum iteration\n", vh0.idx()+1, vh1.idx()+1);
        return NO_ITSC;
    }

    return status;
}

////////////////////////////////////////////////////////////////
/// Constraints
////////////////////////////////////////////////////////////////

//
//   1   
//  / \  
// 2---0 
//  \ /  
//   3   
//
inline bool is_convex(const Vec2 &u0, const Vec2 &u1, const Vec2 &u2, const Vec2 &u3)
{
    const int r0 = orientation(u3, u0, u1);
    const int r1 = orientation(u0, u1, u2);
    const int r2 = orientation(u1, u2, u3);
    const int r3 = orientation(u2, u3, u0);
    return r0 > 0 && r1 > 0 && r2 > 0 && r3 > 0 ||
           r0 < 0 && r1 < 0 && r2 < 0 && r3 < 0 ;
}

//
//   1   
//  / \  
// 2---0 
//  \ /  
//   3   
//
inline bool is_flippable(const TriMesh &mesh, const Hh &hh)
{
    Hh hi = mesh.opposite_halfedge_handle(hh);
    const auto u0 = get_xy(mesh, hh);
    const auto u1 = get_xy(mesh, mesh.next_halfedge_handle(hh));
    const auto u2 = get_xy(mesh, hi);
    const auto u3 = get_xy(mesh, mesh.next_halfedge_handle(hi));
    return is_convex(u0, u1, u2, u3);
}

inline bool is_intersecting(const TriMesh &mesh, Vh vh0, Vh vh1, Hh hh)
{
    const auto u0 = get_xy(mesh, vh0);
    const auto u1 = get_xy(mesh, vh1);
    return intersection_info(mesh, u0, u1, hh) == IN_EDGE;
}

static Hh restore_constraint(TriMesh &mesh, Vh vh0, Vh vh1)
{
    Hh hh_rc = mesh.find_halfedge(vh0, vh1);
    if (hh_rc.is_valid()) return hh_rc;

    std::vector<Hh> hhs {}; // edges intersected with (u0,u1)
    std::vector<Vh> vhs {}; // vertices lying on (u0,u1)

    if (get_intersections(mesh, vh0, vh1, hhs, vhs) == NO_ITSC)
        return Hh {};

    // check self-intersection
    bool has_self_intersection {};

    for (Vh vh : vhs)
    {
        printf("segment (%d, %d) goes thru vertex (%d)\n", vh0.idx()+1, vh1.idx()+1, vh.idx()+1);
        // [TODO] record overlapping vertices
        has_self_intersection = true;
    }

    for (Hh hh : hhs) if (is_sharp(mesh, mesh.edge_handle(hh)))
    {
        printf("segment (%d, %d) intersects with another segment ", vh0.idx()+1, vh1.idx()+1); print_handle(mesh, hh, 1); printf("\n");
        // [TODO] record intersection pairs
        has_self_intersection = true;
    }

    if (has_self_intersection) return Hh {};

    //for (Hh hh : hhs) { printf("segment (%d, %d) intersects with edge ", vh0.idx()+1, vh1.idx()+1); print_handle(mesh, hh, 1); printf("\n"); }

    // keep flipping until the edge is recovered
    std::deque<Hh> frontier(hhs.begin(), hhs.end());

    const int max_num_iter = (int)mesh.n_edges();

    for (int iter = 0; !frontier.empty() && iter < max_num_iter; ++iter)
    {
        Hh hh = frontier.front(); frontier.pop_front();

        //printf("get edge "); print_handle(mesh, hh, 1); printf("\n");

        if (!is_flippable(mesh, hh))
        {
            frontier.push_back(hh);
            continue;
        }

        //printf("flip edge "); print_handle(mesh, hh, 1);

        mesh.flip(mesh.edge_handle(hh));

        //printf(" into "); print_handle(mesh, hh, 1); printf("\n");

        if (is_intersecting(mesh, vh0, vh1, hh))
        {
            frontier.push_back(hh);
        }
    }

    hh_rc = mesh.find_halfedge(vh0, vh1);

    return hh_rc;
}

static int make_delaunay(TriMesh &mesh, Eh eh)
{
    auto delaunifier = make_delaunifier(mesh, EuclideanDelaunay {});

    const int max_num_iter = (int)mesh.n_edges();

    const Eh ehs[4] {
        mesh.edge_handle(mesh.next_halfedge_handle(mesh.halfedge_handle(eh, 0))),
        mesh.edge_handle(mesh.prev_halfedge_handle(mesh.halfedge_handle(eh, 0))),
        mesh.edge_handle(mesh.next_halfedge_handle(mesh.halfedge_handle(eh, 1))),
        mesh.edge_handle(mesh.prev_halfedge_handle(mesh.halfedge_handle(eh, 1))) };

    delaunifier.reset(); delaunifier.to_flip(ehs, 4); int n_flip {};

    for (auto eh = delaunifier.flip(); eh.is_valid() && n_flip < max_num_iter; eh = delaunifier.flip(), ++n_flip) {}

    return n_flip;
}

static int restore_constraints(TriMesh &mesh, const LoopMesh &poly)
{
    for (Eh ep : poly.edges())
    {
        Vh vp0 = poly.from_vertex_handle(poly.halfedge_handle(ep, 0));
        Vh vp1 = poly.to_vertex_handle  (poly.halfedge_handle(ep, 0));
        Vh vh0 = mesh.vertex_handle(vp0.idx());
        Vh vh1 = mesh.vertex_handle(vp1.idx());
        Hh hh = restore_constraint(mesh, vh0, vh1);

        if (!hh.is_valid())
        {
            printf("restoring segment (%d, %d) failed\n", vh0.idx()+1, vh1.idx()+1);
            continue;
        }

        set_sharp(mesh, mesh.edge_handle(hh), true); // mark restored edge

        make_delaunay(mesh, mesh.edge_handle(hh)); // maintain Delaunay after restoration
    }

    return 0;
}

static int remove_exteriors(TriMesh &mesh)
{
    std::unordered_set<Fh> exteriors {};
    std::queue<Fh> frontier {};

    Vh vh_last = *(--mesh.vertices_end());
    Fh fh_last = *(mesh.vf_begin(vh_last));

    frontier.push(fh_last);
    exteriors.insert(fh_last);

    while (!frontier.empty())
    {
        Fh fh = frontier.front(); frontier.pop();

        for (Hh hh : mesh.fh_range(fh))
        {
            Eh eh = mesh.edge_handle(hh);
            if (is_sharp(mesh, eh)) continue;
            if (mesh.is_boundary(eh)) continue;

            Fh fi = mesh.opposite_face_handle(hh);
            if (exteriors.count(fi)) continue;

            frontier.push(fi);
            exteriors.insert(fi);
        }
    }

    if (exteriors.size() == mesh.n_faces())
        return 1;

    for (Fh fh : exteriors)
    {
        mesh.delete_face(fh);
    }

    mesh.garbage_collection();

    return 0;
}

////////////////////////////////////////////////////////////////
/// Wrapping up
////////////////////////////////////////////////////////////////

int triangulate(const LoopMesh &poly, TriMesh &mesh)
{
    // Copy all vertices from polygon. Better do it at
    // beginning so the vertex is ordered accordingly.
    add_points(mesh, poly);

    // Generate extended domain
    set_domain(mesh, poly);

    // Insert points into the domain
    insert_vertices(mesh);

    // Restore constraint edges in the domain
    restore_constraints(mesh, poly);

    // Remove exterior region
    remove_exteriors(mesh);

    return 0;
}
