#include "mesh.hh"
#include "topology.hh"
#include "geometry.hh"

using namespace OpenMesh;

template <class MeshT>
inline Vec2 get_xy(const MeshT &mesh, Vh vh)
{ const auto p = mesh.point(vh); return { p[0], p[1] }; }

template <class MeshT>
inline Vec2 get_xy(const MeshT &mesh, Hh hh)
{ return get_xy(mesh, mesh.to_vertex_handle(hh)); }

//
//   2
//  / \
// 0---1
//  \ /
//   3
//
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

int make_delaunay(TriMesh &mesh)
{
    auto delaunifier = make_delaunifier(mesh, EuclideanDelaunay {});

    const int max_num_edge_flip = (int)mesh.n_faces() * 50;

    delaunifier.reset(); delaunifier.to_flip(); int n_flip {};

    for (auto hdge = delaunifier.flip(); hdge.is_valid() && n_flip < max_num_edge_flip; hdge = delaunifier.flip(), ++n_flip) {}

    //printf("Total flips: %d\n", n_flip);
    return n_flip;
}
