#ifndef OPENMESH_LOOPMESH_HH
#define OPENMESH_LOOPMESH_HH

#include <OpenMesh/Core/Geometry/MathDefs.hh>
#include <OpenMesh/Core/Mesh/Traits.hh>
#include <OpenMesh/Core/Mesh/AttribKernelT.hh>
#include <OpenMesh/Core/Mesh/FinalMeshItemsT.hh>
#include "loop_connectivity.hh"

namespace OpenMesh {

template <class Kernel>
class LoopMeshT : public Kernel
{
public:

    typedef LoopMeshT<Kernel>  This;

    typedef typename Kernel::Scalar              Scalar;
    typedef typename Kernel::Point               Point;
    typedef typename Kernel::Normal              Normal;
    typedef typename Kernel::Color               Color;
    typedef typename Kernel::TexCoord1D          TexCoord1D;
    typedef typename Kernel::TexCoord2D          TexCoord2D;
    typedef typename Kernel::TexCoord3D          TexCoord3D;
    typedef typename Kernel::Vertex              Vertex;
    typedef typename Kernel::Halfedge            Halfedge;
    typedef typename Kernel::Edge                Edge;
    typedef typename Kernel::Face                Face; // LoopMesh exclusive

    //--- handles ---

    typedef typename Kernel::VertexHandle        VertexHandle;
    typedef typename Kernel::HalfedgeHandle      HalfedgeHandle;
    typedef typename Kernel::EdgeHandle          EdgeHandle;
    typedef typename Kernel::FaceHandle          FaceHandle; // LoopMesh exclusive

    //--- iterators ---

    typedef typename Kernel::VertexIter          VertexIter;
    typedef typename Kernel::ConstVertexIter     ConstVertexIter;
    typedef typename Kernel::EdgeIter            EdgeIter;
    typedef typename Kernel::ConstEdgeIter       ConstEdgeIter;
    typedef typename Kernel::FaceIter            FaceIter;      // LoopMesh exclusive
    typedef typename Kernel::ConstFaceIter       ConstFaceIter; // LoopMesh exclusive

    //--- circulators ---

    typedef typename Kernel::VertexVertexIter         VertexVertexIter;
    typedef typename Kernel::VertexOHalfedgeIter      VertexOHalfedgeIter;
    typedef typename Kernel::VertexIHalfedgeIter      VertexIHalfedgeIter;
    typedef typename Kernel::VertexEdgeIter           VertexEdgeIter;
    typedef typename Kernel::ConstVertexVertexIter    ConstVertexVertexIter;
    typedef typename Kernel::ConstVertexOHalfedgeIter ConstVertexOHalfedgeIter;
    typedef typename Kernel::ConstVertexIHalfedgeIter ConstVertexIHalfedgeIter;
    typedef typename Kernel::ConstVertexEdgeIter      ConstVertexEdgeIter;

public:

    inline VertexHandle add_vertex(const Point &p)
    {
        VertexHandle vh = Kernel::new_vertex();
        set_point(vh, p);
        return vh;
    }

    inline VertexHandle add_vertex_dirty(const Point &p)
    {
        VertexHandle vh = Kernel::new_vertex_dirty();
        set_point(vh, p);
        return vh;
    }

    inline VertexHandle split(EdgeHandle eh, const Point &p, bool copy_properties = true)
    {
        VertexHandle vh = add_vertex(p);
        Kernel::split(eh, vh, copy_properties);
        return vh;
    }

    /// Returns the point of _vh
    inline Point calc_centroid(VertexHandle vh) const
    { return point(vh); }

    /// Computes and returns the average of the vertices defining _hh (same as calc_edge_midpoint for edge of halfedge)
    inline Point calc_centroid(HalfedgeHandle hh) const
    { return (point(to_vertex_handle(hh)) + point(from_vertex_handle(hh))) * 0.5; }

    /// Computes and returns the average of the vertices defining _eh (same as calc_edge_midpoint)
    inline Point calc_centroid(EdgeHandle eh) const
    { return calc_centroid(halfedge_handle(eh, 0)); }

    /// Calculates the edge vector as the difference of the
    /// the points defined by to_vertex_handle() and from_vertex_handle()
    inline Normal calc_edge_vector(HalfedgeHandle hh) const
    { return point(to_vertex_handle(hh)) - point(from_vertex_handle(hh)); }

    /// Calculates the edge vector as the vector defined by
    /// the halfedge with id #0 (see below)
    inline Normal calc_edge_vector(EdgeHandle eh) const
    { return calc_edge_vector(halfedge_handle(eh, 0)); }

    inline Scalar calc_edge_sqr_length(HalfedgeHandle hh) const
    { return sqrnorm(calc_edge_vector(hh)); }

    inline Scalar calc_edge_sqr_length(EdgeHandle eh) const
    { return sqrnorm(calc_edge_vector(halfedge_handle(eh, 0))); }

    /// Calculates the length of the edge eh
    inline Scalar calc_edge_length(EdgeHandle _eh) const
    { return calc_edge_length(this->halfedge_handle(_eh,0)); }

    /// Calculates the length of the edge hh
    inline Scalar calc_edge_length(HalfedgeHandle hh) const
    { return (Scalar)sqrt(calc_edge_sqr_length(hh)); }
};

////////////////////////////////////////////////////////////////

/// Helper class to build a PolyMesh-type
template <class Traits>
struct LoopMesh_ArrayKernel_GeneratorT
{
    typedef FinalMeshItemsT<Traits, false>               MeshItems;
    typedef AttribKernelT<MeshItems, LoopConnectivity>   AttribKernel;
    typedef LoopMeshT<AttribKernel>                      Mesh;
};

/// \class LoopMesh_ArrayKernelT
/// \ingroup mesh_types_group
/// Loop mesh based on the ArrayKernel.
/// \see OpenMesh::LoopMeshT
/// \see OpenMesh::ArrayKernel
template <class Traits = DefaultTraits>
struct LoopMesh_ArrayKernelT : public LoopMesh_ArrayKernel_GeneratorT<Traits>::Mesh
{};

} // namespace OpenMesh

#endif