#ifndef OPENMESH_LOOP_CONNECTIVITY_HH
#define OPENMESH_LOOP_CONNECTIVITY_HH

#include <OpenMesh/Core/Mesh/ArrayKernel.hh>
#include "iterator.hh"
#include "loop_circulator.hh"

namespace OpenMesh {

class LoopConnectivity : public ArrayKernel
{
public:

    /// \name Mesh Handles
    //@{
    /// Invalid handle
    static const VertexHandle   InvalidVertexHandle;
    /// Invalid handle
    static const HalfedgeHandle InvalidHalfedgeHandle;
    /// Invalid handle
    static const EdgeHandle     InvalidEdgeHandle;
    //@}

    typedef LoopConnectivity   This;

public: /// Linear iterator

    typedef Iterator2::GenericIteratorT<
            This,
            This::VertexHandle,
            ArrayKernel,
            &ArrayKernel::has_vertex_status,
            &ArrayKernel::n_vertices> VertexIter;

    typedef Iterator2::GenericIteratorT<
            This,
            This::HalfedgeHandle,
            ArrayKernel,
            &ArrayKernel::has_halfedge_status,
            &ArrayKernel::n_halfedges> HalfedgeIter;

    typedef Iterator2::GenericIteratorT<
            This,
            This::EdgeHandle,
            ArrayKernel,
            &ArrayKernel::has_edge_status,
            &ArrayKernel::n_edges> EdgeIter;

    typedef Iterator2::GenericIteratorT<
            This,
            This::FaceHandle,
            ArrayKernel,
            &ArrayKernel::has_face_status,
            &ArrayKernel::n_faces> FaceIter;

    typedef VertexIter      ConstVertexIter;
    typedef HalfedgeIter    ConstHalfedgeIter;
    typedef EdgeIter        ConstEdgeIter;
    typedef FaceIter        ConstFaceIter;

    typedef VertexHandle    VHandle;
    typedef HalfedgeHandle  HHandle;
    typedef EdgeHandle      EHandle;
    typedef FaceHandle      FHandle;

    typedef VertexIter      VIter;
    typedef HalfedgeIter    HIter;
    typedef EdgeIter        EIter;
    typedef FaceIter        FIter;

    typedef ConstVertexIter     CVIter;
    typedef ConstHalfedgeIter   CHIter;
    typedef ConstEdgeIter       CEIter;
    typedef ConstFaceIter       CFIter;

    /// Begin iterator for vertices
    inline VertexIter vertices_begin()
    { return VertexIter(*this, VertexHandle(0), false); }

    /// Const begin iterator for vertices
    inline ConstVertexIter vertices_begin() const
    { return ConstVertexIter(*this, VertexHandle(0), false); }

    /// End iterator for vertices
    inline VertexIter vertices_end()
    { return VertexIter(*this, VertexHandle(int(n_vertices())), false); }

    /// Const end iterator for vertices
    inline ConstVertexIter vertices_end() const
    { return ConstVertexIter(*this, VertexHandle(int(n_vertices())), false); }

    /// Begin iterator for halfedges
    inline HalfedgeIter halfedges_begin()
    { return HalfedgeIter(*this, HalfedgeHandle(0), false); }

    /// Const begin iterator for halfedges
    inline ConstHalfedgeIter halfedges_begin() const
    { return ConstHalfedgeIter(*this, HalfedgeHandle(0), false); }

    /// End iterator for halfedges
    inline HalfedgeIter halfedges_end()
    { return HalfedgeIter(*this, HalfedgeHandle(int(n_halfedges())), false); }

    /// Const end iterator for halfedges
    inline ConstHalfedgeIter halfedges_end() const
    { return ConstHalfedgeIter(*this, HalfedgeHandle(int(n_halfedges())), false); }

    /// Begin iterator for edges
    inline EdgeIter edges_begin()
    { return EdgeIter(*this, EdgeHandle(0), false); }

    /// Const begin iterator for edges
    inline ConstEdgeIter edges_begin() const
    { return ConstEdgeIter(*this, EdgeHandle(0), false); }

    /// End iterator for edges
    inline EdgeIter edges_end()
    { return EdgeIter(*this, EdgeHandle(int(n_edges())), false); }

    /// Const end iterator for edges
    inline ConstEdgeIter edges_end() const
    { return ConstEdgeIter(*this, EdgeHandle(int(n_edges())), false); }

    /// Begin iterator for vertices
    inline VertexIter vertices_sbegin()
    { return VertexIter(*this, VertexHandle(0), true); }

    /// Const begin iterator for vertices
    inline ConstVertexIter vertices_sbegin() const
    { return ConstVertexIter(*this, VertexHandle(0), true); }

    /// Begin iterator for halfedges
    inline HalfedgeIter halfedges_sbegin()
    { return HalfedgeIter(*this, HalfedgeHandle(0), true); }

    /// Const begin iterator for halfedges
    inline ConstHalfedgeIter halfedges_sbegin() const
    { return ConstHalfedgeIter(*this, HalfedgeHandle(0), true); }

    /// Begin iterator for edges
    inline EdgeIter edges_sbegin()
    { return EdgeIter(*this, EdgeHandle(0), true); }

    /// Const begin iterator for edges
    inline ConstEdgeIter edges_sbegin() const
    { return ConstEdgeIter(*this, EdgeHandle(0), true); }

    /// Warning: Face could be meaningless if the underlying
    /// loop is not closed. Use at your own risk.

    /// Begin iterator for faces
    inline FaceIter faces_begin()
    { return FaceIter(*this, FaceHandle(0), false); }

    /// Const begin iterator for faces
    inline ConstFaceIter faces_begin() const
    { return ConstFaceIter(*this, FaceHandle(0), false); }

    /// Begin iterator for faces
    inline FaceIter faces_sbegin()
    { return FaceIter(*this, FaceHandle(0), true); }

    /// Const begin iterator for faces
    inline ConstFaceIter faces_sbegin() const
    { return ConstFaceIter(*this, FaceHandle(0), true); }

    /// End iterator for faces
    inline FaceIter faces_end()
    { return FaceIter(*this, FaceHandle(int(n_faces())), false); }

    /// Const end iterator for faces
    inline ConstFaceIter faces_end() const
    { return ConstFaceIter(*this, FaceHandle(int(n_faces())), false); }

public: /// Linear range

    typedef Iterator2::EntityRange<
            Iterator2::RangeTraitT<
            const LoopConnectivity,
            LoopConnectivity::ConstVertexIter,
            &LoopConnectivity::vertices_begin,
            &LoopConnectivity::vertices_end>> ConstVertexRange;

    typedef Iterator2::EntityRange<
            Iterator2::RangeTraitT<
            const LoopConnectivity,
            LoopConnectivity::ConstVertexIter,
            &LoopConnectivity::vertices_sbegin,
            &LoopConnectivity::vertices_end>> ConstVertexRangeSkipping;

    typedef Iterator2::EntityRange<
            Iterator2::RangeTraitT<
            const LoopConnectivity,
            LoopConnectivity::ConstHalfedgeIter,
            &LoopConnectivity::halfedges_begin,
            &LoopConnectivity::halfedges_end>> ConstHalfedgeRange;

    typedef Iterator2::EntityRange<
            Iterator2::RangeTraitT<
            const LoopConnectivity,
            LoopConnectivity::ConstHalfedgeIter,
            &LoopConnectivity::halfedges_sbegin,
            &LoopConnectivity::halfedges_end>> ConstHalfedgeRangeSkipping;

    typedef Iterator2::EntityRange<
            Iterator2::RangeTraitT<
            const LoopConnectivity,
            LoopConnectivity::ConstEdgeIter,
            &LoopConnectivity::edges_begin,
            &LoopConnectivity::edges_end>> ConstEdgeRange;

    typedef Iterator2::EntityRange<
            Iterator2::RangeTraitT<
            const LoopConnectivity,
            LoopConnectivity::ConstEdgeIter,
            &LoopConnectivity::edges_sbegin,
            &LoopConnectivity::edges_end>> ConstEdgeRangeSkipping;

    typedef Iterator2::EntityRange<
            Iterator2::RangeTraitT<
            const LoopConnectivity,
            LoopConnectivity::ConstFaceIter,
            &LoopConnectivity::faces_begin,
            &LoopConnectivity::faces_end>> ConstFaceRange;

    typedef Iterator2::EntityRange<
            Iterator2::RangeTraitT<
            const LoopConnectivity,
            LoopConnectivity::ConstFaceIter,
            &LoopConnectivity::faces_sbegin,
            &LoopConnectivity::faces_end>> ConstFaceRangeSkipping;

    /// @return The vertices as a range object suitable
    /// for C++11 range based for loops. Will skip deleted vertices.
    inline ConstVertexRangeSkipping vertices() const
    { return ConstVertexRangeSkipping(*this); }

    /// @return The vertices as a range object suitable
    /// for C++11 range based for loops. Will include deleted vertices.
    inline ConstVertexRange all_vertices() const
    { return ConstVertexRange(*this); }

    /// @return The halfedges as a range object suitable
    /// for C++11 range based for loops. Will skip deleted halfedges.
    inline ConstHalfedgeRangeSkipping halfedges() const
    { return ConstHalfedgeRangeSkipping(*this); }

    /// @return The halfedges as a range object suitable
    /// for C++11 range based for loops. Will include deleted halfedges.
    inline ConstHalfedgeRange all_halfedges() const
    { return ConstHalfedgeRange(*this); }

    /// @return The edges as a range object suitable
    /// for C++11 range based for loops. Will skip deleted edges.
    inline ConstEdgeRangeSkipping edges() const
    { return ConstEdgeRangeSkipping(*this); }

    /// @return The edges as a range object suitable
    /// for C++11 range based for loops. Will include deleted edges.
    inline ConstEdgeRange all_edges() const
    { return ConstEdgeRange(*this); }

    /// @return The faces as a range object suitable
    /// for C++11 range based for loops. Will skip deleted faces.
    inline ConstFaceRangeSkipping faces() const
    { return ConstFaceRangeSkipping(*this); }

    /// @return The faces as a range object suitable
    /// for C++11 range based for loops. Will include deleted faces.
    inline ConstFaceRange all_faces() const
    { return ConstFaceRange(*this); }

public: /// Circulator

    struct VertexVertexTraits
    {
        using Mesh = This;
        using CenterEntityHandle = This::VertexHandle;
        using ValueHandle = This::VertexHandle;
        inline static ValueHandle toHandle(const Mesh *mesh, This::HalfedgeHandle hh) { return mesh->to_vertex_handle(hh); }
    };

    struct VertexOHalfedgeTraits
    {
        using Mesh = This;
        using CenterEntityHandle = This::VertexHandle;
        using ValueHandle = This::HalfedgeHandle;
        inline static ValueHandle toHandle(const Mesh*, This::HalfedgeHandle hh) { return hh; }
    };

    struct VertexIHalfedgeTraits
    {
        using Mesh = This;
        using CenterEntityHandle = This::VertexHandle;
        using ValueHandle = This::HalfedgeHandle;
        inline static ValueHandle toHandle(const Mesh *mesh, This::HalfedgeHandle hh) { return mesh->opposite_halfedge_handle(hh); }
    };

    struct VertexEdgeTraits
    {
        using Mesh = This;
        using CenterEntityHandle = This::VertexHandle;
        using ValueHandle = This::EdgeHandle;
        inline static ValueHandle toHandle(const Mesh *mesh, This::HalfedgeHandle hh) { return static_cast<const ArrayKernel*>(mesh)->edge_handle(hh); }
    };

    struct EdgeVertexTraits
    {
        using Mesh = This;
        using CenterEntityHandle = This::EdgeHandle;
        using ValueHandle = This::VertexHandle;
        inline static ValueHandle toHandle(const Mesh *mesh, This::HalfedgeHandle hh) { return static_cast<const ArrayKernel*>(mesh)->from_vertex_handle(hh); }
    };

    struct EdgeHalfedgeTraits
    {
        using Mesh = This;
        using CenterEntityHandle = This::EdgeHandle;
        using ValueHandle = This::HalfedgeHandle;
        inline static ValueHandle toHandle(const Mesh*, This::HalfedgeHandle hh) { return hh; }
    };

    struct FaceHalfedgeTraits
    {
        using Mesh = This;
        using CenterEntityHandle = This::FaceHandle;
        using ValueHandle = This::HalfedgeHandle;
        static ValueHandle toHandle(const Mesh*, This::HalfedgeHandle hh) { return hh; }
    };

    typedef Circulator3::GenericCirculatorT<VertexVertexTraits>       VertexVertexIter;
    typedef Circulator3::GenericCirculatorT<VertexOHalfedgeTraits>    VertexOHalfedgeIter;
    typedef Circulator3::GenericCirculatorT<VertexIHalfedgeTraits>    VertexIHalfedgeIter;
    typedef Circulator3::GenericCirculatorT<VertexEdgeTraits>         VertexEdgeIter;
    typedef Circulator3::GenericCirculatorT<EdgeVertexTraits>         EdgeVertexIter;
    typedef Circulator3::GenericCirculatorT<EdgeHalfedgeTraits>       EdgeHalfedgeIter;
    typedef Circulator3::GenericCirculatorT<FaceHalfedgeTraits>       HalfedgeLoopIter;

    typedef VertexVertexIter            ConstVertexVertexIter;
    typedef VertexOHalfedgeIter         ConstVertexOHalfedgeIter;
    typedef VertexIHalfedgeIter         ConstVertexIHalfedgeIter;
    typedef VertexEdgeIter              ConstVertexEdgeIter;
    typedef EdgeVertexIter              ConstEdgeVertexIter;
    typedef EdgeHalfedgeIter            ConstEdgeHalfedgeIter;
    typedef HalfedgeLoopIter            ConstHalfedgeLoopIter;

    typedef VertexVertexIter            VVIter;
    typedef VertexOHalfedgeIter         VOHIter;
    typedef VertexIHalfedgeIter         VIHIter;
    typedef VertexEdgeIter              VEIter;
    typedef EdgeVertexIter              EVIter;
    typedef EdgeHalfedgeIter            EHIter;

    typedef ConstVertexVertexIter       CVVIter;
    typedef ConstVertexOHalfedgeIter    CVOHIter;
    typedef ConstVertexIHalfedgeIter    CVIHIter;
    typedef ConstVertexEdgeIter         CVEIter;
    typedef ConstEdgeVertexIter         CEVIter;
    typedef ConstEdgeHalfedgeIter       CEHIter;

    inline ConstVertexVertexIter vv_iter(ArrayKernel::VertexHandle vh)
    { return VertexVertexIter(*this, vh); }

    inline VertexIHalfedgeIter vih_iter(ArrayKernel::VertexHandle vh)
    {  return VertexIHalfedgeIter(*this, vh); }

    inline VertexOHalfedgeIter voh_iter(ArrayKernel::VertexHandle vh)
    {  return VertexOHalfedgeIter(*this, vh); }

    inline VertexEdgeIter ve_iter(ArrayKernel::VertexHandle vh)
    {  return VertexEdgeIter(*this, vh); }

    inline ConstVertexVertexIter cvv_iter(ArrayKernel::VertexHandle vh) const
    { return ConstVertexVertexIter(*this, vh); }

    inline ConstVertexIHalfedgeIter cvih_iter(ArrayKernel::VertexHandle vh) const
    { return ConstVertexIHalfedgeIter(*this, vh); }

    inline ConstVertexOHalfedgeIter cvoh_iter(ArrayKernel::VertexHandle vh) const
    { return ConstVertexOHalfedgeIter(*this, vh); }

    inline ConstVertexEdgeIter cve_iter(ArrayKernel::VertexHandle vh) const
    { return ConstVertexEdgeIter(*this, vh); }

    inline EdgeVertexIter ev_iter(ArrayKernel::EdgeHandle eh)
    { return EdgeVertexIter(*this, eh); }

    inline EdgeHalfedgeIter eh_iter(ArrayKernel::EdgeHandle eh)
    { return EdgeHalfedgeIter(*this, eh); }

    inline ConstEdgeVertexIter cev_iter(ArrayKernel::EdgeHandle eh) const
    { return ConstEdgeVertexIter(*this, eh); }

    inline ConstEdgeHalfedgeIter ceh_iter(ArrayKernel::EdgeHandle eh) const
    { return ConstEdgeHalfedgeIter(*this, eh); }

    /// circulator begin

    inline VertexVertexIter vv_begin(VertexHandle vh)
    { return VertexVertexIter(*this, vh); }

    inline VertexIHalfedgeIter vih_begin(VertexHandle vh)
    { return VertexIHalfedgeIter(*this, vh); }

    inline VertexOHalfedgeIter voh_begin(VertexHandle vh)
    { return VertexOHalfedgeIter(*this, vh); }

    inline VertexEdgeIter ve_begin(VertexHandle vh)
    { return VertexEdgeIter(*this, vh); }

    inline ConstVertexVertexIter cvv_begin(VertexHandle vh) const
    { return ConstVertexVertexIter(*this, vh); }

    inline ConstVertexIHalfedgeIter cvih_begin(VertexHandle vh) const
    { return ConstVertexIHalfedgeIter(*this, vh); }

    inline ConstVertexOHalfedgeIter cvoh_begin(VertexHandle vh) const
    { return ConstVertexOHalfedgeIter(*this, vh); }

    inline ConstVertexEdgeIter cve_begin(VertexHandle vh) const
    { return ConstVertexEdgeIter(*this, vh); }

    inline EdgeVertexIter ev_begin(EdgeHandle eh)
    { return EdgeVertexIter(*this, eh); }

    inline EdgeHalfedgeIter eh_begin(EdgeHandle eh)
    { return EdgeHalfedgeIter(*this, eh); }

    inline ConstEdgeVertexIter cev_begin(EdgeHandle eh) const
    { return ConstEdgeVertexIter(*this, eh); }

    inline ConstEdgeHalfedgeIter ceh_begin(EdgeHandle eh) const
    { return ConstEdgeHalfedgeIter(*this, eh); }

    inline HalfedgeLoopIter hl_begin(HalfedgeHandle hh) const
    { return HalfedgeLoopIter(*this, hh); }

    inline ConstHalfedgeLoopIter chl_begin(HalfedgeHandle hh) const
    { return ConstHalfedgeLoopIter(*this, hh); }

    /// circulator end

    inline VertexVertexIter vv_end(VertexHandle vh)
    { return VertexVertexIter(*this, vh, true); }

    inline VertexIHalfedgeIter vih_end(VertexHandle vh)
    { return VertexIHalfedgeIter(*this, vh, true); }

    inline VertexOHalfedgeIter voh_end(VertexHandle vh)
    { return VertexOHalfedgeIter(*this, vh, true); }

    inline VertexEdgeIter ve_end(VertexHandle vh)
    { return VertexEdgeIter(*this, vh, true); }

    inline ConstVertexVertexIter cvv_end(VertexHandle vh) const
    { return ConstVertexVertexIter(*this, vh, true); }

    inline ConstVertexIHalfedgeIter cvih_end(VertexHandle vh) const
    { return ConstVertexIHalfedgeIter(*this, vh, true); }

    inline ConstVertexOHalfedgeIter cvoh_end(VertexHandle vh) const
    { return ConstVertexOHalfedgeIter(*this, vh, true); }

    inline ConstVertexEdgeIter cve_end(VertexHandle vh) const
    { return ConstVertexEdgeIter(*this, vh, true); }

    inline EdgeVertexIter ev_end(EdgeHandle eh)
    { return EdgeVertexIter(*this, eh, true); }

    inline EdgeHalfedgeIter eh_end(EdgeHandle eh)
    { return EdgeHalfedgeIter(*this, eh, true); }

    inline ConstEdgeVertexIter cev_end(EdgeHandle eh) const
    { return ConstEdgeVertexIter(*this, eh, true); }

    inline ConstEdgeHalfedgeIter ceh_end(EdgeHandle eh) const
    { return ConstEdgeHalfedgeIter(*this, eh, true); }

    inline HalfedgeLoopIter hl_end(HalfedgeHandle hh) const
    { return HalfedgeLoopIter(*this, hh, true); }

    inline ConstHalfedgeLoopIter chl_end(HalfedgeHandle hh) const
    { return ConstHalfedgeLoopIter(*this, hh, true); }

public: /// Circulator range

    typedef Circulator3::CirculatorRange<
            Circulator3::CirculatorRangeTraitT<
            LoopConnectivity,
            ConstVertexVertexIter,
            VertexHandle,
            VertexHandle,
            &LoopConnectivity::cvv_begin,
            &LoopConnectivity::cvv_end>> ConstVertexVertexRange;

    typedef Circulator3::CirculatorRange<
            Circulator3::CirculatorRangeTraitT<
            LoopConnectivity,
            ConstVertexIHalfedgeIter,
            VertexHandle,
            HalfedgeHandle,
            &LoopConnectivity::cvih_begin,
            &LoopConnectivity::cvih_end>> ConstVertexIHalfedgeRange;

    typedef Circulator3::CirculatorRange<
            Circulator3::CirculatorRangeTraitT<
            LoopConnectivity,
            ConstVertexOHalfedgeIter,
            VertexHandle,
            HalfedgeHandle,
            &LoopConnectivity::cvoh_begin,
            &LoopConnectivity::cvoh_end>> ConstVertexOHalfedgeRange;

    typedef Circulator3::CirculatorRange<
            Circulator3::CirculatorRangeTraitT<
            LoopConnectivity,
            ConstVertexEdgeIter,
            VertexHandle,
            EdgeHandle,
            &LoopConnectivity::cve_begin,
            &LoopConnectivity::cve_end>> ConstVertexEdgeRange;

    typedef Circulator3::CirculatorRange<
            Circulator3::CirculatorRangeTraitT<
            LoopConnectivity,
            ConstEdgeVertexIter,
            EdgeHandle,
            VertexHandle,
            &LoopConnectivity::cev_begin,
            &LoopConnectivity::cev_end>> ConstEdgeVertexRange;

    typedef Circulator3::CirculatorRange<
            Circulator3::CirculatorRangeTraitT<
            LoopConnectivity,
            ConstEdgeHalfedgeIter,
            EdgeHandle,
            HalfedgeHandle,
            &LoopConnectivity::ceh_begin,
            &LoopConnectivity::ceh_end>> ConstEdgeHalfedgeRange;

    typedef Circulator3::CirculatorRange<
            Circulator3::CirculatorRangeTraitT<
            LoopConnectivity,
            ConstHalfedgeLoopIter,
            HalfedgeHandle,
            HalfedgeHandle,
            &LoopConnectivity::chl_begin,
            &LoopConnectivity::chl_end>> ConstHalfedgeLoopRange;

    inline ConstVertexVertexRange vv_range(VertexHandle vh) const
    { return ConstVertexVertexRange(*this, vh); }

    inline ConstVertexIHalfedgeRange vih_range(VertexHandle vh) const
    { return ConstVertexIHalfedgeRange(*this, vh); }

    inline ConstVertexOHalfedgeRange voh_range(VertexHandle vh) const
    { return ConstVertexOHalfedgeRange(*this, vh); }

    inline ConstVertexEdgeRange ve_range(VertexHandle vh) const
    { return ConstVertexEdgeRange(*this, vh); }

    inline ConstEdgeVertexRange ev_range(EdgeHandle eh) const
    { return ConstEdgeVertexRange(*this, eh); }

    inline ConstEdgeHalfedgeRange eh_range(EdgeHandle eh) const
    { return ConstEdgeHalfedgeRange(*this, eh); }

    inline ConstHalfedgeLoopRange hl_range(HalfedgeHandle hh) const
    { return ConstHalfedgeLoopRange(*this, hh); }

public:

    LoopConnectivity() {}

    virtual ~LoopConnectivity() {}

    inline bool is_endian(HalfedgeHandle hh) const
    { return next_halfedge_handle(hh) == opposite_halfedge_handle(hh); }

    inline bool is_boundary(HalfedgeHandle hh) const
    { return !face_handle(hh).is_valid(); }

    /// Vertex valence
    uint valence(VertexHandle vh) const;

    /// Find any halfedge from vh0 to vh1. Returns invalid handle if not found.
    HalfedgeHandle find_halfedge(VertexHandle vh0, VertexHandle vh1) const;

    /// Create an edge without assigning connection if assign_random_connection is false,
    /// or assign an arbitrary connection if true.
    HalfedgeHandle add_edge(VertexHandle vh0, VertexHandle vh1, bool assign_random_connection);

    /// Create an edge and assign connection.
    HalfedgeHandle add_edge(HalfedgeHandle hh0, HalfedgeHandle hh1);

    void delete_vertex(VertexHandle vh, bool delete_isolated_vertices = true);

    void delete_edge(EdgeHandle eh, bool delete_isolated_vertices = true);

    void split(EdgeHandle eh, VertexHandle vh, bool copy_properties = true);
 
    void collapse(HalfedgeHandle hh, bool delete_isolated_vertices = true);

    void merge_vertices(HalfedgeHandle hh0, HalfedgeHandle hh1);
};

} // namespace OpenMesh

#endif