#include "loop_connectivity.hh"

namespace OpenMesh {

const LoopConnectivity::VertexHandle    LoopConnectivity::InvalidVertexHandle;

const LoopConnectivity::HalfedgeHandle  LoopConnectivity::InvalidHalfedgeHandle;

const LoopConnectivity::EdgeHandle      LoopConnectivity::InvalidEdgeHandle;

uint OpenMesh::LoopConnectivity::valence(VertexHandle vh) const
{
    uint count {};
    for (auto vi : vv_range(vh)) ++count;
    return count;
}

HalfedgeHandle OpenMesh::LoopConnectivity::find_halfedge(VertexHandle vh0, VertexHandle vh1) const
{
    assert(vh0.is_valid() && vh1.is_valid());

    for (auto hh : voh_range(vh0))
        if (to_vertex_handle(hh) == vh1)
            return hh;

    return InvalidHalfedgeHandle;
}

HalfedgeHandle OpenMesh::LoopConnectivity::add_edge(VertexHandle vh0, VertexHandle vh1, bool assign_random_connection)
{
    HalfedgeHandle hh0 = ArrayKernel::new_edge(vh0, vh1);
    HalfedgeHandle hh1 = opposite_halfedge_handle(hh0);

    //    --- h0 -->
    // v0 ---------- v1
    //    <-- h1 ---

    // the edge forms a loop
    set_next_halfedge_handle(hh0, hh1);
    set_next_halfedge_handle(hh1, hh0);

    if (is_isolated(vh0)) set_halfedge_handle(vh0, hh0);
    if (is_isolated(vh1)) set_halfedge_handle(vh1, hh1);

    // I don't want to do this, but since loop-mesh has been decoupled from
    // graph-mesh and adjacent halfedge list is not available anymore, one
    // cannot get all incident halfedges of a vertex without connection.
    // Unable to get all incident halfedge will cause hazard like dangling
    // vertex handle in halfedge after deletion of a vertex.
    // Meanwhile, assigning random connection will undermine original topology
    // that was constructed early, and cause trouble in other manipulations
    // in which a new edge is created like split.
    // So better to use at initialization only, like importing data.
    if (assign_random_connection)
    {
        HalfedgeHandle hh_1st = halfedge_handle(vh0);
        HalfedgeHandle hh_nth = ccw_rotated_halfedge_handle(hh_1st);
        set_next_halfedge_handle(opposite_halfedge_handle(hh_nth), hh0);
        set_next_halfedge_handle(hh1, hh_1st);

        hh_1st = halfedge_handle(vh1);
        hh_nth = ccw_rotated_halfedge_handle(hh_1st);
        set_next_halfedge_handle(opposite_halfedge_handle(hh_nth), hh1);
        set_next_halfedge_handle(hh0, hh_1st);
    }

    return hh0;
}

HalfedgeHandle OpenMesh::LoopConnectivity::add_edge(HalfedgeHandle hh0, HalfedgeHandle hh1)
{
    VertexHandle vh0 = to_vertex_handle(hh0);
    VertexHandle vh1 = to_vertex_handle(hh1);

    assert(vh0 != vh1);

    HalfedgeHandle hh0_next = next_halfedge_handle(hh0);
    HalfedgeHandle hh1_next = next_halfedge_handle(hh1);

    HalfedgeHandle hh = add_edge(vh0, vh1, false);
    HalfedgeHandle hh_opp = opposite_halfedge_handle(hh);

    set_next_halfedge_handle(hh0, hh);
    set_next_halfedge_handle(hh, hh1_next);
    set_next_halfedge_handle(hh1, hh_opp);
    set_next_halfedge_handle(hh_opp, hh0_next);

    return hh;
}

void OpenMesh::LoopConnectivity::delete_edge(EdgeHandle eh, bool delete_isolated_vertices)
{
    HalfedgeHandle hh0 = halfedge_handle(eh, 0);
    HalfedgeHandle hh1 = halfedge_handle(eh, 1);

    VertexHandle vh0 = from_vertex_handle(hh0);
    VertexHandle vh1 = from_vertex_handle(hh1);

    HalfedgeHandle hh0_next = next_halfedge_handle(hh0); // or hh1
    HalfedgeHandle hh0_prev = prev_halfedge_handle(hh0); // or hh1
    HalfedgeHandle hh1_next = next_halfedge_handle(hh1); // or hh0
    HalfedgeHandle hh1_prev = prev_halfedge_handle(hh1); // or hh0

    if (halfedge_handle(vh0) == hh0) // need to update vh0's outgoing halfedge
    {
        if (hh1_next == hh0)  set_isolated(vh0); // vh0 is an endian vertex
        else set_halfedge_handle(vh0, hh1_next);
    }

    if (halfedge_handle(vh1) == hh1) // need to update vh1's outgoing halfedge
    {
        if (hh0_next == hh1)  set_isolated(vh1); // vh1 is an endian vertex
        else set_halfedge_handle(vh1, hh0_next);
    }

    //if (hh1_next != hh0) // v0 is not an endian vertex in eh
    set_next_halfedge_handle(hh0_prev, hh1_next);

    //if (hh0_next != hh1) // v1 is not an endian vertex in eh
    set_next_halfedge_handle(hh1_prev, hh0_next);

    status(eh).set_deleted(true);

    if (has_halfedge_status())
    {
        status(hh0).set_deleted(true);
        status(hh1).set_deleted(true);
    }

    if (delete_isolated_vertices)
    {
        if (is_isolated(vh0)) status(vh0).set_deleted(true);
        if (is_isolated(vh1)) status(vh1).set_deleted(true);
    }
}

void OpenMesh::LoopConnectivity::delete_vertex(VertexHandle vh, bool delete_isolated_vertices)
{
    // make every incoming halfedge of vh an endian on its edge
    for (HalfedgeHandle hh : vih_range(vh))
        set_next_halfedge_handle(hh, opposite_halfedge_handle(hh));

    // remove incident edges
    std::vector<HalfedgeHandle> hhs;
    hhs.reserve(16);

    for (auto hh : voh_range(vh))
        hhs.push_back(hh);

    for (auto hh : hhs)
        delete_edge(edge_handle(hh), delete_isolated_vertices);

    status(vh).set_deleted(true);
}

void OpenMesh::LoopConnectivity::split(EdgeHandle eh, VertexHandle vh, bool copy_properties)
{
    HalfedgeHandle hh0 = halfedge_handle(eh, 0);
    HalfedgeHandle hh1 = halfedge_handle(eh, 1);

    VertexHandle vh0 = from_vertex_handle(hh0);
    VertexHandle vh1 = from_vertex_handle(hh1);

    assert(is_isolated(vh));
    assert(vh0 != vh && vh1 != vh);

    HalfedgeHandle hh0_next = next_halfedge_handle(hh0); // or hh1
    HalfedgeHandle hh0_prev = prev_halfedge_handle(hh0); // or hh1
    HalfedgeHandle hh1_next = next_halfedge_handle(hh1); // or hh0
    HalfedgeHandle hh1_prev = prev_halfedge_handle(hh1); // or hh0

    bool is_end0 = (hh1_next == hh0);
    bool is_end1 = (hh0_next == hh1);

    //       --- h0 -->
    // v0 ----------------- v1
    //       <-- h1 ---

    set_vertex_handle(hh0, vh); // set hh0 point to the new vertex

    HalfedgeHandle hh2 = add_edge(vh, vh1, false);
    HalfedgeHandle hh3 = opposite_halfedge_handle(hh2);

    set_halfedge_handle(vh, hh2); // assign new vertex with new edge

    if (halfedge_handle(vh1) == hh1) // adjust vh1's halfedge if necessary
        set_halfedge_handle(vh1, hh3);

    //    - h0 ->   - h2 ->
    // v0 ------- v ------- v1
    //    <- h1 -   <- h3 -

    set_next_halfedge_handle(hh0, hh2);
    set_next_halfedge_handle(hh3, hh1);

    if (!is_end0)
    {
        set_next_halfedge_handle(hh0_prev, hh0);
        set_next_halfedge_handle(hh1, hh1_next);
    }
    else
    {
        set_next_halfedge_handle(hh1, hh0);
    }

    if (!is_end1)
    {
        set_next_halfedge_handle(hh2, hh0_next);
        set_next_halfedge_handle(hh1_prev, hh3);
    }
    else
    {
        set_next_halfedge_handle(hh2, hh3);
    }

    EdgeHandle ei = edge_handle(hh2);

    if (copy_properties)
        copy_all_properties(eh, ei, true);
}

void OpenMesh::LoopConnectivity::collapse(HalfedgeHandle hh0, bool delete_isolated_vertices)
{
    HalfedgeHandle hh1 = opposite_halfedge_handle(hh0);

    VertexHandle vh0 = from_vertex_handle(hh0);
    VertexHandle vh1 = from_vertex_handle(hh1);

    HalfedgeHandle hh0_next = next_halfedge_handle(hh0); // or hh1
    HalfedgeHandle hh0_prev = prev_halfedge_handle(hh0); // or hh1
    HalfedgeHandle hh1_next = next_halfedge_handle(hh1); // or hh0
    HalfedgeHandle hh1_prev = prev_halfedge_handle(hh1); // or hh0

    bool is_end0 = (hh1_next == hh0);
    bool is_end1 = (hh0_next == hh1);

    // \ h0p          h0n / 
    //  \                /  
    //   \   -- h0 ->   /   
    //    v0 -------- v1    
    //   /   <- h1 --   \   
    //  /                \  
    // / h1n          h1p \ 

    // make all vh0's incoming halfedges point to vh1
    for (auto hh : vih_range(vh0)) if (hh != hh1)
        set_vertex_handle(hh, vh1);

    // adjust vh1's outgoing halfedges if necessary
    if (halfedge_handle(vh1) == hh1)
        set_halfedge_handle(vh1, hh0_next);

    delete_edge(edge_handle(hh0), false);

    status(vh0).set_deleted(true);

    // \ h0p      h0n / 
    //  \            /  
    //   \          /   
    //        v1        
    //   /          \   
    //  /            \  
    // / h1n      h1p \ 

    if (!is_end0 && is_end1)
    {
        set_next_halfedge_handle(hh0_prev, hh1_next);
    }
    else if (is_end0 && !is_end1)
    {
        set_next_halfedge_handle(hh1_prev, hh0_next);
    }
    else if (!is_end0 && !is_end1)
    {
        set_next_halfedge_handle(hh0_prev, hh0_next);
        set_next_halfedge_handle(hh1_prev, hh1_next);
    }

    if (delete_isolated_vertices)
        if (is_isolated(vh1)) status(vh1).set_deleted(true);
}

void OpenMesh::LoopConnectivity::merge_vertices(HalfedgeHandle hh0, HalfedgeHandle hh1)
{
    HalfedgeHandle hh2 = next_halfedge_handle(hh0);
    HalfedgeHandle hh3 = next_halfedge_handle(hh1);

    VertexHandle vh0 = to_vertex_handle(hh0);
    VertexHandle vh1 = to_vertex_handle(hh1);

    // -- h0 -> v0 -- h2 ->
    //
    // <- h3 -- v1 <- h1 --

    if (vh0 == vh1)
    {
        set_next_halfedge_handle(hh0, hh3);
        set_next_halfedge_handle(hh1, hh2);
        return;
    }

    // iff vh0 and vh1 are connected, collapse vh0 into vh1
    for (auto hh : vih_range(vh0)) if (from_vertex_handle(hh) == vh1)
    {
        collapse(hh, true);
        return;
    }

    // make all vh0's incoming halfedges point to vh1
    for (auto hh : vih_range(vh0))
        set_vertex_handle(hh, vh1);

    status(vh0).set_deleted(true);

    // -- h0 ->    -- h2 ->
    //          v1
    // <- h3 --    <- h1 --

    set_next_halfedge_handle(hh0, hh3);
    set_next_halfedge_handle(hh1, hh2);
}

} // namespace OpenMesh
