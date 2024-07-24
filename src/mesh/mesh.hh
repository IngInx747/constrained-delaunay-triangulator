#ifndef MESH_DEFINITION_HH
#define MESH_DEFINITION_HH

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Casts.hh>
#include "loop_mesh.hh"
#include "property.hh"
#include "vector_n.hh"

using Hh = OpenMesh::HalfedgeHandle;
using Vh = OpenMesh::VertexHandle;
using Fh = OpenMesh::FaceHandle;
using Eh = OpenMesh::EdgeHandle;
using Mh = OpenMesh::MeshHandle;

struct MeshTraits : public OpenMesh::DefaultTraitsDouble
{
    // Default types
    typedef double TexCoord1D;
    typedef Vec2   TexCoord2D;
    typedef Vec3   TexCoord3D;

    // Default attributes
    VertexAttributes   (OpenMesh::Attributes::Status);
    FaceAttributes     (OpenMesh::Attributes::Status);
    EdgeAttributes     (OpenMesh::Attributes::Status);
    HalfedgeAttributes (OpenMesh::Attributes::Status);

    // Customized attributes
    VertexTraits   {};
    FaceTraits     {};
    EdgeTraits     {};
    HalfedgeTraits {};
};

using PolyMesh = OpenMesh::PolyMesh_ArrayKernelT<MeshTraits>;

using TriMesh = OpenMesh::TriMesh_ArrayKernelT<MeshTraits>;

using LoopMesh = OpenMesh::LoopMesh_ArrayKernelT<MeshTraits>;

template <class MeshT>
inline bool is_marked(const MeshT &mesh, const Vh &vh) { return mesh.status(vh).selected(); }
template <class MeshT>
inline void set_marked(MeshT &mesh, const Vh &vh, const bool val) { mesh.status(vh).set_selected(val); }

template <class MeshT>
inline bool is_marked(const MeshT &mesh, const Fh &fh) { return mesh.status(fh).selected(); }
template <class MeshT>
inline void set_marked(MeshT &mesh, const Fh &fh, const bool val) { mesh.status(fh).set_selected(val); }

template <class MeshT>
inline bool is_marked(const MeshT &mesh, const Eh &eh) { return mesh.status(eh).selected(); }
template <class MeshT>
inline void set_marked(MeshT &mesh, const Eh &eh, const bool val) { mesh.status(eh).set_selected(val); }

template <class MeshT>
inline bool is_marked(const MeshT &mesh, const Hh &hh) { return mesh.status(hh).selected(); }
template <class MeshT>
inline void set_marked(MeshT &mesh, const Hh &hh, const bool val) { mesh.status(hh).set_selected(val); }

template <class MeshT>
inline bool is_sharp(const MeshT &mesh, const Vh &vh) { return mesh.status(vh).feature(); }
template <class MeshT>
inline void set_sharp(MeshT &mesh, const Vh &vh, const bool val) { mesh.status(vh).set_feature(val); }

template <class MeshT>
inline bool is_sharp(const MeshT &mesh, const Eh &eh) { return mesh.status(eh).feature(); }
template <class MeshT>
inline void set_sharp(MeshT &mesh, const Eh &eh, const bool val) { mesh.status(eh).set_feature(val); }

template <class MeshT>
inline bool is_fixed(const MeshT &mesh, const Eh &eh) { return mesh.status(eh).locked(); }
template <class MeshT>
inline void set_fixed(MeshT &mesh, const Eh &eh, const bool val) { mesh.status(eh).set_locked(val); }

template <class MeshT>
inline bool is_fixed(const MeshT &mesh, const Vh &vh) { return mesh.status(vh).locked(); }
template <class MeshT>
inline void set_fixed(MeshT &mesh, const Vh &vh, const bool val) { mesh.status(vh).set_locked(val); }

template <class MeshT>
inline bool is_sync(const MeshT &mesh, const Hh &hh)
{
    return mesh.halfedge_handle(mesh.edge_handle(hh), 0) == hh;
}

template <class MeshT>
inline void print_handle(const MeshT &mesh, const Vh &vh, const int offset = 0)
{
    printf("(%d)", vh.idx() + offset);
}

template <class MeshT>
inline void print_handle(const MeshT &mesh, const Hh &hh, const int offset = 0)
{
    printf("(%d, %d)",
        mesh.from_vertex_handle(hh).idx() + offset,
        mesh.to_vertex_handle(hh).idx()   + offset);
}

template <class MeshT>
inline void print_handle(const MeshT &mesh, const Eh &eh, const int offset = 0)
{
    printf("(%d, %d)",
        mesh.from_vertex_handle(mesh.halfedge_handle(eh, 0)).idx() + offset,
        mesh.to_vertex_handle(mesh.halfedge_handle(eh, 0)).idx()   + offset);
}

inline void print_handle(const TriMesh &mesh, const Fh &fh, const int offset = 0)
{
    printf("(%d, %d, %d)",
        mesh.from_vertex_handle(mesh.halfedge_handle(fh)).idx() + offset,
        mesh.to_vertex_handle(mesh.halfedge_handle(fh)).idx()   + offset,
        mesh.to_vertex_handle(mesh.next_halfedge_handle(mesh.halfedge_handle(fh))).idx() + offset);
}

inline void print_handle(const PolyMesh &mesh, const Fh &fh, const int offset = 0)
{
    int nh {};
    for (auto hdge : mesh.fh_range(fh))
        printf("%s%d", (nh++) ? ", " : "(",
            hdge.to().idx() + offset);
    printf(")");
}

const char *var_v_index();
const char *var_f_index();
const char *var_e_index();
const char *var_h_index();

const char *var_v_label();
const char *var_e_label();
const char *var_h_label();
const char *var_f_label();

const char *var_m_name();
const char *var_m_path();

#endif