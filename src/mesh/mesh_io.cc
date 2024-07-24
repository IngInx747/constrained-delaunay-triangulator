#include <string>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "mesh.hh"

using namespace OpenMesh;

enum IO_ERR_TYPE
{
    NO_ERROR,
    INTERNAL,
    CANNOT_OPEN,
    UNSUPPORTED_FORMAT,
};

static const char *__io_err_msg[] = {
    "",
    "internal error",
    "cannot open file",
    "unsupported format"
};

inline std::string get_file_extension(const char *filename)
{
    std::string s { filename };
    auto found = s.find_last_of(".");
    return s.substr(found + 1);
}

int save_obj(
    const double *vs, const int nv,
    const int    *fs, const int nf,
    const int    *es, const int ne,
    const char   *filename,
    const std::streamsize prec = 17i64);

int save_obj(
    const double *vs, const int nv,
    const int    *fs, const int nf,
    const int    *es, const int ne,
    const char   *filename,
    const std::streamsize prec)
{
    std::ofstream out(filename, std::ios::out);
    if (!out) return out.bad();

    out << std::fixed << std::setprecision(prec);

    for (int i = 0; i < nv; ++i)
        out << "v "
            << vs[i*3 + 0] << " "
            << vs[i*3 + 1] << " "
            << vs[i*3 + 2] << " "
            << "\n";

    for (int i = 0; i < nf; ++i)
        out << "f "
            << fs[i*3 + 0] << " "
            << fs[i*3 + 1] << " "
            << fs[i*3 + 2] << " "
            << "\n";

    for (int i = 0; i < ne; ++i)
        out << "l "
            << es[i*2 + 0] << " "
            << es[i*2 + 1] << " "
            << "\n";

    return 0;
}

int save_nodes(const double *vs, const int nv, const char* filename, const std::streamsize prec)
{
    std::ofstream out(filename, std::ios::out);
    if (!out) return out.bad();

    out << std::fixed << std::setprecision(prec);

    out << nv << " 2\n";

    for (int i = 0; i < nv; ++i)
        out << i + 1 << " "
            << vs[i*2 + 0] << " "
            << vs[i*2 + 1] << " "
            << "\n";

    return 0;
}

int save_poly(const double *vs, const int nv, const int *es, const int ne, const char* filename, const std::streamsize prec)
{
    std::ofstream out(filename, std::ios::out);
    if (!out) return out.bad();

    out << std::fixed << std::setprecision(prec);

    out << nv << " 2 0 1\n";

    for (int i = 0; i < nv; ++i)
        out << i + 1 << " "
            << vs[i*2 + 0] << " "
            << vs[i*2 + 1] << " "
            << -1 << "\n";

    out << ne << " 1\n";

    for (int i = 0; i < ne; ++i)
        out << i + 1 << " "
            << es[i*2 + 0] << " "
            << es[i*2 + 1] << " "
            << -1 << "\n";

    out << "0\n";

    return 0;
}

template <class MeshT>
inline int read_mesh_builtin(MeshT &mesh, const char *filename)
{
    int err { IO_ERR_TYPE::NO_ERROR };

    IO::Options opt;

    opt += IO::Options::VertexTexCoord;
    mesh.request_vertex_texcoords2D();

    opt += IO::Options::VertexNormal;
    mesh.request_vertex_normals();

    opt += IO::Options::VertexColor;
    mesh.request_vertex_colors();

    opt += IO::Options::FaceNormal;
    mesh.request_face_normals();

    opt += IO::Options::FaceColor;
    mesh.request_face_colors();

    if (!IO::read_mesh(mesh, filename, opt)) err = IO_ERR_TYPE::INTERNAL;

    if (!opt.vertex_has_texcoord()) mesh.release_vertex_texcoords2D();

    if (!opt.vertex_has_normal()) mesh.release_vertex_normals();

    if (!opt.vertex_has_color()) mesh.release_vertex_colors();

    if (!opt.face_has_normal()) mesh.release_face_normals();

    if (!opt.face_has_color()) mesh.release_face_colors();

    return err;
}

template <class MeshT>
inline int save_mesh_builtin(const MeshT &mesh, const char *filename)
{
    int err { IO_ERR_TYPE::NO_ERROR };

    IO::Options opt;

    if (mesh.has_vertex_texcoords2D()) opt += IO::Options::VertexTexCoord;

    if (mesh.has_vertex_normals()) opt += IO::Options::VertexNormal;

    if (mesh.has_vertex_colors()) opt += IO::Options::VertexColor;

    if (mesh.has_face_normals()) opt += IO::Options::FaceNormal;

    if (mesh.has_face_colors()) opt += IO::Options::FaceColor;

    if (mesh.has_halfedge_texcoords2D()) opt += IO::Options::FaceTexCoord;

    if (!IO::write_mesh(mesh, filename, opt, 17i64)) err = IO_ERR_TYPE::INTERNAL;

    return err;
}

template <class MeshT>
inline int read_mesh_detri2(MeshT &mesh, const char *filename)
{
    constexpr size_t kInf = std::numeric_limits<std::streamsize>::max();

    std::ifstream in(filename, std::ios::in);
    if (!in) return IO_ERR_TYPE::CANNOT_OPEN;

    std::vector<std::tuple<int, int, int>> face_list;
    std::vector<std::tuple<int, int>> edge_list;
    int ver {}, dim {}, nv {}, nf {}, ne {};

    while (in)
    {
        std::string buf {};
        in >> buf;

        if (buf.compare("MeshVersionFormatted") == 0) // version
        {
            in >> ver;
        }
        else if (buf.compare("Dimension") == 0) // dimension
        {
            in >> dim;
        }
        else if (buf.compare("Vertices") == 0) // vertices
        {
            in >> nv;
            in.ignore(kInf, '\n');
            for (int i = 0; i < nv; ++i)
            {
                TriMesh::Point p { 0,0,0 };
                for (int j = 0; j < dim; ++j)
                    in >> p[j];
                mesh.new_vertex(p);
                // optional: process vertex tag
                in.ignore(kInf, '\n');
            }
        }
        else if (buf.compare("Triangles") == 0) // faces
        {
            in >> nf;
            in.ignore(kInf, '\n');
            face_list.reserve(nf);
            for (int i = 0; i < nf; ++i)
            {
                int ids[3];
                in >> ids[0] >> ids[1] >> ids[2];
                face_list.emplace_back(ids[0]-1, ids[1]-1, ids[2]-1);
                // optional: process facet tag
                in.ignore(kInf, '\n');
            }
        }
        else if (buf.compare("Edges") == 0) // edge
        {
            in >> ne;
            in.ignore(kInf, '\n');
            edge_list.reserve(ne);
            for (int i = 0; i < ne; ++i)
            {
                int ids[2];
                in >> ids[0] >> ids[1];
                edge_list.emplace_back(ids[0]-1, ids[1]-1);
                // optional: process edge tag
                in.ignore(kInf, '\n');
            }
        }
        else if (buf.compare("End") == 0)
        {}
    }

    for (const auto &face : face_list)
    {
        int i0 = std::get<0>(face);
        int i1 = std::get<1>(face);
        int i2 = std::get<2>(face);
        mesh.add_face({
            mesh.vertex_handle(i0),
            mesh.vertex_handle(i1),
            mesh.vertex_handle(i2)
        });
    }

    for (const auto &edge : edge_list)
    {
        int i0 = std::get<0>(edge);
        int i1 = std::get<1>(edge);
        auto hdge = mesh.find_halfedge(mesh.vertex_handle(i0), mesh.vertex_handle(i1));
        if (hdge.is_valid()) set_sharp(mesh, hdge.edge(), true);
    }

    return IO_ERR_TYPE::NO_ERROR;
}

template <class MeshT>
inline int save_mesh_detri2(MeshT &mesh, const char *filename)
{
    return IO_ERR_TYPE::UNSUPPORTED_FORMAT;
}

template <>
static int save_mesh_detri2(const TriMesh &mesh, const char *filename)
{
    std::ofstream out(filename, std::ios::out);
    if (!out) return IO_ERR_TYPE::CANNOT_OPEN;

    const_cast<TriMesh&>(mesh).garbage_collection(true, true, true);

    out << std::fixed << std::setprecision(17);

    out << "MeshVersionFormatted 1\n\nDimension\n3\n\n";

    out << "Vertices\n" << mesh.n_vertices() << "\n";

    if (!hasProperty<Vh, int>(mesh, var_v_label()))
    {
        for (auto vert : mesh.vertices())
        {
            const auto &p = mesh.point(vert);
            out << p[0] << "  "
                << p[1] << "  "
                << p[2] << "  "
                << "0\n"; // arbitrary tag
        }
    }
    else
    {
        auto v_tag = getProperty<Vh, int>(mesh, var_v_label());

        for (auto vert : mesh.vertices())
        {
            const auto &p = mesh.point(vert);
            out << p[0] << "  "
                << p[1] << "  "
                << p[2] << "  "
                << v_tag[vert] + 1 << "\n";
        }
    }

    out << "Triangles\n" << mesh.n_faces() << "\n";

    if (!hasProperty<Fh, int>(mesh, var_f_label()))
    {
        for (auto face : mesh.faces())
        {
            for (auto vert : face.vertices())
                out << vert.idx() + 1 << " ";
            out << "1\n"; // arbitrary positive tag
        }
    }
    else
    {
        auto f_tag = getProperty<Fh, int>(mesh, var_f_label());

        for (auto face : mesh.faces())
        {
            for (auto vert : face.vertices())
                out << vert.idx() + 1 << " ";
            out << f_tag[face] + 1 << "\n";
        }
    }

    int ne {};
    for (auto edge : mesh.edges()) if (is_sharp(mesh, edge)) ++ne;

    out << "Edges\n" << ne << "\n";

    if (!hasProperty<Eh, int>(mesh, var_e_label()))
    {
        for (auto edge : mesh.edges())
        {
            if (is_sharp(mesh, edge))
            {
                auto vert0 = edge.v0(), vert1 = edge.v1();
                out << vert0.idx() + 1 << " "
                    << vert1.idx() + 1 << " -1\n";
            }
        }
    }
    else
    {
        auto e_tag = getProperty<Eh, int>(mesh, var_e_label());

        for (auto edge : mesh.edges())
        {
            if (is_sharp(mesh, edge))
            {
                auto vert0 = edge.v0(), vert1 = edge.v1();
                out << vert0.idx() + 1 << " "
                    << vert1.idx() + 1 << " "
                    << e_tag[edge] + 1 << "\n";
            }
        }
    }

    out << "End\n";

    return IO_ERR_TYPE::NO_ERROR;
}

template <class MeshT>
int read_mesh(MeshT &mesh, const char *filename)
{
    int err { IO_ERR_TYPE::NO_ERROR };

    const auto ext = get_file_extension(filename);

    if (ext.compare("obj") == 0 || ext.compare("off") == 0)
    {
        err = read_mesh_builtin(mesh, filename);
    }
    else if (ext.compare("mesh") == 0)
    {
        err = read_mesh_detri2(mesh, filename);
    }
    else
    {
        err = IO_ERR_TYPE::UNSUPPORTED_FORMAT;
    }

    if (err) fprintf(stderr, "read_mesh error %d: %s\n", err, __io_err_msg[err]);

    return err;
}

template <class MeshT>
int save_mesh(const MeshT &mesh, const char *filename)
{
    int err { IO_ERR_TYPE::NO_ERROR };

    const auto ext = get_file_extension(filename);

    if (ext.compare("obj") == 0 || ext.compare("off") == 0)
    {
        err = save_mesh_builtin(mesh, filename);
    }
    else if (ext.compare("mesh") == 0)
    {
        err = save_mesh_detri2(mesh, filename);
    }
    else
    {
        err = IO_ERR_TYPE::UNSUPPORTED_FORMAT;
    }

    if (err) fprintf(stderr, "save_mesh error %d: %s\n", err, __io_err_msg[err]);

    return err;
}

template
int read_mesh(TriMesh&, const char*);

template
int read_mesh(PolyMesh&, const char*);

template
int save_mesh(const TriMesh&, const char*);

template
int save_mesh(const PolyMesh&, const char*);

template <class MeshT>
inline int read_poly_detri2(MeshT &mesh, const char *filename)
{
    return IO_ERR_TYPE::UNSUPPORTED_FORMAT;
}

//template <>
//struct std::hash<Int2>
//{
//    inline size_t operator()(const Int2 &i2) const noexcept
//    {
//        size_t h1 = hash<int>{}(i2[0]);
//        size_t h2 = hash<int>{}(i2[1]);
//        return h1 == h2 ? h1 : h1 ^ h2;
//    }
//};

//inline bool operator==(const Int2 &lhs, const Int2 &rhs)
//{ return lhs[0] == rhs[0] && lhs[1] == rhs[1] || lhs[0] == rhs[1] && lhs[1] == rhs[0]; }

template <>
inline int read_poly_detri2(LoopMesh &mesh, const char *filename)
{
    constexpr size_t kInf = std::numeric_limits<std::streamsize>::max();

    std::ifstream in(filename, std::ios::in);
    if (!in) return IO_ERR_TYPE::CANNOT_OPEN;

    int ver {}, dim {}, nv {}, ne {};
    std::vector<Int2> edge_list;

    in >> nv;
    in >> dim;
    in.ignore(kInf, '\n');

    for (int i = 0; i < nv; ++i)
    {
        int vid {}; in >> vid;
        TriMesh::Point p { 0,0,0 };
        for (int j = 0; j < dim; ++j)
            in >> p[j];
        Vh vh = mesh.add_vertex(p);
        // optional: process vertex tag
        in.ignore(kInf, '\n');
    }

    in >> ne;
    in.ignore(kInf, '\n');

    struct HashInt2
    {
        inline size_t operator()(const Int2 &lhs) const
        {
            size_t h1 = std::hash<int>{}(lhs[0]);
            size_t h2 = std::hash<int>{}(lhs[1]);
            return h1 == h2 ? h1 : h1 ^ h2;
        }
    };

    struct EqualToInt2
    {
        inline bool operator()(const Int2 &lhs, const Int2 &rhs) const
        {
            return
                lhs[0] == rhs[0] && lhs[1] == rhs[1] ||
                lhs[0] == rhs[1] && lhs[1] == rhs[0];
        }
    };

    std::unordered_set<Int2, HashInt2, EqualToInt2> vvs {};

    for (int i = 0; i < ne; ++i)
    {
        int eid {}; in >> eid;
        Int2 vv; in >> vv[0] >> vv[1];
        if (vvs.count(vv)) continue;
        Vh vh0 = mesh.vertex_handle(vv[0] - 1);
        Vh vh1 = mesh.vertex_handle(vv[1] - 1);
        Hh hh0 = mesh.add_edge(vh0, vh1, true);
        // optional: process edge tag
        in.ignore(kInf, '\n');
        vvs.insert(vv);
    }

    return IO_ERR_TYPE::NO_ERROR;
}

template <class MeshT>
inline int save_poly_obj(const MeshT &mesh, const char *filename)
{
    std::vector<int> es;

    for (auto eh : mesh.edges())
    {
        es.push_back(mesh.from_vertex_handle(mesh.halfedge_handle(eh, 0)).idx() + 1);
        es.push_back(mesh.to_vertex_handle  (mesh.halfedge_handle(eh, 0)).idx() + 1);
    }

    return save_obj(
        mesh.points()->data(), (int)mesh.n_vertices(),
        nullptr, 0,
        es.data(), (int)mesh.n_edges(),
        filename);
}

template <class MeshT>
static int save_poly_detri2(const MeshT &mesh, const char *filename)
{
    std::ofstream out(filename, std::ios::out);
    if (!out) return IO_ERR_TYPE::CANNOT_OPEN;

    const_cast<MeshT&>(mesh).garbage_collection(true, true, false);

    out << std::fixed << std::setprecision(17);

    out << mesh.n_vertices() << " 3 0 1\n";

    if (!hasProperty<Vh, int>(mesh, var_v_label()))
    {
        for (auto vh : mesh.vertices())
        {
            const auto &p = mesh.point(vh);
            out << vh.idx() + 1 << "  "
                << p[0] << "  "
                << p[1] << "  "
                << p[2] << "  "
                << "-1\n"; // arbitrary tag
        }
    }
    else
    {
        auto v_tag = getProperty<Vh, int>(mesh, var_v_label());

        for (auto vh : mesh.vertices())
        {
            const auto &p = mesh.point(vh);
            out << vh.idx() + 1 << "  "
                << p[0] << "  "
                << p[1] << "  "
                << p[2] << "  "
                << v_tag[vh] + 1 << "\n";
        }
    }

    out << mesh.n_edges() << " 1\n";

    if (!hasProperty<Eh, int>(mesh, var_e_label()))
    {
        for (auto eh : mesh.edges())
        {
            Hh hh = mesh.halfedge_handle(eh, 0);
            out << eh.idx() + 1 << " "
                << mesh.from_vertex_handle(hh).idx() + 1 << " "
                << mesh.to_vertex_handle  (hh).idx() + 1 << " "
                << -1 << "\n";
        }
    }
    else
    {
        auto e_tag = getProperty<Eh, int>(mesh, var_e_label());

        for (auto eh : mesh.edges())
        {
            Hh hh = mesh.halfedge_handle(eh, 0);
            out << eh.idx() + 1 << " "
                << mesh.from_vertex_handle(hh).idx() + 1 << " "
                << mesh.to_vertex_handle  (hh).idx() + 1 << " "
                << e_tag[eh] + 1 << "\n";
        }
    }

    out << "0\n";

    return IO_ERR_TYPE::NO_ERROR;
}

template <class MeshT>
static int save_poly_mesh_detri2(const MeshT &mesh, const char *filename)
{
    std::ofstream out(filename, std::ios::out);
    if (!out) return IO_ERR_TYPE::CANNOT_OPEN;

    const_cast<MeshT&>(mesh).garbage_collection(true, true, false);

    out << std::fixed << std::setprecision(17);

    out << "MeshVersionFormatted 1\n\nDimension\n3\n\n";

    out << "Vertices\n" << mesh.n_vertices() << "\n";

    if (!hasProperty<Vh, int>(mesh, var_v_label()))
    {
        for (auto vert : mesh.vertices())
        {
            const auto &p = mesh.point(vert);
            out << p[0] << "  "
                << p[1] << "  "
                << p[2] << "  "
                << "0\n"; // arbitrary tag
        }
    }
    else
    {
        auto v_tag = getProperty<Vh, int>(mesh, var_v_label());

        for (auto vert : mesh.vertices())
        {
            const auto &p = mesh.point(vert);
            out << p[0] << "  "
                << p[1] << "  "
                << p[2] << "  "
                << v_tag[vert] + 1 << "\n";
        }
    }

    out << "Triangles\n" << 0 << "\n";

    out << "Edges\n" << mesh.n_edges() << "\n";

    if (!hasProperty<Eh, int>(mesh, var_e_label()))
    {
        for (auto eh : mesh.edges())
        {
            Hh hh = mesh.halfedge_handle(eh, 0);
            out << mesh.from_vertex_handle(hh).idx() + 1 << " "
                << mesh.to_vertex_handle  (hh).idx() + 1 << " "
                << "-1\n";
        }
    }
    else
    {
        auto e_tag = getProperty<Eh, int>(mesh, var_e_label());

        for (auto eh : mesh.edges())
        {
            Hh hh = mesh.halfedge_handle(eh, 0);
            out << mesh.from_vertex_handle(hh).idx() + 1 << " "
                << mesh.to_vertex_handle  (hh).idx() + 1 << " "
                << e_tag[eh] + 1 << "\n";
        }
    }

    out << "End\n";

    return IO_ERR_TYPE::NO_ERROR;
}

template <class MeshT>
int read_poly(MeshT &mesh, const char *filename)
{
    int err { IO_ERR_TYPE::NO_ERROR };

    const auto ext = get_file_extension(filename);

    if (ext.compare("poly") == 0)
    {
        err = read_poly_detri2(mesh, filename);
    }
    else
    {
        err = IO_ERR_TYPE::UNSUPPORTED_FORMAT;
    }

    if (err) fprintf(stderr, "read_mesh error %d: %s\n", err, __io_err_msg[err]);

    return err;
}

template <class MeshT>
int save_poly(const MeshT &mesh, const char *filename)
{
    int err { IO_ERR_TYPE::NO_ERROR };

    const auto ext = get_file_extension(filename);

    if (ext.compare("obj") == 0)
    {
        err = save_poly_obj(mesh, filename);
    }
    else if (ext.compare("poly") == 0)
    {
        err = save_poly_detri2(mesh, filename);
    }
    else if (ext.compare("mesh") == 0)
    {
        err = save_poly_mesh_detri2(mesh, filename);
    }
    else
    {
        err = IO_ERR_TYPE::UNSUPPORTED_FORMAT;
    }

    if (err) fprintf(stderr, "save_mesh error %d: %s\n", err, __io_err_msg[err]);

    return err;
}

template
int read_poly(LoopMesh&, const char*);

template
int save_poly(const PolyMesh&, const char*);

template
int save_poly(const TriMesh&, const char*);

template
int save_poly(const LoopMesh&, const char*);

template <class MeshT>
inline int save_marked_face_centroids(const MeshT &mesh, const char *filename, const double offset)
{
    std::vector<Vec3> ps {};
    int np {};

    for (const auto face : mesh.faces()) if (is_marked(mesh, face))
    {
        const auto p = mesh.calc_face_centroid(face);
        const auto n = mesh.calc_normal(face) * offset;
        ps.push_back(p + n);
        ++np;
    }

    return save_obj(
        (const double*)ps.data(), (const int)ps.size(),
        nullptr,                  0,
        nullptr,                  0,
        filename);
}

template
int save_marked_face_centroids(const TriMesh&, const char*, const double);

template
int save_marked_face_centroids(const PolyMesh&, const char*, const double);

template <class MeshT>
int save_marked_vertices(const MeshT &mesh, const char *filename, const double offset)
{
    std::vector<Vec3> ps {};
    int np {};

    for (const auto vert : mesh.vertices()) if (is_marked(mesh, vert))
    {
        const auto p = mesh.point(vert);
        const auto n = mesh.calc_normal(vert) * offset;
        ps.push_back(p + n);
        ++np;
    }

    return save_obj(
        (const double*)ps.data(), (const int)ps.size(),
        nullptr,                  0,
        nullptr,                  0,
        filename);
}

template
int save_marked_vertices(const TriMesh&, const char*, const double);

template
int save_marked_vertices(const PolyMesh&, const char*, const double);
