#ifndef MESH_IO_HH
#define MESH_IO_HH

#include "mesh.hh"

template <class MeshT>
int read_mesh(MeshT&, const char*);

template <class MeshT>
int save_mesh(const MeshT&, const char*);

//template <class MeshT>
//int read_poly(MeshT&, const char*);

template <class MeshT>
int save_poly(const MeshT&, const char*);

int read_node(
    std::vector<Vec2> &vs,
    const char *filename);

int read_poly(
    std::vector<Vec2> &vs,
    std::vector<Int2> &es,
    const char *filename);

int save_obj(
    const double *vs, const int nv,
    const int    *fs, const int nf,
    const int    *es, const int ne,
    const char   *filename,
    const std::streamsize prec = 17i64);

int save_node(
    const double* vs, const int nv,
    const char* filename,
    const std::streamsize prec = 17i64);

int save_poly(
    const double *vs, const int nv,
    const int    *es, const int ne,
    const char* filename,
    const std::streamsize prec = 17i64);

#endif