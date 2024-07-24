#include "mesher.hh"
#include "triangulator.hh"
#include "delaunay.hh"
#include "mesh_io.hh"

using namespace OpenMesh;

int main(const int argc, const char **argv)
{
    std::string filename, prefix, path;
    if (argc < 2) { printf("No file provided.\n"); return 1; }

    filename.append(argv[1]);
    prefix = filename.substr(0, filename.find_last_of("."));
    path = filename.substr(0, filename.find_last_of("/\\"));

    LoopMesh poly;
    if (read_poly(poly, filename.c_str()) != 0)
    { printf("Cannot load polygon.\n"); return 1; }

    getOrMakeProperty<Mh, std::string>(poly, var_m_name())() = prefix;
    getOrMakeProperty<Mh, std::string>(poly, var_m_path())() = path;

    poly.delete_isolated_vertices();
    poly.garbage_collection();

    TriMesh mesh;
    getOrMakeProperty<Mh, std::string>(mesh, var_m_name())() = prefix;
    getOrMakeProperty<Mh, std::string>(mesh, var_m_path())() = path;

    triangulate(poly, mesh);
    make_delaunay(mesh);

    save_mesh(mesh, (prefix + ".out.mesh").c_str());

    return 0;
}