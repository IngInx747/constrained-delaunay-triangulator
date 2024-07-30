#include "triangulator.hh"
#include "mesh_io.hh"

using namespace OpenMesh;

int main(const int argc, const char **argv)
{
    std::string filename, prefix, path;
    if (argc < 2) { printf("No file provided.\n"); return 1; }

    filename.append(argv[1]);
    prefix = filename.substr(0, filename.find_last_of("."));
    path = filename.substr(0, filename.find_last_of("/\\"));

    std::vector<Vec2> vs {};
    std::vector<Int2> es {};
    if (read_poly(vs, es, filename.c_str()) != 0)
    { printf("Cannot load poly file.\n"); return 1; }

    TriMesh mesh;
    int err = triangulate(vs, es, mesh);
    save_mesh(mesh, (prefix + ".out.mesh").c_str());

    return err;
}