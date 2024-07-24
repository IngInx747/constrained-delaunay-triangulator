#ifndef VECTOR_N_HH
#define VECTOR_N_HH

#include <OpenMesh/Core/Geometry/VectorT.hh>

template <typename T, size_t N>
using VecN = OpenMesh::VectorT<T, N>;

using Vec2 = VecN<double, 2>;
using Vec3 = VecN<double, 3>;
using Vec4 = VecN<double, 4>;
using Vec5 = VecN<double, 5>;
using Vec6 = VecN<double, 6>;
using Vec7 = VecN<double, 7>;
using Vec8 = VecN<double, 8>;

using Int2 = VecN<int, 2>;
using Int3 = VecN<int, 3>;
using Int4 = VecN<int, 4>;
using Int5 = VecN<int, 5>;
using Int6 = VecN<int, 6>;
using Int7 = VecN<int, 7>;
using Int8 = VecN<int, 8>;

#endif