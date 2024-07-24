#ifndef OPENMESH_LOOP_CIRCULATOR_HH
#define OPENMESH_LOOP_CIRCULATOR_HH

#include <cassert>
#include <cstddef>
#include <iterator>
#include <OpenMesh/Core/Mesh/SmartRange.hh>

namespace OpenMesh {

namespace Circulator3 {

template <class Mesh, class CenterEntityHandle>
struct GenericCirculator_CenterEntityFnsT
{
    static void increment(const Mesh *mesh, typename Mesh::HalfedgeHandle &hh, const typename Mesh::HalfedgeHandle &hh0, int &lap_counter);
    static void decrement(const Mesh *mesh, typename Mesh::HalfedgeHandle &hh, const typename Mesh::HalfedgeHandle &hh0, int &lap_counter);
};

template <class Mesh>
struct GenericCirculator_CenterEntityFnsT<Mesh, typename Mesh::VertexHandle>
{
    inline static void increment(const Mesh *mesh, typename Mesh::HalfedgeHandle &hh, const typename Mesh::HalfedgeHandle &hh0, int &lap_counter)
    {
        hh = mesh->cw_rotated_halfedge_handle(hh);
        if (hh == hh0) ++lap_counter;
    }

    inline static void decrement(const Mesh *mesh, typename Mesh::HalfedgeHandle &hh, const typename Mesh::HalfedgeHandle &hh0, int &lap_counter)
    {
        if (hh == hh0) --lap_counter;
        hh = mesh->ccw_rotated_halfedge_handle(hh);
    }
};

template <class Mesh>
struct GenericCirculator_CenterEntityFnsT<Mesh, typename Mesh::EdgeHandle>
{
    inline static void increment(const Mesh *mesh, typename Mesh::HalfedgeHandle &hh, const typename Mesh::HalfedgeHandle &hh0, int &lap_counter)
    {
        hh = mesh->opposite_halfedge_handle(hh);
        if (hh == hh0) ++lap_counter;
    }

    inline static void decrement(const Mesh *mesh, typename Mesh::HalfedgeHandle &hh, const typename Mesh::HalfedgeHandle &hh0, int &lap_counter)
    {
        if (hh == hh0) --lap_counter;
        hh = mesh->opposite_halfedge_handle(hh);
    }
};

template<class Mesh>
struct GenericCirculator_CenterEntityFnsT<Mesh, typename Mesh::FaceHandle>
{
    inline static void increment(const Mesh *mesh, typename Mesh::HalfedgeHandle &hh, const typename Mesh::HalfedgeHandle &hh0, int &lap_counter)
    {
        hh = mesh->next_halfedge_handle(hh);
        if (hh == hh0) ++lap_counter;
    }

    inline static void decrement(const Mesh *mesh, typename Mesh::HalfedgeHandle &hh, const typename Mesh::HalfedgeHandle &hh0, int &lap_counter)
    {
        if (hh == hh0) --lap_counter;
        hh = mesh->prev_halfedge_handle(hh);
    }
};

////////////////////////////////////////////////////////////////

template <class Mesh, class CenterEntityHandle, class ValueHandle>
struct GenericCirculator_ValueHandleFnsT
{
    inline static bool is_valid(const typename Mesh::HalfedgeHandle &hh, const int lap_counter)
    { return hh.is_valid() && (lap_counter == 0); }

    inline static void init(const Mesh*, typename Mesh::HalfedgeHandle&, typename Mesh::HalfedgeHandle&, int&)
    {};

    inline static void increment(const Mesh *mesh, typename Mesh::HalfedgeHandle &hh, const typename Mesh::HalfedgeHandle &hh0, int &lap_counter)
    { GenericCirculator_CenterEntityFnsT<Mesh, CenterEntityHandle>::increment(mesh, hh, hh0, lap_counter); }

    inline static void decrement(const Mesh *mesh, typename Mesh::HalfedgeHandle &hh, const typename Mesh::HalfedgeHandle &hh0, int &lap_counter)
    { GenericCirculator_CenterEntityFnsT<Mesh, CenterEntityHandle>::decrement(mesh, hh, hh0, lap_counter); }
};

////////////////////////////////////////////////////////////////

template<class Mesh>
class GenericCirculatorBaseT
{
public:

    typedef const Mesh *mesh_ptr;
    typedef const Mesh &mesh_ref;

    template <typename> friend class CirculatorRange;

public:

    GenericCirculatorBaseT() : mesh(nullptr), lap_counter(0) {}

    GenericCirculatorBaseT(mesh_ref mesh, typename Mesh::HalfedgeHandle hh, bool end = false) :
    mesh(&mesh), hh(hh), hh0(hh), lap_counter(static_cast<int>(end && hh.is_valid())) {}

    GenericCirculatorBaseT(const GenericCirculatorBaseT &rhs) :
    mesh(rhs.mesh), hh(rhs.hh), hh0(rhs.hh0), lap_counter(rhs.lap_counter) {}

    inline typename Mesh::HalfedgeHandle toHalfedgeHandle() const
    { return hh; }

    inline typename Mesh::HalfedgeHandle toOppositeHalfedgeHandle() const
    { return mesh->opposite_halfedge_handle(hh); }

    inline typename Mesh::EdgeHandle toEdgeHandle() const
    { return mesh->edge_handle(hh); }

    inline typename Mesh::VertexHandle toVertexHandle() const
    { return mesh->to_vertex_handle(hh); }

    inline typename Mesh::FaceHandle toFaceHandle() const
    { return mesh->face_handle(hh); }

    inline GenericCirculatorBaseT &operator=(const GenericCirculatorBaseT &rhs)
    {
        mesh = rhs.mesh;
        hh0 = rhs.hh0;
        hh = rhs.hh;
        lap_counter = rhs.lap_counter;
        return *this;
    }

    inline bool operator==(const GenericCirculatorBaseT &rhs) const
    { return mesh == rhs.mesh && hh0 == rhs.hh0 && hh == rhs.hh && lap_counter == rhs.lap_counter; }

    inline bool operator!=(const GenericCirculatorBaseT &rhs) const
    { return !operator==(rhs); }

protected:

    mesh_ptr mesh;

    typename Mesh::HalfedgeHandle hh; // current handle

    typename Mesh::HalfedgeHandle hh0; // starting handle

    int lap_counter;
};

////////////////////////////////////////////////////////////////

template <typename GenericCirculatorT_TraitsT>
class GenericCirculatorT : protected GenericCirculatorBaseT<typename GenericCirculatorT_TraitsT::Mesh>
{
public:

    using Mesh = typename GenericCirculatorT_TraitsT::Mesh;
    using value_type = typename GenericCirculatorT_TraitsT::ValueHandle;
    using CenterEntityHandle = typename GenericCirculatorT_TraitsT::CenterEntityHandle;

    typedef std::ptrdiff_t difference_type;
    typedef const value_type& reference;
    typedef const value_type* pointer;
    typedef std::bidirectional_iterator_tag iterator_category;

    typedef typename GenericCirculatorBaseT<Mesh>::mesh_ptr mesh_ptr;
    typedef typename GenericCirculatorBaseT<Mesh>::mesh_ref mesh_ref;
    typedef GenericCirculator_ValueHandleFnsT<Mesh, CenterEntityHandle, value_type> GenericCirculator_ValueHandleFns;

    template <typename> friend class CirculatorRange;

public:

    GenericCirculatorT() {}

    GenericCirculatorT(mesh_ref mesh, CenterEntityHandle xh0, bool end = false) : GenericCirculatorBaseT<Mesh>(mesh, mesh.halfedge_handle(xh0), end)
    { GenericCirculator_ValueHandleFns::init(this->mesh, this->hh, this->hh0, this->lap_counter); }

    GenericCirculatorT(mesh_ref mesh, typename Mesh::HalfedgeHandle hh, bool end = false) : GenericCirculatorBaseT<Mesh>(mesh, hh, end)
    { GenericCirculator_ValueHandleFns::init(this->mesh, this->hh, this->hh0, this->lap_counter); }

    GenericCirculatorT(const GenericCirculatorT &rhs) : GenericCirculatorBaseT<Mesh>(rhs) {}

    inline GenericCirculatorT& operator++()
    {
        assert(this->mesh);
        GenericCirculator_ValueHandleFns::increment(this->mesh, this->hh, this->hh0, this->lap_counter);
        return *this;
    }

    inline GenericCirculatorT& operator--()
    {
        assert(this->mesh_);
        GenericCirculator_ValueHandleFns::decrement(this->mesh, this->hh, this->hh0, this->lap_counter);
        return *this;
    }

    /// Post-increment
    inline GenericCirculatorT operator++(int)
    {
        assert(this->mesh);
        GenericCirculatorT cpy(*this);
        ++(*this);
        return cpy;
    }

    /// Post-decrement
    inline GenericCirculatorT operator--(int)
    {
        assert(this->mesh);
        GenericCirculatorT cpy(*this);
        --(*this);
        return cpy;
    }

    /// Standard dereferencing operator.
    inline value_type operator*() const
    {
#ifndef NDEBUG
        assert(this->hh.is_valid());
        value_type res = GenericCirculatorT_TraitsT::toHandle(this->mesh, this->hh);
        assert(res.is_valid());
        return res;
#else
        return GenericCirculatorT_TraitsT::toHandle(this->mesh, this->hh);
#endif
    }

    inline GenericCirculatorT &operator=(const GenericCirculatorT &rhs)
    {
        GenericCirculatorBaseT<Mesh>::operator=(rhs);
        return *this;
    };

    inline bool operator==(const GenericCirculatorT &rhs) const
    { return GenericCirculatorBaseT<Mesh>::operator==(rhs); }

    inline bool operator!=(const GenericCirculatorT &rhs) const
    { return GenericCirculatorBaseT<Mesh>::operator!=(rhs); }

    inline bool is_valid() const
    { return GenericCirculator_ValueHandleFns::is_valid(this->hh, this->lap_counter); }
};

} // namespace Circulator3

////////////////////////////////////////////////////////////////

namespace Circulator3 {

template <
    typename CONTAINER_T,
    typename ITER_T,
    typename CENTER_ENTITY_T,
    typename TO_ENTITY_T,
    ITER_T (CONTAINER_T::*begin_fn)(CENTER_ENTITY_T) const,
    ITER_T (CONTAINER_T::*end_fn)  (CENTER_ENTITY_T) const>
struct CirculatorRangeTraitT
{
    using CONTAINER_TYPE = CONTAINER_T;
    using ITER_TYPE = ITER_T;
    using CENTER_ENTITY_TYPE = CENTER_ENTITY_T;
    using TO_ENTITYE_TYPE = TO_ENTITY_T;

    static ITER_TYPE begin(const CONTAINER_TYPE& container, CENTER_ENTITY_TYPE xh) { return (container.*begin_fn)(xh); }
    static ITER_TYPE end  (const CONTAINER_TYPE& container, CENTER_ENTITY_TYPE xh) { return (container.*end_fn)  (xh); }

    static ITER_TYPE begin(const CONTAINER_TYPE& container, HalfedgeHandle hh, int) { return ITER_TYPE(container, hh); }
    static ITER_TYPE end  (const CONTAINER_TYPE& container, HalfedgeHandle hh, int) { return ITER_TYPE(container, hh, true); }
};

template <typename CirculatorRangeTraitT>
class CirculatorRange : public SmartRangeT<CirculatorRange<CirculatorRangeTraitT>, typename CirculatorRangeTraitT::TO_ENTITYE_TYPE>
{
public:

    typedef typename CirculatorRangeTraitT::ITER_TYPE ITER_TYPE;
    typedef typename CirculatorRangeTraitT::CENTER_ENTITY_TYPE CENTER_ENTITY_TYPE;
    typedef typename CirculatorRangeTraitT::CONTAINER_TYPE CONTAINER_TYPE;
    typedef ITER_TYPE iterator;
    typedef ITER_TYPE const_iterator;

    CirculatorRange(const CONTAINER_TYPE &container, CENTER_ENTITY_TYPE center) : container(container), hh()
    {
        auto it = CirculatorRangeTraitT::begin(container, center);
        hh = it.hh;
    }

    CirculatorRange(const CONTAINER_TYPE &container, HalfedgeHandle hh, int) : container(container), hh(hh) {}

    ITER_TYPE begin() const { return CirculatorRangeTraitT::begin(container, hh, 1); }
    ITER_TYPE end()   const { return CirculatorRangeTraitT::end  (container, hh, 1); }

private:

    const CONTAINER_TYPE &container;
    HalfedgeHandle hh;
};

} // namespace Circulator3

} // namespace OpenMesh

#endif