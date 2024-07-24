#ifndef OPENMESH_GENERIC_ITERATOR_HH
#define OPENMESH_GENERIC_ITERATOR_HH

#include <cassert>
#include <cstddef>
#include <iterator>
#include <OpenMesh/Core/Mesh/Status.hh>
#include <OpenMesh/Core/Mesh/SmartRange.hh>

namespace OpenMesh {

namespace Iterator2 {

template <class Mesh> class ConstVertexIterT;
template <class Mesh> class VertexIterT;
template <class Mesh> class ConstHalfedgeIterT;
template <class Mesh> class HalfedgeIterT;
template <class Mesh> class ConstEdgeIterT;
template <class Mesh> class EdgeIterT;

template <class Mesh, class ValueHandle, class MemberOwner, bool (MemberOwner::*PrimitiveStatusMember)() const, size_t (MemberOwner::*PrimitiveCountMember)() const>
class GenericIteratorT
{
public:

    typedef ValueHandle                     value_handle;
    typedef value_handle                    value_type;
    typedef std::bidirectional_iterator_tag iterator_category;
    typedef std::ptrdiff_t                  difference_type;
    typedef const Mesh*                     mesh_ptr;
    typedef const Mesh&                     mesh_ref;
    typedef const value_handle&             reference;
    typedef const value_handle*             pointer;

public:

    /// Default constructor.
    GenericIteratorT() : m(nullptr), h(), skip_bits(0)
    {}

    /// Construct with mesh and a target handle.
    GenericIteratorT(mesh_ref m, value_handle h, bool skip = false) : m(&m), h(h), skip_bits(0)
    { if (skip) { enable_skipping(); } }

    /// Standard dereferencing operator.
    inline reference operator*() const
    { return h; }

    /// Standard pointer operator.
    inline pointer operator->() const
    { return &h; }

    /// Are two iterators equal? Only valid if they refer to the same mesh!
    inline bool operator==(const GenericIteratorT& rhs) const
    { return ((m == rhs.m) && (h == rhs.h)); }

    /// Not equal?
    inline bool operator!=(const GenericIteratorT& rhs) const
    { return !operator==(rhs); }

    /// Standard pre-increment operator
    inline GenericIteratorT& operator++()
    {
        h.__increment();
        if (skip_bits)
            skip_fwd();
        return *this;
    }

    /// Standard post-increment operator
    inline GenericIteratorT operator++(int)
    {
        GenericIteratorT cpy(*this);
        ++(*this);
        return cpy;
    }

    /// Standard pre-decrement operator
    inline GenericIteratorT& operator--()
    {
        h.__decrement();
        if (skip_bits)
            skip_bwd();
        return *this;
    }

    /// Standard post-decrement operator
    inline GenericIteratorT operator--(int)
    {
        GenericIteratorT cpy(*this);
        --(*this);
        return cpy;
    }

    /// Turn on skipping: automatically skip deleted/hidden elements
    inline void enable_skipping()
    {
        if (m && (m->*PrimitiveStatusMember)())
        {
            Attributes::StatusInfo status;
            status.set_deleted(true);
            status.set_hidden(true);
            skip_bits = status.bits();
            skip_fwd();
        }
        else
        {
            skip_bits = 0;
        }
    }

    /// Turn on skipping: automatically skip deleted/hidden elements
    inline void disable_skipping()
    {
        skip_bits = 0;
    }

protected:

    void skip_fwd()
    {
        assert(m && skip_bits);
        while ((h.idx() < (signed) (m->*PrimitiveCountMember)()) && (m->status(h).bits() & skip_bits))
            h.__increment();
    }

    void skip_bwd()
    {
        assert(m && skip_bits);
        while ((h.idx() >= 0) && (m->status(h).bits() & skip_bits))
            h.__decrement();
    }

protected:

    mesh_ptr m;
    value_handle h;
    unsigned int skip_bits;
};

} // namespace Iterator2

////////////////////////////////////////////////////////////////

namespace Iterator2 {

template <typename CONTAINER_T, typename ITER_T, ITER_T (CONTAINER_T::*begin_fn)() const, ITER_T (CONTAINER_T::*end_fn)() const>
struct RangeTraitT
{
    using CONTAINER_TYPE = CONTAINER_T;
    using ITER_TYPE = ITER_T;
    static ITER_TYPE begin(const CONTAINER_TYPE& _container) { return (_container.*begin_fn)(); }
    static ITER_TYPE end  (const CONTAINER_TYPE& _container) { return (_container.*end_fn)  (); }
};

template <typename RangeTraitT>
class EntityRange : public SmartRangeT<EntityRange<RangeTraitT>, typename RangeTraitT::ITER_TYPE::value_handle>
{
public:

    typedef typename RangeTraitT::ITER_TYPE iterator;
    typedef typename RangeTraitT::ITER_TYPE const_iterator;

    explicit EntityRange(typename RangeTraitT::CONTAINER_TYPE &container) : container(container) {}
    typename RangeTraitT::ITER_TYPE begin() const { return RangeTraitT::begin(container); }
    typename RangeTraitT::ITER_TYPE end  () const { return RangeTraitT::end  (container); }

private:

    typename RangeTraitT::CONTAINER_TYPE &container;
};

} // namespace Iterator2

} // namespace OpenMesh

#endif