#pragma once

/**
 * @file JIArrayExpr.h
 * @brief Expression templates for JIArray element-wise arithmetic.
 *
 * The binary/unary/scalar arithmetic operators (`+ - * /`, unary `-`) return
 * lazy expression nodes instead of materialised arrays.  A whole expression
 * such as `d = a + b + c` is then evaluated in a **single** allocation-free
 * pass when it is assigned to (or used to construct) a JIArray — eliminating
 * the per-operator temporaries that made the eager operators a hot-path trap
 * (see docs/EVALUATION.md, F4).
 *
 * Design contract (mirrors Eigen/xtensor):
 *   - Nodes evaluate **linearly** over flat storage (`eval(i)`, i in [0,nn)),
 *     exactly like the previous eager operators.  This is layout-agnostic:
 *     operands share an identical shape, so column-/row-major need no special
 *     casing here.
 *   - Nodes are lightweight temporaries; array leaves capture the operand's
 *     data/shape **pointers** (valid for the lifetime of the full expression).
 *   - Because operators no longer return JIArray, **do not capture an
 *     expression with `auto`** — write the concrete array type
 *     (`zdouble1 r = a + b;`) so it materialises.  `auto r = a + b;` binds the
 *     node (holding operand references); using it after the operands die is a
 *     dangling read.
 *
 * Everything is `JIARRAY_HD` so the header compiles under NVCC on the host
 * side; the device container (JICudaArray) carries its own arithmetic.
 */

#include "pch.h"
#include <type_traits>
#include <utility>

namespace dnegri::jiarray {

// Forward declaration (JIArray.h includes this header before the class body).
template <class T, size_t RANK, class SEQ>
class JIArray;

template <typename T>
struct is_jiarray;

// ============================================================================
// Expression-node base + trait
// ============================================================================

/// CRTP base tagging every expression node; enables is_jiexpr detection.
template <class Derived>
struct JIArrayExprBase {
    JIARRAY_HD const Derived& derived() const noexcept {
        return static_cast<const Derived&>(*this);
    }
};

/// True when X is a JIArray expression node.
template <class X>
inline constexpr bool is_jiexpr_v =
    std::is_base_of_v<JIArrayExprBase<std::decay_t<X>>, std::decay_t<X>>;

/// True when X is a JIArray or an expression node (an "array-like" operand).
template <class X>
inline constexpr bool is_arraylike_v =
    is_jiarray<std::decay_t<X>>::value || is_jiexpr_v<X>;

/// True when X may participate as an operand (array-like or a raw scalar).
template <class X>
inline constexpr bool is_operand_v =
    is_arraylike_v<X> || std::is_arithmetic_v<std::decay_t<X>>;

/// Enable an operator only when at least one side is array-like and neither
/// side is something unrelated (so we never hijack `int+int`, string+string…).
template <class L, class R>
inline constexpr bool enable_expr_op =
    (is_arraylike_v<L> || is_arraylike_v<R>) && is_operand_v<L> && is_operand_v<R>;

// ============================================================================
// Leaf nodes
// ============================================================================

/**
 * @brief Expression leaf wrapping a JIArray operand (captures data + shape).
 */
template <class T, size_t RANK>
struct JIArrayLeaf : JIArrayExprBase<JIArrayLeaf<T, RANK>> {
    using value_type              = T;
    static constexpr size_t rank  = RANK;
    static constexpr bool   scalar = false;

    const T*   mm;
    int        nn;
    const int* rankSize_;
    const int* offset_;

    template <class SEQ>
    JIARRAY_HD explicit JIArrayLeaf(const JIArray<T, RANK, SEQ>& a) noexcept
        : mm(a.data()), nn(a.getSize()), rankSize_(a.getRankSize()), offset_(a.getOffset()) {}

    JIARRAY_HD T          eval(int i)   const noexcept { return mm[i]; }
    JIARRAY_HD int        size()        const noexcept { return nn; }
    JIARRAY_HD const int* rankSize()    const noexcept { return rankSize_; }
    JIARRAY_HD const int* offset()      const noexcept { return offset_; }
};

/**
 * @brief Expression leaf wrapping a broadcast scalar (no shape).
 */
template <class S>
struct JIArrayScalar : JIArrayExprBase<JIArrayScalar<S>> {
    using value_type              = S;
    static constexpr size_t rank  = 0;
    static constexpr bool   scalar = true;

    S v;
    JIARRAY_HD explicit JIArrayScalar(S v_) noexcept : v(v_) {}
    JIARRAY_HD S eval(int) const noexcept { return v; }
};

// ============================================================================
// Operation functors
// ============================================================================

struct JIExprAdd { template <class A, class B> JIARRAY_HD static auto apply(A a, B b) { return a + b; } };
struct JIExprSub { template <class A, class B> JIARRAY_HD static auto apply(A a, B b) { return a - b; } };
struct JIExprMul { template <class A, class B> JIARRAY_HD static auto apply(A a, B b) { return a * b; } };
struct JIExprDiv { template <class A, class B> JIARRAY_HD static auto apply(A a, B b) { return a / b; } };
struct JIExprNeg { template <class A>          JIARRAY_HD static auto apply(A a)      { return -a;    } };

// ============================================================================
// Composite nodes
// ============================================================================

/**
 * @brief Binary expression node. Shape is taken from the non-scalar operand.
 */
template <class Op, class L, class R>
struct JIArrayBinaryExpr : JIArrayExprBase<JIArrayBinaryExpr<Op, L, R>> {
    static constexpr bool   l_scalar = L::scalar;
    static constexpr size_t rank     = l_scalar ? R::rank : L::rank;
    using value_type = decltype(Op::apply(std::declval<typename L::value_type>(),
                                          std::declval<typename R::value_type>()));

    L l;
    R r;
    JIARRAY_HD JIArrayBinaryExpr(L l_, R r_) noexcept : l(l_), r(r_) {}

    JIARRAY_HD value_type eval(int i) const { return Op::apply(l.eval(i), r.eval(i)); }
    JIARRAY_HD int        size()      const {
        if constexpr (!l_scalar) return l.size(); else return r.size();
    }
    JIARRAY_HD const int* rankSize() const {
        if constexpr (!l_scalar) return l.rankSize(); else return r.rankSize();
    }
    JIARRAY_HD const int* offset() const {
        if constexpr (!l_scalar) return l.offset(); else return r.offset();
    }
    static constexpr bool scalar = false;
};

/**
 * @brief Unary expression node (currently negation). Shape from its operand.
 */
template <class Op, class E>
struct JIArrayUnaryExpr : JIArrayExprBase<JIArrayUnaryExpr<Op, E>> {
    static constexpr size_t rank = E::rank;
    using value_type             = decltype(Op::apply(std::declval<typename E::value_type>()));

    E e;
    JIARRAY_HD explicit JIArrayUnaryExpr(E e_) noexcept : e(e_) {}

    JIARRAY_HD value_type eval(int i) const { return Op::apply(e.eval(i)); }
    JIARRAY_HD int        size()      const { return e.size(); }
    JIARRAY_HD const int* rankSize()  const { return e.rankSize(); }
    JIARRAY_HD const int* offset()    const { return e.offset(); }
    static constexpr bool scalar = false;
};

// ============================================================================
// to_expr: wrap any operand into an expression node
// ============================================================================

template <class T, size_t RANK, class SEQ>
JIARRAY_HD inline JIArrayLeaf<T, RANK> to_expr(const JIArray<T, RANK, SEQ>& a) noexcept {
    return JIArrayLeaf<T, RANK>(a);
}

template <class E, std::enable_if_t<is_jiexpr_v<E>, int> = 0>
JIARRAY_HD inline E to_expr(const E& e) noexcept {
    return e;   // already an expression node — copy through (nodes are light)
}

template <class S, std::enable_if_t<std::is_arithmetic_v<std::decay_t<S>>, int> = 0>
JIARRAY_HD inline JIArrayScalar<std::decay_t<S>> to_expr(S s) noexcept {
    return JIArrayScalar<std::decay_t<S>>(s);
}

// ============================================================================
// Operators (lazy). Enabled only when at least one operand is array-like.
// ============================================================================

#define JIARRAY_DEFINE_BINOP(sym, opfunctor)                                                       \
    template <class L, class R, std::enable_if_t<enable_expr_op<L, R>, int> = 0>                    \
    JIARRAY_HD inline auto operator sym(const L& l, const R& r) {                                   \
        auto le = to_expr(l);                                                                      \
        auto re = to_expr(r);                                                                      \
        return JIArrayBinaryExpr<opfunctor, decltype(le), decltype(re)>(le, re);                   \
    }

JIARRAY_DEFINE_BINOP(+, JIExprAdd)
JIARRAY_DEFINE_BINOP(-, JIExprSub)
JIARRAY_DEFINE_BINOP(*, JIExprMul)
JIARRAY_DEFINE_BINOP(/, JIExprDiv)

#undef JIARRAY_DEFINE_BINOP

template <class E, std::enable_if_t<is_arraylike_v<E>, int> = 0>
JIARRAY_HD inline auto operator-(const E& e) {
    auto ee = to_expr(e);
    return JIArrayUnaryExpr<JIExprNeg, decltype(ee)>(ee);
}

} // namespace dnegri::jiarray
