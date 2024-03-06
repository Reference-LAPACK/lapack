#ifndef LAPACK_CPP_ENUMP_HPP
#define LAPACK_CPP_ENUMP_HPP

#include <cassert>

namespace lapack_cpp {

enum class Layout : char { RowMajor = 'R', ColMajor = 'C' };

enum class Op : char { NoTrans = 'N', Trans = 'T', ConjTrans = 'C' };

inline constexpr Op char2op(char t)
{
    switch (t) {
        case 'N':
            return Op::NoTrans;
        case 'T':
            return Op::Trans;
        case 'C':
            return Op::ConjTrans;
        default:
            assert(false);
            return Op::NoTrans;
    }
}

enum class Side : char { Left = 'L', Right = 'R' };

inline constexpr Side char2side(char t)
{
    switch (t) {
        case 'L':
            return Side::Left;
        case 'R':
            return Side::Right;
        default:
            assert(false);
            return Side::Left;
    }
}

enum class Uplo : char { Upper = 'U', Lower = 'L' };

inline constexpr Uplo char2uplo(char t)
{
    switch (t) {
        case 'U':
            return Uplo::Upper;
        case 'R':
            return Uplo::Lower;
        default:
            assert(false);
            return Uplo::Upper;
    }
}

enum class Diag : char { Unit = 'U', NonUnit = 'N' };

inline constexpr Diag char2diag(char t)
{
    switch (t) {
        case 'U':
            return Diag::Unit;
        case 'N':
            return Diag::NonUnit;
        default:
            assert(false);
            return Diag::NonUnit;
    }
}

enum class Direction : char { Forward = 'F', Backward = 'B' };

inline constexpr Direction char2direct(char t)
{
    switch (t) {
        case 'F':
            return Direction::Forward;
        case 'N':
            return Direction::Backward;
        default:
            assert(false);
            return Direction::Forward;
    }
}

enum class IncDec : char { Increasing = 'I', Decreasing = 'D' };

inline constexpr IncDec char2incdec(char t)
{
    switch (t) {
        case 'I':
            return IncDec::Increasing;
        case 'D':
            return IncDec::Decreasing;
        default:
            assert(false);
            return IncDec::Increasing;
    }
}

enum class CompQ : char { No = 'N', Update = 'V', Initialize = 'I' };

inline constexpr CompQ char2compq(char t)
{
    switch (t) {
        case 'N':
            return CompQ::No;
        case 'V':
            return CompQ::Update;
        case 'I':
            return CompQ::Initialize;
        default:
            assert(false);
            return CompQ::No;
    }
}

}  // namespace lapack_cpp

#endif  // LAPACK_CPP_ENUMP_HPP