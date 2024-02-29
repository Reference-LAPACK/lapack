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
    }
}

}  // namespace lapack_cpp

#endif  // LAPACK_CPP_ENUMP_HPP