/**
 * @file
 * @brief More verbose assertions than STD's assert function
 * 
 * Implements ASSERT(|_EQ|_NE|_LS|_LE|_GT|_GE).
 * In case an assertion fails, an error message including the code location
 * and parameters is provided and the program is terminated using an abort() call.
 * 
 * Even for production code, the asserts are not removed.
 * 
 * @author Manuel Penschuck
 * @copyright
 * Copyright (C) 2017 Manuel Penschuck
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * @copyright
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * @copyright
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once

//FIX-ME: Read assertion blog again!

/**
 * @ingroup utils
 * @defgroup utils-assert Assertions
 * 
 * @addtogroup utils-assert
 * @{
 */

#ifndef NDEBUG
    #define ASSERT_STREAM std::cout
    #define ASSERT_FUNC __PRETTY_FUNCTION__
    #include <iostream>
    #include <cstdlib>

    #define ASSERTIONS

    //! Requires condition X to be satisfied.
    #define ASSERT(X) {\
        if (!(X)) {\
            ASSERT_STREAM <<  "Assertion (" << #X << ") at " << __FILE__ << ":" << __LINE__ << "\n   in " << ASSERT_FUNC << " failed" << std::endl;\
            abort();\
        }}

    #define ASSERT_BINOP(OP1, OP, OP2) {\
        {auto x1 = (OP1); auto x2 = (OP2);\
        if (!(x1 OP x2)) {\
            ASSERT_STREAM <<  "Assertion (" << #OP1 << " " << #OP << " " << #OP2  << ") at " << __FILE__ << ":" << __LINE__ << "\n   in " << ASSERT_FUNC << " failed.\n   Actuals [" << x1 << " " << #OP << " " << x2 << "]" << std::endl;\
            abort();\
        }}}

    /**
    * @brief Report an error with code location and terminate execution using abort().
    * 
    * @code
    * ASSERT_FAIL("Description of error with some parameter X=" << X << " and Y=" << Y);
    * @endcode
    */
    #define ASSERT_FAIL(X) {\
        ASSERT_STREAM <<  "Error at " << __FILE__ << ":" << __LINE__ << "\n   in " << ASSERT_FUNC << " failed: " << std::endl << X << std::endl; \
        abort();\
    }
#else
    #define ASSERT(X) {}
    #define ASSERT_BINOP(OP1, OP2, OP) {}
#endif    

/**
 * Helper function to mark intentionally unused variables in order to prevent compiler warnings
 */
template <typename U>
inline void ASSERT_UNUSED(const U&) {}

#define ASSERT_EQ(OP1, OP2) ASSERT_BINOP((OP1), ==, (OP2)) //!< Requires OP1 == OP2
#define ASSERT_NE(OP1, OP2) ASSERT_BINOP((OP1), !=, (OP2)) //!< Requires OP1 != OP2
#define ASSERT_LS(OP1, OP2) ASSERT_BINOP((OP1), < , (OP2)) //!< Requires OP1 <  OP2
#define ASSERT_LE(OP1, OP2) ASSERT_BINOP((OP1), <=, (OP2)) //!< Requires OP1 <= OP2
#define ASSERT_GT(OP1, OP2) ASSERT_BINOP((OP1), > , (OP2)) //!< Requires OP1 >  OP2
#define ASSERT_GE(OP1, OP2) ASSERT_BINOP((OP1), >=, (OP2)) //!< Requires OP1 >= OP2

/** @} @} */
