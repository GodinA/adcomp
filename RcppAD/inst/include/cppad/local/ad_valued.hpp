/* $Id: ad_valued.hpp 2506 2012-10-24 19:36:49Z bradbell $ */
# ifndef CPPAD_AD_VALUED_INCLUDED
# define CPPAD_AD_VALUED_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    GNU General Public License Version 3.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

/*
$begin ADValued$$
$spell
$$

$index operation, AD valued$$
$index function, AD valued$$

$section AD Valued Operations and Functions$$

$childtable%
	cppad/local/arithmetic.hpp%
	cppad/local/std_math_ad.hpp%
	cppad/local/math_other.hpp%
	cppad/local/cond_exp.hpp%
	cppad/local/discrete.hpp%
	cppad/local/user_atomic.hpp
%$$

$end
*/

// include MathOther.h after CondExp.h because some MathOther.h routines use 
// CondExp.h and CondExp.h is not sufficently declared in Declare.h

# include <cppad/local/arithmetic.hpp>
# include <cppad/local/std_math_ad.hpp>
# include <cppad/local/cond_exp.hpp>
# include <cppad/local/math_other.hpp>
# include <cppad/local/discrete.hpp>
# include <cppad/local/user_atomic.hpp>

# endif
