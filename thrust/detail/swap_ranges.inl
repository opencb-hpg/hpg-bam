/*
 *  Copyright 2008-2012 NVIDIA Corporation
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */


/*! \file swap_ranges.inl
 *  \brief Inline file for swap_ranges.h.
 */

#include <thrust/swap.h>
#include <thrust/iterator/iterator_traits.h>
#include <thrust/system/detail/generic/select_system.h>
#include <thrust/system/detail/generic/swap_ranges.h>
#include <thrust/system/detail/adl/swap_ranges.h>

namespace thrust
{

template<typename ForwardIterator1,
         typename ForwardIterator2>
  ForwardIterator2 swap_ranges(ForwardIterator1 first1,
                               ForwardIterator1 last1,
                               ForwardIterator2 first2)
{
  using thrust::system::detail::generic::select_system;
  using thrust::system::detail::generic::swap_ranges;

  typedef typename thrust::iterator_system<ForwardIterator1>::type system1;
  typedef typename thrust::iterator_system<ForwardIterator2>::type system2;

  return swap_ranges(select_system(system1(),system2()), first1, last1, first2);
} // end swap_ranges()

} // end namespace thrust
