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


/*! \file sequence.inl
 *  \brief Inline file for sequence.h.
 */

#include <thrust/detail/config.h>
#include <thrust/sequence.h>
#include <thrust/iterator/iterator_traits.h>
#include <thrust/system/detail/generic/select_system.h>
#include <thrust/system/detail/generic/sequence.h>
#include <thrust/system/detail/adl/sequence.h>

namespace thrust
{


template<typename ForwardIterator>
  void sequence(ForwardIterator first,
                ForwardIterator last)
{
  using thrust::system::detail::generic::select_system;
  using thrust::system::detail::generic::sequence;

  typedef typename thrust::iterator_system<ForwardIterator>::type system;

  return sequence(select_system(system()), first, last);
} // end sequence()


template<typename ForwardIterator, typename T>
  void sequence(ForwardIterator first,
                ForwardIterator last,
                T init)
{
  using thrust::system::detail::generic::select_system;
  using thrust::system::detail::generic::sequence;

  typedef typename thrust::iterator_system<ForwardIterator>::type system;

  return sequence(select_system(system()), first, last, init);
} // end sequence()


template<typename ForwardIterator, typename T>
  void sequence(ForwardIterator first,
                ForwardIterator last,
                T init,
                T step)
{
  using thrust::system::detail::generic::select_system;
  using thrust::system::detail::generic::sequence;

  typedef typename thrust::iterator_system<ForwardIterator>::type system;

  return sequence(select_system(system()), first, last, init, step);
} // end sequence()


} // end namespace thrust
