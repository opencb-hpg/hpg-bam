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


/*! \file gather.inl
 *  \brief Inline file for gather.h.
 */

#include <thrust/gather.h>
#include <thrust/iterator/iterator_traits.h>
#include <thrust/system/detail/generic/select_system.h>
#include <thrust/system/detail/generic/gather.h>
#include <thrust/system/detail/adl/gather.h>

namespace thrust
{

template<typename InputIterator,
         typename RandomAccessIterator,
         typename OutputIterator>
  OutputIterator gather(InputIterator        map_first,
                        InputIterator        map_last,
                        RandomAccessIterator input_first,
                        OutputIterator       result)
{
  using thrust::system::detail::generic::select_system;
  using thrust::system::detail::generic::gather;

  typedef typename thrust::iterator_system<InputIterator>::type        system1; 
  typedef typename thrust::iterator_system<RandomAccessIterator>::type system2; 
  typedef typename thrust::iterator_system<OutputIterator>::type       system3; 

  return gather(select_system(system1(),system2(),system3()), map_first, map_last, input_first, result);
} // end gather()


template<typename InputIterator1,
         typename InputIterator2,
         typename RandomAccessIterator,
         typename OutputIterator>
  OutputIterator gather_if(InputIterator1       map_first,
                           InputIterator1       map_last,
                           InputIterator2       stencil,
                           RandomAccessIterator input_first,
                           OutputIterator       result)
{
  using thrust::system::detail::generic::select_system;
  using thrust::system::detail::generic::gather_if;

  typedef typename thrust::iterator_system<InputIterator1>::type       system1; 
  typedef typename thrust::iterator_system<InputIterator2>::type       system2; 
  typedef typename thrust::iterator_system<RandomAccessIterator>::type system3; 
  typedef typename thrust::iterator_system<OutputIterator>::type       system4; 

  return gather_if(select_system(system1(),system2(),system3(),system4()), map_first, map_last, stencil, input_first, result);
} // end gather_if()


template<typename InputIterator1,
         typename InputIterator2,
         typename RandomAccessIterator,
         typename OutputIterator,
         typename Predicate>
  OutputIterator gather_if(InputIterator1       map_first,
                           InputIterator1       map_last,
                           InputIterator2       stencil,
                           RandomAccessIterator input_first,
                           OutputIterator       result,
                           Predicate            pred)
{
  using thrust::system::detail::generic::select_system;
  using thrust::system::detail::generic::gather_if;

  typedef typename thrust::iterator_system<InputIterator1>::type       system1; 
  typedef typename thrust::iterator_system<InputIterator2>::type       system2; 
  typedef typename thrust::iterator_system<RandomAccessIterator>::type system3; 
  typedef typename thrust::iterator_system<OutputIterator>::type       system4; 

  return gather_if(select_system(system1(),system2(),system3(),system4()), map_first, map_last, stencil, input_first, result, pred);
} // end gather_if()

} // end namespace thrust

