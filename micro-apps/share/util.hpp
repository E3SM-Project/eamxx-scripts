#ifndef INCLUDE_UTIL
#define INCLUDE_UTIL

#include "types.hpp"
#include "scream_assert.hpp"

#include <cstdio>
#include <cstring>
#include <memory>
#include <map>

#ifndef KOKKOS_ENABLE_CUDA
# include <cmath>
# include <algorithm>
#endif

/*
 * STL-like utilities including some CUDA-compatible replacements
 * for the STL.
 */

namespace util {

struct FILECloser { void operator() (FILE* fh) { fclose(fh); } };
using FILEPtr = std::unique_ptr<FILE, FILECloser>;

template<typename T>
void write (const T* v, size_t sz, const FILEPtr& fid) {
  size_t nwrite = fwrite(v, sizeof(T), sz, fid.get());
  micro_require_msg(nwrite == sz, "write: nwrite = " << nwrite << " sz = " << sz);
}

template<typename T>
void read (T* v, size_t sz, const FILEPtr& fid) {
  size_t nread = fread(v, sizeof(T), sz, fid.get());
  micro_require_msg(nread == sz, "read: nread = " << nread << " sz = " << sz);
}

bool eq(const std::string& a, const char* const b1, const char* const b2 = 0);

template <typename Scalar> struct is_single_precision {};
template <> struct is_single_precision<float> { enum : bool { value = true }; };
template <> struct is_single_precision<double> { enum : bool { value = false }; };

#ifdef KOKKOS_ENABLE_CUDA
// Replacements for namespace std functions that don't run on the GPU.
template <typename T> KOKKOS_INLINE_FUNCTION
const T& min (const T& a, const T& b) { return a < b ? a : b; }
template <typename T> KOKKOS_INLINE_FUNCTION
const T& max (const T& a, const T& b) { return a > b ? a : b; }
KOKKOS_INLINE_FUNCTION bool isfinite (const Real& a) {
  return a == a && a != INFINITY && a != -INFINITY;
}
template <typename T> KOKKOS_INLINE_FUNCTION
const T* max_element (const T* const begin, const T* const end) {
  const T* me = begin;
  for (const T* it = begin + 1; it < end; ++it)
    if ( ! (*it < *me)) // use operator<
      me = it;
  return me;
}
KOKKOS_INLINE_FUNCTION
size_t strlen(const char* str)
{
  micro_kassert(str != NULL);
  const char *char_ptr;
  for (char_ptr = str; ; ++char_ptr)  {
    if (*char_ptr == '\0') return char_ptr - str;
  }
}
KOKKOS_INLINE_FUNCTION
void strcpy(char* dst, const char* src)
{
  micro_kassert(dst != NULL && src != NULL);
  while(*dst++ = *src++);
}
KOKKOS_INLINE_FUNCTION
int strcmp(const char* first, const char* second)
{
  while(*first && (*first == *second))
  {
    first++;
    second++;
  }
  return *(const unsigned char*)first - *(const unsigned char*)second;
}
template<class T>
KOKKOS_INLINE_FUNCTION
const T* upper_bound(const T* first, const T* last, const T& value)
{
  const T* it;
  int count, step;
  count = last - first;

  while (count > 0) {
    it = first;
    step = count / 2;
    it += step;
    if (value >= *it) {
      first = ++it;
      count -= step + 1;
    }
    else {
      count = step;
    }
  }
  return first;
}
#else
using std::min;
using std::max;
using std::isfinite;
using std::max_element;
using std::strlen;
using std::strcpy;
using std::strcmp;
using std::upper_bound;
#endif

template <typename Integer> KOKKOS_INLINE_FUNCTION
void set_min_max (const Integer& lim0, const Integer& lim1,
                  Integer& min, Integer& max) {
  min = util::min(lim0, lim1);
  max = util::max(lim0, lim1);
}

template <typename Integer, typename Integer1> KOKKOS_INLINE_FUNCTION
void set_min_max (const Integer& lim0, const Integer& lim1,
                  Integer& min, Integer& max, const Integer1& vector_size) {
  min = util::min(lim0, lim1) / vector_size;
  max = util::max(lim0, lim1) / vector_size;
}

template <typename T> KOKKOS_INLINE_FUNCTION T reldif (const T& a, const T& b) {
  return std::abs((b - a)/a);
}

template <typename T>
const T* flatten(const vector_2d_t<T>& data)
{
  const size_t dim0 = data.size();
  const size_t dim1 = data[0].size();
  T* result = new T[dim0 * dim1];

  for (size_t i = 0; i < dim0; ++i) {
    std::copy(data[i].begin(), data[i].end(), result + dim1 * i);
  }

  return result;
}

} // namespace util

#endif
