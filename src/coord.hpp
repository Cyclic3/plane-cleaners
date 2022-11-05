#pragma once

#include <array>
#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <ostream>

template<typename F, F ptr>
struct c_deleter {
  inline auto operator()(auto* obj) { return ptr(obj); }
};
#define C_DELETER(F) c_deleter<decltype(&F), &F>

/// A simple class that handles operations on coordinates
template<size_t dims, typename CoordT = double>
class coord_t {
  static_assert(dims >= 1);
private:
  std::array<CoordT, dims> _arr;

public:
  constexpr CoordT& operator[](size_t idx) noexcept { return _arr[idx]; }
  constexpr CoordT const& operator[](size_t idx) const noexcept { return _arr[idx]; }
  constexpr coord_t operator+(coord_t const& other) const noexcept {
    coord_t ret;
    std::transform(_arr.begin(), _arr.end(), other.begin(), ret._arr.begin(), std::plus<>{});
    return ret;
  }
  constexpr coord_t& operator+=(coord_t const& other) noexcept {
    std::transform(_arr.begin(), _arr.end(), other.begin(), _arr.begin(), std::plus<>{});
    return *this;
  }
  constexpr coord_t operator-(coord_t const& other) const noexcept {
    coord_t ret;
    std::transform(_arr.begin(), _arr.end(), other.begin(), ret._arr.begin(), std::minus<>{});
    return ret;
  }
  constexpr coord_t& operator-=(coord_t const& other) noexcept {
    std::transform(_arr.begin(), _arr.end(), other.begin(), _arr.begin(), std::minus<>{});
    return *this;
  }
  constexpr coord_t operator/(CoordT const& scalar) const noexcept {
    coord_t ret;
    std::transform(_arr.begin(), _arr.end(), ret._arr.begin(), [&scalar](auto& i) { return i / scalar; });
    return ret;
  }
  constexpr coord_t operator*(CoordT const& scalar) const noexcept {
    coord_t ret;
    std::transform(_arr.begin(), _arr.end(), ret._arr.begin(), [&scalar](auto& i) { return i * scalar; });
    return ret;
  }

  constexpr std::partial_ordering operator<=>(coord_t const& other) const noexcept = default;

  constexpr auto begin() noexcept { return _arr.begin(); }
  constexpr auto begin() const noexcept { return _arr.begin(); }
  constexpr auto cbegin() const noexcept { return _arr.cbegin(); }
  constexpr auto end() noexcept { return _arr.end(); }
  constexpr auto end() const noexcept { return _arr.end(); }
  constexpr auto cend() const noexcept { return _arr.cend(); }

public:
  template<typename... Args, typename = std::enable_if_t<sizeof...(Args) == dims && (std::is_constructible_v<CoordT, Args> && ...)>>
  constexpr coord_t(Args... args) noexcept : _arr{static_cast<CoordT>(std::forward<Args>(args))...} {}

  constexpr coord_t() noexcept = default;
};

template<size_t dims, typename CoordT = double>
inline std::ostream& operator<<(std::ostream& os, coord_t<dims, CoordT> const& coord) {
  auto iter = coord.begin();
  os << '(' << *iter;
  for (++iter; iter != coord.end(); ++iter)
    os << ", " << *iter;
  os << ')';
  return os;
}

namespace std {
  template<size_t dims, typename CoordT>
  constexpr CoordT abs(coord_t<dims, CoordT> const& x) noexcept {
    return std::sqrt(std::accumulate(x.begin(), x.end(), CoordT{}, [](CoordT acc, CoordT next) { return acc + (next*next); }));
  }
}

static_assert(std::ranges::contiguous_range<coord_t<2>>);
