#include "coord.hpp"

#include <pcg_random.hpp>

#include <iostream>
#include <condition_variable>
#include <queue>
#include <random>
#include <thread>
#include <type_traits>
#include <latch>

/// Describes the points in some subset of the plane (usually a square)
class chunk {
private:
  std::vector<coord_t<2>> points;

public:
  constexpr size_t n_points() const noexcept { return points.size(); }
  constexpr bool empty() const noexcept { return points.empty(); }
  using find_nearest_ret_t = std::pair<decltype(points)::const_iterator, double>;
  // I tried optimising this with a tree, but I couldn't find a number of points per chunk where trees were faster than a big vector,
  // let alone faster than a small one
  inline find_nearest_ret_t find_nearest(coord_t<2> target) const {
    // If there are no points, then there is no nearest point
    if (points.empty())
      throw std::runtime_error{"find_nearest for empty chunk"};

    // Set the first point as the current best candidate (as we know there is at least one point)
    auto iter = points.begin();
    find_nearest_ret_t ret = {iter, std::abs(*iter - target)};

    // Look through the remanining points, replacing `ret` if we find a closer one
    for (++iter; iter != points.end(); ++iter)
      if (auto dist = std::abs(*iter - target); dist < ret.second)
        ret = {iter, dist};

    return ret;
  }
  inline void remove(find_nearest_ret_t::first_type point) {
    points.erase(point);
  }

  static coord_t<2> const& get_point(find_nearest_ret_t::first_type const& x) { return *x; }
  constexpr auto begin() noexcept { return points.begin(); }
  constexpr auto begin() const noexcept { return points.begin(); }
  constexpr auto cbegin() const noexcept { return points.cbegin(); }
  constexpr auto end() noexcept { return points.end(); }
  constexpr auto end() const noexcept { return points.end(); }
  constexpr auto cend() const noexcept { return points.cend(); }

public:
  chunk(std::vector<coord_t<2>> points) : points{std::move(points)} {}
};

/// A class that generates affine chunks on many threads
class chunk_generator {
private:
  /// Describes the state of the rng for each worker
  struct rng_state {
    std::poisson_distribution<size_t> n_points_dist;
    std::uniform_real_distribution<double> pos_dist;

    std::vector<coord_t<2>> generate(auto& rng) {
      std::vector<coord_t<2>> ret(n_points_dist(rng));

      for (auto& i : ret) {
        i[0] = pos_dist(rng);
        i[1] = pos_dist(rng);
      }

      return ret;
    }
  };

private:
  //
  size_t single_len = 8192;
  size_t backlog_size;

  // A condvar that is triggered when we need another batch for the backlog
  std::condition_variable used_backlog;
  std::mutex used_backlog_mutex;

  // This will still have some issues with slower ones being pushed to the front, but I have no good way of dealing with that, and it should *hopefully* not accumulate
  std::vector<std::vector<std::vector<coord_t<2>>>> backlog;
  std::mutex backlog_mutex;

  // The current batch
  std::vector<std::vector<coord_t<2>>> working_set;

  // The RNG used when we are waiting for a new batch
  std::random_device fallback_rng;
  rng_state fallback_generator;

  // Putting this last means that it is the first to be destroyed, joining threads before we destroy their dependencies
  std::vector<std::jthread> workers;

private:
  void worker(std::stop_token stop, rng_state state) {
    std::vector<std::vector<coord_t<2>>> new_batch;

    std::random_device base_rng;
    pcg32_fast batch_rng;
    std::random_device rng;
    std::uniform_int_distribution<decltype(batch_rng)::result_type> seed_dist;
    while (true) {
      // Wait for either a stop or a batch request
      {
        std::unique_lock lock{used_backlog_mutex};
        used_backlog.wait(lock, [&stop, this]() { return backlog.size() < backlog_size || stop.stop_requested(); });
      }

      if (stop.stop_requested())
        break;

      // Reseed the rng to ensure reasonable randomness
      //
      // We could just use std::random_device everywhere, but that makes zooming impractically slow
      batch_rng.seed(seed_dist(base_rng));

      new_batch.reserve(single_len);
      for (size_t i = 0; i < single_len; ++i)
        new_batch.push_back(state.generate(batch_rng));

      // Shove our new batch on the end
      {
        std::unique_lock lock{backlog_mutex};
        backlog.emplace_back(std::move(new_batch));
      }
    }

  }

private:
  std::vector<coord_t<2>> obtain_affine_chunk() {
    if (working_set.size() == 0) {
      std::unique_lock lock{backlog_mutex};
      if (working_set.size() != 0)
        ;
      else if (backlog.size() == 0) {
        lock.unlock();
        std::cout << "STALL" << std::endl;
        return fallback_generator.generate(fallback_rng);
      }
      else {
        working_set = std::move(backlog.back());
        backlog.pop_back();
        used_backlog.notify_one();
      }
    }

    std::vector<coord_t<2>> ret = std::move(working_set.back());
    working_set.pop_back();
    return ret;
  }

public:
  chunk create_chunk(coord_t<2> const& centre) {
    // Grab an affine chunk...
    auto dat = obtain_affine_chunk();

    // ... and translate it to our new centre
    for (auto& i : dat)
      i += centre;

    return {std::move(dat)};
  }

public:
  chunk_generator(double radius, size_t backlog_size_ = 16, int n_workers = -1) : backlog_size{backlog_size_} {
    // Create the fallback generator, that we will copy for each worker
    std::poisson_distribution<size_t> n_points_dist{radius*radius*4};
    std::uniform_real_distribution<double> pos_dist{-radius, radius};
    fallback_generator = {n_points_dist, pos_dist};

    // If the user hasn't limited our worker count, assume they want us to completely lock up their cpu
    if (n_workers < 0)
      n_workers = std::thread::hardware_concurrency();

    // Create our worker threads
    for (int i = 0; i < n_workers; ++i)
      workers.emplace_back([this](std::stop_token tok) { worker(tok, fallback_generator); });
  }

  ~chunk_generator() {
    for (auto& i : workers)
      i.request_stop();
    used_backlog.notify_all();
  }
};
