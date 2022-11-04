#include "chunk.hpp"

#include <chrono>
#include <iostream>
#include <map>
#include <mutex>
#include <shared_mutex>

class space {
private:
  using chunk_id_t = coord_t<2, int64_t>;
  struct chunk_info;
  using chunk_map = std::map<chunk_id_t, chunk_info>;
  struct chunk_info {
    /// The actual points in the chunk
    chunk data;
    /// Either `active_chunks.end()` or an iterator for the next chunk to the right (+ {1, 0})
    chunk_map::iterator cached_right;
    /// Either `active_chunks.end()` or an iterator for the next chunk up (+ {0, 1})
    chunk_map::iterator cached_up;
    /// either `active_chunks.end()` or an iterator for the next chunk (- {1, 1})
    ///
    /// This allows us to easily iterate a L_\infty search for a point, as the bottom right point of an L_\infty circle with 1 greater radius is 1 down and 1 to the right
    chunk_map::iterator cached_next;
  };

private:
  /// The side length of each chunk
  //
  // 6 seems to be near some local maximum on my machine, but this has not been rigorously selected
  static constexpr double chunk_size = 6;
  /// The set of chunks that currently reside in memory
  //
  // TODO: implement swapping/memory efficient storage for LRU points
  chunk_map active_chunks;
  chunk_generator chunk_gen{chunk_size / 2};

public:
  /// Find the chunk that contains a point
  constexpr chunk_id_t calc_chunk_id(coord_t<2> const& x) const noexcept {
    return {std::round(x[0] / chunk_size), std::round(x[1] / chunk_size)};
  }
  // XXX: These chunks will overlap a tiny bit due to float imprecision, but that should (hopefully) not have too much of an effect on the result
  constexpr coord_t<2> get_chunk_centre(chunk_id_t chunk) const noexcept {
    return {chunk[0] * chunk_size, chunk[1] * chunk_size};
  }
private:
  chunk_map::iterator load_chunk(chunk_id_t chunk_id) {
    // Lower bound makes insert faster, roughly doubling the speed of this function
    auto iter = active_chunks.lower_bound(chunk_id);
    // If our chunk already exists, return it
    if (iter != active_chunks.end() && iter->first == chunk_id)
      return iter;

    // Otherwise, try to find the nearby chunks to cache
    auto right = active_chunks.find(chunk_id + chunk_id_t{1, 0});
    auto up = active_chunks.find(chunk_id + chunk_id_t{0, 1});
    auto next = active_chunks.find(chunk_id - chunk_id_t{1, 1});

    return active_chunks.emplace_hint(iter, chunk_id, chunk_info{chunk_gen.create_chunk(get_chunk_centre(chunk_id)), right, up, next});
  }
  chunk_map::iterator load_right(chunk_map::iterator const& chunk) {
    if (chunk->second.cached_right != active_chunks.end())
      return chunk->second.cached_right;
    return chunk->second.cached_right = load_chunk(chunk->first + chunk_id_t{1, 0});
  }
  chunk_map::iterator load_up(chunk_map::iterator const& chunk) {
    if (chunk->second.cached_up != active_chunks.end())
      return chunk->second.cached_up;
    return chunk->second.cached_up = load_chunk(chunk->first + chunk_id_t{0, 1});
  }
  chunk_map::iterator load_next(chunk_map::iterator const& chunk) {
    if (chunk->second.cached_next != active_chunks.end())
      return chunk->second.cached_next;
    return chunk->second.cached_next = load_chunk(chunk->first - chunk_id_t{1, 1});
  }

public:
  // A convienience wrapper for the result of a point search
  struct find_nearest_ret_t {
    std::pair<chunk_map::iterator, chunk::find_nearest_ret_t::first_type> iter;
    double distance;

    using start_t = chunk_map::iterator;
    constexpr coord_t<2> const& get_point() const noexcept { return chunk::get_point(iter.second); }
    constexpr start_t const& get_start() const noexcept { return iter.first; }
  };
  find_nearest_ret_t find_nearest(coord_t<2> target) {
    auto chunk_id = calc_chunk_id(target);

    return find_nearest(target, load_chunk(chunk_id));
  }

  find_nearest_ret_t find_nearest(coord_t<2> target, chunk_map::iterator start_point) {
    // Start with some hypothetical point infinitely far away
    find_nearest_ret_t best{{}, std::numeric_limits<double>::infinity()};
    double l_inf_radius_upper_bound = std::numeric_limits<double>::infinity();
    auto target_offset = target - get_chunk_centre(start_point->first);
    double l_inf_centre_dist = std::max(std::abs(target_offset[0]), std::abs(target_offset[1]));

    int64_t l_inf_radius = 1;
    // DRY function for checking the points in a chunk
    auto handle_one = [&best, &target, &l_inf_radius_upper_bound, &l_inf_centre_dist](chunk_map::iterator const& chunk_iter) {
      if (chunk_iter->second.data.empty())
        return;

      // Calculating a lower bound for the distance to points in this chunk actually takes more time that it saves

      // Check all the points in the chunk
      if (auto closest_in_chunk = chunk_iter->second.data.find_nearest(target); closest_in_chunk.second < best.distance) {
        best = {{chunk_iter, closest_in_chunk.first}, closest_in_chunk.second};
        // Checking to see if we're on the last iteration slows us down for some reason
        //
        // So too does breaking this out into the main loop, and explicitly checking if we need another loop

        // This bound is justified in the README
        l_inf_radius_upper_bound = (l_inf_centre_dist + best.distance) / chunk_size + 0.5;
      }
    };

    // Check our centre point
    handle_one(start_point);

    // Iterate through remaining points
    chunk_map::iterator last_start = start_point;
    for (; l_inf_radius <= l_inf_radius_upper_bound; ++l_inf_radius) {
      // Step out one
      last_start = load_next(last_start);

      chunk_map::iterator pos = last_start;
      // Load bottom side
      for (int64_t x = -l_inf_radius; x < l_inf_radius; ++x, pos = load_right(pos))
        handle_one(pos);
      // Handle row end
      handle_one(pos);
      // Load right side
      for (int64_t y = -l_inf_radius; y < l_inf_radius; ++y)
        handle_one(pos = load_up(pos));
      // Load left side
      pos = last_start;
      for (int64_t y = -l_inf_radius; y < l_inf_radius; ++y)
        handle_one(pos = load_up(pos));
      // Load remaining top side
      for (int64_t x = -l_inf_radius; x < l_inf_radius; ++x)
        handle_one(pos = load_right(pos));
    }

    // Return our best candidate
    return best;
  }

  /// Iterates through the given rectangle, but single threaded
  ///
  /// XXX: assumes bottom_right is neither above nor to the right of top_right
  using examine_arg_t = chunk_map::iterator;
  template<typename F>
  void examine_rectangle_st(coord_t<2> bottom_left, coord_t<2> top_right, F func) {
    chunk_id_t start_chunk_id = calc_chunk_id(bottom_left);
    chunk_id_t finish_chunk_id = calc_chunk_id(top_right);

    auto x_chunk = load_chunk(start_chunk_id);
    for (int64_t x = start_chunk_id[0]; x < finish_chunk_id[0]; ++x, x_chunk = load_right(x_chunk)) {
      auto y_chunk = x_chunk;
      func(y_chunk);
      for (int64_t y = start_chunk_id[1]; y < finish_chunk_id[1]; ++y)
        func(y_chunk = load_up(y_chunk));
    }
    func(x_chunk);
    for (int64_t y = start_chunk_id[1]; y < finish_chunk_id[1]; ++y)
      func(x_chunk = load_up(x_chunk));
  }
  /// Iterates through the given rectangle, but multithreaded
  ///
  /// XXX: assumes bottom_right is neither above nor to the right of top_right
  template<typename F>
  void examine_rectangle(coord_t<2> bottom_left, coord_t<2> top_right, F func) {
    chunk_id_t start_chunk_id = calc_chunk_id(bottom_left);
    chunk_id_t finish_chunk_id = calc_chunk_id(top_right);

    // This mutex protects the map from being added to with load_chunk concurrently,
    // but allows many threads to read through it concurrently
    std::shared_mutex map_mutex;


    // This mutex protects the simple scheduling logic, so that the threads can just take the next one they need
    std::mutex scheduler_mutex;
    auto next_x_chunk = load_chunk(start_chunk_id);
    int64_t next_column = start_chunk_id[0];

    int n_threads = std::jthread::hardware_concurrency();
    std::vector<std::jthread> threads;

    for (int thread_no = 0; thread_no < n_threads; ++thread_no) {
      threads.emplace_back([&map_mutex, &next_column, &finish_chunk_id, &scheduler_mutex, &next_x_chunk, &start_chunk_id, &func, this] {
        // The condition is too complex to reasonably fit into the while condition
        while (true) {
          // These will be initialised by the locked region
          chunk_map::iterator iter;
          {
            // Lock the scheduler state
            std::unique_lock x_lock{scheduler_mutex};
            // If we have run out of columns, stop processing
            if (next_column > finish_chunk_id[0])
              return;

            // Otherwise, we take the initial chunk
            iter = next_x_chunk;

            // Now we prepare the scheduler state for the next thread, because we're nice
            ++next_column;
            // If we're not at the end, then get the next chunk to the right
            if (next_column <= finish_chunk_id[0]) {
              // If we haven't cached the right chunk, we need to write-lock the map and load it
              if (next_x_chunk->second.cached_right == active_chunks.end()) {
                std::unique_lock lock{map_mutex};
                next_x_chunk = load_right(next_x_chunk);
              }
              // Otherwise, we can just set it without worrying about map state
              else
                next_x_chunk = next_x_chunk->second.cached_right;
            }
          }
          std::shared_lock map_lock{map_mutex};

          // Run the function on the first element
          func(iter);
          for (int64_t y = start_chunk_id[1]; y < finish_chunk_id[1]; ++y) {
            // If we haven't cached the up chunk, we need to write-lock the map and load it
            if (iter->second.cached_up == active_chunks.end()) {
              // Stop holding a read lock, as we want a write lock, and so to do otherwise would cause a deadlock
              map_lock.unlock();
              {
                // Write lock the map
                std::unique_lock lock{map_mutex};
                iter = load_up(iter);
              }
              // Reaquire our read lock
              map_lock.lock();
            }
            // Otherwise, we can just set it without worrying about map state
            else
              iter = iter->second.cached_up;
            // Run the function on the next element
            func(iter);
          }
        }
      });
    }
  }

  /// XXX: all functions on this class will now become invalid, except for using it as a start for a search
  void remove(decltype(find_nearest_ret_t::iter) point) {
    point.first->second.data.remove(point.second);
  }
};
