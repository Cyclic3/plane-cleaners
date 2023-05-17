#include "space.hpp"

#include <SDL2/SDL.h>
// why does sdl do this on Windows?
#undef main

#include <iostream>
#include <mutex>
#include <thread>
#include <latch>
#include <optional>

struct sdl_exception : std::runtime_error {
  sdl_exception() : std::runtime_error{SDL_GetError()} {}
};

class thread_status {
  std::atomic<bool> _stopped = false;
  std::atomic<bool> _started = false;

public:
  bool has_started() const { return _started; }
  bool has_stopped() const { return _stopped; }
  void lock() { _started = true; }
  void unlock() { _stopped = true; }
};

class graphics {
private:
  std::unique_ptr<SDL_Renderer, C_DELETER(SDL_DestroyRenderer)> renderer;
  std::unique_ptr<SDL_Window, C_DELETER(SDL_DestroyWindow)> window;
  thread_status worker_status;
  int cleaner_size = 1;
  std::mutex renderer_mutex;
  std::atomic<uint64_t> iter_count = 0;

  std::jthread worker_thread;

  std::mutex param_mutex;
  struct {
    /// Current zoom level
    double scale = 1;
    /// Current middle of the screen
    coord_t<2> centre = {0,0};
    /// Current screen dimensions
    coord_t<2, int> dims_int;
    /// The number of iterations we do in a single step
    uint64_t iter_step = 1;
  } _params;

  struct {
    /// Controls how sensitive the scroll wheel is
    double zoom_power = 1.2;
    /// Controls how sensitive the keyboard navigation is
    double rel_step = 10;
    /// Controls how sensitive dragging is. The value of 1 means that dragging moves precisely with the mouse
    double mouse_rel_step = 1;
  } tweaks;

private:
  /// Handles user input
  void worker_body(std::stop_token tok, std::latch& ready_latch) {
    // Windows needs this to be in the same thread as the the user input logic
    SDL_Init(SDL_INIT_VIDEO);
    {
      SDL_Renderer* renderer_ptr;
      SDL_Window* window_ptr;
      // We don't need to lock yet, as no other threads are running
      if (SDL_CreateWindowAndRenderer(_params.dims_int[0], _params.dims_int[1], SDL_WINDOW_RESIZABLE, &window_ptr, &renderer_ptr))
        throw sdl_exception{};
      if (!renderer_ptr)
        throw sdl_exception{};
      renderer.reset(renderer_ptr);
      window.reset(window_ptr);
    }

    ready_latch.count_down();

    std::unique_lock lock{worker_status};
    SDL_Event event;

    bool mouse_down = false;

    while (!tok.stop_requested()) {
      SDL_WaitEvent(&event);
      switch (event.type) {
        // Scroll zoom
        case SDL_MOUSEWHEEL: {
          if (event.wheel.y < 0) {
            std::unique_lock lock{param_mutex};
            _params.scale *= tweaks.zoom_power;
          }
          else if (event.wheel.y > 0) {
            std::unique_lock lock{param_mutex};
            _params.scale /= tweaks.zoom_power;
          }
        } break;

         // Drag panning
        case SDL_MOUSEBUTTONDOWN: {
          if (event.button.button == SDL_BUTTON_LEFT)
            mouse_down = true;
        } break;
        case SDL_MOUSEBUTTONUP: {
          if (event.button.button == SDL_BUTTON_LEFT)
            mouse_down = false;
        } break;
        case SDL_MOUSEMOTION: {
          if (!mouse_down)
            break;

          std::unique_lock lock{param_mutex};
          _params.centre -= coord_t<2>{event.motion.xrel, event.motion.yrel} * (tweaks.mouse_rel_step * _params.scale);
        } break;

        // Keyboard navigation
        case SDL_KEYDOWN: {
          auto key = SDL_GetKeyFromScancode(event.key.keysym.scancode);
          switch (key) {
            case 'i': {
              std::unique_lock lock{param_mutex};
              _params.scale /= tweaks.zoom_power;
            } break;

            case 'o':{
              std::unique_lock lock{param_mutex};
              _params.scale *= tweaks.zoom_power;
            } break;

            case 'a':
            case SDLK_LEFT:
            {
              std::unique_lock lock{param_mutex};
              _params.centre -= coord_t<2>{tweaks.rel_step * _params.scale, 0};
            } break;

            case 'd':
            case SDLK_RIGHT:
            {
              std::unique_lock lock{param_mutex};
              _params.centre += coord_t<2>{tweaks.rel_step * _params.scale, 0};
            } break;

            case 'w':
            case SDLK_UP:
            {
              std::unique_lock lock{param_mutex};
              _params.centre -= coord_t<2>{0, tweaks.rel_step * _params.scale};
            } break;

            case 's':
            case SDLK_DOWN:
            {
              std::unique_lock lock{param_mutex};
              _params.centre += coord_t<2>{0, tweaks.rel_step * _params.scale};
            } break;

            case 'p':
            {
              std::unique_lock lock{param_mutex};
              std::unique_lock renderer_lock{renderer_mutex};
              std::vector<uint32_t> out(_params.dims_int[0] * _params.dims_int[1]);
              std::unique_ptr<SDL_Surface, C_DELETER(SDL_FreeSurface)> surface(SDL_CreateRGBSurfaceWithFormat(0, _params.dims_int[0], _params.dims_int[1], 24, SDL_PIXELFORMAT_RGB24));
              if (!surface)
                throw sdl_exception{};
              if (SDL_RenderReadPixels(renderer.get(), nullptr, SDL_PIXELFORMAT_RGB24, surface->pixels, surface->pitch))
                throw sdl_exception{};
              std::string fname = "plane-cleaners_" + std::to_string(iter_count) + ".bmp";
              if (SDL_SaveBMP(surface.get(), fname.c_str()))
                throw sdl_exception{};
            } break;

            case '+':
            case '=':
            {
              std::unique_lock lock{param_mutex};
              _params.iter_step = std::max(_params.iter_step, _params.iter_step*2);
            } break;

            case '-':
            {
              std::unique_lock lock{param_mutex};
              _params.iter_step = std::max<uint64_t>(1, _params.iter_step/2);
            } break;

            case '0':
            {
              std::unique_lock lock{param_mutex};
              _params.iter_step = 1;
            } break;

            default: {}
          }
        } break;

        // Window resizing
        case SDL_WINDOWEVENT: {
          if (event.window.event != SDL_WINDOWEVENT_RESIZED)
            break;
          std::unique_lock lock{param_mutex};
          _params.dims_int = {event.window.data1, event.window.data2};
        } break;

        // Closing
        case SDL_QUIT: {
          return;
        } break;

        default: {}
      }
    }
  }

  /// Walks through the given rectangle, generating black/white pixels for each point
  ///
  // TODO: fix this being pretty slow
  std::vector<uint32_t> get_pixels(space& s, coord_t<2> start, coord_t<2> finish, double scale, coord_t<2, int> const& dims_int) {
    std::vector<uint32_t> out(dims_int[0] * dims_int[1], 0xffffff);
    s.examine_rectangle(start, finish, [&start, &out, &scale, &dims_int](space::examine_arg_t const& c) {
      for (auto& point : c->second.data) {
        auto pos = (point - start) / scale;
        coord_t<2, int64_t> pos_int{pos[0], pos[1]};
        if (pos_int[0] < 0 || pos_int[1] < 0 || pos_int[0] >= dims_int[0] || pos_int[1] >= dims_int[1])
          return;
        out[pos_int[0] + pos_int[1] * dims_int[0]] = 0;
      }
    });
    return out;
  }


public:
  void set_iter_count(uint64_t iters) {
    iter_count = iters;
  }
  // A safe way of reading the parameters
  decltype(_params) get_params() {
    std::unique_lock lock{param_mutex};
    return _params;
  }

  void draw_space(space& s, std::vector<coord_t<2>> const& cleaners) {
    auto params = get_params();

    // Calculate the window position and dimensions
    coord_t<2> dims{params.dims_int[0], params.dims_int[1]};
    coord_t<2> start = params.centre - ((dims / 2) * params.scale);
    coord_t<2> finish = params.centre + ((dims / 2) * params.scale);

    // Load the pixels
    auto pix = get_pixels(s, start, finish, params.scale, params.dims_int);

    std::unique_lock renderer_lock{renderer_mutex};

    // Draw the pixels
    std::unique_ptr<SDL_Surface, C_DELETER(SDL_FreeSurface)> surface{SDL_CreateRGBSurfaceFrom(pix.data(), params.dims_int[0], params.dims_int[1], 32, params.dims_int[0] * sizeof(uint32_t), 0x000000ff, 0x0000ff00, 0x00ff0000, 0x00000000)};
    if (!surface)
      throw sdl_exception{};
    std::unique_ptr<SDL_Texture, C_DELETER(SDL_DestroyTexture)> texture{SDL_CreateTextureFromSurface(renderer.get(), surface.get())};
    if (!texture)
      throw sdl_exception{};
    if (SDL_RenderCopy(renderer.get(), texture.get(), nullptr, nullptr))
      throw sdl_exception{};

    // Draw our cleaners on top, as red boxes rather than pixels so that they can be seen
    if (SDL_SetRenderDrawColor(renderer.get(), 255, 0, 0, SDL_ALPHA_OPAQUE))
      throw sdl_exception{};
    for (auto& cleaner : cleaners) {
      auto pix_pos = (cleaner - start) / params.scale;
      SDL_Rect rect{static_cast<int>(pix_pos[0] - cleaner_size), static_cast<int>(pix_pos[1] - cleaner_size), 1 + cleaner_size * 2, 1 + cleaner_size * 2};
      if (SDL_RenderDrawRect(renderer.get(), &rect))
        throw sdl_exception{};
    }

    // Finally, actually draw the result
    SDL_RenderPresent(renderer.get());
  }

  /// Returns true iff the user requested an exit
  bool should_loop() {
    return !worker_status.has_stopped();
  }

public:
  graphics(coord_t<2, int> dims_int_) {
    _params.dims_int = dims_int_;
    std::latch ready_latch(1);
    worker_thread = std::jthread{[this, &ready_latch](std::stop_token tok) { worker_body(tok, ready_latch); }};
    ready_latch.wait();
  }
};

class cleaner_manager {
private:
  // The maximum value of `current_time` before we reset the clocks, and move back `when_next` for each cleaner;
  static constexpr double time_bound = 1e9;

private:
  struct cleaner_time_info {
    size_t cleaner_idx;
    double when_next;
  };
  static constexpr auto cleaner_time_cmp = [](auto& a, auto& b) { return a.when_next > b.when_next; };
  struct cleaner_info {
    coord_t<2> pos = {0, 0};
    std::optional<space::find_nearest_ret_t::start_t> start;
  };

private:
  size_t n_cleaners;
  std::vector<cleaner_info> cleaners;
  std::vector<cleaner_time_info> next_cleaner_heap;
  double current_time;

  // These are opaque, so it does not matter if these change unexpectedly,
  // as long as it isn't by two threads simultaneously
  mutable pcg32_fast cleaner_time_rng;
  mutable std::exponential_distribution<double> cleaner_time_dist;

private:
  double gen_time() const {
    return current_time + cleaner_time_dist(cleaner_time_rng);
  }
  /// XXX: this mutates the underlying state, and can change the `when_next` on literally every single cleaner, as well as reset current_time to 0
  size_t get_next_cleaner() {
    std::pop_heap(next_cleaner_heap.begin(), next_cleaner_heap.end(), cleaner_time_cmp);
    // Copy the cleaner, because we are going to mutate it
    cleaner_time_info next_cleaner = next_cleaner_heap.back();
    // Time travel
    current_time = next_cleaner.when_next;
    // If we are too far in the future, rewind all the clocks
    if (current_time > time_bound) {
      for (auto& i : next_cleaner_heap)
        i.when_next -= current_time;
      current_time = 0;
    }
    // Give a new time to the next cleaner, and insert it into the heap
    next_cleaner_heap.back().when_next = gen_time();
    std::push_heap(next_cleaner_heap.begin(), next_cleaner_heap.end(), cleaner_time_cmp);
    // Return the index of the next cleaner
    return next_cleaner.cleaner_idx;
  }
public:
  void step_one(space& s) {
    auto& cleaner = cleaners.at(get_next_cleaner());

    auto nearest = cleaner.start ? s.find_nearest(cleaner.pos, *cleaner.start) : s.find_nearest(cleaner.pos);

    cleaner.pos = nearest.get_point();
    cleaner.start = nearest.get_start();
    // Remove after read to prevent early invalidation
    s.remove(nearest.iter);
  }
  std::vector<coord_t<2>> get_cleaners_pos() const {
    std::vector<coord_t<2>> ret;
    ret.reserve(n_cleaners);
    for (auto& i: cleaners)
      ret.push_back(i.pos);
    return ret;
  }

public:
  cleaner_manager(size_t _n_cleaners) : n_cleaners{_n_cleaners}, cleaners(_n_cleaners), current_time{0} {
    if (n_cleaners <= 0)
      throw std::invalid_argument{"Asked for zero cleaners?"};

    // Seed the rng used for generating cleaner times
    std::random_device good_rng;
    std::uniform_int_distribution<decltype(cleaner_time_rng)::result_type> cleaner_time_rng_seed_dist;
    cleaner_time_rng.seed(cleaner_time_rng_seed_dist(good_rng));

    // Generate the cleaner times
    for (size_t cleaner_idx = 0; cleaner_idx < n_cleaners; ++cleaner_idx)
      next_cleaner_heap.emplace_back(cleaner_time_info{.cleaner_idx = cleaner_idx, .when_next = gen_time()});
    // Order the cleaners into a heap
    std::make_heap(next_cleaner_heap.begin(), next_cleaner_heap.end(), cleaner_time_cmp);
  }
};

int main(int argc, char** argv) {
  space s;
  graphics g({480, 480});

  // Default to 1 cleaner
  size_t n_cleaners = (argc > 1 ? std::stoull(argv[1]) : 2);
  cleaner_manager cleaners{n_cleaners};

  auto start_time = std::chrono::system_clock::now();
  uint64_t count = 0;
  while (g.should_loop()) {
    auto step_start = std::chrono::system_clock::now();
    auto params = g.get_params();
    for (uint64_t this_count = 0; this_count < params.iter_step; ++this_count, ++count)
      cleaners.step_one(s);

    auto step_end = std::chrono::system_clock::now();
    auto total_time = std::chrono::duration_cast<std::chrono::duration<double>>(step_end - start_time).count();
    auto step_time = std::chrono::duration_cast<std::chrono::duration<double>>(step_end - step_start).count();
    std::cout << count << " iters (" << (params.iter_step/n_cleaners) << " avg per worker per step) in " << total_time << "s (step has " << (params.iter_step / (step_time * 1000)) << "kiter/s)" << std::endl;

    g.set_iter_count(count*n_cleaners);
    g.draw_space(s, cleaners.get_cleaners_pos());
  }

  // Dump out all the cleaner final positions
  for (auto& cleaner: cleaners.get_cleaners_pos()) {
    std::cout << cleaner << ": " << std::abs(cleaner) << std::endl;
  }

  return 0;
}
