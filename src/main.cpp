#include "space.hpp"

#include <SDL2/SDL.h>

#include <iostream>
#include <mutex>
#include <thread>
#include <latch>

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
    coord_t<2, int64_t> dims_int;
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
  void worker_body(std::stop_token tok) {
    std::unique_lock lock{worker_status};
    SDL_Event event;

    bool mouse_down = false;

    while (!tok.stop_requested()) {
      std::this_thread::sleep_for(std::chrono::milliseconds{1});
      SDL_PollEvent(&event);
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
  std::vector<uint32_t> get_pixels(space& s, coord_t<2> start, coord_t<2> finish, double scale, coord_t<2, int64_t> const& dims_int) {
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

  void draw_space(space& s, std::vector<coord_t<2>>& cleaners) {
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
  graphics(coord_t<2, int64_t> dims_int_) {
    _params.dims_int = dims_int_;
    SDL_Init(SDL_INIT_VIDEO);
    {
      SDL_Renderer* renderer_ptr;
      SDL_Window* window_ptr;
      if (SDL_CreateWindowAndRenderer(dims_int_[0], dims_int_[1], SDL_WINDOW_RESIZABLE, &window_ptr, &renderer_ptr))
        throw sdl_exception{};
      if (!renderer_ptr)
        throw sdl_exception{};
      renderer.reset(renderer_ptr);
      window.reset(window_ptr);
    }
    worker_thread = std::jthread{&graphics::worker_body, this};
  }
};

int main(int argc, char** argv) {
  space s;
  graphics g({480, 480});

  // Default to 1 cleaner
  size_t n_cleaners = (argc > 1 ? std::stoull(argv[1]) : 1);

  // Load the cleaners at the origin
  std::vector<coord_t<2>> cleaners(n_cleaners, {0,0});
  std::vector<space::find_nearest_ret_t::start_t> cleaner_starts(n_cleaners);

  // Do initial iteration on cleaners to find their start position
  for (size_t cleaner_no = 0; cleaner_no < n_cleaners; ++cleaner_no) {
    auto& cleaner = cleaners[cleaner_no];
    auto& cleaner_start = cleaner_starts[cleaner_no];

    auto nearest = s.find_nearest(cleaner);

    cleaner = nearest.get_point();
    cleaner_start = nearest.get_start();
    // Remove after read to prevent early invalidation
    s.remove(nearest.iter);
  }

  auto start_time = std::chrono::system_clock::now();
  uint64_t count = 0;
  while (g.should_loop()) {
    auto step_start = std::chrono::system_clock::now();
    auto params = g.get_params();
    for (uint64_t this_count = 0; this_count < params.iter_step; ++this_count, ++count){
      for (size_t cleaner_no = 0; cleaner_no < n_cleaners; ++cleaner_no) {
        auto& cleaner = cleaners[cleaner_no];
        auto& cleaner_start = cleaner_starts[cleaner_no];

        auto nearest = s.find_nearest(cleaner, cleaner_start);

        cleaner = nearest.get_point();
        cleaner_start = nearest.get_start();
        // Remove after read to prevent early invalidation
        s.remove(nearest.iter);
      }
    }

    auto step_end = std::chrono::system_clock::now();
    auto total_time = std::chrono::duration_cast<std::chrono::duration<double>>(step_end - start_time).count();
    auto step_time = std::chrono::duration_cast<std::chrono::duration<double>>(step_end - step_start).count();
    std::cout << count*n_cleaners << " iters (" << params.iter_step << " per worker per step) in " << total_time << "s (step has " << (n_cleaners * params.iter_step / (step_time * 1000)) << "kiter/s)" << std::endl;

    g.set_iter_count(count*n_cleaners);
    g.draw_space(s, cleaners);
  }

  // Dump out all the cleaner final positions
  for (auto& cleaner: cleaners) {
    std::cout << cleaner << ": " << std::abs(cleaner) << std::endl;
  }

  return 0;
}
