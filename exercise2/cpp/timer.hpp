#pragma once
#ifndef TIMER_HPP
#define TIMER_HPP
#include <chrono>
#include <iostream>

// Timer class for measuring elapsed time
template <typename T = std::chrono::milliseconds> struct Timer {
  using Clock = std::chrono::high_resolution_clock;
  using Duration = T;

  Clock::time_point start_time_wall;

  Timer() {}
  ~Timer() {
    // Print
  }
  void reset() { /* reset start time */ }

  T elapsed() const {
    // Return during time
  }

  friend std::ostream &operator<<(std::ostream &os, const Timer &timer) {
    // TODO:
  }

  template <typename U> static constexpr std::string_view time_unit_name() {
    if constexpr (std::is_same_v<U, std::chrono::seconds>)
      return "s"; // Support More units
    return "?";
  }
}; // Timer
#endif // TIMER_HPP