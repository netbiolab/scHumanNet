// Taken from RcppParallel.h and then modified slightly to rename header guards
// and namespaces to avoid any potential clashes. Parallel is licensed under
// GPLv2

#ifndef PERPENDICULAR
#define PERPENDICULAR

#include <thread>
#include <utility>
#include <vector>

namespace Perpendicular {

using IndexRange = std::pair<std::size_t, std::size_t>;

template <typename Worker>
auto worker_thread(Worker &worker, const IndexRange &range) -> void {
  try {
    worker(range.first, range.second);
  } catch (...) {
  }
}

// Function to calculate the ranges for a given input
inline auto split_input_range(const IndexRange &range, std::size_t thread_no,
                              std::size_t grain_size)
    -> std::vector<IndexRange> {
  // determine max number of threads
  if (thread_no <= 0) {
    thread_no = SYS_THREADS_DEF;
  }

  // compute grain_size (including enforcing requested minimum)
  std::size_t length = range.second - range.first;
  if (thread_no == 1)
    grain_size = length;
  else if ((length % thread_no) == 0)  // perfect division
    grain_size = (std::max)(length / thread_no, grain_size);
  else  // imperfect division, divide by threads - 1
    grain_size = (std::max)(length / (thread_no - 1), grain_size);

  // allocate ranges
  std::vector<IndexRange> ranges;
  std::size_t begin = range.first;
  while (begin < range.second) {
    std::size_t end = (std::min)(begin + grain_size, range.second);
    ranges.emplace_back(std::make_pair(begin, end));
    begin = end;
  }

  return ranges;
}

// Execute the Worker over the IndexRange in parallel
template <typename Worker>
inline void parallel_for(std::size_t begin, std::size_t end, Worker &worker,
                         std::size_t thread_no, std::size_t grain_size = 1) {
  // split the work
  IndexRange input_range(begin, end);
  std::vector<IndexRange> ranges =
      split_input_range(input_range, thread_no, grain_size);

  std::vector<std::thread> threads;
  for (auto &range : ranges) {
    threads.push_back(
        std::thread(&worker_thread<Worker>, std::ref(worker), range));
  }

  for (auto &thread : threads) {
    thread.join();
  }
}

}  // namespace Perpendicular

#endif  // PERPENDICULAR
