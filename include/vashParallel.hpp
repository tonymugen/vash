/*
 * Copyright (c) 2026 Anthony J. Greenberg
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
 * IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

/// Parallel-execution configuration
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2026 Anthony J. Greenberg
 * \version 0.6
 *
 * Single point of control for the parallel STL backend. `VASH_HAVE_TBB` is
 * defined by CMake when `std::execution::par` is implemented with Intel [Threading Building Blocks (TBB)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onetbb.html),
 * the only standard library configuration that lets us cap the worker-thread count. When it is absent the
 * algorithms still run, but serially (`std::execution::seq`) and the thread ceiling is a no-op.
 *
 * Use `BayesicSpace::parallelPolicy` as the first argument to standard
 * algorithms, and instantiate a `BayesicSpace::ThreadCeiling` for the duration
 * of a parallel region to limit concurrency.
 */

#pragma once

#include <cstddef>
#include <execution>

#if defined(VASH_HAVE_TBB)
#include <optional>
#include <thread>
#include <tbb/global_control.h>
#endif

namespace BayesicSpace {
#if defined(VASH_HAVE_TBB)
	/** \brief Execution policy for parallel algorithms.
	 *
	 * Parallel because a capping-capable TBB backend is present.
	 */
	inline constexpr std::execution::parallel_policy parallelPolicy{std::execution::par};

	/** \brief RAII worker-thread ceiling.
	 *
	 * While an instance is alive, TBB (and therefore `std::execution::par`) uses
	 * at most the requested number of worker threads process-wide. The limit is
	 * lifted when the object is destroyed. A request of 0, or one at or above
	 * `std::thread::hardware_concurrency()`, leaves the TBB default in place
	 * (effectively uncapped) rather than installing a redundant control.
	 *
	 * The ceiling is global, not per-call, so construct one in the outermost
	 * scope that owns the user-requested thread count.
	 */
	class ThreadCeiling {
	public:
		/** \brief Constructor
		 *
		 * \param[in] maxThreads maximum number of worker threads; 0 (or >= the hardware maximum) means no limit
		 */
		explicit ThreadCeiling(size_t maxThreads) {
			const unsigned int hwConcurrency{std::thread::hardware_concurrency()};
			if ( maxThreads > 0 && (hwConcurrency == 0 || maxThreads < static_cast<size_t>(hwConcurrency)) ) {
				control_.emplace(tbb::global_control::max_allowed_parallelism, maxThreads);
			}
		}
	private:
		/** \brief Active limit, engaged only when a sub-maximal ceiling is requested */
		std::optional<tbb::global_control> control_;
	};
#else
	/** \brief Execution policy for parallel algorithms.
	 *
	 * Sequential because no capping-capable backend was found at build time; this
	 * keeps results correct without an enforceable thread ceiling.
	 */
	inline constexpr std::execution::sequenced_policy parallelPolicy{std::execution::seq};

	/** \brief No-op worker-thread ceiling for builds without a TBB backend. */
	class ThreadCeiling {
	public:
		/** \brief Constructor (ignores the requested count)
		 *
		 * \param[in] maxThreads ignored
		 */
		explicit ThreadCeiling(size_t /* maxThreads */) noexcept {}
	};
#endif
} // namespace BayesicSpace
