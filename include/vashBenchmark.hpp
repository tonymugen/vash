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

/// In-source phase timers for benchmarking (benchmarking branch only)
/** \file
 * \author Anthony J. Greenberg
 * \version 0.1
 *
 * Lightweight wall-clock phase markers used to locate serial bottlenecks inside otherwise
 * threaded methods (e.g. `GenoTableHash::ldInGroups` and its serial `makeLDgroups` phase).
 *
 * The markers expand to real `std::chrono` measurements only when `VASH_BENCHMARK` is defined
 * (the benchmark build sets it); in every other build they expand to a discarded no-op so the
 * production/test code generation is byte-for-byte unchanged. Measurements are written to
 * `std::cerr` with a `[vash-bench]` prefix so they are visible immediately without touching the
 * logger. Use `VASH_BENCH_TP(name)` to stamp a starting time point and `VASH_BENCH_LAP(label, name)`
 * to report the elapsed milliseconds since that stamp and re-arm it for the next phase.
 */
#ifndef VASH_BENCHMARK_HPP
#define VASH_BENCHMARK_HPP

#ifdef VASH_BENCHMARK

#include <chrono>
#include <iostream>

// Stamp a monotonic starting time point in a local variable. A function cannot introduce a
// caller-scoped variable, so this is deliberately a macro; `name` is a declarator, not an
// expression, so it must not be parenthesized.
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage,bugprone-macro-parentheses)
#define VASH_BENCH_TP(name) auto name = std::chrono::steady_clock::now()

// Report milliseconds since `name` was stamped, then re-arm `name` for the next phase. Kept a
// macro so the disabled build (below) can drop `label` unevaluated, avoiding its string work.
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define VASH_BENCH_LAP(label, name)                                                       \
	do {                                                                                  \
		const auto vashBenchNow_ = std::chrono::steady_clock::now();                      \
		std::cerr << "[vash-bench] " << (label) << ": "                                   \
			<< std::chrono::duration<double, std::milli>(vashBenchNow_ - (name)).count()  \
			<< " ms\n";                                                                   \
		(name) = vashBenchNow_;                                                           \
	} while (false)

#else

// No-ops: `name` is never declared, so nothing downstream references it.
// Expand to a discarded expression (not do-while) so the statement-terminating
// semicolon is well-formed without tripping the avoid-do-while lint.
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define VASH_BENCH_TP(name)         ( (void)0 )
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define VASH_BENCH_LAP(label, name) ( (void)0 )

#endif // VASH_BENCHMARK

#endif // VASH_BENCHMARK_HPP
