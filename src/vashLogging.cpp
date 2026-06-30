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

/// Log events from vash program execution
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2026 Anthony J. Greenberg
 * \version 0.5
 *
 * Message logging class implementation.
 *
 */

#include <ctime>
#include <iomanip>  // for put_time
#include <chrono>
#include <string>
#include <sstream>

#include "vashLogging.hpp"

using namespace BayesicSpace;

VashLog::VashLog(const LogFileNameWithMessage &lfNameWithMessage) : logFile_(lfNameWithMessage.logFileName, std::ios::trunc | std::ios::out), toSave_{true} {
	std::stringstream logStream;
	const time_t startTime = std::time(nullptr);
	struct tm buf{};
	logStream << std::put_time(localtime_r(&startTime, &buf), "%b %e %Y %H:%M %Z");
	messages_ = lfNameWithMessage.initialMessage + " started on " + logStream.str() + "\n";
	logStream.clear();
	startTime_ = std::chrono::steady_clock::now();
};

VashLog::~VashLog() {
	if (toSave_) {
		// might fail due to external factors (e.g. full disk)
		logFile_ << messages_;
	}
}

VashLog::VashLog(VashLog &&toMove) noexcept {
	*this = std::move(toMove);
}

VashLog& VashLog::operator=(VashLog &&toMove) noexcept {
	if (this != &toMove) {
		messages_      = std::move(toMove.messages_);
		logFile_       = std::move(toMove.logFile_);
		toSave_        = toMove.toSave_;
		startTime_     = toMove.startTime_;
		toMove.toSave_ = false;
	}
	return *this;
}

void VashLog::add(std::string entry) {
	constexpr auto twoDigitThreshold{10}; // pad the seconds remainder to two digits below this
	const auto elapsed{std::chrono::steady_clock::now() - startTime_};
	const auto minutes{std::chrono::duration_cast<std::chrono::minutes>(elapsed)};
	const auto seconds{std::chrono::duration_cast<std::chrono::seconds>(elapsed - minutes)};
	messages_ += "[" + std::to_string( minutes.count() ) + ':'
		+ (seconds.count() < twoDigitThreshold ? "0" : "") + std::to_string( seconds.count() )
		+ "] " + std::move(entry) + "\n";
}
