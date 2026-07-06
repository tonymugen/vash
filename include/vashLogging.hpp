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
 * \version 0.6
 *
 * Definitions and interface documentation for the message logging class.
 *
 */

#pragma once

#include <fstream>
#include <string>
#include <chrono>

namespace BayesicSpace {
	/** \brief Log file name and start message */
	struct LogFileNameWithMessage {
		std::string logFileName;
		std::string initialMessage;
	};
	/** \brief Execution message logging
	 *
	 * Time stamps are added to the logged messages as elapsed minutes:seconds.
	 */
	class VashLog {
	public:
		/** \brief Default constructor */
		VashLog() = default;
		/** \brief Constructor with log file name
		 * 
		 * \param[in] lfNameWithMessage log file name with start message
		 */
		VashLog(const LogFileNameWithMessage &lfNameWithMessage);
		/** \brief Copy constructor
		 *
		 * \param[in] toCopy object to copy
		 */
		VashLog(const VashLog &toCopy) = delete;
		/** \brief Copy assignment operator
		 *
		 * \param[in] toCopy object to copy
		 * \return `VashLog` object
		 */
		VashLog& operator=(const VashLog &toCopy) = delete;
		/** \brief Move constructor
		 *
		 * Leaves the moved-from object with saving disabled.
		 *
		 * \param[in] toMove object to move
		 */
		VashLog(VashLog &&toMove) noexcept;
		/** \brief Move assignment operator
		 *
		 * Leaves the moved-from object with saving disabled.
		 *
		 * \param[in] toMove object to move
		 * \return `VashLog` object
		 */
		VashLog& operator=(VashLog &&toMove) noexcept;
		/** \brief Destructor
		 * 
		 * Possibly saves the log to a file.
		 */
		~VashLog();

		/** \brief Adds an entry to the log
		 *
		 * Elapsed time from logging start is prepended as a bracketed [m:ss] stamp.
		 * New line is added automatically to the end.
		 * The operation is not thread-safe, call only from the orchestrating thread.
		 *
		 * \param[in] entry log entry
		 */
		void add(std::string entry);
	private:
		/** \brief Message content */
		std::string messages_;
		/** \brief Log file */
		std::fstream logFile_;
		/** \brief Log saving flag */
		bool toSave_{false};
		/** \brief Start time */
		std::chrono::steady_clock::time_point startTime_;
	};
}
