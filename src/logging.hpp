#ifndef LOGGING_HPP
#define LOGGING_HPP

#include "version.hpp"
#include <ctime>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace call_chets {

inline void printHeader(const std::string &toolName,
                        const std::string &description,
                        const std::map<std::string, std::string> &files,
                        const std::map<std::string, std::string> &parameters) {
  std::time_t now = std::time(nullptr);
  char timestr[100];
  std::strftime(timestr, sizeof(timestr), "%d/%m/%Y - %H:%M:%S",
                std::localtime(&now));

  std::cerr << "\n[" << toolName << "] " << description
            << "\n  * Author        : Frederik Lassen, University of Oxford"
            << "\n  * Contact       : frederik.heymann@gmail.com"
            << "\n  * Version       : " << getFullVersion()
            << "\n  * Run date      : " << timestr << "\n"
            << std::endl;

  if (!files.empty()) {
    std::cerr << "Files:" << std::endl;
    for (const auto &pair : files) {
      if (!pair.second.empty()) {
        std::cerr << "  * " << pair.first;
        // Padding for alignment could be improved but simple is fine
        std::string padding(14 - pair.first.length(), ' ');
        std::cerr << padding << ": [" << pair.second << "]" << std::endl;
      }
    }
    std::cerr << std::endl;
  }

  if (!parameters.empty()) {
    std::cerr << "Parameters:" << std::endl;
    for (const auto &pair : parameters) {
      std::cerr << "  * " << pair.first;
      std::string padding(14 - pair.first.length(), ' ');
      std::cerr << padding << ": " << pair.second << std::endl;
    }
    // Don't add blank line here - let caller control it
  }
}

inline void log(const std::string &message) {
  std::cerr << "[LOG] " << message << std::endl;
}

inline void logProgress(size_t current, size_t total, size_t step = 1000) {
  if (total == 0)
    return;
  if (current % step == 0) {
    float percent = (float)current / total * 100.0;
    std::cerr << "\r[PROGRESS] Processing... " << percent << "%" << std::flush;
  }
}

} // namespace call_chets

#endif // LOGGING_HPP
