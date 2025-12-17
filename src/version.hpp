#ifndef VERSION_HPP
#define VERSION_HPP

#include <string>

// Default values if macros are not defined
#ifndef GIT_COMMIT
#define GIT_COMMIT "unknown"
#endif

#ifndef GIT_DATE
#define GIT_DATE "unknown"
#endif

#ifndef VERSION
#define VERSION "0.4.0"
#endif

inline std::string getVersion() { return VERSION; }

inline std::string getGitCommit() { return GIT_COMMIT; }

inline std::string getGitDate() { return GIT_DATE; }

inline std::string getFullVersion() {
  return std::string(VERSION) + " / commit = " + GIT_COMMIT +
         " / release = " + GIT_DATE;
}

#endif // VERSION_HPP
