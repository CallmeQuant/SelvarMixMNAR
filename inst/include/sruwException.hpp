#ifndef SRUW_EXCEPTION_H
#define SRUW_EXCEPTION_H

#include <stdexcept>
#include <string>

class SRUWStatusException : public std::runtime_error {
  std::string errorClass_;

public:
  SRUWStatusException(const std::string& message,
                      const std::string& errorClass)
      : std::runtime_error(message), errorClass_(errorClass) {}

  const std::string& errorClass() const { return errorClass_; }
};

#endif
