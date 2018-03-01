#pragma once
#include <sstream>

template <typename T>
std::string strFromNumber(T n) {
  std::string result;
  std::ostringstream convert;
  convert << n;
  result = convert.str();
  return result;
}
