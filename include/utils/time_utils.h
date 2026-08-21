#pragma once

#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

namespace june {

// =============================================================================
// Time Utilities
// =============================================================================

// Parse "HH:MM" string to minutes since midnight
inline int parseTimeToMinutes(const std::string& time_str) {
  int hours, minutes;
  char colon;
  std::istringstream ss(time_str);

  if (!(ss >> hours >> colon >> minutes) || colon != ':') {
    throw std::runtime_error("Invalid time format: " + time_str +
                             " (expected HH:MM)");
  }

  if (hours < 0 || hours > 23 || minutes < 0 || minutes > 59) {
    throw std::runtime_error("Invalid time values: " + time_str);
  }

  return hours * 60 + minutes;
}

// Calculate duration in hours between two times
// Handles wraparound (e.g., 21:00 to 08:00 = 11 hours)
inline double calculateDuration(const std::string& start,
                                const std::string& end) {
  int start_min = parseTimeToMinutes(start);
  int end_min = parseTimeToMinutes(end);

  int duration_min;
  if (end_min >= start_min) {
    duration_min = end_min - start_min;
  } else {
    // Wraparound case (crosses midnight)
    duration_min = (24 * 60 - start_min) + end_min;
  }

  return duration_min / 60.0;  // Convert to hours
}

// Parse date string "YYYY-MM-DD" to std::tm
inline std::tm parseDate(const std::string& date_str) {
  std::tm tm = {};
  std::istringstream ss(date_str);
  ss >> std::get_time(&tm, "%Y-%m-%d");

  if (ss.fail()) {
    throw std::runtime_error("Invalid date format: " + date_str +
                             " (expected YYYY-MM-DD)");
  }

  return tm;
}

// Convert a calendar date to a Julian Day Number. Works for any year,
// including pre-1970 dates that mktime cannot handle.
inline long long toJulianDay(int year, int month, int day) {
  long long a = (14 - month) / 12;
  long long y = year + 4800 - a;
  long long m = month + 12 * a - 3;
  return day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045;
}

inline long long tmToJulianDay(const std::tm& date) {
  return toJulianDay(date.tm_year + 1900, date.tm_mon + 1, date.tm_mday);
}

// Convert a Julian Day Number back to a std::tm (year/month/day only)
inline std::tm julianDayToTm(long long jd) {
  long long a = jd + 32044;
  long long b = (4 * a + 3) / 146097;
  long long c = a - (146097 * b) / 4;
  long long d = (4 * c + 3) / 1461;
  long long e = c - (1461 * d) / 4;
  long long m = (5 * e + 2) / 153;
  std::tm result = {};
  result.tm_mday = static_cast<int>(e - (153 * m + 2) / 5 + 1);
  result.tm_mon = static_cast<int>(m + 3 - 12 * (m / 10)) - 1;
  result.tm_year = static_cast<int>(100 * b + d - 4800 + m / 10) - 1900;
  return result;
}

inline bool isPreEpoch(const std::tm& date) {
  return date.tm_year + 1900 < 1970;
}

// Get day of week (0 = Monday, 6 = Sunday)
inline int getDayOfWeek(const std::tm& date) {
  if (isPreEpoch(date)) {
    // JD 0 is a Monday, so JD % 7 gives 0=Monday ... 6=Sunday
    return static_cast<int>(tmToJulianDay(date) % 7);
  }
  std::tm temp = date;
  std::mktime(&temp);  // Normalize the tm structure
  // tm_wday: 0 = Sunday, 1 = Monday, ..., 6 = Saturday
  // Convert to: 0 = Monday, ..., 6 = Sunday
  return (temp.tm_wday + 6) % 7;
}

// Check if day is weekend
inline bool isWeekend(int day_of_week) {
  return day_of_week == 5 || day_of_week == 6;  // Saturday or Sunday
}

// Add days to a date
inline std::tm addDays(const std::tm& date, int days) {
  if (isPreEpoch(date) ||
      isPreEpoch(julianDayToTm(tmToJulianDay(date) + days))) {
    return julianDayToTm(tmToJulianDay(date) + days);
  }
  std::tm result = date;
  result.tm_mday += days;
  std::mktime(&result);  // Normalize
  return result;
}

// Format date as "YYYY-MM-DD"
inline std::string formatDate(const std::tm& date) {
  if (isPreEpoch(date)) {
    // Sized for the widest output three ints can produce, so an out-of-range
    // tm prints as visibly wrong rather than being silently truncated
    char buffer[36];
    std::snprintf(buffer, sizeof(buffer), "%04d-%02d-%02d", date.tm_year + 1900,
                  date.tm_mon + 1, date.tm_mday);
    return std::string(buffer);
  }
  char buffer[11];
  std::strftime(buffer, sizeof(buffer), "%Y-%m-%d", &date);
  return std::string(buffer);
}

// Calculate number of days between two dates
inline int daysBetween(const std::tm& start, const std::tm& end) {
  if (isPreEpoch(start) || isPreEpoch(end)) {
    return static_cast<int>(tmToJulianDay(end) - tmToJulianDay(start));
  }
  std::tm start_copy = start;
  std::tm end_copy = end;
  std::time_t start_time = std::mktime(&start_copy);
  std::time_t end_time = std::mktime(&end_copy);
  double diff_seconds = std::difftime(end_time, start_time);
  return static_cast<int>(diff_seconds / (60 * 60 * 24));
}

class Timer {
 public:
  static std::string timestamp() {
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << "[" << std::put_time(std::localtime(&in_time_t), "%H:%M:%S") << "] ";
    return ss.str();
  }
};

}  // namespace june
