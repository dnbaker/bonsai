#ifndef _LOG_UTIL_H__
#define _LOG_UTIL_H__
#define __STDC_FORMAT_MACROS
#include <cinttypes>
#include <cstdlib>
#include <cstdio>
#include <cstdarg>

#define _FUNCTION_MACRO_ __PRETTY_FUNCTION__
#define LOG_INFO(...) log_info(__func__, ##__VA_ARGS__)
#define LOG_WARNING(...) log_warning(_FUNCTION_MACRO_, ##__VA_ARGS__)
#define LOG_EXIT(...) log_error(_FUNCTION_MACRO_, __LINE__, ##__VA_ARGS__)
#if !NDEBUG
#    define LOG_DEBUG(...) log_debug(_FUNCTION_MACRO_, __LINE__, ##__VA_ARGS__)
#else
#    define LOG_DEBUG(...)
#endif
#define LOG_ASSERT(condition) log_assert(_FUNCTION_MACRO_, __LINE__, ((u64)(condition)), (#condition))

static inline int log_debug(const char *func, int line, const char *fmt, ...) {
    std::va_list args;
    va_start(args, fmt);
    int ret(std::fprintf(stderr, "[D:%s:%d] ", func, line));
    ret += std::vfprintf(stderr, fmt, args);
    va_end(args);
    return ret;
}

static inline int log_warning(const char *func, const char *fmt, ...) {
    std::va_list args;
    va_start(args, fmt);
    int ret(std::fprintf(stderr, "[W:%s] ", func));
    ret += std::vfprintf(stderr, fmt, args);
    va_end(args);
    return ret;
}

static inline int log_info(const char *func, const char *fmt, ...) {
    std::va_list args;
    va_start(args, fmt);
    int ret(std::fprintf(stderr, "[%s] ", func));
    ret += std::vfprintf(stderr, fmt, args);
    va_end(args);
    return ret;
}

static inline void log_error(const char *func, int line, const char *fmt, ...) {
    std::va_list args;
    va_start(args, fmt);
    std::fprintf(stderr, "[E:%s:%d] ", func, line);
    std::vfprintf(stderr, fmt, args);
    va_end(args);
    std::exit(EXIT_FAILURE);
}

static inline void log_assert(const char *func, int line, int assertion, const char *assert_str) {
    if(!assertion) {
        std::fprintf(stderr, "[E:%s:%d] Assertion '%s' failed.\n",
                func, line, assert_str);
        std::exit(EXIT_FAILURE);
    }
}

#endif /* #ifndef _LOG_UTIL_H__ */
