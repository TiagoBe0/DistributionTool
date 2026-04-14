#pragma once
#include <cmath>
#include <cstdio>
#include <exception>
#include <functional>
#include <string>
#include <vector>

// ── Minimal self-registering test harness (no external dependencies) ─────────
//
// Usage:
//   TEST_CASE("description") { REQUIRE(x == y); APPROX_EQ(a, b, 1e-9); }
//   int main() { return Test::run(); }

namespace Test {

struct Failure {
    std::string msg;
    const char* file;
    int         line;
};

struct Case {
    const char*           name;
    std::function<void()> fn;
};

inline std::vector<Case>& registry() {
    static std::vector<Case> v;
    return v;
}

struct Registrar {
    Registrar(const char* name, std::function<void()> fn) {
        registry().push_back({name, fn});
    }
};

inline int run() {
    int passed = 0, failed = 0;
    for (const auto& c : registry()) {
        try {
            c.fn();
            std::printf("  [ OK ] %s\n", c.name);
            ++passed;
        } catch (const Failure& f) {
            std::printf("  [FAIL] %s\n         %s:%d  %s\n",
                        c.name, f.file, f.line, f.msg.c_str());
            ++failed;
        } catch (const std::exception& e) {
            std::printf("  [FAIL] %s\n         exception: %s\n",
                        c.name, e.what());
            ++failed;
        }
    }
    std::printf("\n%d passed, %d failed\n", passed, failed);
    return (failed > 0) ? 1 : 0;
}

} // namespace Test

// ── Macro helpers ─────────────────────────────────────────────────────────────
// Use __LINE__ so all three references inside TEST_CASE resolve to the same name.
// Two TEST_CASEs in different files can share a line number without conflict
// because the functions/registrars are file-static (internal linkage).
#define _TC_CAT_(a, b) a##b
#define _TC_CAT(a, b)  _TC_CAT_(a, b)

#define TEST_CASE(name)                                                         \
    static void _TC_CAT(_tfn_, __LINE__)();                                     \
    static ::Test::Registrar _TC_CAT(_treg_, __LINE__)(name,                   \
                                     _TC_CAT(_tfn_, __LINE__));                 \
    static void _TC_CAT(_tfn_, __LINE__)()

#define REQUIRE(expr)                                                  \
    do {                                                               \
        if (!(expr)) {                                                 \
            char _b[512];                                              \
            std::snprintf(_b, sizeof(_b), "REQUIRE(%s)", #expr);      \
            throw ::Test::Failure{_b, __FILE__, __LINE__};            \
        }                                                              \
    } while (0)

#define CHECK_THROWS(expr)                                             \
    do {                                                               \
        bool _threw = false;                                           \
        try { (void)(expr); } catch (...) { _threw = true; }          \
        if (!_threw) {                                                 \
            char _b[512];                                              \
            std::snprintf(_b, sizeof(_b),                             \
                          "expected exception from: %s", #expr);      \
            throw ::Test::Failure{_b, __FILE__, __LINE__};            \
        }                                                              \
    } while (0)

// Relative epsilon: |a-b| <= eps * max(1, |a|, |b|)
inline bool _approx(double a, double b, double eps) {
    using std::abs;
    return abs(a - b) <= eps * std::max({1.0, abs(a), abs(b)});
}

#define APPROX_EQ(a, b, eps)                                           \
    do {                                                               \
        double _a=(a), _b=(b), _e=(eps);                              \
        if (!_approx(_a, _b, _e)) {                                   \
            char _buf[512];                                            \
            std::snprintf(_buf, sizeof(_buf),                         \
                "APPROX_EQ(%s, %s): got %.10g vs %.10g (eps=%.3g)",   \
                #a, #b, _a, _b, _e);                                  \
            throw ::Test::Failure{_buf, __FILE__, __LINE__};          \
        }                                                              \
    } while (0)
