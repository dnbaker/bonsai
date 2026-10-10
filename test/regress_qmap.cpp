// Sliding-window minimizer queue: the monotone QueueMap must reproduce the
// lexicographic minimum (score, then element) of every window exactly,
// including across the ring-buffer wrap, with duplicate elements and with
// context-dependent scores (the same element carrying different scores,
// where the previous tree-based implementation miscounted entries).
#include "test/catch.hpp"

#include "qmap.h"
#include "hll/include/circularqueue/cq.h"
#include "util.h"
#include "aesctr/wy.h"

#include <deque>
#include <vector>

using namespace bns;

namespace {

// Brute-force reference: recompute the window minimum on every slide.
template<typename T, typename ScoreType>
struct BruteWindow {
    std::deque<ElScore<T, ScoreType>> window;
    size_t wsz;
    explicit BruteWindow(size_t w): wsz(w) {}
    T next_value(const T el, const ScoreType score) {
        window.emplace_back(el, score);
        if(window.size() > wsz) window.pop_front();
        if(window.size() < wsz) return std::numeric_limits<T>::max();
        auto best(window.front());
        for(const auto &p: window) if(p < best) best = p;
        return best.el_;
    }
};

} // namespace

TEST_CASE("QueueMap matches brute-force window minimum across ring wrap") {
    wy::WyRand<u64> rng(9182);
    for(const unsigned wsz: {1u, 2u, 3u, 17u, 60u, 128u}) {
        QueueMap<u64, u64> q(wsz);
        BruteWindow<u64, u64> ref(wsz);
        u64 state(rng());
        for(unsigned i(0); i < 5000; ++i) {
            state = state * 6364136223846793005ULL + 1442695040888963407ULL;
            const u64 k(state >> 50); // small range: heavy ties and duplicates
            REQUIRE(q.next_value(k, k) == ref.next_value(k, k));
        }
    }
}

TEST_CASE("QueueMap handles context-dependent scores for the same element") {
    // The same k-mer with different scores in the window must be selected by
    // (score, element), not collapse by element identity.
    QueueMap<u64, u64> q(3);
    BruteWindow<u64, u64> ref(3);
    const u64 k(42);
    for(const u64 score: {9u, 1u, 5u, 5u, 1u, 9u}) {
        REQUIRE(q.next_value(k, score) == ref.next_value(k, score));
    }
}

TEST_CASE("circ::deque pop_back wraps past index zero") {
    // Regression: pop_back decremented stop_ without masking it, so popping
    // at stop_ == 0 underflowed the unsigned index. Monotone window
    // maintenance is the first pop_back user. The pushes below wrap stop_
    // past zero; the pop_back then lands on stop_ == 0.
    circ::deque<elscore_t, u32> d(2); // capacity 4
    d.emplace_back(elscore_t(1, 1));
    d.emplace_back(elscore_t(2, 2));
    d.emplace_back(elscore_t(3, 3));
    d.pop_front();
    d.pop_front();   // start advances to 2; live: {3}
    d.emplace_back(elscore_t(4, 4)); // stop_ wraps to 0; live: {3, 4}
    d.pop_back();    // old code: stop_ underflows here
    d.emplace_back(elscore_t(5, 5));
    REQUIRE(d.size() == 2);
    REQUIRE(d.front().el_ == 3);
    REQUIRE(d.back().el_ == 5);
}
