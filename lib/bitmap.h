#ifndef COUNTER_H__
#define COUNTER_H__

namespace emp {


class Counter {

    std::unordered_map<T, std::vector<bool>> map_;
    size_t n_;

    void add(T elem) {
        auto match(map_.find(elem));
        if(match == map_.end()) map_.emplace(elem, 1);
        else ++match->second;
        ++n_;
    }
};

} // namespace emp

#endif  // COUNTER_H__
