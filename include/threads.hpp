#pragma once

#include <mutex>
#include <string>
#include "policies.hpp"

template <template <class...> class Container, class T>
class MultiThreadSafeQueue
  : private LinearContainerPolicies<Container<T>>
    , private DataProducerPolicy<T> {
private:
    using self_ = MultiThreadSafeQueue<Container, T>;

    using container_policies = LinearContainerPolicies<Container<T>>;
    using container_policies::Reserve;
    using container_policies::Push;
    using container_policies::Capacity;
    using producer_policies = DataProducerPolicy<T>;

    using source_type = typename producer_policies::source_type;
    using producer_policies::Produce;

public:
    using container_type = typename container_policies::container_type;
    using size_type = typename container_policies::size_type;
    using value_type = typename container_policies::value_type;

private:
    source_type& source_;
    container_type data_;
    std::mutex mx_;
    size_type size_;
    size_type capacity_;
public:
    explicit MultiThreadSafeQueue(source_type& s, size_type cap)
      : source_(s)
        , size_(0)
        , capacity_(cap) {
        Reserve(data_, capacity_);
    }

    template <class... Args>
    void Fill(Args&& ... args);

    container_type Pop();

    template <class... Args>
    container_type FillAndPop(Args&& ... args);
};

template <template <class...> class Container, class T>
template <class... Args>
void MultiThreadSafeQueue<Container, T>::Fill(Args&& ... args) {
    std::lock_guard<std::mutex> lock(mx_);
    while(size_ < capacity_) {
        auto r = Produce(source_, std::forward<Args>(args)...);
        if(!r.second) return;
        Push(data_, std::move(r.first));
        ++size_;
    }
};

template <template <class...> class Container, class T>
auto MultiThreadSafeQueue<Container, T>::Pop() -> container_type {
    std::lock_guard<std::mutex> lock(mx_);
    container_type newdata;
    newdata.reserve(Capacity(data_));
    data_.swap(newdata);
    size_ = 0;
    return newdata;
};


template <template <class...> class Container, class T>
template <class... Args>
auto MultiThreadSafeQueue<Container, T>::FillAndPop(Args&& ... args) -> container_type {
    std::lock_guard<std::mutex> lock(mx_);
    while(size_ < capacity_) {
        auto r = Produce(source_, std::forward<Args>(args)...);
        if(!r.second) break;
        Push(data_, std::move(r.first));
        ++size_;
    }
    container_type newdata;
    newdata.reserve(Capacity(data_));
    data_.swap(newdata);
    size_ = 0;
    return newdata;
};
