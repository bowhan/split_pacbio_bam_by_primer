#pragma once

#include <vector>
#include <list>
#include <deque>
#include <pbbam/BamReader.h>
#include <pbbam/BamRecord.h>

template <class T>
struct LinearContainerPolicies;

template <class... Args>
struct LinearContainerPolicies<std::vector<Args...> > {

    using container_type = std::vector<Args...>;
    using value_type = typename container_type::value_type;
    using reference = typename container_type::reference;
    using iterator = typename container_type::iterator;
    using const_iterator = typename container_type::const_iterator;
    using size_type = typename container_type::size_type;

    static void Reserve(container_type& c, size_type s) {
        c.reserve(s);
    }

    static size_type Capacity(const container_type& c) {
        return c.capacity();
    }

    static value_type& At(const container_type& c, size_type idx) {
        return c[idx];
    }

//    static void Push(container_type& c, const value_type& val) {
//        c.push_back(val);
//    }

    static void Push(container_type& c, value_type&& val) {
        c.push_back(std::forward<value_type>(val));
    }

    template <class... TT>
    static void Push(container_type& c, TT&& ... args) {
        c.emplace_bacK(std::forward<TT>(args)...);
    }

};

// TODO: add policies for other containers

template <class Format>
struct DataProducerPolicy;

template < >
struct DataProducerPolicy<PacBio::BAM::BamRecord> {
    using source_type = PacBio::BAM::BamReader;
    using data_type = PacBio::BAM::BamRecord;

    template <class... Args>
    std::pair<data_type, bool> Produce(source_type& s, Args&& ... args) {
        data_type d(std::forward<Args>(args)...);
        auto success = s.GetNext(d);
        return std::make_pair(d, success);
    }

};