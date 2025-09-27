//
// Created by Christos on 9/26/2025.
//

#ifndef CONVEX_HULL_3D_CUSTOMUNORDEREDSET_H
#define CONVEX_HULL_3D_CUSTOMUNORDEREDSET_H

#endif //CONVEX_HULL_3D_CUSTOMUNORDEREDSET_H
#pragma once
#include <vector>
#include <list>
#include <functional>

template<typename Key,
        typename Hash = std::hash<Key>,
        typename KeyEqual = std::equal_to<Key>>
class CustomUnorderedSet {
private:
    // --- Internal Data ---
    std::vector<std::list<Key>> buckets; // bucket array
    size_t num_elements;
    size_t bucket_count;// number of elements stored
    float max_load_factor;               // threshold for rehashing
    Hash hasher;                         // hash functor
    KeyEqual key_equal;                  // equality functor



    size_t getBucketIndex(const Key &key) const {
        return hasher(key) % bucket_count;
    }

    void rehash(size_t new_bucket_count) {
        bucket_count = new_bucket_count;

        std::vector<std::list<Key>> new_buckets(new_bucket_count);
        for (auto &bucket: buckets) {
            for (auto &key: bucket) {
                int bucket_index = getBucketIndex(key);
                new_buckets[bucket_index].push_back(key);
            }
        }
        buckets = std::move(new_buckets);

    }


public:


    // --- Constructors & Destructor ---
    explicit CustomUnorderedSet(size_t bucket_count = 16,
                            const Hash &hash = Hash(),
                            const KeyEqual &equal = KeyEqual())
            : buckets(bucket_count),
              num_elements(0),
              bucket_count(bucket_count),
              max_load_factor(1.0f),
              hasher(hash),
              key_equal(equal) {}

    ~CustomUnorderedSet() = default;


    bool insert(const Key &key) {
        int bucket_index = getBucketIndex(key);
        for (auto &b_key: buckets[bucket_index]) {
            if (key_equal(key, b_key)) {
                buckets[bucket_index].remove(b_key);
                num_elements--; }
        }
        buckets[bucket_index].push_back(key);
        num_elements++;
        float load = static_cast<float>(num_elements) / bucket_count;

        if (load > max_load_factor) rehash(2 * bucket_count);

        return true;
    }

    void clear() {
        for (auto &bucket: buckets) bucket.clear();
        num_elements = 0;
    }

    bool erase(const Key &key) {
        int bucket_index = getBucketIndex(key);
        for (auto &b_key: buckets[bucket_index]) {
            if (key_equal(key, b_key)) {
                buckets[bucket_index].remove(b_key);
                num_elements--;
                return true;
            }
        }
        return false;
    }

    bool contains(const Key &key) const {
        int bucket_index = getBucketIndex(key);
        for (auto &b_key: buckets[bucket_index]) {
            if (key_equal(key, b_key)) return true;
        }
        return false;
    }

    Key &operator[](const Key &key) {
        size_t bucket_index = getBucketIndex(key);
        for (auto &b_key: buckets[bucket_index]) {
            if (key_equal(key, b_key)) {
                return b_key;
            }
        }
        throw "Key not found";
    }

};