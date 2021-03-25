#ifndef __MEM_MANAGER_HPP__
#define __MEM_MANAGER_HPP__

#include "chromosome.hpp"
#include "util.hpp"
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <typeindex>
#include <typeinfo>
#include <vector>

/*! The memory manager
 */

namespace MemoryManager
{

typedef std::pair<bool, void *> space_entry;

std::mutex memory_map_mutex;
std::map<std::string, std::vector<space_entry> *> memory_map;

template <class T> std::vector<T> &getMemorySpace(int length)
{
    // Must be thread safe
    memory_map_mutex.lock();

    std::string type_name(typeid(T).name());

    // Allocate vector if does not exist
    if (memory_map.find(type_name) == memory_map.end())
    {
        std::pair<std::string, std::vector<space_entry> *> key_value;

        key_value.first  = type_name;
        key_value.second = new std::vector<space_entry>();

        memory_map.insert(key_value);
    }

    // Searches for best already allocated memory space (may be made faster by
    // ordering vectors by reverse size)
    space_entry *best_match = nullptr;

    for (auto &allocation : *memory_map.at(type_name))
    {
        std::vector<T> *second = (std::vector<T> *)allocation.second;

        if (allocation.first && second->size() >= length &&
            second->size() < length * 1.3)
        { // Check if available and has at least size length and not many more
          // items than required
            if (!best_match ||
                ((std::vector<T> *)best_match->second)->size() > second->size())
                best_match = &allocation;
        }
    }

    // Allocates new vector if no fitting spaces were found
    if (!best_match)
    {
        auto &v = *memory_map.at(type_name);

        auto vec = new std::vector<T>(length, 0);

        space_entry pair;
        pair.first  = false;
        pair.second = vec;

        v.push_back(pair);

        best_match = &v.back();
    }

    memory_map_mutex.unlock();

    return *((std::vector<T> *)best_match->second);
}

template <class T> void freeMemorySpace(std::vector<T> &v)
{
    std::string type_name(typeid(T).name());
    if (memory_map.find(type_name) == memory_map.end()) return;

    for (auto &allocation : *memory_map.at(std::string(typeid(T).name())))
    {
        if (&v == (std::vector<T> *)allocation.second)
        {
            memory_map_mutex.lock();
            allocation.first = true;
            memory_map_mutex.unlock();
            break;
        }
    }
}
} // namespace MemoryManager

#endif