#pragma once
#ifndef ACTIVE_MANAGER_HPP
#define ACTIVE_MANAGER_HPP

#include <Vc/Allocator>

#include <vector>
#include <algorithm>
#include <unordered_map>

#include "Assert.hpp"
#include "Definitions.hpp"
#include "Point.hpp"

class ActiveManager {
public:
    std::vector<Node> req_ids;

    template<typename T>
    using VVec = std::vector<T, Vc::Allocator<T> >;

    VVec<Coord_v> req_phi;
    VVec<Coord_v> req_poin_x;
    VVec<Coord_v> req_poin_y;
    VVec<Coord_v> req_poin_r;

    VVec<Coord_m> req_old;

    static_assert(sizeof(Coord_v) == sizeof(Coord_b) * Coord_v::size(),
        "Unexpected padding"
    );

    static constexpr auto Packing = Coord_v::size();

    bool requestAlreadyPending(const Request &req) const;

    void addRequest(const Request &req, bool old = false) {
        _starts.push_back(req);
        std::push_heap(_starts.begin(), _starts.end(), _start_comp);
    }

    template <bool Endgame, typename InsertCallback>
    void update(Coord_b deleteLimit, Coord_b insertLimit, InsertCallback cb = [] (const Request&) {}) {
        if (_verbose) {
            std::cout << "Perform update with del=" << deleteLimit << ", ins=" << insertLimit << "\n";
            std::cout << *this << std::endl;
        }

        checkInvariants();



#ifndef NDEBUG
        // in debug mode we count the number of insertions/deletions
        // and compare them later to the actual size change.
        int64_t expected_size_after = _size;
        {
            for (const auto &r : _starts) {
                expected_size_after += (r.range.first <= insertLimit) && (r.range.second > deleteLimit);
            }

            for (const auto &m : _stops) {
                expected_size_after -= (m.first <= deleteLimit);
            }

            if (_verbose) {
                std::cout << "Current size: " << _size << " Expected size after: " << expected_size_after << std::endl;
            }

            ASSERT_GE(expected_size_after, 0);
        }
#endif

        auto insert = [&] (const unsigned int& pos, const Request& req) {
            if (_verbose)
                std::cout << "Insert " << req << " at pos " << pos << std::endl;

            // insert into map (and reverse map)
            auto res = _map.emplace(req.id, pos);
            if (!res.second) {
                // TODO: Is it worth preventing double insertions?
#ifndef NDEBUG
                expected_size_after--;
#endif
                auto stop = std::find_if(_stops.begin(), _stops.end(), [&req] (const StopMessage& msg) {
                    return msg.second == req.id;
                });

                ASSERT(stop != _stops.end());

                if (stop->first < req.range.second) {
                    std::cout << "WARNING: Increase request size" << std::endl;
                    stop->first = req.range.second;
                    std::make_heap(_stops.begin(), _stops.end(), _stop_comp);
                }


                if (_verbose)
                    std::cout << " ... skipped due to duplicate" << std::endl;
                return false;
            }

            // copy data
            req_ids[pos] = req.id;
            _req_phi_data[pos] = req.phi;
            _req_poin_x_data[pos] = req.poinX;
            _req_poin_y_data[pos] = req.poinY;
            #ifdef LOG_TRANSFORM
            _req_poin_r_data[pos] = req.poinLogInvLen;
            #else
            _req_poin_r_data[pos] = req.poinInvLen;
            #endif

            if (Endgame) {
                unsigned int i = pos / Packing;
                unsigned int j = pos % Packing;
                req_old[i][j] = req.old();
            }

            _stops.emplace_back(req.range.second, req.id);
            std::push_heap(_stops .begin(), _stops .end(), _stop_comp);

            return true;
        };

        auto shift_remove = [&] (const unsigned int target, const unsigned source) {
            if (_verbose)
                std::cout << "Shift remove " << source << " -> " << target << std::endl;

            if (target != source) {
                req_ids[target] = req_ids[source];
                _map[req_ids[target]] = target;

                _req_phi_data[target] = _req_phi_data[source];
                _req_poin_x_data[target] = _req_poin_x_data[source];
                _req_poin_y_data[target] = _req_poin_y_data[source];
                _req_poin_r_data[target] = _req_poin_r_data[source];

                if (Endgame) {
                    unsigned int ti = target / Packing;
                    unsigned int tj = target % Packing;
                    unsigned int si = source / Packing;
                    unsigned int sj = source % Packing;

                    req_old[ti][tj] = req_old[si][sj];
                }
            }

            // make it impossible to connect to this point
            _req_poin_r_data[source] = std::numeric_limits<Coord_b>::max();
#ifndef NDEBUG
            req_ids[source] = std::numeric_limits<Node>::max();
#endif
        };

        auto popObsoleteInserts = [&] (bool removeFirst) {
            while (!_starts.empty() && (removeFirst ||  _starts.front().range.second <= deleteLimit)) {
                cb(_starts.front());
                std::pop_heap(_starts.begin(), _starts.end(), _start_comp);
                _starts.pop_back();
                removeFirst = false;
            }
        };

        auto tryInsert = [&] (unsigned int pos) {
            while(!_starts.empty() && _starts.front().range.first <= insertLimit) {
                bool success = insert(pos, _starts.front());

                if (!success) {
                    // insertion was not sucessful -- undo all changes
                    // and silently drop request
                    std::pop_heap(_starts.begin(), _starts.end(), _start_comp);
                    _starts.pop_back();
                }

                popObsoleteInserts(success);

                if (success)
                    return true;
            }

            return false;
        };

        popObsoleteInserts(false);

        while (!_stops.empty() && _stops.front().first <= deleteLimit) {
            if (_verbose)
                std::cout << "Remove req id " << _stops.front().second << " at position " << _stops.front().first << std::endl;

            // fetch and remove old position
            const auto it = _map.find(_stops.front().second);
            ASSERT (it != _map.cend());

            auto oldPos = it->second;
            _map.erase(it);

            // TODO: Handle overlaps in replaceRemove
            if (!tryInsert(oldPos)) {
                ASSERT_GT(_size, 0);
                shift_remove(oldPos, --_size);
            }

            std::pop_heap(_stops.begin(), _stops.end(), _stop_comp);
            _stops.pop_back();
        }


        while(!_starts.empty() && _starts.front().range.first <= insertLimit) {
            _size++;
            _update_size();

            bool success = tryInsert(_size - 1);
            if (!success) {
                // insertion was not sucessful -- undo all changes
                // and silently drop request
                _size--;
            }
        }

        _end = (_size + Packing - 1) / Packing;

        if (_verbose) {
            std::cout << "Performed update with del=" << deleteLimit << ", ins=" << insertLimit << std::endl;
            for(const auto& i : _map) {
                std::cout << "End: Active id:" << i.first << " @ " << i.second << std::endl;
            }
        }

        for(unsigned int i=_size; i < _end * Packing; ++i) {
            _req_poin_r_data[i] = std::numeric_limits<Coord_b>::max();
#ifndef NDEBUG
            req_ids[i] = std::numeric_limits<Node>::max();
#endif
        }

#ifndef NDEBUG
        ASSERT_EQ(_size, expected_size_after);
#endif

        checkInvariants();
    }

    bool empty() const {
        return !_size && _starts.empty();
    }

    Coord_b nextInsertionPending() const {
        if (_starts.empty())
            return std::numeric_limits<Coord_b>::max();

        return _starts.front().range.first;
    }

    void copyFrom(const ActiveManager& src,
                  bool markOld = false,
                  Coord_b thresh = std::numeric_limits<Coord_b>::max()
    );

    const size_t& end() const {
        return _end;
    }

    Coord_b maxRange() const {
        if (unlikely(_stops.empty()))
            return std::numeric_limits<Coord_b>::min();
        return _stops.back().first;
    }

    Coord_b maxReplayRange() const;

    friend std::ostream& operator<<(std::ostream& stream, const ActiveManager& o);

    // in case a stop message was "behind" the last band, i.e. phi-stop > 2PI
    // we move it towards the front
    void fixRange(Coord_b maxRange);

#ifndef NDEBUG
    void checkInvariants() const;
#else
    void checkInvariants() const {}
#endif


protected:
    static constexpr bool _verbose {VERBOSITY(true)};


    using StopMessage = std::pair<Coord_b, Node>;

    std::greater<StopMessage> _stop_comp;
    std::vector<StopMessage> _stops;
    std::greater<Request> _start_comp;
    std::vector<Request> _starts;

    std::unordered_map<Node, unsigned int> _map;

    size_t _size{0};
    size_t _end{0};

    Coord_b _last_delete {0.0};
    Coord_b _last_insert {0.0};


    Coord_b* _req_phi_data;
    Coord_b* _req_poin_x_data;
    Coord_b* _req_poin_y_data;
    Coord_b* _req_poin_r_data;

    void _update_size() {
        if (_size < req_ids.size())
            return;

        const size_t vsize = (2*_size + Packing - 1) / Packing;
        req_phi.resize(vsize);
        req_poin_x.resize(vsize);
        req_poin_y.resize(vsize);
        req_poin_r.resize(vsize);
        req_old.resize(vsize);

        req_ids.resize(vsize * Packing);

        _req_phi_data    = reinterpret_cast<Coord_b*>(req_phi.data());
        _req_poin_x_data = reinterpret_cast<Coord_b*>(req_poin_x.data());
        _req_poin_y_data = reinterpret_cast<Coord_b*>(req_poin_y.data());
        _req_poin_r_data = reinterpret_cast<Coord_b*>(req_poin_r.data());
    }
};
#endif