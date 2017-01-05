#include "ActiveManager.hpp"
#include "Assert.hpp"

void ActiveManager::copyFrom(const ActiveManager& src, bool markOld, Coord_b thresh){
    if (markOld) {
        for (auto msg : src._starts) {
            if (msg.range.first < thresh) {
                msg.setOld();
                _starts.push_back(msg);
            }
        }
    } else {
        for (const auto& msg : src._starts) {
            if (msg.range.first < thresh) {
                _starts.push_back(msg);
            }
        }
    }

    for(const auto msg : src._stops) {
        if (msg.first < thresh)
            _stops.push_back(msg);
    }

    std::make_heap(_starts.begin(), _starts.end(), _start_comp);
    std::make_heap(_stops .begin(), _stops .end(), _stop_comp);


    if (src._size) {
        ASSERT_EQ(_size, 0);

        _map = src._map;

        _size = src._size;
        _end = src._end;
        _update_size();

        std::copy(src.req_phi.cbegin(),    src.req_phi.cbegin() + _end,    req_phi.begin());
        std::copy(src.req_poin_x.cbegin(), src.req_poin_x.cbegin() + _end, req_poin_x.begin());
        std::copy(src.req_poin_y.cbegin(), src.req_poin_y.cbegin() + _end, req_poin_y.begin());
        std::copy(src.req_poin_r.cbegin(), src.req_poin_r.cbegin() + _end, req_poin_r.begin());
        std::copy(src.req_ids.cbegin(),    src.req_ids.cbegin() + _size,   req_ids.begin());

        if (markOld) {
            Vc::Mask<Coord_b> tm(true);
            for(unsigned int i=0; i<_end; i++)
                req_old[i] = tm;
        }
    }
}

void ActiveManager::fixRange(Coord_b maxRange) {
    {
        bool fixed = false;
        for (StopMessage &msg : _stops) {
            if (msg.first >= maxRange) {
                msg.first -= 2 * M_PI;
                fixed = true;
            }
        }

        if (fixed) {
            std::make_heap(_stops.begin(), _stops.end(), _stop_comp);
        }
    }

    {
        bool fixed = false;
        for (auto &req : _starts) {
            if (req.phi >= maxRange) {
                req.phi -= 2 * M_PI;
            }

            if (req.range.first >= maxRange) {
                req.range.first -= 2 * M_PI;
                fixed = true;
            }

            if (req.range.second >= maxRange) {
                req.range.second -= 2 * M_PI;
            }
        }

        if (fixed) {
            std::make_heap(_starts.begin(), _starts.end(), _start_comp);
        }
    }


    checkInvariants();
}

#ifndef NDEBUG
void ActiveManager::checkInvariants() const {
    // maintain heap structure on both streams
    ASSERT(std::is_heap(_stops.cbegin(), _stops.cend(), _stop_comp));
    ASSERT(std::is_heap(_starts.cbegin(), _starts.cend(), _start_comp));

    // every active entry has a stop/book-keeping entry
    ASSERT_EQ(_stops.size(), _size);
    ASSERT_GE(_map.size(), _size);

    // properly invalidated
    for(unsigned int i=0; i < Packing*_end; ++i) {
        if (i < _size) {
            ASSERT_LS(_req_poin_r_data[i], std::numeric_limits<Coord_b>::max());
        } else {
            ASSERT_EQ(_req_poin_r_data[i], std::numeric_limits<Coord_b>::max());
        }
    }

    // 1:1 mapping
    {
        std::vector<bool> assigned(_size, false);
        for(const auto& it : _map) {
            ASSERT(!assigned.at(it.second));
            assigned.at(it.second) = true;
        }
    }
}
#endif

std::ostream& operator<<(std::ostream& os, const ActiveManager& am) {
    auto start(am._starts);
    auto stop(am._stops);

    std::sort_heap(start.begin(), start.end(), am._start_comp);
    std::sort_heap(stop.begin(), stop.end(), am._stop_comp);

    for(auto it = start.crbegin(); it != start.crend(); ++it) {
        std::cout << "Start: " << *it << std::endl;
    }

    for(auto it = stop.crbegin(); it != stop.crend(); ++it) {
        std::cout << "Stop: id:" << it->second << " due @ " << it->first << "\n";
    }

    os << "Active: [";

    std::string tmp;
    for(auto it = am._map.cbegin(); it != am._map.cend(); ++it) {
        tmp = "id:" + std::to_string(it->first) + "@" + std::to_string(it->second)
              + (tmp.empty() ? "" : ", ") + tmp;
    }
    os << tmp << "]" << std::endl;

    return os;
}

bool ActiveManager::requestAlreadyPending(const Request &req) const {
    // pending to be inserted
    for(const auto &r : _starts) {
        if (r.id == req.id && (
                (r.range.first <= req.range.first && r.range.second >= req.range.first) ||
                (r.range.first <= req.range.first && r.range.second >= req.range.first))) {
           return true;
        }
    }

    // already inserted
    const auto it = std::find_if(_stops.cbegin(), _stops.cend(), [&req] (const StopMessage& sm) {
        return sm.second == req.id;
    });
    if (it != _stops.cend() && it->first >= req.range.first) {
        return true;
    }

    return false;
}

Coord_b ActiveManager::maxReplayRange() const {
    unsigned int i=0;

    auto stops(_stops);
    std::sort_heap(stops.begin(), stops.end(), _stop_comp);

    for(auto it=stops.cbegin(); it != stops.cend(); ++it) {
        const auto pos = _map.find(it->second)->second;
        if (req_old[pos / Packing][pos % Packing]) {
            return it->first;
        }
    }

    return std::numeric_limits<Coord_b>::min();
}