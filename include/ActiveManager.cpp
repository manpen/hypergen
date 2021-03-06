/**
 * @file
 * @brief Implementation of ActiveManager
 *
 * @author Manuel Penschuck
 * @copyright
 * Copyright (C) 2017 Manuel Penschuck
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * @copyright
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * @copyright
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "ActiveManager.hpp"
#include "Assert.hpp"

ActiveManager::ActiveManager(Node avg_degree) {
    _size = avg_degree;
    //_update_size();
    _size = 0;
}

ActiveManager::~ActiveManager() {

}

template<typename T>
static void deallocate_vector(T& v) {
    T empty;
    v.swap(empty);
}


void ActiveManager::clear(bool keep_queues) {
    deallocate_vector(_req_phi);
    deallocate_vector(_req_poin_x);
    deallocate_vector(_req_poin_y);
    deallocate_vector(_req_poin_r);
    deallocate_vector(_req_ids);

    if (!keep_queues) {
        deallocate_vector(_stops);
        deallocate_vector(_starts);
    }

    deallocate_vector(_map);

    _size = 0;
    _end = 0;
}

void ActiveManager::copyFrom(const ActiveManager& src, bool markOld, Coord_b thresh){
    for (const auto& msg : src._starts) {
        if (msg.range.first < thresh) {
            _starts.push_back(msg);
            if (markOld)
                _starts.back().setOld();
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

        std::copy(src._req_phi.cbegin(),    src._req_phi.cbegin() + _end,    _req_phi.begin());
        std::copy(src._req_poin_x.cbegin(), src._req_poin_x.cbegin() + _end, _req_poin_x.begin());
        std::copy(src._req_poin_y.cbegin(), src._req_poin_y.cbegin() + _end, _req_poin_y.begin());
        std::copy(src._req_poin_r.cbegin(), src._req_poin_r.cbegin() + _end, _req_poin_r.begin());
        std::copy(src._req_ids.cbegin(),    src._req_ids.cbegin() + _size,   _req_ids.begin());

        if (markOld) {
            for(unsigned int i=0; i < _size; ++i)
                _req_ids[i] |= Point::OLD_MASK;
        }
    }
}



void ActiveManager::copyFromBelow(const ActiveManager &src, const Geometry& geometry,
                                  const SinhCosh &newRad, const CoordInter phi) {
    auto addRequest = [&](const CoordInter range,
                          const Request &req) {
        if (range.first > phi.second)
            return;

        if (range.second < phi.first)
            return;

        _starts.emplace_back(req, CoordInter{
            std::max<Coord>(range.first, phi.first),
            std::min<Coord>(range.second, phi.second)
        });
    };


    for(const  Request& old : src._starts) {
        Request req(old, geometry, newRad);

        if (req.range.first >= 0.0 && req.range.second <= 2.0 * M_PI) {
            addRequest(req.range, req);
        } else if (req.range.second - req.range.first >= 2 * M_PI - 1e-3) {
            addRequest({0, 2 * M_PI}, req);
        } else if (req.range.first <= 0) {
            addRequest({0, req.range.second}, req);
            addRequest({2 * M_PI + req.range.first, 2 * M_PI}, req);
        } else {
            ASSERT_GT(req.range.second, 2 * M_PI);
            addRequest({0, req.range.second - 2 * M_PI}, req);
            addRequest({req.range.first, 2 * M_PI}, req);
        }
    }

    std::make_heap(_starts.begin(), _starts.end(), _start_comp);
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
                req.range.second = maxRange;
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

    // no old in aux structures
    for(const auto& it : _map)
        ASSERT(!(it.first & Point::OLD_MASK));

    for(const auto& msg : _stops)
        ASSERT(!(msg.second & Point::OLD_MASK));
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

        if (_req_ids[pos] & Point::OLD_MASK) {
            return it->first;
        }
    }

    return std::numeric_limits<Coord_b>::min();
}