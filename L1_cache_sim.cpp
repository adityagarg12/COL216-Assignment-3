#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <deque>
#include <iomanip>
#include <string>
#include <bitset>
#include <algorithm>
#include <cassert>

using namespace std;

// Forward declaration of Cache class
class Cache;

class BusManager {
    public:
        vector<Cache*> cores;
    
        void register_core(Cache* c) {
            cores.push_back(c);
        }
    
        // For Read miss
        bool broadcast_BusRd(int from_core_id, unsigned int addr) {
            bool shared = false;
            for (std::vector<Cache*>::iterator it = cores.begin(); it != cores.end(); ++it) {
                Cache* core = *it;
                if (core->core_id != from_core_id) {
                    shared = shared || core->snoop(addr, 'R', from_core_id);
                }
            }
            
            return shared;
        }
    
        // For Write miss
        bool broadcast_BusRdX(int from_core_id, unsigned int addr) {
            bool shared = false;
            for (std::vector<Cache*>::iterator it = cores.begin(); it != cores.end(); ++it) {
                Cache* core = *it;
                if (core->core_id != from_core_id) {
                    shared = shared || core->snoop(addr, 'W', from_core_id);
                }
            }

            return shared;
        }

        void broadcast_Invalidate(int from_core_id, unsigned int addr) {
            for (std::vector<Cache*>::iterator it = cores.begin(); it != cores.end(); ++it) {
                Cache* core = *it;
                if (core->core_id != from_core_id) {
                    core->snoop(addr, 'I', from_core_id);
                }
            }

        }
    };
    

enum MESIState { INVALID, SHARED, EXCLUSIVE, MODIFIED };

struct CacheLine {
    bool valid = false;
    bool dirty = false;
    MESIState state = INVALID;
    unsigned int tag = 0;
    int lru_counter = 0;
    vector<char> block_data; // Placeholder for block data
};

struct CacheSet {
    vector<CacheLine> lines;
};

class Cache {
public:
    Cache(int s_bits, int e_assoc, int b_bits, int id, BusManager* bus_ptr);
    void access(char op, unsigned int address, int cycle);
    bool snoop(unsigned int address, char op, int from_core);
    void print_stats();
    void update_lru(CacheSet& set, int accessed_idx);

    int core_id;
    int accesses = 0;
    int reads = 0, writes = 0;
    int misses = 0, evictions = 0, writebacks = 0;
    int idle_cycles = 0, total_cycles = 0;
    int invalidations =0;

private:
    int s, E, b;
    int num_sets;
    BusManager* bus;
    vector<CacheSet> sets;

    unsigned int get_tag(unsigned int addr);
    unsigned int get_index(unsigned int addr);
    int find_line(CacheSet& set, unsigned int tag);
};

Cache::Cache(int s_bits, int e_assoc, int b_bits, int id, BusManager* bus_ptr)
    : s(s_bits), E(e_assoc), b(b_bits), core_id(id), bus(bus_ptr) {
    num_sets = 1 << s; // 2^s sets
    sets.resize(num_sets);
    for (int i = 0; i < num_sets; ++i) {
        sets[i].lines.resize(E);
        for (int j = 0; j < E; ++j) {
            sets[i].lines[j].block_data.resize(1 << b); // Block size 2^b bytes
        }
    }
}

unsigned int Cache::get_tag(unsigned int addr) {
    return addr >> (s + b);
}

unsigned int Cache::get_index(unsigned int addr) {
    return (addr >> b) & ((1 << s) - 1);
}

int Cache::find_line(CacheSet& set, unsigned int tag) {
    for (int i = 0; i < E; ++i) {
        if (set.lines[i].valid && set.lines[i].tag == tag) {
            return i;
        }
    }
    return -1;
}

void Cache::update_lru(CacheSet& set, int accessed_idx) {
    for (int i = 0; i < E; ++i) {
        if (i == accessed_idx) set.lines[i].lru_counter = 0;
        else set.lines[i].lru_counter++;
    }
}

void Cache::access(char op, unsigned int address, int cycle) {
    unsigned int index = get_index(address);
    unsigned int tag = get_tag(address);
    CacheSet& set = sets[index];
    accesses++;
    total_cycles++;
    if (op == 'R') reads++;
    else writes++;

    int line_idx = find_line(set, tag);

    if (line_idx != -1 && set.lines[line_idx].state != INVALID) {
        // Cache hit
        update_lru(set, line_idx);
        // MESI protocol logic here for hits
        if (op == 'W') {
            if (set.lines[line_idx].state == SHARED) {
                bus->broadcast_Invalidate(core_id, address); // Broadcast BusRdX for write
                set.lines[line_idx].state = MODIFIED;
            }
            else if (set.lines[line_idx].state == EXCLUSIVE) {
                set.lines[line_idx].state = MODIFIED;
            }
            set.lines[line_idx].dirty = true;
        }
        else if (op == 'R') {
            //just return value
        }

    } else {
        // Cache miss
        misses++;

        // Find replacement
        int victim = -1;
        if (line_idx == -1) {
            
            int max_lru = -1;
            for (int i = 0; i < E; ++i) {
                if (!set.lines[i].valid) {
                    victim = i;
                    break;
                }
                if (set.lines[i].lru_counter > max_lru) {
                    max_lru = set.lines[i].lru_counter;
                    victim = i;
                }
            }
            
            if (set.lines[victim].valid) {
                evictions++;
                // Evict the line
                if (set.lines[victim].dirty) {
                    writebacks++;
                    idle_cycles += 100; // eviction writeback cost
                }
            }
            
        set.lines[victim].valid = true;
        set.lines[victim].tag = tag;
        set.lines[victim].dirty = (op == 'W');
       // set.lines[victim].state = (op == 'R') ? EXCLUSIVE : MODIFIED;
        

        }
        else if (set.lines[line_idx].state == INVALID){
            victim = line_idx;
        }

        update_lru(set, victim);

        if (op == 'R'){
            // Broadcast BusRd to other cores
            if (bus->broadcast_BusRd(core_id, address)) {
                // If another core has the block in SHARED state, we can get it from there
                
                idle_cycles += ((1 << b) / 4) * 2; // Update idle_cycles by number of words transferred * 2 cycles
                
                set.lines[victim].state = SHARED;
                set.lines[victim].valid = true;
                set.lines[victim].tag = tag;
                set.lines[victim].dirty = false;

            } else {
                // If no other core has it, we can load it into EXCLUSIVE state
                
                idle_cycles += 100; // simulate memory access time

                set.lines[victim].valid = true;
                set.lines[victim].tag = tag;
                set.lines[victim].dirty = false;
                set.lines[victim].state = EXCLUSIVE;
            }

        } else if (op == 'W') {
            // Broadcast BusRdX to invalidate other caches
            if(bus->broadcast_BusRdX(core_id, address)){               
                // never 
            }
            else {

                idle_cycles+=100; // Value read from memory to local cache
                
                set.lines[victim].valid = true;
                set.lines[victim].tag = tag;
                set.lines[victim].dirty = true;
                set.lines[victim].state = MODIFIED;

            }
            
        }

        

        


        // Simulate broadcast of miss to other cores (snooping)
        // Needs global coordination: currently placeholder
    }
}

bool Cache::snoop(unsigned int address, char op, int from_core) {
    unsigned int index = get_index(address);
    unsigned int tag = get_tag(address);
    CacheSet& set = sets[index];
    int line_idx = find_line(set, tag);

    if (line_idx == -1 || set.lines[line_idx].state == INVALID) return false;

    CacheLine& line = set.lines[line_idx];

    if (op == 'R') {
        if(line.state == EXCLUSIVE) {
            line.state = SHARED;

            return true;
        }
        else if (line.state == SHARED){
            return true;
        }
        else if (line.state == MODIFIED){
            line.state = SHARED;
            writebacks++;
            idle_cycles += 100; // simulate writeback
            line.dirty = false;
        }

    } else if (op == 'W') {
        if (line.state == SHARED || line.state == EXCLUSIVE){
            line.state = INVALID;
            invalidations++;
            return false;
        }
        else if (line.state == MODIFIED){
            line.state = INVALID;
            invalidations++;
            writebacks++;
            idle_cycles +=100 ; //writeback
            return false;
        }
    }
    else if (op == 'I'){
        line.state = INVALID;
        invalidations++;
    }

    return true;
}

void Cache::print_stats() {
    cout << "Core " << core_id << ":\n";
    cout << "Reads: " << reads << ", Writes: " << writes << ", Misses: " << misses << endl;
    cout << "Evictions: " << evictions << ", Writebacks: " << writebacks << endl;
    cout << "Idle cycles: " << idle_cycles << ", Total cycles: " << total_cycles << endl;
    cout << fixed << setprecision(2);
    cout << "Miss rate: " << (misses * 100.0 / (reads + writes)) << "%\n";
    cout << "------------------------\n";
}

// Add global bus coordinator and per-core trace handling

int main(int argc, char* argv[]) {
    // Placeholder: parse args, load traces, simulate per-core accesses
    cout << "L1 Cache Simulator starting...\n";

    return 0;
}
