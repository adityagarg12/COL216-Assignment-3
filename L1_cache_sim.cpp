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

enum MESIState { INVALID, SHARED, EXCLUSIVE, MODIFIED };

struct CacheLine {
    bool valid = false;
    bool dirty = false;
    MESIState state = INVALID;
    unsigned int tag = 0;
    int lru_counter = 0;
};

struct CacheSet {
    vector<CacheLine> lines;
};

class Cache {
public:
    Cache(int s_bits, int e_assoc, int b_bits, int id);
    void access(char op, unsigned int address, int cycle);
    void snoop(unsigned int address, char op, int from_core);
    void print_stats();
    void update_lru(CacheSet& set, int accessed_idx);

    int core_id;
    int reads = 0, writes = 0;
    int misses = 0, evictions = 0, writebacks = 0;
    int idle_cycles = 0, total_cycles = 0;

private:
    int s, E, b;
    int num_sets;
    vector<CacheSet> sets;

    unsigned int get_tag(unsigned int addr);
    unsigned int get_index(unsigned int addr);
    int find_line(CacheSet& set, unsigned int tag);
};

Cache::Cache(int s_bits, int e_assoc, int b_bits, int id)
    : s(s_bits), E(e_assoc), b(b_bits), core_id(id) {
    num_sets = 1 << s;
    sets.resize(num_sets);
    for (int i = 0; i < num_sets; ++i) {
        sets[i].lines.resize(E);
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

    total_cycles++;
    if (op == 'R') reads++;
    else writes++;

    int line_idx = find_line(set, tag);

    if (line_idx != -1 && set.lines[line_idx].state != INVALID) {
        // Cache hit
        update_lru(set, line_idx);
        // MESI protocol logic here for hits
    } else {
        // Cache miss
        misses++;
        idle_cycles += 100;

        // Find replacement
        int victim = -1, max_lru = -1;
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

        if (set.lines[victim].valid && set.lines[victim].dirty) {
            writebacks++;
            idle_cycles += 100; // eviction writeback cost
        }

        set.lines[victim].valid = true;
        set.lines[victim].tag = tag;
        set.lines[victim].dirty = (op == 'W');
        set.lines[victim].state = (op == 'R') ? EXCLUSIVE : MODIFIED;
        update_lru(set, victim);

        // Simulate broadcast of miss to other cores (snooping)
        // Needs global coordination: currently placeholder
    }
}

void Cache::snoop(unsigned int address, char op, int from_core) {
    // To be implemented: snooping logic for MESI
    // Example: If this cache has a block in M, must writeback and downgrade
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
