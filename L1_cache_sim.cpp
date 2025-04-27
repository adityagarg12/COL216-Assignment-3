#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <iomanip>
#include <string>
#include <cassert>
#include <getopt.h>
#include <algorithm>
#include <cstring>
#include <map>

using namespace std;

enum MESIState { INVALID, SHARED, EXCLUSIVE, MODIFIED };

string mesi_to_string(MESIState state) {
    switch (state) {
        case INVALID: return "I";
        case SHARED: return "S";
        case EXCLUSIVE: return "E";
        case MODIFIED: return "M";
        default: return "?";
    }
}

// Forward declaration of Cache class
class Cache;

class BusManager {
public:
    vector<Cache*> cores;
    // int invalidation_count = 0;
    unsigned long long data_traffic_bytes = 0;
    int bus_transactions = 0;
    void register_core(Cache* c) {
        cores.push_back(c);
    }
    
    bool broadcast_BusRd(int from_core_id, unsigned int addr, int block_size_bytes, int& cycles_of_operation);
    bool broadcast_BusRdX(int from_core_id, unsigned int addr, int block_size_bytes, int& cycles_of_operation);
    void broadcast_Invalidate(int from_core_id, unsigned int addr, int& cycles_of_operation);
    
    // int get_invalidations() const { return invalidation_count; }
    unsigned long long get_data_traffic() const { return data_traffic_bytes; }
};

struct CacheLine {
    bool valid = false;
    bool dirty = false;
    MESIState state = INVALID;
    unsigned int tag = 0;
    
    CacheLine(int block_size) {}
};

struct CacheSet {
    vector<CacheLine> lines;
    
    CacheSet(int associativity, int block_size) {
        lines.reserve(associativity);
        for (int i = 0; i < associativity; i++) {
            lines.emplace_back(block_size);
        }
    }
};

class Cache {
public:
    Cache(int s_bits, int e_assoc, int b_bits, int id, BusManager* bus_ptr);
    int access(char op, unsigned int address);
    bool snoop(unsigned int address, char op, int from_core, int block_size_bytes, int& cycles_of_operation);
    void print_stats(ostream& out);
    void update_lru(CacheSet& set, int accessed_idx);
    
    double get_cache_utilization(int& valid_lines, int& total_lines) const {
        valid_lines = 0;
        total_lines = (1 << s) * E;
        
        for (const auto& set : sets) {
            for (const auto& line : set.lines) {
                if (line.valid && line.state != INVALID) valid_lines++;
            }
        }
        
        return (valid_lines * 100.0 / total_lines);
    }
    
    int M_to_S = 0, M_to_I = 0, E_to_S = 0, E_to_M = 0, E_to_I = 0, S_to_M = 0, S_to_I = 0;

    int core_id;
    int accesses = 0;
    int reads = 0, writes = 0;
    int misses = 0, read_misses = 0, write_misses = 0;
    int evictions = 0, writebacks = 0;
    int idle_cycles = 0, total_cycles = 0;
    int invalidation_count = 0; 
    int data_traffic_bytes = 0;
private:
    int s, E, b;
    int num_sets;
    int block_size;
    BusManager* bus;
    vector<CacheSet> sets;

    unsigned int get_tag(unsigned int addr);
    unsigned int get_index(unsigned int addr);
    unsigned int get_block_offset(unsigned int addr);
    int find_line(CacheSet& set, unsigned int tag);
    int get_victim_line(CacheSet& set);
};

Cache::Cache(int s_bits, int e_assoc, int b_bits, int id, BusManager* bus_ptr)
    : core_id(id), s(s_bits), E(e_assoc), b(b_bits), bus(bus_ptr) {
    num_sets = 1 << s;
    block_size = 1 << b;
    
    sets.reserve(num_sets);
    for (int i = 0; i < num_sets; ++i) {
        sets.emplace_back(E, block_size);
    }
}

unsigned int Cache::get_tag(unsigned int addr) {
    return addr >> (s + b);
}

unsigned int Cache::get_index(unsigned int addr) {
    return (addr >> b) & ((1 << s) - 1);
}

unsigned int Cache::get_block_offset(unsigned int addr) {
    return addr & ((1 << b) - 1);
}

int Cache::find_line(CacheSet& set, unsigned int tag) {
    for (int i = 0; i < E; ++i) {
        if (set.lines[i].valid && set.lines[i].state != INVALID && set.lines[i].tag == tag) {
            return i;
        }
    }
    return -1;
}

int Cache::get_victim_line(CacheSet& set) {
    // First check for any invalid lines
    for (int i = 0; i < E; ++i) {
        if (!set.lines[i].valid || set.lines[i].state == INVALID) {
            return i;
        }
    }
    
    // Otherwise, use LRU policy - last element is LRU
    return E - 1;
}

void Cache::update_lru(CacheSet& set, int accessed_idx) {
    // Move the accessed line to the front (MRU position)
    CacheLine temp = set.lines[accessed_idx];
    for (int i = accessed_idx; i > 0; --i) {
        set.lines[i] = set.lines[i-1];
    }
    set.lines[0] = temp;
}

int Cache::access(char op, unsigned int address) {
    unsigned int index = get_index(address);
    unsigned int tag = get_tag(address);
    CacheSet& set = sets[index];
    // total_cycles++; // Increment total cycles for each access
    int cycles_for_operation = 1; 
    
    accesses++;
    if (op == 'R') reads++;
    else writes++;

    int line_idx = find_line(set, tag);

    if (line_idx != -1) {
        // Cache hit
        CacheLine& line = set.lines[line_idx];
        
        if (op == 'W') {
            if (line.state == SHARED) {
                S_to_M++;
                bus->broadcast_Invalidate(core_id, address,cycles_for_operation);
                line.state = MODIFIED;
                line.dirty = true;
            } else if (line.state == EXCLUSIVE) {
                E_to_M++;
                line.state = MODIFIED;
                line.dirty = true;
            }
            // If already MODIFIED, stay in MODIFIED
        }
        
        update_lru(set, line_idx);
        cycles_for_operation = 1; // Cache hit takes 1 cycle
    } else {
        // Cache miss
        misses++;
        if (op == 'R') read_misses++;
        else write_misses++;
        
        int victim = get_victim_line(set);
        CacheLine& line = set.lines[victim];
        
        if (line.valid && line.state != INVALID) {
            evictions++;
            if (line.dirty) {
                writebacks++;
                cycles_for_operation += 100; // Memory write takes 100 cycles
                idle_cycles += 100;
            }
        }
        
        line.valid = true;
        line.tag = tag;
        
        if (op == 'R') {
            // Check if other cores have this cache line
            bool shared = bus->broadcast_BusRd(core_id, address, block_size,cycles_for_operation);
            data_traffic_bytes += block_size; // Data traffic for BusRd
            if (shared) {
                line.state = SHARED;
                line.dirty = false;
                // Fetch from another cache (words_per_block * 2 cycles)
                int words_per_block = block_size / 4; // 4 bytes per word
                cycles_for_operation += words_per_block * 2;
                idle_cycles += words_per_block * 2;
            } else {
                line.state = EXCLUSIVE;
                line.dirty = false;
                cycles_for_operation += 100; // Memory read takes 100 cycles
                idle_cycles += 100;
            }
        } else { // Write
            bus->broadcast_BusRdX(core_id, address, block_size,cycles_for_operation);
            data_traffic_bytes += block_size; // Data traffic for BusRdX
            line.state = MODIFIED;
            line.dirty = true;
            cycles_for_operation += 100; // Memory read takes 100 cycles
            idle_cycles += 100;
        }
        
        update_lru(set, victim);
    }
    
    total_cycles += cycles_for_operation;
    return cycles_for_operation;
}

bool Cache::snoop(unsigned int address, char op, int from_core, int block_size_bytes, int& cycles_of_operation) {
    unsigned int index = get_index(address);
    unsigned int tag = get_tag(address);
    CacheSet& set = sets[index];
    
    int line_idx = find_line(set, tag);
    
    if (line_idx == -1 || set.lines[line_idx].state == INVALID)
        return false;
    
    CacheLine& line = set.lines[line_idx];
    bool shared = false;
    
    if (op == 'R') { // BusRd
        if (line.state == EXCLUSIVE) {
            E_to_S++;
            line.state = SHARED;
            shared = true;
        } else if (line.state == SHARED) {
            // Stay in SHARED
            shared = true;
        } else if (line.state == MODIFIED) {
            M_to_S++;
            line.state = SHARED;
            line.dirty = false; // Write back to memory
            writebacks++; // Need to write back data to memory
            cycles_of_operation += 100; // Memory write takes 100 cycles
            idle_cycles+=100;
            data_traffic_bytes+=block_size_bytes;
            shared = true;
        }
    } else if (op == 'W') { // BusRdX
        if (line.state == SHARED) {
            S_to_I++;
            invalidation_count++;
            shared = true;
            line.state = INVALID;
            line.valid = false;
        } else if (line.state == EXCLUSIVE) {
            E_to_I++;
            shared = true;
            invalidation_count++;
            line.state = INVALID;
            line.valid = false;
        } else if (line.state == MODIFIED) {
            M_to_I++;
            invalidation_count++;
            if (line.dirty) {
                writebacks++; // Need to write back data to memory
                cycles_of_operation += 100; // Memory write takes 100 cycles
                idle_cycles+=100;
                data_traffic_bytes+=block_size;
            }
            shared = true;
            line.state = INVALID;
            line.valid = false;
            line.dirty = false;
        }
    } else if (op == 'I') { // Invalidate operation
        if (line.state == SHARED) {
            S_to_I++;
            shared = true;
            line.state = INVALID;
            invalidation_count++;
            line.valid = false;
        } else if (line.state == EXCLUSIVE) {
            E_to_I++;
            shared = true;
            line.state = INVALID;
            line.valid = false;
            invalidation_count++;
        } else if (line.state == MODIFIED) {
            M_to_I++;
            invalidation_count++;
            shared = true;
            if (line.dirty) {
                writebacks++; // Need to write back data to memory
                cycles_of_operation += 100; // Memory write takes 100 cycles
                idle_cycles+=100;
                data_traffic_bytes+=block_size;
            }
            line.state = INVALID;
            line.valid = false;
            line.dirty = false;
        }
    }
    
    return shared;
}

void Cache::print_stats(ostream& out) {
    out << "Core " << core_id << " Statistics:\n";
    out << "  Reads: " << reads << ", Writes: " << writes << " (Total: " << accesses << ")\n";
    out << "  Read Misses: " << read_misses << ", Write Misses: " << write_misses << " (Total: " << misses << ")\n";
    out << "  Total execution cycles: " << total_cycles << "\n";
    out << "  Idle cycles: " << idle_cycles << " (" << fixed << setprecision(2) << (total_cycles ? idle_cycles * 100.0 / total_cycles : 0.0) << "%)\n";
    out << "  Miss rate: " << fixed << setprecision(2) << (accesses ? misses * 100.0 / accesses : 0.0) << "%\n";
    out << "  Hit rate: " << fixed << setprecision(2) << (accesses ? (accesses - misses) * 100.0 / accesses : 0.0) << "%\n";
    out << "  Evictions: " << evictions << "\n";
    out << "  Writebacks: " << writebacks << "\n";
    out << "  Bus Invalidations: " << invalidation_count << "\n";  
    out << "  Data traffic (Bytes): " << data_traffic_bytes << " bytes\n";  
    // out << "  MESI transitions: M→S: " << M_to_S << ", M→I: " << M_to_I << ", E→S: " << E_to_S << 
    //        ", E→M: " << E_to_M << ", E→I: " << E_to_I << ", S→M: " << S_to_M << ", S→I: " << S_to_I << "\n";
    out << "------------------------\n";
}

bool BusManager::broadcast_BusRd(int from_core_id, unsigned int addr, int block_size_bytes, int& cycles_of_operation) {
    bool shared = false;
    bus_transactions++;
    for (Cache* core : cores) {
        if (core->core_id != from_core_id) {
            if (core->snoop(addr, 'R', from_core_id, block_size_bytes, cycles_of_operation)) {
                shared = true;
                data_traffic_bytes += block_size_bytes;
            }
        }
    }
    
    return shared;
}

bool BusManager::broadcast_BusRdX(int from_core_id, unsigned int addr, int block_size_bytes, int& cycles_of_operation) {
    bool shared = false;
    bus_transactions++;
    for (Cache* core : cores) {
        if (core->core_id != from_core_id) {
            if (core->snoop(addr, 'W', from_core_id, block_size_bytes, cycles_of_operation)) {
                shared = true;
                // data_traffic_bytes += block_size_bytes;
            }
        }
        data_traffic_bytes += block_size_bytes; // Always add data traffic for BusRdX
    }
    
    return shared;
}

void BusManager::broadcast_Invalidate(int from_core_id, unsigned int addr, int& cycles_of_operation) {
    bus_transactions++;
    for (Cache* core : cores) {
        if (core->core_id != from_core_id) {
            if (core->snoop(addr, 'I', from_core_id, 0, cycles_of_operation)) {
                // If the core had a valid copy
                
            }
        }
    }
}

void print_help() {
    cout << "Usage: ./L1simulate [OPTIONS]\n";
    cout << "  -t <tracefile>:  name of parallel application (e.g., app1) whose 4 traces are to be used\n";
    cout << "  -s <s>:          number of set index bits (number of sets = 2^s)\n";
    cout << "  -E <E>:          associativity (number of cache lines per set)\n";
    cout << "  -b <b>:          number of block bits (block size = 2^b bytes)\n";
    cout << "  -o <outfile>:    logs output in file for plotting etc.\n";
    cout << "  -h:              prints this help\n";
    cout << "\nExample: ./L1simulate -t app1 -s 5 -E 2 -b 5 -o output.txt\n";
}

struct Operation {
    char op;
    unsigned int address;
};

struct Event {
    int cycle;
    int core_id;
    size_t trace_idx;
    bool is_completion;
    
    bool operator<(const Event& other) const {
        // Use > to create a min-heap (earliest cycle first)
        if (cycle != other.cycle)
            return cycle > other.cycle;
        
        // If cycles are equal, prioritize completions over starts
        if (is_completion != other.is_completion)
            return is_completion < other.is_completion;
            
        // If both are completions or both are starts, use core_id to break ties
        return core_id > other.core_id;
    }
};

int main(int argc, char* argv[]) {
    char *tracefile_prefix = nullptr;
    int s_bits = 0;
    int E_assoc = 0;
    int b_bits = 0;
    char *output_file = nullptr;
    bool output_to_file = false;
    
    int opt;
    while ((opt = getopt(argc, argv, "t:s:E:b:o:h")) != -1) {
        switch (opt) {
            case 't':
                tracefile_prefix = optarg;
                break;
            case 's':
                s_bits = atoi(optarg);
                break;
            case 'E':
                E_assoc = atoi(optarg);
                break;
            case 'b':
                b_bits = atoi(optarg);
                break;
            case 'o':
                output_file = optarg;
                output_to_file = true;
                break;
            case 'h':
                print_help();
                return 0;
            default:
                cerr << "Invalid option\n";
                print_help();
                return 1;
        }
    }
    
    if (!tracefile_prefix || s_bits <= 0 || E_assoc <= 0 || b_bits <= 0) {
        cerr << "Error: Missing or invalid required arguments\n";
        print_help();
        return 1;
    }
    
    ofstream outfile;
    if (output_to_file) {
        outfile.open(output_file);
        if (!outfile) {
            cerr << "Error opening output file: " << output_file << endl;
            return 1;
        }
    }
    
    ostream& out = output_to_file ? outfile : cout;
    
    out << "Simulation Parameters:\n";
    out << "Trace Prefix: " << tracefile_prefix << "\n";
    out << "Set Index Bits: " << s_bits << "\n";
    out << "Associativity: " << E_assoc << "\n";
    out << "Block Bits: " << b_bits << "\n";
    out << "Block Size (Bytes): " << (1 << b_bits) << "\n";
    out << "Number of Sets: " << (1 << s_bits) << "\n";
    out << "Cache Size (KB per core): " << ((1 << s_bits) * E_assoc * (1 << b_bits)) / 1024 << " KB\n";
    out << "MESI Protocol: Enabled\n";
    out << "Write Policy: Write-back, Write-allocate\n";
    out << "Replacement Policy: LRU\n";
    out << "Bus: Central snooping bus\n";
    out << "------------------------\n";
    
    BusManager bus_manager;
    vector<Cache*> cores;
    
    for (int i = 0; i < 4; ++i) {
        cores.push_back(new Cache(s_bits, E_assoc, b_bits, i, &bus_manager));
        bus_manager.register_core(cores[i]);
    }
    
    // Preprocess trace files
    vector<vector<Operation>> traces(4);
    for (int i = 0; i < 4; ++i) {
        stringstream ss;
        ss << tracefile_prefix << "_proc" << i << ".trace";
        string filename = ss.str();
        
        ifstream trace_file(filename);
        if (!trace_file) {
            cerr << "Error: Could not open trace file: " << filename << endl;
            return 1;
        }
        
        string line;
        while (getline(trace_file, line)) {
            if (line.empty() || line[0] != 'R' && line[0] != 'W') continue;
            
            char op = line[0];
            unsigned int addr;
            
            // Extract the address part - handle both formats
            size_t pos = line.find("0x");
            if (pos != string::npos) {
                string addr_str = line.substr(pos);
                addr = stoul(addr_str, nullptr, 0);
            } else {
                istringstream iss(line);
                iss >> op >> hex >> addr;
            }
            
            traces[i].push_back({op, addr});
        }
        
        trace_file.close();
    }
    
    // Event-driven simulation
    priority_queue<Event> event_queue;
    vector<size_t> trace_indices(4, 0);
    vector<bool> blocked(4, false);
    int global_cycles = 0;
    
    // Initialize event queue with the first operation from each core
    for (int i = 0; i < 4; ++i) {
        if (!traces[i].empty()) {
            event_queue.push({0, i, trace_indices[i], false});
        }
    }
    
    while (!event_queue.empty()) {
        Event e = event_queue.top();
        event_queue.pop();
        
        global_cycles = max(global_cycles, e.cycle);
        
        if (e.is_completion) {
            // Operation completed, core is no longer blocked
            blocked[e.core_id] = false;
            
            // Schedule the next operation from this core
            if (trace_indices[e.core_id] < traces[e.core_id].size()) {
                event_queue.push({global_cycles, e.core_id, trace_indices[e.core_id], false});
            }
        } else {
            // Starting a new operation
            if (!blocked[e.core_id]) {
                Operation& op = traces[e.core_id][e.trace_idx];
                int cycle_count = cores[e.core_id]->access(op.op, op.address);
                
                trace_indices[e.core_id]++;
                
                if (cycle_count > 1) {
                    // Operation takes multiple cycles, block the core
                    blocked[e.core_id] = true;
                    event_queue.push({global_cycles + cycle_count, e.core_id, e.trace_idx, true});
                } else if (trace_indices[e.core_id] < traces[e.core_id].size()) {
                    // Operation took 1 cycle (cache hit), schedule next operation
                    event_queue.push({global_cycles + 1, e.core_id, trace_indices[e.core_id], false});
                }
            }
        }
    }
    
    out << "\nSimulation Results:\n";
    for (int i = 0; i < 4; ++i) {
        cores[i]->print_stats(out);
    }
    
    out << "Overall Bus Summary:\n";
    out << "  Total Bus Transactions: " << bus_manager.bus_transactions << "\n";
    int total_data_traffic = 0;
    for (int i = 0; i < 4; ++i) {
        total_data_traffic += cores[i]->data_traffic_bytes;
    }
    out << "  Total Bus Traffic (Bytes): " << total_data_traffic << "\n";
    out << "------------------------\n";
    
    int max_cycles = 0;
    for (int i = 0; i < 4; ++i) {
        max_cycles = max(max_cycles, cores[i]->total_cycles);
    }
    out << "Maximum execution time: " << max_cycles << " cycles\n";
    out << "Total global cycles: " << global_cycles << " cycles\n";
    
    out << "\nCache Utilization:\n";
    for (int i = 0; i < 4; ++i) {
        int valid_lines = 0;
        int total_lines = 0;
        double utilization = cores[i]->get_cache_utilization(valid_lines, total_lines);
        
        out << "  Core " << i << ": " << fixed << setprecision(2) 
            << utilization << "% (" 
            << valid_lines << "/" << total_lines << " lines)\n";
    }
    
    for (Cache* core : cores) {
        delete core;
    }
    
    if (output_to_file) {
        outfile.close();
        cout << "Results written to " << output_file << endl;
    }
    
    return 0;
}