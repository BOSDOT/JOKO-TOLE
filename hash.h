/*******************************************************************************                                                                    
						     Joko Tole Chess Engine
						   Copyright (C) 2016 BOSDOT						
             a lot inspired by chessprogramming.wikispaces.com and 
		    many others engine stockfish, Gull, glaurung, fruit, etc    
		   as long as no code is not fail at build engine can be used 
--------------------------------------------------------------------------------					 					  					  
 ******************************************************************************/
 
#ifndef HASH_H
#define HASH_H

#include "board.h"
#include "general.h"

const uint8_t PV_NODE = 0;
const uint8_t CUT_NODE = 1;
const uint8_t ALL_NODE = 2;
const uint8_t NO_NODE_INFO = 3;


// Pack the information stored in a hash entry into a single 64-bit integer
uint64_t packHashData(int depth, Move m, int score, uint8_t nodeType, uint8_t age);

// Functions for unpacking hash data
inline int getHashDepth(uint64_t data) {
    return (int8_t) ((data >> 48) & 0xFF);
}

inline Move getHashMove(uint64_t data) {
    return (data >> 16) & 0xFFFF;
}

inline int getHashscore(uint64_t data) {
    return (int16_t) (data & 0xFFFF);
}

inline uint8_t getHashAge(uint64_t data) {
    return (data >> 40) & 0xFF;
}

inline uint8_t getHashNodeType(uint64_t data) {
    return (data >> 32) & 0x3;
}

struct HashEntry {
    uint64_t zobristKey;
    uint64_t data;

    HashEntry() {
        clearEntry();
    }

    void setEntry(Board &b, uint64_t _data) {
        zobristKey = b.getZobristKey() ^ _data;
        data = _data;
    }

    void clearEntry() {
        zobristKey = 0;
        data = 0;
    }

    ~HashEntry() {}
};

// This contains each of the hash table entries
class HashNode {
public:
    HashEntry slot1;
    HashEntry slot2;

    HashNode() {}
    ~HashNode() {}
};

class Hash {
private:
    HashNode *table;
    uint64_t size;

// Prevent direct copying and assignment
    Hash(const Hash &other);
    Hash& operator=(const Hash &other);

public:
    uint64_t keys;

    Hash(uint64_t MB);
    ~Hash();

    void add(Board &b, uint64_t data, int depth, uint8_t age);
    uint64_t get(Board &b);
    uint64_t getSize();
    void setSize(uint64_t MB);
    void clear();
};

#endif
