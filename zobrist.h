/*******************************************************************************                                                                    
						     Joko Tole Chess Engine
						   Copyright (C) 2016 BOSDOT						
             a lot inspired by chessprogramming.wikispaces.com and 
		    many others engine stockfish, Gull, glaurung, fruit, etc    
		   as long as no code is not fail at build engine can be used 
--------------------------------------------------------------------------------					 					  					  
 ******************************************************************************/
 
#ifndef ZOBRIST_H
#define ZOBRIST_H

#include "board.h"
#include "general.h"

const int EVAL_HASH_OFFSET = (1 << 20);

struct EvalHashEntry {
    uint64_t zobristKey;
    uint64_t score;

    EvalHashEntry() {
        clearEntry();
    }

    void setEntry(Board &b, int _score) {
        score = (uint64_t) (_score + EVAL_HASH_OFFSET);
        zobristKey = ((uint64_t) (b.getZobristKey() >> 32)) ^ score;
    }

    void clearEntry() {
        zobristKey = 0;
        score = 0;
    }

    ~EvalHashEntry() {}
};

class EvalHash {
private:
    EvalHashEntry *table;
    uint64_t size;
    EvalHash(const EvalHash &other);
    EvalHash& operator=(const EvalHash &other);

public:
    uint64_t keys;

    EvalHash(uint64_t MB);
    ~EvalHash();

    void add(Board &b, int score);
    int get(Board &b);
    void setSize(uint64_t MB);
    void clear();
};

#endif
