/*******************************************************************************                                                                    
 Joko Tole 1.0.0 beta Chess Engine Copyright (C) 2016 BOSDOT
 a lot inspired by chessprogramming.wikispaces.com and many 
 others engine stockfish, Gull, glaurung, fruit, etc
 as long as no code is not fail at build engine can be used 
--------------------------------------------------------------------------------					 					  					  
 ******************************************************************************/
 
#include "zobrist.h"

// Adds key and move
void EvalHash::add(Board &b, int score) {

    uint64_t h = (uint64_t) (b.getZobristKey() & 0xFFFFFFFF);
    uint64_t index = h % size;
    table[index].setEntry(b, score);
}

// Get the hash
int EvalHash::get(Board &b) {

    uint64_t h = (uint64_t) (b.getZobristKey() & 0xFFFFFFFF);
    uint64_t index = h % size;

    if((table[index].zobristKey ^ table[index].score) == (uint64_t) (b.getZobristKey() >> 32))
        return table[index].score;

    return 0;
}

void EvalHash::setSize(uint64_t MB) {
    delete[] table;
    uint64_t arrSize = MB << 20;
    arrSize /= sizeof(EvalHashEntry);
    size = arrSize;
    table = new EvalHashEntry[size];
    keys = 0;
}

EvalHash::EvalHash(uint64_t MB) {
    uint64_t arrSize = MB << 20;
    arrSize /= sizeof(EvalHashEntry);
    table = new EvalHashEntry[arrSize];
    size = arrSize;
    keys = 0;
}

EvalHash::~EvalHash() {
    delete[] table;
}

void EvalHash::clear() {
    delete[] table;
    table = new EvalHashEntry[size];
    keys = 0;
}
