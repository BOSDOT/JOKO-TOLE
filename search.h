/*******************************************************************************                                                                    
 Joko Tole 1.0.0 beta Chess Engine Copyright (C) 2016 BOSDOT
 a lot inspired by chessprogramming.wikispaces.com and many 
 others engine stockfish, Gull, glaurung, fruit, etc
 as long as no code is not fail at build engine can be used 
--------------------------------------------------------------------------------					 					  					  
 ******************************************************************************/
 
#ifndef SEARCH_H
#define SEARCH_H

#include "board.h"
#include "general.h"

struct TwoFoldStack {
public:
    uint64_t keys[128];
    unsigned int length;

    TwoFoldStack() {
        length = 0;
    }
    ~TwoFoldStack() {}

    unsigned int size() { return length; }

    void push(uint64_t pos) {
        keys[length] = pos;
        length++;
    }

    void pop() { length--; }

    void clear() {
        length = 0;
    }

    bool find(uint64_t pos) {
        for (unsigned int i = length; i > 0; i--) {
            if (keys[i-1] == pos)
                return true;
        }
        return false;
    }
};

void getBestMove(Board *b, int mode, int value, Move *bestMove);
void clearTables();
void setHashSize(uint64_t MB);
void setEvalCacheSize(uint64_t MB);
uint64_t getNodes();
void setMultiPV(unsigned int n);
void setNumThreads(int n);

int getBestMoveForSort(Board *b, MoveList &legalMoves, int depth, int threadID);

// Search modes
const int TIME = 1;
const int DEPTH = 2;
// const int NODES = 3;
const int MOVETIME = 4;

// Time constants
const uint64_t ONE_SECOND = 1000;
const uint64_t MAX_TIME = (1ULL << 63) - 1;

// Search parameters
const int EASYMOVE_MARGIN = 150;
const int MAX_POS_score = 120;

// Time management constants
const int MOVE_HORIZON = 40; // expect this many moves left in the game
const int BUFFER_TIME = 100; // try to leave this much time in case of an emergency
const double TIME_FACTOR = 0.7; // timeFactor = log b / (b - 1) where b is branch factor
const double MAX_TIME_FACTOR = 2.5; // do not spend more than this multiple of time over the limit

#endif
