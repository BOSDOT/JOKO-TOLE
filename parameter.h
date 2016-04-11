/*******************************************************************************                                                                    
Joko Tole UCI Chess Engine Copyright(C) 2016 Bosdot (Javanese)
influenzed by chessprogramming.wikispaces.com and many others 
chess engine with open source code.
--------------------------------------------------------------------------------					 					  					  
 ******************************************************************************/
 
#ifndef PARAMETER_H
#define PARAMETER_H

#include "general.h"

struct SearchParameters {
    int ply;
    int nullMoveCount;
    int selectiveDepth;
    ChessTime startTime;
    uint64_t timeLimit;
    Move killers[MAX_DEPTH][2];
    int historyTable[2][6][64];
    uint8_t rootMoveNumber;

    SearchParameters() {
        reset();
        resetHistoryTable();
    }
	
    void reset() {
        ply = 0;
        nullMoveCount = 0;
        for (int i = 0; i < MAX_DEPTH; i++) {
            killers[i][0] = NULL_MOVE;
            killers[i][1] = NULL_MOVE;
        }
        //resetHistoryTable();
    }
    
    void resetHistoryTable() {
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 6; j++) {
                for (int k = 0; k < 64; k++)
                    historyTable[i][j][k] = 0;
            }
        }
    }

    void ageHistoryTable(int depth, bool isEndOfSearch) {
        int posHistoryScale, negHistoryScale;
        if (isEndOfSearch) {
            posHistoryScale = depth * depth;
            negHistoryScale = depth;
        }
        else {
            posHistoryScale = depth;
            negHistoryScale = std::max(1, ((int) std::sqrt(depth)) / 2);
        }

        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 6; j++) {
                for (int k = 0; k < 64; k++) {
                    if (historyTable[i][j][k] > 0)
                        historyTable[i][j][k] /= posHistoryScale;
                    else
                        historyTable[i][j][k] /= negHistoryScale;
                }
            }
        }
    }
};

#endif
