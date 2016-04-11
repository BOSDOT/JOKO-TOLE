/*******************************************************************************                                                                    
Joko Tole UCI Chess Engine Copyright(C) 2016 Bosdot (Javanese)
influenzed by chessprogramming.wikispaces.com and many others 
chess engine with open source code.
--------------------------------------------------------------------------------					 					  					  
 ******************************************************************************/

#ifndef MOVER_H
#define MOVER_H

#include <cmath>
#include "board.h"
#include "general.h"
#include "parameter.h"

enum MoveGenStage {
    STAGE_NONE, STAGE_HASH_MOVE, STAGE_IID_MOVE, STAGE_CAPTURES, STAGE_QUIETS
};

struct mover {
	Board *b;
	int color;
	int depth;
    int threadID;
	bool isPVNode;
	bool isInCheck;
	SearchParameters *parameter;
    MoveGenStage mgStage;
    Move hashed;
	MoveList legalMoves;
	scoreList scores;
    unsigned int quietStart;
	unsigned int index;

	mover(Board *_b, int _color, int _depth, int _threadID, bool _isPVNode, bool _isInCheck,
		SearchParameters *_parameter, Move _hashed, MoveList _legalMoves);

	// Node is reducible if not PV node and not in check
	bool nodeIsReducible();
    bool doIID();

	void generateMoves();
	Move nextMove();
    void reduceBadHistories(Move bestMove);

private:
    void scoreCaptures(bool isIIDMove);
    void scoreQuiets();
    void scoreIIDMove();
    void findQuietStart();
};

#endif
