/*******************************************************************************                                                                    
 Joko Tole 1.0.0 beta Chess Engine Copyright (C) 2016 BOSDOT
 a lot inspired by chessprogramming.wikispaces.com and many 
 others engine stockfish, Gull, glaurung, fruit, etc
 as long as no code is not fail at build engine can be used 
--------------------------------------------------------------------------------					 					  					  
 ******************************************************************************/
 
#ifndef EVALUATION_H
#define EVALUATION_H

#define V(mg, eg) ((score) ((eg << 16) | mg))
#include "general.h"

typedef uint64_t score;

// Evaluation constants and encoding.
const int SEE_PIECE_VALS[6] = {105, 380, 405, 720, 1140, MATE_score};
const int KNOWN_WIN = PAWN_VALUE_EG * 100;

// Retrieves the final evaluation score to return from the packed eval value
int decEval(score encodedValue, int egFactor) {
    int valueMg = (int) (encodedValue & 0xFFFF) - 0x8000;
    int valueEg = (int) (encodedValue >> 16) - 0x8000;
    return (valueMg * (EG_FACTOR_RES - egFactor) + valueEg * egFactor) / EG_FACTOR_RES;
}
// PENALTY
const score SEMIOPEN_OWN_PENALTY = V(9, 0);
const score SEMIOPEN_OPP_PENALTY = V(4, 1);
const score CENTRAL_ISOLATED_PENALTY = V(4, 4);
const score bishop_PAWN_COLOR_PENALTY = V(2, 2);
const score OPEN_PENALTY = V(4, 1);
const score ISOLATED_PENALTY = V(11, 13);
const score KNIGHT_C3_CLOSED_PENALTY = V(15, 0);
const score BLOCKADED_PASSER_PENALTY = V(8, 7);
const score BACKWARD_PENALTY = V(12, 13);

const score EVAL_ZERO = 0x80008000;
const int bishop_PAIR_VALUE = 58;
const int TEMPO_VALUE = 10;
const score CASTLING_RIGHTS_VALUE[3] = {V(0, 0), V(24, 0), V(40, 0)};
const score PAWN_SHIELD_VALUE = V(12, 0);
const score P_PAWN_SHIELD_BONUS = V(8, 0);
const score KNIGHT_PAWN_BONUS = V(1, 1);
const score KNIGHT_OUTPOST_BONUS = V(8, 1);
const score OUTPOST_PAWN_DEF_BONUS = V(3, 1);
const score ROOK_OPEN_FILE_BONUS = V(10, 10);

// Doubled pawns
const score ISOLATED_DOUBLED_PENALTY = V(8, 8);
const score DOUBLED_PENALTY_SCALE[9] = {0, 0, 3, 2, 1, 1, 1, 1, 1};

// Pawn structure, Passed pawns
const score PASSER_BONUS[8] = {
	V(0, 0), V(5, 15), V(5, 15), V(10, 25),
	V(25, 40), V(60, 65), V(100, 100), V(0, 0)};

const score PASSER_FILE_BONUS[8] = {
	V(11, 8), V(7, 4), V(3, 1), V(1, 0),
    V(0, 0),  V(2, 1), V(5, 4), V(10, 8)};
	
const score DOUBLED_PENALTY[7] = {
	V(0, 0), V(0, 0), V(11, 16), V(33, 48),
    V(80, 120), V(150, 230), V(250, 350)};
								  
#endif
