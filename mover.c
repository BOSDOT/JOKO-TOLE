/*******************************************************************************                                                                    
Joko Tole UCI Chess Engine Copyright(C) 2016 Bosdot (Javanese)
influenzed by chessprogramming.wikispaces.com and many others 
chess engine with open source code.
--------------------------------------------------------------------------------					 					  					  
 ******************************************************************************/

#include "search.h"
#include "mover.h"

const int score_IID_MOVE = (1 << 20);
const int score_WINNING_CAPTURE = (1 << 18);
const int score_QUEEN_PROMO = (1 << 17);
const int score_EVEN_CAPTURE = (1 << 16);
const int score_LOSING_CAPTURE = -(1 << 30) - (1 << 28);
const int score_QUIET_MOVE = -(1 << 30);
const int score_LOSING_QUIET = -(1 << 30) - (1 << 29) - (1 << 28);
const int QUIET_SEE_DEPTH = 2;

mover::mover(Board *_b, int _color, int _depth, int _threadID, bool _isPVNode, bool _isInCheck,
	SearchParameters *_parameter, Move _hashed, MoveList _legalMoves) {
	b = _b;
	color = _color;
	depth = _depth;
    threadID = _threadID;
	isPVNode = _isPVNode;
	isInCheck = _isInCheck;
	parameter = _parameter;
    mgStage = STAGE_NONE;
    quietStart = 0;
    index = 0;
    hashed = _hashed;
    legalMoves = _legalMoves;
}

bool mover::nodeIsReducible() {
	return !isPVNode && !isInCheck;
}

// Returns true if there are still moves remaining, false if we have
void mover::generateMoves() {
    switch (mgStage) {
        // The hash move, if any, is handled separately from the rest of the list
        case STAGE_NONE:
            if (hashed != NULL_MOVE) {
                mgStage = STAGE_HASH_MOVE;

                // Remove the hash move from the list, since it has already been tried
                for (unsigned int i = 0; i < legalMoves.size(); i++) {
                    if (legalMoves.get(i) == hashed) {
                        legalMoves.remove(i);
                        break;
                    }
                }

                break;
            }
       
        // where the quiet moves start in the list, and then do IID or score captures.
        case STAGE_HASH_MOVE:
            findQuietStart();
            if (hashed == NULL_MOVE && doIID()) {
                mgStage = STAGE_IID_MOVE;
                scoreIIDMove();
            }
            else {
                mgStage = STAGE_CAPTURES;
                scoreCaptures(false);
            }
            break;

        // After searching the IID move, we score captures
        case STAGE_IID_MOVE:
            mgStage = STAGE_CAPTURES;
            scoreCaptures(true);
            break;

        // After winnning captures, we score quiets
        case STAGE_CAPTURES:
            mgStage = STAGE_QUIETS;
            scoreQuiets();
            break;

        // We are done
        case STAGE_QUIETS:
            break;
    }
}

// Sort captures using SEE and MVV/LVA
void mover::scoreCaptures(bool isIIDMove) {
    for (unsigned int i = isIIDMove; i < quietStart; i++) {
        Move m = legalMoves.get(i);

        // We want the best move first for PV nodes
        if (isPVNode) {
            int see = b->getSEEForMove(color, m);

            if (see > 0)
                scores.add(score_WINNING_CAPTURE + see);
            else if (see == 0)
                scores.add(score_EVEN_CAPTURE);
            else
                // If we are doing SEE on quiets, score losing captures lower
                scores.add(((depth >= QUIET_SEE_DEPTH) ? 0 : score_LOSING_CAPTURE) + see);
        }

        // Otherwise, MVV/LVA for cheaper cutoffs might help
        else {
            int exchange = b->getExchangescore(color, m);
            if (exchange > 0)
                scores.add(score_WINNING_CAPTURE + b->getMVVLVAscore(color, m));

            else if (exchange == 0)
                scores.add(score_EVEN_CAPTURE + b->getMVVLVAscore(color, m));
            else {
                int see = b->getSEEForMove(color, m);
                if (see > 0)
                    scores.add(score_WINNING_CAPTURE + b->getMVVLVAscore(color, m));
                else if (see == 0)
                    scores.add(score_EVEN_CAPTURE + b->getMVVLVAscore(color, m));
                else
                    scores.add(((depth >= QUIET_SEE_DEPTH) ? 0 : score_LOSING_CAPTURE) + see);
            }
        }
    }
}

void mover::scoreQuiets() {
    for (unsigned int i = quietStart; i < legalMoves.size(); i++) {
        Move m = legalMoves.get(i);

        if (m == parameter->killers[parameter->ply][0])
            scores.add(score_EVEN_CAPTURE - 1);

        else if (getPromotion(m) == QUEENS)
            scores.add(score_QUEEN_PROMO);
        else {
            int startSq = getStartSq(m);
            int endSq = getEndSq(m);
            int pieceID = b->getPieceOnSquare(color, startSq);

            if (depth >= QUIET_SEE_DEPTH) {
                int see = b->getSEEForMove(color, m);
                scores.add(((see < 0) ? score_LOSING_QUIET : score_QUIET_MOVE)
                    + parameter->historyTable[color][pieceID][endSq]);
            }
            else {
                scores.add(score_QUIET_MOVE
                    + parameter->historyTable[color][pieceID][endSq]);
            }
        }
    }
}

bool mover::doIID() {
    return depth >= (isPVNode ? 4 : 6);
}

void mover::scoreIIDMove() {
    // Sort the moves with what we have so far
    for (Move m = nextMove(); m != NULL_MOVE;
              m = nextMove());
    index = 0;

    int iidDepth = isPVNode ? depth-2 : (depth - 5) / 2;
    int bestIndex = getBestMoveForSort(b, legalMoves, iidDepth, threadID);

    // Mate check to prevent crashes
    if (bestIndex == -1) {
        legalMoves.clear();
        mgStage = STAGE_QUIETS;
    }

    else {
        scores.add(score_IID_MOVE);

        if (isCapture(legalMoves.get(bestIndex))) {
            legalMoves.swap(0, bestIndex);
        }
        else {
            legalMoves.swap(quietStart, bestIndex);
            legalMoves.swap(0, quietStart);
            quietStart++;
        }
    }
}

// Retrieves the next move with the highest score, starting from index using a
Move mover::nextMove() {

    if (mgStage == STAGE_HASH_MOVE)
        return hashed;
    if (mgStage == STAGE_IID_MOVE) {
        generateMoves();
        index++;
        return legalMoves.get(0);
    }

    while (index >= scores.size()) {
        if (mgStage == STAGE_QUIETS)
            return NULL_MOVE;
        else {
            generateMoves();
        }
    }

    // Find the index of the next best move
    int bestIndex = index;
    int bestscore = scores.get(index);
    for (unsigned int i = index + 1; i < scores.size(); i++) {
        if (scores.get(i) > bestscore) {
            bestIndex = i;
            bestscore = scores.get(bestIndex);
        }
    }

    // Swap the best move to the correct position
    legalMoves.swap(bestIndex, index);
    scores.swap(bestIndex, index);

    if (mgStage == STAGE_CAPTURES && bestscore < score_WINNING_CAPTURE)
        generateMoves();

    return legalMoves.get(index++);
}

// When a PV or cut move is found, the histories of all
void mover::reduceBadHistories(Move bestMove) {

    if (index <= 0)
        return;
    for (unsigned int i = 0; i < index-1; i++) {
        if (legalMoves.get(i) == bestMove)
            break;
        if (isCapture(legalMoves.get(i)))
            continue;
        parameter->historyTable[color][b->getPieceOnSquare(color, getStartSq(legalMoves.get(i)))][getEndSq(legalMoves.get(i))] -= depth;
    }
}

void mover::findQuietStart() {
    for (unsigned int i = 0; i < legalMoves.size(); i++) {
        if (!isCapture(legalMoves.get(i))) {
            quietStart = i;
            return;
        }
    }

    quietStart = legalMoves.size();
}
