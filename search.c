/*******************************************************************************                                                                    
Joko Tole UCI Chess Engine Copyright(C) 2016 Bosdot (Javanese)
influenzed by chessprogramming.wikispaces.com and many others 
chess engine with open source code.
--------------------------------------------------------------------------------					 					  					  
 ******************************************************************************/
 
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <thread>
#include "zobrist.h"
#include "search.h"
#include "mover.h"
#include "parameter.h"
#include "board.h"

using namespace std;

struct SearchStatistics {
    uint64_t nodes;
    uint64_t hashProbes, hashHits, hashscoreCuts;
    uint64_t hashMoveAttempts, hashMoveCuts;
    uint64_t failHighs, firstFailHighs;
    uint64_t qsNodes;
    uint64_t qsFailHighs, qsFirstFailHighs;
    uint64_t evalCacheProbes, evalCacheHits;

    SearchStatistics() {
        reset();
    }

    void reset() {
        nodes = 0;
        hashProbes = hashHits = hashscoreCuts = 0;
        hashMoveAttempts = hashMoveCuts = 0;
        failHighs = firstFailHighs = 0;
        qsNodes = 0;
        qsFailHighs = qsFirstFailHighs = 0;
        evalCacheProbes = evalCacheHits = 0;
    }
};

struct SearchPV {
    int pvLength;
    Move pv[MAX_DEPTH+1];

    SearchPV() {
        pvLength = 0;
    }
};

struct EasyMove {
    Move prevBest;
    Move streakBest;
    int pvStreak;

    void reset() {
        prevBest = NULL_MOVE;
        streakBest = NULL_MOVE;
        pvStreak = 0;
    }
};

// Futility pruning margins indexed by depth. If static eval is at least this
const int FUTILITY_MARGIN[5] = { 0, 285, 372,402, 352 };

// Reverse futility pruning margins indexed by depth. If static eval is at least
const int REVERSE_FUTILITY_MARGIN[5] = {2, 117,  277, 483, 930 };
	    
// Futility and reductions lookup tables, initialized at startup
int FutilityMoveCounts[2][16];  // [improving][depth]

// Razor margins indexed by depth. If static eval is far below alpha, use a
const int RAZOR_MARGIN[4] = {0, 357, 520, 610 };
    
// Move count pruning
const unsigned int LMP_MOVE_COUNTS[6] = {0, 9, 12, 18, 36, 52 };
    
//variables
static Hash transpositionTable(DEFAULT_HASH_SIZE);
static EvalHash evalCache(DEFAULT_HASH_SIZE);
static SearchParameters parameterArray[MAX_THREADS];
static SearchStatistics searchStatsArray[MAX_THREADS];
static EasyMove easyMoveInfo;

TwoFoldStack twoFoldPositions[MAX_THREADS];
extern bool isStop;
volatile bool stopSignal;
static volatile int threadsRunning;
mutex threadsRunningMutex;
unsigned int multiPV;
int numThreads;

void getBestMoveAtDepth(Board *b, MoveList *legalMoves, int depth, int alpha,
int beta, int *bestMoveIndex, int *bestscore, unsigned int startMove,
int threadID, SearchPV *pvLine);
int PVS(Board &b, int depth, int alpha, int beta, int threadID, SearchPV *pvLine);
int quiescence(Board &b, int plies, int alpha, int beta, int threadID);
int checkQuiescence(Board &b, int plies, int alpha, int beta, int threadID);

// help
int scoreMate(bool isInCheck, int plies);
int adjustHashscore(int score, int plies);

// Other utility functions
Move nextMove(MoveList &moves, scoreList &scores, unsigned int index);
void changePV(Move best, SearchPV *parent, SearchPV *child);
string retrievePV(SearchPV *pvLine);
int getSelectiveDepth();
double getPercentage(uint64_t numerator, uint64_t denominator);
void printStatistics();

// Best Move
void getBestMove(Board *b, int mode, int value, Move *bestMove) {
    for (int i = 0; i < numThreads; i++) {
        parameterArray[i].reset();
        searchStatsArray[i].reset();
        parameterArray[i].rootMoveNumber = (uint8_t) (b->getMoveNumber());
        parameterArray[i].selectiveDepth = 0;
    }

    int color = b->getPlayerToMove();
    MoveList legalMoves = b->getAllLegalMoves(color);

    // Special case if we are given a mate/stalemate position
    if (legalMoves.size() <= 0) {
        *bestMove = NULL_MOVE;
        isStop = true;
        cout << "bestmove none" << endl;
        return;
    }

    *bestMove = legalMoves.get(0);
    
// Set up timing
    parameterArray[0].timeLimit = (mode == TIME) ? (uint64_t) (MAX_TIME_FACTOR * value)
                                                    : (mode == MOVETIME) ? value
                                                                         : MAX_TIME;
    parameterArray[0].startTime = ChessClock::now();
    double timeSoFar = getTimeElapsed(parameterArray[0].startTime);

// legal move
    if (legalMoves.size() == 1 && mode == TIME) {
        parameterArray[0].timeLimit = min(parameterArray[0].timeLimit / 32, ONE_SECOND);
    }

    int bestscore, bestMoveIndex;
    int rootDepth = 1;
    Move prevBest = NULL_MOVE;
    int pvStreak = 0;
    stopSignal = false;

// Iterative deepening loop
    do {
        SearchPV pvLine;

// PVs
        for (unsigned int multiPVNum = 1;
                          multiPVNum <= multiPV
                            && multiPVNum <= legalMoves.size();
                          multiPVNum++) {
            int aspAlpha = -MATE_score;
            int aspBeta = MATE_score;

// Aspiration windows
            if (rootDepth >= 10 && multiPV == 1 && abs(bestscore) < 2 * QUEEN_VALUE) {
                aspAlpha = bestscore - 16;
                aspBeta = bestscore + 16;
            }

            int deltaAlpha = 20;
            int deltaBeta = 20;

// Aspiration looping
            while (!isStop) {
                for (int i = 0; i < numThreads; i++)
                    parameterArray[i].reset();
                pvLine.pvLength = 0;

                if (rootDepth >= 8 && numThreads > 1) {
                    thread *threadPool = new thread[numThreads-1];
                    threadsRunning = numThreads;

// Dummy variables since we don't care about these results
                    int *dummyBestIndex = new int[numThreads-1];
                    int *dummyBestscore = new int[numThreads-1];
                    SearchPV *dummyPVLine = new SearchPV[numThreads-1];

// Secondary threads
                    for (int i = 1; i < numThreads; i++) {
                        twoFoldPositions[i] = twoFoldPositions[0];
                        threadPool[i-1] = thread(getBestMoveAtDepth, b,
                            &legalMoves, rootDepth + (i % 2), aspAlpha, aspBeta,
                            dummyBestIndex+i-1, dummyBestscore+i-1,
                            multiPVNum-1, i, dummyPVLine+i-1);                    
                        threadPool[i-1].detach();
                    }

// Start the primary result thread
                    getBestMoveAtDepth(b, &legalMoves, rootDepth, aspAlpha, aspBeta,
                        &bestMoveIndex, &bestscore, multiPVNum-1, 0, &pvLine);

                    stopSignal = true;
// Wait for all other threads to finish
                    while (threadsRunning > 0)
                        std::this_thread::yield();
                    stopSignal = false;

                    delete[] threadPool;
                    delete[] dummyBestIndex;
                    delete[] dummyBestscore;
                    delete[] dummyPVLine;
                }
// Otherwise, just search with one thread
                else {
                    getBestMoveAtDepth(b, &legalMoves, rootDepth, aspAlpha, aspBeta,
                        &bestMoveIndex, &bestscore, multiPVNum-1, 0, &pvLine);
                }

                timeSoFar = getTimeElapsed(parameterArray[0].startTime);
// Calculate values for printing
                uint64_t nps = (uint64_t) ((double) getNodes() / timeSoFar);
                string pvStr = retrievePV(&pvLine);

// Handle fail highs and fail lows
                if (bestMoveIndex == -1 && !isStop) {
                    cout << "info depth " << rootDepth;
                    cout << " seldepth " << getSelectiveDepth();
                    cout << " score";
                    cout << " cp " << bestscore * 100 / PAWN_VALUE_EG << " upperbound";
                    cout << " time " << (int) (timeSoFar * ONE_SECOND)
                         << " nodes " << getNodes() << " nps " << nps * 3
                         << " hashfull " << 1000 * transpositionTable.keys
                                                 / transpositionTable.getSize()
                         << " pv " << pvStr << endl;

                    aspAlpha = bestscore - deltaAlpha;
                    deltaAlpha *= 2;
                    if (aspAlpha < -2 * QUEEN_VALUE)
                        aspAlpha = -MATE_score;
                }
// Fail high: best score is at least beta
                else if (bestscore >= aspBeta) {
                    cout << "info depth " << rootDepth;
                    cout << " seldepth " << getSelectiveDepth();
                    cout << " score";
                    cout << " cp " << bestscore * 100 / PAWN_VALUE_EG << " lowerbound";
                    cout << " time " << (int) (timeSoFar * ONE_SECOND)
                         << " nodes " << getNodes() << " nps " << nps * 3 
                         << " hashfull " << 5000 * transpositionTable.keys
                                                 / transpositionTable.getSize()
                         << " pv " << pvStr << endl;
                    aspBeta = bestscore + deltaBeta;
                    deltaBeta *= 2;
                    if (aspBeta > 2 * QUEEN_VALUE)
                        aspBeta = MATE_score;
                    legalMoves.swap(multiPVNum-1, bestMoveIndex);
                    *bestMove = legalMoves.get(0);
                }
                else break;
            }
            timeSoFar = getTimeElapsed(parameterArray[0].startTime);
            if (bestMoveIndex == -1)
                break;
// Swap PV to be searched first next iteration
            legalMoves.swap(multiPVNum-1, bestMoveIndex);
            *bestMove = legalMoves.get(0);
            uint64_t nps = (uint64_t) ((double) getNodes() / timeSoFar);
            string pvStr = retrievePV(&pvLine);
            
// Output info using UCI
            cout << "info depth " << rootDepth;
            cout << " seldepth " << getSelectiveDepth();
            if (multiPV > 1)
                cout << " multipv " << multiPVNum;
                cout << " score";			
				
// score in mate
            if (bestscore >= MATE_score - MAX_DEPTH)
                cout << " mate " << (MATE_score - bestscore) / 2 + 1;
            else if (bestscore <= -MATE_score + MAX_DEPTH)
                cout << " mate " << (-MATE_score - bestscore) / 2;
            else
                cout << " cp " << bestscore * 100 / PAWN_VALUE_EG;

            cout << " time " << (int) (timeSoFar * ONE_SECOND)
                 << " nodes " << getNodes() << " nps " << nps * 3 
                 << " hashfull " << 5000 * transpositionTable.keys
                                         / transpositionTable.getSize()
                 << " pv " << pvStr << endl;

// Heuristic table
            for (int i = 0; i < numThreads; i++)
                parameterArray[i].ageHistoryTable(rootDepth, false);
        }

// Record candidate easymoves
        if (multiPV == 1 && pvLine.pvLength >= 3) {
            if (pvLine.pv[2] == easyMoveInfo.streakBest) {
                easyMoveInfo.pvStreak++;
            }
            else {
                easyMoveInfo.streakBest = pvLine.pv[2];
                easyMoveInfo.pvStreak = 1;
            }
        }

        if (*bestMove == prevBest) {
            pvStreak++;
        }
        else {
            prevBest = *bestMove;
            pvStreak = 1;
        }

// Easymove confirmation
        if (mode == TIME && multiPV == 1
         && timeSoFar * ONE_SECOND > value / 16
         && timeSoFar * ONE_SECOND < value / 2
         && abs(bestscore) < MATE_score - MAX_DEPTH) {
            if ((*bestMove == easyMoveInfo.prevBest && pvStreak >= 7)
                || pvStreak >= 10) {
                int secondBestMove;
                int secondBestscore;
                int easymoveWindow = bestscore - EASYMOVE_MARGIN;

                getBestMoveAtDepth(b, &legalMoves, rootDepth-5, easymoveWindow - 1, easymoveWindow,//-MATE_score, MATE_score,
                    &secondBestMove, &secondBestscore, 1, 0, &pvLine);

                if (secondBestscore < easymoveWindow)
                    break;
                else
                    pvStreak = -128;
            }
        }

        rootDepth++;
    }
// Iterative deepening loop
    while (!isStop
        && ((mode == TIME && (timeSoFar * ONE_SECOND < value * TIME_FACTOR)
                          && (rootDepth <= MAX_DEPTH))
         || (mode == MOVETIME && timeSoFar < value)
         || (mode == DEPTH && rootDepth <= value)));

    if (easyMoveInfo.pvStreak >= 8) {
        easyMoveInfo.prevBest = easyMoveInfo.streakBest;
    }
    else
        easyMoveInfo.reset();
    
    printStatistics();

    for (int i = 0; i < numThreads; i++)
        parameterArray[i].ageHistoryTable(rootDepth, true);

    transpositionTable.keys = 0;
    
    // Output best move to UCI interface
    isStop = true;
    cout << "bestmove " << moveToString(*bestMove) << endl;
    return;
}

void getBestMoveAtDepth(Board *b, MoveList *legalMoves, int depth, int alpha,
        int beta, int *bestMoveIndex, int *bestscore, unsigned int startMove,
        int threadID, SearchPV *pvLine) {
    SearchParameters *parameter = &(parameterArray[threadID]);
    SearchStatistics *searchStats = &(searchStatsArray[threadID]);
    SearchPV line;
    int color = b->getPlayerToMove();
    int tempMove = -1;
    int score = -MATE_score;
    *bestscore = -INFTY;
    

    twoFoldPositions[threadID].push(b->getZobristKey());
    std::this_thread::yield();

    for (unsigned int i = startMove; i < legalMoves->size(); i++) {
        
        double timeSoFar = getTimeElapsed(parameterArray[0].startTime);
        if (threadID == 0 && timeSoFar * ONE_SECOND > 5 * ONE_SECOND)
            cout << "info depth " << depth << " currmove " << moveToString(legalMoves->get(i))
                 << " currmovenumber " << i+1 << " nodes " << getNodes() << endl;

        Board copy = b->staticCopy();
        copy.doMove(legalMoves->get(i), color);
        searchStats->nodes++;
        
        if (i != 0) {
            parameter->ply++;
            score = -PVS(copy, depth-1, -alpha-1, -alpha, threadID, &line);
            parameter->ply--;
            if (alpha < score && score < beta) {
                line.pvLength = 0;
                parameter->ply++;
                score = -PVS(copy, depth-1, -beta, -alpha, threadID, &line);
                parameter->ply--;
            }
        }
        else {
            parameter->ply++;
            score = -PVS(copy, depth-1, -beta, -alpha, threadID, &line);
            parameter->ply--;
        }

// Stop condition. If stopping, return search
        if (isStop || stopSignal)
            break;

        if (score > *bestscore) {
            *bestscore = score;
            if (score > alpha) {
                alpha = score;
                tempMove = (int) i;
                changePV(legalMoves->get(i), pvLine, &line);
            }
// To get a PV if failing low
            else if (i == 0)
                changePV(legalMoves->get(i), pvLine, &line);
        }

        if (score >= beta)
            break;
    }

    twoFoldPositions[threadID].pop();

    *bestMoveIndex = tempMove;

// This thread is finished running.
    threadsRunningMutex.lock();
    threadsRunning--;
    threadsRunningMutex.unlock();
}

// Get best move for hash move is not available.
int getBestMoveForSort(Board *b, MoveList &legalMoves, int depth, int threadID) {
    SearchParameters *parameter = &(parameterArray[threadID]);
    SearchStatistics *searchStats = &(searchStatsArray[threadID]);
    SearchPV line;
    int color = b->getPlayerToMove();
    int tempMove = -1;
    int score = -MATE_score;
    int alpha = -MATE_score;
    int beta = MATE_score;

// Push current position to two fold stack
    twoFoldPositions[threadID].push(b->getZobristKey());
    
    for (unsigned int i = 0; i < legalMoves.size(); i++) {
        Board copy = b->staticCopy();
        if(!copy.doPseudoLegalMove(legalMoves.get(i), color))
            continue;
        searchStats->nodes++;
        
        if (i != 0) {
            parameter->ply++;
            score = -PVS(copy, depth-1, -alpha-1, -alpha, threadID, &line);
            parameter->ply--;
            if (alpha < score && score < beta) {
                parameter->ply++;
                score = -PVS(copy, depth-1, -beta, -alpha, threadID, &line);
                parameter->ply--;
            }
        }
        else {
            parameter->ply++;
            score = -PVS(copy, depth-1, -beta, -alpha, threadID, &line);
            parameter->ply--;
        }

// Stop condition as fast
        if (isStop || stopSignal)
            return i;
        
        if (score > alpha) {
            alpha = score;
            tempMove = i;
        }
    }

    twoFoldPositions[threadID].pop();

    return tempMove;
}

// The standard implementation of a fail-soft PVS search.
int PVS(Board &b, int depth, int alpha, int beta, int threadID, SearchPV *pvLine) {
    SearchParameters *parameter = &(parameterArray[threadID]);
    SearchStatistics *searchStats = &(searchStatsArray[threadID]);

    // Static board evaluation is done there.
    if (depth <= 0 || parameter->ply >= MAX_DEPTH) {

        if (parameter->ply > parameter->selectiveDepth)
            parameter->selectiveDepth = parameter->ply;
        pvLine->pvLength = 0;
        return quiescence(b, 0, alpha, beta, threadID);
    }

    if (b.isDraw())
        return 0;
    if (twoFoldPositions[threadID].find(b.getZobristKey()))
        return 0;

// Mate distance pruning
    int matingscore = MATE_score - parameter->ply;
    if (matingscore < beta) {
        beta = matingscore;
        if (alpha >= matingscore)
            return alpha;
    }

    int matedscore = -MATE_score + parameter->ply;
    if (matedscore > alpha) {
        alpha = matedscore;
        if (beta <= matedscore)
            return beta;
    }
       
    int prevAlpha = alpha;
    int color = b.getPlayerToMove();
    bool isPVNode = (beta - alpha != 1);

// Transposition table probe
    Move hashed = NULL_MOVE;
    int hashscore = -INFTY;
    int hashDepth = 0;
    uint8_t nodeType = NO_NODE_INFO;
    searchStats->hashProbes++;

    uint64_t hashEntry = transpositionTable.get(b);
    if (hashEntry != 0) {
        searchStats->hashHits++;
        hashscore = getHashscore(hashEntry);
        nodeType = getHashNodeType(hashEntry);
        hashDepth = getHashDepth(hashEntry);

// Adjust the hash score to mate distance from root if necessary
        if (hashscore >= MATE_score - MAX_DEPTH)
            hashscore -= parameter->ply;
        else if (hashscore <= -MATE_score + MAX_DEPTH)
            hashscore += parameter->ply;

 // Return hash score failing soft if hash depth >= current depth and:
        if (!isPVNode && hashDepth >= depth) {
            if ((nodeType == ALL_NODE && hashscore <= alpha)
             || (nodeType == CUT_NODE && hashscore >= beta)
             || (nodeType == PV_NODE)) {
                searchStats->hashscoreCuts++;
                return hashscore;
            }
        }

// Otherwise, store the hash move for later use
        hashed = getHashMove(hashEntry);
        // Don't use hash moves from too low of a depth
        if ((hashDepth < 1 && depth >= 7)
         || (isPVNode && hashDepth < depth - 3))
            hashed = NULL_MOVE;
    }

// Keeps PV  up to root
    SearchPV line;
    bool isInCheck = b.isInCheck(color);
	
// A static evaluation, used to activate null move pruning and futility
    int staticEval = INFTY;
    if (!isInCheck) {
        searchStats->evalCacheProbes++;
		
// Probe the eval cache for a saved calculation
        int ehe = evalCache.get(b);
        if (ehe != 0) {
            searchStats->evalCacheHits++;
            staticEval = ehe - EVAL_HASH_OFFSET;
        }
        else {
            staticEval = (color == WHITE) ? b.evaluate()
                                          : -b.evaluate();
            evalCache.add(b, staticEval);
        }
    }

// TT score as a better "static" eval, if available.
    if (hashscore != -INFTY) {
        if ((nodeType == ALL_NODE && hashscore < staticEval)
         || (nodeType == CUT_NODE && hashscore > staticEval)
         || (nodeType == PV_NODE))
            staticEval = hashscore;
    }
    
// Reverse futility pruning
    if (!isPVNode && !isInCheck
     && (depth <= 4 && staticEval - REVERSE_FUTILITY_MARGIN[depth] >= beta)
     && b.getNonPawnMaterial(color))
        return staticEval - REVERSE_FUTILITY_MARGIN[depth];

// Razoring
    if (!isPVNode && !isInCheck && abs(alpha) < 2 * QUEEN_VALUE
     && depth <= 3 && staticEval <= alpha - RAZOR_MARGIN[depth]) {
        if (depth == 1)
            return quiescence(b, 0, alpha, beta, threadID);

        int value = quiescence(b, 0, alpha, beta, threadID);
        // Fail hard here to be safe
        if (value <= alpha)
            return alpha;
    }

// Null move pruning
    if (!isPVNode && !isInCheck
     && depth >= 2 && staticEval >= beta
     && parameter->nullMoveCount < 2
     && b.getNonPawnMaterial(color)) {
        int reduction;
		
// Reduce
        reduction = 1 + (int) ((depth + 1.5) / 4.5 + (staticEval - beta) / 118.0);

        uint16_t epCaptureFile = b.getEPCaptureFile();
        b.doNullMove();
        parameter->nullMoveCount++;
        parameter->ply++;
        int nullscore = -PVS(b, depth-1-reduction, -beta, -alpha, threadID, &line);
        parameter->ply--;

// Undo the null move
        b.undoNullMove(epCaptureFile);
        parameter->nullMoveCount = 0;

        if (nullscore >= beta) {
            return nullscore;
        }
    }

// Move generation
    PieceMoveList pml = b.getPieceMoveList<PML_LEGAL_MOVES>(color);
    MoveList legalMoves = isInCheck ? b.getPseudoLegalCheckEscapes(color, pml)
                                    : b.getAllPseudoLegalMoves(color, pml);
    mover moveSorter(&b, color, depth, threadID, isPVNode, isInCheck,
        parameter, hashed, legalMoves);
    moveSorter.generateMoves();

    Move toHash = NULL_MOVE;
	
// separate counter only incremented when valid move is searched
    unsigned int movesSearched = 0;
    int bestscore = -INFTY;
    int score = -INFTY;

//Search loop
    for (Move m = moveSorter.nextMove(); m != NULL_MOVE;
              m = moveSorter.nextMove()) {

        double timeSoFar = getTimeElapsed(parameterArray[0].startTime);
        if (timeSoFar * ONE_SECOND > parameterArray[0].timeLimit)
            isStop = stopSignal = true;
        // Stop condition to help break out as quickly as possible
        if (isStop || stopSignal)
            return INFTY;

        bool moveIsPrunable = moveSorter.nodeIsReducible()
                           && !isCapture(m)
                           && !isPromotion(m)
                           && m != hashed
                           && abs(alpha) < 2 * QUEEN_VALUE
                           && !b.isCheckMove(color, m);

        int startSq = getStartSq(m);
        int endSq = getEndSq(m);
        int pieceID = b.getPieceOnSquare(color, startSq);

        // Futility pruning
        if (moveIsPrunable
         && depth <= 4 && staticEval <= alpha - FUTILITY_MARGIN[depth]) {
            if (bestscore < staticEval + FUTILITY_MARGIN[depth])
                bestscore = staticEval + FUTILITY_MARGIN[depth];
            continue;
        }

        if(moveIsPrunable
        && depth == 1 && staticEval <= alpha
        && ((!isCapture(m) && b.getSEEForMove(color, m) < 0)
         || (isCapture(m) && b.getExchangescore(color, m) < 0 && b.getSEEForMove(color, m) < -200))) {
            score = alpha;
            continue;
        }

        // As used in Fruit/Stockfish:
        // https://chessprogramming.wikispaces.com/Futility+Pruning#MoveCountBasedPruning
        if (moveIsPrunable
         && depth <= 5 && movesSearched > LMP_MOVE_COUNTS[depth]
         && alpha <= prevAlpha
         && m != parameter->killers[parameter->ply][0]
         && m != parameter->killers[parameter->ply][1]) {
            int historyValue = parameter->historyTable[color][pieceID][endSq];
            if (depth < 3 || historyValue < 0) {
                if (bestscore < alpha)
                    bestscore = alpha;
                continue;
            }
        }

        Board copy = b.staticCopy();

        if (m == hashed) {
            if (!copy.doHashMove(m, color)) {
                hashed = NULL_MOVE;
                moveSorter.hashed = NULL_MOVE;
                moveSorter.generateMoves();
                continue;
            }
            moveSorter.generateMoves();
        }
        else if (!copy.doPseudoLegalMove(m, color))
            continue;
        searchStats->nodes++;

        int reduction = 0;

        if (moveSorter.nodeIsReducible()
         && depth >= 3 && movesSearched > 2 && alpha <= prevAlpha
         && !isCapture(m) && !isPromotion(m)
         && m != parameter->killers[parameter->ply][0]
         && m != parameter->killers[parameter->ply][1]
         && !copy.isInCheck(color^1)) {
            // Increase reduction with higher depth and later moves
            reduction = 1 + (int) (sqrt((depth-3) * movesSearched) / 4.0);
            // Reduce more for moves with poor history
            int historyValue = parameter->historyTable[color][pieceID][endSq];
            if (historyValue < 0)
                reduction++;

            // Do not let search descend directly into q-search
            reduction = min(reduction, depth - 2);
            reduction = min(reduction, 1 + (int) (movesSearched - 3) / 2);
        }

        int extension = 0;
        // Check extensions
        if (reduction == 0
         && copy.isInCheck(color^1)
         && (isCapture(m) || b.getSEEForMove(color, m) >= 0)) {
            extension++;
        }

        // Extension for transition into pawn endgame
        if (depth >= 3 && reduction == 0) {
            uint64_t nonPawns = b.getNonPawnMaterial(color^1);
            if (INDEX_TO_BIT[getEndSq(m)] == nonPawns) {
                extension += 1;
                if (!b.getNonPawnMaterial(color))
                    extension += 2;
            }
        }

        // Record two fold stack 
        twoFoldPositions[threadID].push(b.getZobristKey());

        // Singular extensions
        if (depth >= 6 && reduction == 0 && extension == 0
         && m == hashed
         && abs(hashscore) < 2 * QUEEN_VALUE
         && ((hashscore >= beta && (nodeType == CUT_NODE || nodeType == PV_NODE)
                                && hashDepth >= depth - 4)
          || (isPVNode && nodeType == PV_NODE && hashDepth >= depth - 2))) {

            bool isSingular = true;

            // Do a reduced depth search fail low check
            for (unsigned int i = 0; i < legalMoves.size(); i++) {
                Move seMove = legalMoves.get(i);
                Board seCopy = b.staticCopy();
                // Search every move except the hash move
                if (seMove == hashed)
                    continue;
                if (!seCopy.doPseudoLegalMove(seMove, color))
                    continue;

                // The window is lowered more for PV nodes and for higher depths
                int SEWindow = isPVNode ? hashscore - 50 - 2 * depth
                                        : alpha - 10 - depth;
                // Do a reduced search for fail-low confirmation
                int SEDepth = isPVNode ? 2 * depth / 3 - 1
                                       : depth / 2 - 1;

                parameter->ply++;
                score = -PVS(seCopy, SEDepth, -SEWindow - 1, -SEWindow, threadID, &line);
                parameter->ply--;

                // If a move did not fail low, no singular extension
                if (score > SEWindow) {
                    isSingular = false;
                    break;
                }
            }

            // If all moves other than the hash move failed low, we extend for
            if (isSingular)
                extension++;
        }

        // Reset the PV line just in case
        line.pvLength = 0;

        // Null-window search, with re-search if applicable
        if (movesSearched != 0) {
            parameter->ply++;
            score = -PVS(copy, depth-1-reduction+extension, -alpha-1, -alpha, threadID, &line);
            parameter->ply--;

            // LMR re-search if the reduced search did not fail low
            if (reduction > 0 && score > alpha) {
                line.pvLength = 0;
                parameter->ply++;
                score = -PVS(copy, depth-1+extension, -alpha-1, -alpha, threadID, &line);
                parameter->ply--;
            }

            // Re-search for a scout window at PV nodes
            else if (alpha < score && score < beta) {
                line.pvLength = 0;
                parameter->ply++;
                score = -PVS(copy, depth-1+extension, -beta, -alpha, threadID, &line);
                parameter->ply--;
            }
        }

        // The first move is always searched at a normal depth
        else {
            parameter->ply++;
            score = -PVS(copy, depth-1+extension, -beta, -alpha, threadID, &line);
            parameter->ply--;
        }

        // Pop the position in case we return early from this search
        twoFoldPositions[threadID].pop();

        // Stop condition to help break out as quickly as possible
        if (isStop || stopSignal)
            return INFTY;
        
        // Beta cutoff
        if (score >= beta) {
            searchStats->failHighs++;
            if (movesSearched == 0)
                searchStats->firstFailHighs++;
            if (hashed != NULL_MOVE && nodeType != ALL_NODE) {
                searchStats->hashMoveAttempts++;
                if (m == hashed)
                    searchStats->hashMoveCuts++;
            }

            // Hash the cut move and score
            uint64_t hashData = packHashData(depth, m,
                adjustHashscore(score, parameter->ply), CUT_NODE,
                parameter->rootMoveNumber);
            transpositionTable.add(b, hashData, depth, parameter->rootMoveNumber);

            // Record killer if applicable
            if (!isCapture(m)) {

                if (m != parameter->killers[parameter->ply][0]) {
                    parameter->killers[parameter->ply][1] =
                        parameter->killers[parameter->ply][0];
                    parameter->killers[parameter->ply][0] = m;
                }

                parameter->historyTable[color][pieceID][endSq]
                    += depth * depth;
                moveSorter.reduceBadHistories(m);
            }

            changePV(m, pvLine, &line);

            return score;
        }

        if (score > bestscore) {
            bestscore = score;
            if (score > alpha) {
                alpha = score;
                toHash = m;
                changePV(m, pvLine, &line);
            }
        }

        movesSearched++;
    }
    // End main search loop
    if (bestscore == -INFTY && movesSearched == 0)
        return scoreMate(moveSorter.isInCheck, parameter->ply);
    
    // Exact scores indicate a principal variation
    if (prevAlpha < alpha && alpha < beta) {
        if (hashed != NULL_MOVE && nodeType != ALL_NODE) {
            searchStats->hashMoveAttempts++;
            if (toHash == hashed)
                searchStats->hashMoveCuts++;
        }

        uint64_t hashData = packHashData(depth, toHash,
            adjustHashscore(alpha, parameter->ply), PV_NODE,
            parameter->rootMoveNumber);
        transpositionTable.add(b, hashData, depth, parameter->rootMoveNumber);

        // Update the history table
        if (!isCapture(toHash)) {
            parameter->historyTable[color][b.getPieceOnSquare(color, getStartSq(toHash))][getEndSq(toHash)]
                += depth * depth;
            moveSorter.reduceBadHistories(toHash);
        }
    }

    // Record all-nodes
    else if (alpha <= prevAlpha) {
        
		// If we would have done IID
        if (!isPVNode && moveSorter.doIID()) {
            uint64_t hashData = packHashData(depth,
                (hashed == NULL_MOVE) ? moveSorter.legalMoves.get(0) : hashed,
                adjustHashscore(bestscore, parameter->ply), ALL_NODE,
                parameter->rootMoveNumber);
            transpositionTable.add(b, hashData, depth, parameter->rootMoveNumber);
        }

        else {
            uint64_t hashData = packHashData(depth, NULL_MOVE,
                adjustHashscore(bestscore, parameter->ply), ALL_NODE,
                parameter->rootMoveNumber);
            transpositionTable.add(b, hashData, depth, parameter->rootMoveNumber);
        }
    }

    return bestscore;
}

int quiescence(Board &b, int plies, int alpha, int beta, int threadID) {
    SearchParameters *parameter = &(parameterArray[threadID]);
    SearchStatistics *searchStats = &(searchStatsArray[threadID]);
    int color = b.getPlayerToMove();

    // legal check 
    if (b.isInCheck(color))
        return checkQuiescence(b, plies, alpha, beta, threadID);

    if (b.isInsufficientMaterial())
        return 0;
    if (plies <= 2 && twoFoldPositions[threadID].find(b.getZobristKey()))
        return 0;

    // search hash table probe
    uint64_t hashEntry = transpositionTable.get(b);
    if (hashEntry != 0) {
        int hashscore = getHashscore(hashEntry);

        // Adjust the hash score
        if (hashscore >= MATE_score - MAX_DEPTH)
            hashscore -= parameter->ply + plies;
        else if (hashscore <= -MATE_score + MAX_DEPTH)
            hashscore += parameter->ply + plies;

        uint8_t nodeType = getHashNodeType(hashEntry);

        if (getHashDepth(hashEntry) >= -plies) {
            // Check for the correct node type and bounds
            if ((nodeType == ALL_NODE && hashscore <= alpha)
             || (nodeType == CUT_NODE && hashscore >= beta)
             || (nodeType == PV_NODE))
                return hashscore;
        }
    }

    // Stand pat: if our current position is already way too good or way too bad
    int standPat;
    // Probe the eval cache for a saved calculation
    searchStats->evalCacheProbes++;
    int ehe = evalCache.get(b);
    if (ehe != 0) {
        searchStats->evalCacheHits++;
        standPat = ehe - EVAL_HASH_OFFSET;
    }
    else {
        standPat = (color == WHITE) ? b.evaluate() : -b.evaluate();
        evalCache.add(b, standPat);
    }
    
    // The stand pat cutoff
    if (standPat >= beta || standPat < alpha - MAX_POS_score - QUEEN_VALUE)
        return standPat;

    if (alpha < standPat)
        alpha = standPat;

    // Generate captures and order by MVV/LVA
    PieceMoveList pml = b.getPieceMoveList<PML_LEGAL_MOVES>(color);
    MoveList legalCaptures = b.getPseudoLegalCaptures(color, pml, false);
    scoreList scores;
    for (unsigned int i = 0; i < legalCaptures.size(); i++) {
        scores.add(b.getMVVLVAscore(color, legalCaptures.get(i)));
    }
    
    int bestscore = -INFTY;
    int score = -INFTY;
    unsigned int i = 0;
    unsigned int j = 0; 
    for (Move m = nextMove(legalCaptures, scores, i); m != NULL_MOVE;
              m = nextMove(legalCaptures, scores, ++i)) {
        // Delta prune
        if (standPat + b.valueOfPiece(b.getPieceOnSquare(color^1, getEndSq(m))) < alpha - MAX_POS_score)
            continue;
        // Static exchange evaluation pruning
        if (b.getExchangescore(color, m) < 0 && b.getSEEForMove(color, m) < -MAX_POS_score)
            continue;
        
        Board copy = b.staticCopy();
        if (!copy.doPseudoLegalMove(m, color))
            continue;
        
        searchStats->nodes++;
        searchStats->qsNodes++;
        score = -quiescence(copy, plies+1, -beta, -alpha, threadID);
        
        if (score >= beta) {
            searchStats->qsFailHighs++;
            if (j == 0)
                searchStats->qsFirstFailHighs++;

            uint64_t hashData = packHashData(-plies, m,
                adjustHashscore(score, parameter->ply + plies), CUT_NODE,
                parameter->rootMoveNumber);
            transpositionTable.add(b, hashData, -plies, parameter->rootMoveNumber);

            return score;
        }

        if (score > bestscore) {
            bestscore = score;
            if (score > alpha)
                alpha = score;
        }

        j++;
    }

    // Generate and search promotions
    MoveList legalPromotions = b.getPseudoLegalPromotions(color);
    for (unsigned int i = 0; i < legalPromotions.size(); i++) {
        Move m = legalPromotions.get(i);

        // Static exchange evaluation pruning
        if (b.getSEEForMove(color, m) < 0)
            continue;

        Board copy = b.staticCopy();
        if (!copy.doPseudoLegalMove(m, color))
            continue;
        
        searchStats->nodes++;
        searchStats->qsNodes++;
        score = -quiescence(copy, plies+1, -beta, -alpha, threadID);
        
        if (score >= beta) {
            searchStats->qsFailHighs++;
            if (j == 0)
                searchStats->qsFirstFailHighs++;

            uint64_t hashData = packHashData(-plies, m,
                adjustHashscore(score, parameter->ply + plies), CUT_NODE,
                parameter->rootMoveNumber);
            transpositionTable.add(b, hashData, -plies, parameter->rootMoveNumber);

            return score;
        }

        if (score > bestscore) {
            bestscore = score;
            if (score > alpha)
                alpha = score;
        }

        j++;
    }

    // Checks: only on the first two plies of q-search
    if (plies <= 1) {
        MoveList legalMoves = b.getPseudoLegalChecks(color);

        for (unsigned int i = 0; i < legalMoves.size(); i++) {
            Move m = legalMoves.get(i);

            // Static exchange evaluation pruning
            if (b.getSEEForMove(color, m) < 0)
                continue;

            Board copy = b.staticCopy();
            if (!copy.doPseudoLegalMove(m, color))
                continue;
            
            searchStats->nodes++;
            searchStats->qsNodes++;
            twoFoldPositions[threadID].push(b.getZobristKey());

            int score = -checkQuiescence(copy, plies+1, -beta, -alpha, threadID);
            
            twoFoldPositions[threadID].pop();

            if (score >= beta) {
                searchStats->qsFailHighs++;
                if (j == 0)
                    searchStats->qsFirstFailHighs++;

                uint64_t hashData = packHashData(-plies, m,
                    adjustHashscore(score, parameter->ply + plies), CUT_NODE,
                    parameter->rootMoveNumber);
                transpositionTable.add(b, hashData, -plies, parameter->rootMoveNumber);

                return score;
            }

            if (score > bestscore) {
                bestscore = score;
                if (score > alpha)
                    alpha = score;
            }

            j++;
        }
    }

    // Fail low and hard if there are no captures
    if (bestscore == -INFTY)
        bestscore = alpha;

    return bestscore;
}

int checkQuiescence(Board &b, int plies, int alpha, int beta, int threadID) {
    if (twoFoldPositions[threadID].find(b.getZobristKey()))
        return 0;

    SearchParameters *parameter = &(parameterArray[threadID]);
    SearchStatistics *searchStats = &(searchStatsArray[threadID]);
    int color = b.getPlayerToMove();
    PieceMoveList pml = b.getPieceMoveList<PML_LEGAL_MOVES>(color);
    MoveList legalMoves = b.getPseudoLegalCheckEscapes(color, pml);

    int bestscore = -INFTY;
    int score = -INFTY;
    unsigned int j = 0; // separate counter only incremented when valid move is searched
    for (unsigned int i = 0; i < legalMoves.size(); i++) {
        Move m = legalMoves.get(i);

        if (bestscore > -INFTY && b.getSEEForMove(color, m) < 0)
            continue;

        Board copy = b.staticCopy();
        if (!copy.doPseudoLegalMove(m, color))
            continue;
        
        searchStats->nodes++;
        searchStats->qsNodes++;
        twoFoldPositions[threadID].push(b.getZobristKey());

        score = -quiescence(copy, plies+1, -beta, -alpha, threadID);
        
        twoFoldPositions[threadID].pop();

        if (score >= beta) {
            searchStats->qsFailHighs++;
            if (j == 0)
                searchStats->qsFirstFailHighs++;
            return score;
        }

        if (score > bestscore) {
            bestscore = score;
            if (score > alpha)
                alpha = score;
        }

        j++;
    }

    if (bestscore == -INFTY) {
        return (-MATE_score + parameter->ply + plies);
    }
    
    return bestscore;
}

// Used to get a score when we have realized that we have no legal moves.
int scoreMate(bool isInCheck, int plies) {
    // If we are in check, then it is a checkmate
    if (isInCheck)
        // Adjust score so that quicker mates are better
        return (-MATE_score + plies);

    // Else, it is a stalemate
    else return 0;
}

// Adjust a mate score to accurately reflect distance to mate from the
int adjustHashscore(int score, int plies) {
    if (score >= MATE_score - MAX_DEPTH)
        return score + plies;
    if (score <= -MATE_score + MAX_DEPTH)
        return score - plies;
    return score;
}

// Used in uci.c
void clearTables() {
    transpositionTable.clear();
    evalCache.clear();
    for (int i = 0; i < MAX_THREADS; i++)
        parameterArray[i].resetHistoryTable();
}

void setHashSize(uint64_t MB) {
    transpositionTable.setSize(MB);
}

void setEvalCacheSize(uint64_t MB) {
    evalCache.setSize(MB);
}

uint64_t getNodes() {
    uint64_t total = 0;
    for (int i = 0; i < numThreads; i++) {
        total += searchStatsArray[i].nodes;
    }
    return total;
}

void setMultiPV(unsigned int n) {
    multiPV = n;
}

void setNumThreads(int n) {
    numThreads = n;
}

// Retrieves the next move with the highest score, starting from index using a
Move nextMove(MoveList &moves, scoreList &scores, unsigned int index) {
    if (index >= moves.size())
        return NULL_MOVE;
    // Find the index of the next best move
    int bestIndex = index;
    int bestscore = scores.get(index);
    for (unsigned int i = index + 1; i < moves.size(); i++) {
        if (scores.get(i) > bestscore) {
            bestIndex = i;
            bestscore = scores.get(bestIndex);
        }
    }

	// Swap the best move to the correct position
    moves.swap(bestIndex, index);
    scores.swap(bestIndex, index);
    return moves.get(index);
}

// New PV line when alpha is raised
void changePV(Move best, SearchPV *parent, SearchPV *child) {
    parent->pv[0] = best;
    for (int i = 0; i < child->pvLength; i++) {
        parent->pv[i+1] = child->pv[i];
    }
    parent->pvLength = child->pvLength + 1;
}

// Recover PV for output
string retrievePV(SearchPV *pvLine) {
    string pvStr = moveToString(pvLine->pv[0]);
    for (int i = 1; i < pvLine->pvLength; i++) {
        pvStr += " " + moveToString(pvLine->pv[i]);
    }

    return pvStr;
}

// The selective depth in a parallel search is the max
int getSelectiveDepth() {
    int max = 0;
    for (int i = 0; i < numThreads; i++)
        if (parameterArray[i].selectiveDepth > max)
            max = parameterArray[i].selectiveDepth;
    return max;
}

// Formats a fraction into a percentage value (0 to 100)
double getPercentage(uint64_t numerator, uint64_t denominator) {
    if (denominator == 0)
        return 0;
    uint64_t tenThousandths = (numerator * 10000) / denominator;
    double percent = ((double) tenThousandths) / 100.0;
    return percent;
}

// Prints the statistics gathered during search
void printStatistics() {
    // Aggregate statistics over all threads
    SearchStatistics searchStats;
    for (int i = 0; i < numThreads; i++) {
        searchStats.nodes +=            searchStatsArray[i].nodes;
        searchStats.hashProbes +=       searchStatsArray[i].hashProbes;
        searchStats.hashHits +=         searchStatsArray[i].hashHits;
        searchStats.hashscoreCuts +=    searchStatsArray[i].hashscoreCuts;
        searchStats.hashMoveAttempts += searchStatsArray[i].hashMoveAttempts;
        searchStats.hashMoveCuts +=     searchStatsArray[i].hashMoveCuts;
        searchStats.failHighs +=        searchStatsArray[i].failHighs;
        searchStats.firstFailHighs +=   searchStatsArray[i].firstFailHighs;
        searchStats.qsNodes +=          searchStatsArray[i].qsNodes;
        searchStats.qsFailHighs +=      searchStatsArray[i].qsFailHighs;
        searchStats.qsFirstFailHighs += searchStatsArray[i].qsFirstFailHighs;
        searchStats.evalCacheProbes +=  searchStatsArray[i].evalCacheProbes;
        searchStats.evalCacheHits +=    searchStatsArray[i].evalCacheHits;
    }

    cerr << setw(22) << "Hash hit rate: " << getPercentage(searchStats.hashHits, searchStats.hashProbes)
         << '%' << " of " << searchStats.hashProbes << " probes" << endl;
    cerr << setw(22) << "Hash score cut rate: " << getPercentage(searchStats.hashscoreCuts, searchStats.hashHits)
         << '%' << " of " << searchStats.hashHits << " hash hits" << endl;
    cerr << setw(22) << "Hash move cut rate: " << getPercentage(searchStats.hashMoveCuts, searchStats.hashMoveAttempts)
         << '%' << " of " << searchStats.hashMoveAttempts << " hash moves" << endl;
    cerr << setw(22) << "First fail high rate: " << getPercentage(searchStats.firstFailHighs, searchStats.failHighs)
         << '%' << " of " << searchStats.failHighs << " fail highs" << endl;
    cerr << setw(22) << "QS Nodes: " << searchStats.qsNodes << " ("
         << getPercentage(searchStats.qsNodes, searchStats.nodes) << '%' << ")" << endl;
    cerr << setw(22) << "QS FFH rate: " << getPercentage(searchStats.qsFirstFailHighs, searchStats.qsFailHighs)
         << '%' << " of " << searchStats.qsFailHighs << " qs fail highs" << endl;
    cerr << setw(22) << "Eval cache hit rate: " << getPercentage(searchStats.evalCacheHits, searchStats.evalCacheProbes)
         << '%' << " of " << searchStats.evalCacheProbes << " probes" << endl;
}
