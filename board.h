/*******************************************************************************                                                                    
Joko Tole Chess Engine Copyright (C) 2016 BOSDOT						
experimet code from chessprogramming.wikispaces.com and many others engine 
stockfish, Gull, glaurung, fruit, etc    
as long as no code is not fail at build/ compile engine can be used
--------------------------------------------------------------------------------					 					  					  
 ******************************************************************************/
 
#ifndef BOARD_H
#define BOARD_H

#include "general.h"

typedef enum
    {
	A1, B1, C1, D1, E1, F1, G1, H1, A2, B2, C2, D2, E2, F2, G2, H2,
	A3, B3, C3, D3, E3, F3, G3, H3, A4, B4, C4, D4, E4, F4, G4, H4,
	A5, B5, C5, D5, E5, F5, G5, H5, A6, B6, C6, D6, E6, F6, G6, H6,
	A7, B7, C7, D7, E7, F7, G7, H7, A8, B8, C8, D8, E8, F8, G8, H8
    } EnumSquares;

	
const uint64_t AFILE = 0x0101010101010101;
const uint64_t HFILE = 0x8080808080808080;
const uint64_t FILES[8] = {
    0x0101010101010101, 0x0202020202020202, 0x0404040404040404,
    0x0808080808080808, 0x1010101010101010, 0x2020202020202020,
    0x4040404040404040, 0x8080808080808080
};
const uint64_t RANKS[8] = {
    0x00000000000000FF, 0x000000000000FF00, 0x0000000000FF0000,
    0x00000000FF000000, 0x000000FF00000000, 0x0000FF0000000000,
    0x00FF000000000000, 0xFF00000000000000
};
const uint64_t DIAGONAL = 0x8040201008040201;
const uint64_t ANTIDIAGONAL = 0x0102040810204080;
const uint64_t LIGHT = 0x55AA55AA55AA55AA;
const uint64_t DARK = 0xAA55AA55AA55AA55;
const uint64_t NOTA = 0xFEFEFEFEFEFEFEFE;
const uint64_t NOTH = 0x7F7F7F7F7F7F7F7F;

// Converts square number to 64-bit integer
const uint64_t INDEX_TO_BIT[64] = {
0x0000000000000001, 0x0000000000000002, 0x0000000000000004, 0x0000000000000008,
0x0000000000000010, 0x0000000000000020, 0x0000000000000040, 0x0000000000000080,
0x0000000000000100, 0x0000000000000200, 0x0000000000000400, 0x0000000000000800,
0x0000000000001000, 0x0000000000002000, 0x0000000000004000, 0x0000000000008000,
0x0000000000010000, 0x0000000000020000, 0x0000000000040000, 0x0000000000080000,
0x0000000000100000, 0x0000000000200000, 0x0000000000400000, 0x0000000000800000,
0x0000000001000000, 0x0000000002000000, 0x0000000004000000, 0x0000000008000000,
0x0000000010000000, 0x0000000020000000, 0x0000000040000000, 0x0000000080000000,
0x0000000100000000, 0x0000000200000000, 0x0000000400000000, 0x0000000800000000,
0x0000001000000000, 0x0000002000000000, 0x0000004000000000, 0x0000008000000000,
0x0000010000000000, 0x0000020000000000, 0x0000040000000000, 0x0000080000000000,
0x0000100000000000, 0x0000200000000000, 0x0000400000000000, 0x0000800000000000,
0x0001000000000000, 0x0002000000000000, 0x0004000000000000, 0x0008000000000000,
0x0010000000000000, 0x0020000000000000, 0x0040000000000000, 0x0080000000000000,
0x0100000000000000, 0x0200000000000000, 0x0400000000000000, 0x0800000000000000,
0x1000000000000000, 0x2000000000000000, 0x4000000000000000, 0x8000000000000000
};

const uint8_t WHITEKSIDE = 0x1;
const uint8_t WHITEQSIDE = 0x2;
const uint8_t BLACKKSIDE = 0x4;
const uint8_t BLACKQSIDE = 0x8;
const uint8_t WHITECASTLE = 0x3;
const uint8_t BLACKCASTLE = 0xC;
const uint64_t WHITE_KSIDE_PASSTHROUGH_SQS = INDEX_TO_BIT[5] | INDEX_TO_BIT[6];
const uint64_t WHITE_QSIDE_PASSTHROUGH_SQS = INDEX_TO_BIT[1] | INDEX_TO_BIT[2] | INDEX_TO_BIT[3];
const uint64_t BLACK_KSIDE_PASSTHROUGH_SQS = INDEX_TO_BIT[61] | INDEX_TO_BIT[62];
const uint64_t BLACK_QSIDE_PASSTHROUGH_SQS = INDEX_TO_BIT[57] | INDEX_TO_BIT[58] | INDEX_TO_BIT[59];
const uint16_t NO_EP_POSSIBLE = 0x8;
const bool PML_LEGAL_MOVES = true;
const bool PML_PSEUDO_MOBILITY = false;
const bool MOVEGEN_CAPTURES = true;
const bool MOVEGEN_QUIETS = false;

struct PieceMoveInfo {
    int pieceID;
    int startSq;
    uint64_t legal;

    PieceMoveInfo() {}
    PieceMoveInfo(int _pieceID, int _startSq, uint64_t _legal) {
        pieceID = _pieceID;
        startSq = _startSq;
        legal = _legal;
    }
};

typedef SearchArrayList<PieceMoveInfo> PieceMoveList;
void initMagicTables(uint64_t seed);
void initZobristTable();
void initInBetweenTable();

class Board {
public:

    // Getter methods
    bool getWhiteCanKCastle();
    bool getBlackCanKCastle();
    bool getWhiteCanQCastle();
    bool getBlackCanQCastle();
    uint16_t getEPCaptureFile();
    uint8_t getFiftyMoveCounter();
    uint16_t getMoveNumber();
    int getPlayerToMove();
    uint64_t getPieces(int color, int piece);
    uint64_t getAllPieces(int color);
    int *getMailbox();
    uint64_t getZobristKey();
    void initZobristKey(int *mailbox);
	
    Board();
    Board(int *mailboxBoard, bool _whiteCanKCastle, bool _blackCanKCastle,
            bool _whiteCanQCastle, bool _blackCanQCastle,
            uint16_t _epCaptureFile, int _fiftyMoveCounter, int _moveNumber,
            int _playerToMove);
    ~Board();
    Board staticCopy();

    void doMove(Move m, int color);
    bool doPseudoLegalMove(Move m, int color);
    bool doHashMove(Move m, int color);
    void doNullMove();
    void undoNullMove(uint16_t _epCaptureFile);

template <bool isMoveGen> PieceMoveList getPieceMoveList(int color);
    MoveList getAllLegalMoves(int color);
    MoveList getAllPseudoLegalMoves(int color, PieceMoveList &pml);
    MoveList getPseudoLegalQuiets(int color, PieceMoveList &pml);
    MoveList getPseudoLegalCaptures(int color, PieceMoveList &pml, bool includePromotions);
    MoveList getPseudoLegalPromotions(int color);
    MoveList getPseudoLegalChecks(int color);
    MoveList getPseudoLegalCheckEscapes(int color, PieceMoveList &pml);

    // Get a bitboard of all xray-ers attacking a square if a blocker has been moved or removed
    uint64_t getXRayPieceMap(int color, int sq, int blockerColor,
    uint64_t blockerStart, uint64_t blockerEnd);
    uint64_t getAttackMap(int color, int sq);
    uint64_t getRookXRays(int sq, uint64_t occ, uint64_t blockers);
    uint64_t getbishopXRays(int sq, uint64_t occ, uint64_t blockers);
    uint64_t getPinnedMap(int color);
	
    bool isCheckMove(int color, Move m);
    bool isInCheck(int color);
    bool isWInMate();
    bool isBInMate();
    bool isStalemate(int sideToMove);
    bool isDraw();
    bool isInsufficientMaterial();

    // Evaluation
    int evaluate();
    int getPseudoMobility(int color, PieceMoveList &pml, int egFactor);
    int getKingSafety(PieceMoveList &pmlWhite, PieceMoveList &pmlBlack,
        uint64_t wKingSqs, uint64_t bKingSqs, int egFactor);
    int checkEndgameCases();
    int scoreSimpleKnownWin(int winningColor);
    int getEGFactor();
    int getMaterial(int color);
    int getMaterialEG(int color);
	int getPieceOnSquare(int color, int sq);

    uint64_t getNonPawnMaterial(int color);
    int getManhattanDistance(int sq1, int sq2);
    uint64_t getLeastValuableAttacker(uint64_t attackers, int color, int &piece);
    int getSEE(int color, int sq);
    int getSEEForMove(int color, Move m);
    int valueOfPiece(int piece);
    int getMVVLVAscore(int color, Move m);
    int getExchangescore(int color, Move m);

private:

    uint64_t allPieces[2];
    uint64_t pieces[2][6];
    uint64_t zobristKey;
    uint16_t epCaptureFile;
    uint16_t playerToMove;
    uint16_t moveNumber;
    uint8_t castlingRights;
    uint8_t fiftyMoveCounter;

    Board(Board *b);

    void addPawnMovesToList(MoveList &quiets, int color);
    void addPawnCapturesToList(MoveList &captures, int color, uint64_t otherPieces, bool includePromotions);
    template <bool isCapture> void addMovesToList(MoveList &moves, int stSq,
    uint64_t allEndSqs, uint64_t otherPieces = 0);
    template <bool isCapture> void addPromotionsToList(MoveList &moves,
    int stSq, int endSq);
    void addCastlesToList(MoveList &moves, int color);

    uint64_t getWPawnSingleMoves(uint64_t pawns);
    uint64_t getBPawnSingleMoves(uint64_t pawns);
    uint64_t getWPawnDoubleMoves(uint64_t pawns);
    uint64_t getBPawnDoubleMoves(uint64_t pawns);
    uint64_t getWPawnLeftCaptures(uint64_t pawns);
    uint64_t getBPawnLeftCaptures(uint64_t pawns);
    uint64_t getWPawnRightCaptures(uint64_t pawns);
    uint64_t getBPawnRightCaptures(uint64_t pawns);
    uint64_t getWPawnCaptures(uint64_t pawns);
    uint64_t getBPawnCaptures(uint64_t pawns);
    uint64_t getKnightSquares(int single);
    uint64_t getbishopSquares(int single, uint64_t occ);
    uint64_t getRookSquares(int single, uint64_t occ);
    uint64_t getQueenSquares(int single, uint64_t occ);
    uint64_t getKingSquares(int single);
    uint64_t getOccupancy();
    int epVictimSquare(int victimColor, uint16_t file);
    uint64_t rankAttacks(uint64_t occ, int single);
    uint64_t fileAttacks(uint64_t occ, int single);
    uint64_t diagAttacks(uint64_t occ, int single);
    uint64_t antiDiagAttacks(uint64_t occ, int single);
			
};
uint64_t perft(Board &b, int color, int depth, uint64_t &captures);

const int DEFAULT_HASH_SIZE = 24;
const int MIN_HASH_SIZE = 1;
const int MAX_HASH_SIZE = 8172;
const int DEFAULT_MULTI_PV = 1;
const int MIN_MULTI_PV = 1;
const int MAX_MULTI_PV = 256;
const int DEFAULT_THREADS = 1;
const int MIN_THREADS = 1;
const int MAX_THREADS = 32;

#endif


