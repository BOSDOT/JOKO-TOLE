/*******************************************************************************
Joko Tole UCI Chess Engine Copyright(C) 2016 Bosdot (Javanese)
influenzed by chessprogramming.wikispaces.com and many others 
chess engine with open source code.
--------------------------------------------------------------------------------
 ******************************************************************************/
 
#include <algorithm>
#include <iostream>
#include <random>
#include "board.h"
#include "tables.h"
#include "evaluation.h"
#define MSEED 161803398

const int CenterDistance[64] = { 
  3, 3, 3, 3, 3, 3, 3, 3,
  3, 2, 2, 2, 2, 2, 2, 3,
  3, 2, 1, 1, 1, 1, 2, 3,
  3, 2, 1, 0, 0, 1, 2, 3,
  3, 2, 1, 0, 0, 1, 2, 3,
  3, 2, 1, 1, 1, 1, 2, 3,
  3, 2, 2, 2, 2, 2, 2, 3,
  3, 3, 3, 3, 3, 3, 3, 3
};

// NORTH SOUTH WEST EAST
// SOUTH EAST, SOUTHWEST NORTH EAST, NORTH WEST
// DummbFill7

uint64_t SOUTHOccl(uint64_t gen, uint64_t pro) {
   uint64_t flood = 0;
   while (gen) {
      flood |= gen;
      gen = (gen >> 8) & pro;
   }
   return flood;
}  

uint64_t SOUTHAttacks(uint64_t rooks, uint64_t kosong) {
   uint64_t flood = rooks;
   flood |= rooks = (rooks >> 8) & kosong;
   flood |= rooks = (rooks >> 8) & kosong;
   flood |= rooks = (rooks >> 8) & kosong;
   flood |= rooks = (rooks >> 8) & kosong;
   flood |= rooks = (rooks >> 8) & kosong;
   flood |=         (rooks >> 8) & kosong;
   return            flood >> 8;
}

uint64_t NORTHAttacks(uint64_t rooks, uint64_t kosong) {
   uint64_t flood = rooks;
   flood |= rooks = (rooks << 8) & kosong;
   flood |= rooks = (rooks << 8) & kosong;
   flood |= rooks = (rooks << 8) & kosong;
   flood |= rooks = (rooks << 8) & kosong;
   flood |= rooks = (rooks << 8) & kosong;
   flood |=         (rooks << 8) & kosong;
   return            flood << 8;
}

// Horizontal H-file to A-file & vice versa.
uint64_t EASTAttacks(uint64_t rooks, uint64_t kosong) {
   const uint64_t notA = 0xfefefefefefefefe;
   uint64_t flood = rooks;
   kosong &= notA;
   flood |= rooks = (rooks << 1) & kosong;
   flood |= rooks = (rooks << 1) & kosong;
   flood |= rooks = (rooks << 1) & kosong;
   flood |= rooks = (rooks << 1) & kosong;
   flood |= rooks = (rooks << 1) & kosong;
   flood |=         (rooks << 1) & kosong;
   return           (flood << 1) & notA ;
}
 
uint64_t WESTAttacks(uint64_t rooks, uint64_t kosong) {
   const uint64_t notH = 0x7f7f7f7f7f7f7f7f;
   uint64_t flood = rooks;
   kosong &= notH;
   flood |= rooks = (rooks >> 1) & kosong;
   flood |= rooks = (rooks >> 1) & kosong;
   flood |= rooks = (rooks >> 1) & kosong;
   flood |= rooks = (rooks >> 1) & kosong;
   flood |= rooks = (rooks >> 1) & kosong;
   flood |=         (rooks >> 1) & kosong;
   return           (flood >> 1) & notH ;
}

uint64_t SOUTHEASTaAttacks(uint64_t bishops, uint64_t kosong) {
   const uint64_t notA = 0xfefefefefefefefe;
   uint64_t flood = bishops;
   kosong &= notA;
   flood |= bishops = (bishops >> 7) & kosong;
   flood |= bishops = (bishops >> 7) & kosong;
   flood |= bishops = (bishops >> 7) & kosong;
   flood |= bishops = (bishops >> 7) & kosong;
   flood |= bishops = (bishops >> 7) & kosong;
   flood |=           (bishops >> 7) & kosong;
   return               (flood >> 7) & notA ;
}

uint64_t SOUTHWESTAttacks(uint64_t bishops, uint64_t kosong) {
   const uint64_t notH = 0x7f7f7f7f7f7f7f7f;
   uint64_t flood = bishops;
   kosong &= notH;
   flood |= bishops = (bishops >> 9) & kosong;
   flood |= bishops = (bishops >> 9) & kosong;
   flood |= bishops = (bishops >> 9) & kosong;
   flood |= bishops = (bishops >> 9) & kosong;
   flood |= bishops = (bishops >> 9) & kosong;
   flood |=           (bishops >> 9) & kosong;
   return               (flood >> 9) & notH ;
}

uint64_t NORTHEASTAttacks(uint64_t bishops, uint64_t kosong) {
   const uint64_t notA = 0xfefefefefefefefe;
   uint64_t flood = bishops;
   kosong &= notA;
   flood |= bishops = (bishops << 9) & kosong;
   flood |= bishops = (bishops << 9) & kosong;
   flood |= bishops = (bishops << 9) & kosong;
   flood |= bishops = (bishops << 9) & kosong;
   flood |= bishops = (bishops << 9) & kosong;
   flood |=           (bishops << 9) & kosong;
   return               (flood << 9) & notA ;
}

uint64_t NORTHWESTAttacks(uint64_t bishops, uint64_t kosong) {
   const uint64_t notA = 0x7f7f7f7f7f7f7f7f;
   uint64_t flood = bishops;
   kosong &= notA;
   flood |= bishops = (bishops << 7) & kosong;
   flood |= bishops = (bishops << 7) & kosong;
   flood |= bishops = (bishops << 7) & kosong;
   flood |= bishops = (bishops << 7) & kosong;
   flood |= bishops = (bishops << 7) & kosong;
   flood |=           (bishops << 7) & kosong;
   return               (flood << 7) & notA ;
}

static uint64_t mseed = 0, mstate = 0;
const int NORTH_SOUTH_FILL = 8;
const int EAST_WEST_FILL = 1;
const int SOUTHEAST_FILL = 7;
const int SOUTHWEST_FILL = 9;
const int NORTHEAST_FILL = 9;
const int NORTHWEST_FILL = 7;
const uint64_t notA = 0xfefefefefefefefe; // ~0x0101010101010101
const uint64_t notH = 0x7f7f7f7f7f7f7f7f; // ~0x8080808080808080

uint64_t SOUTHOne (uint64_t b) {return  b >> 8;}
uint64_t NORTHOne (uint64_t b) {return  b << 8;}
uint64_t rotateLeft (uint64_t x, int s) {return _rotl64(x, s);}
uint64_t rotateRight(uint64_t x, int s) {return _rotr64(x, s);}

int shift[8] = {9, 1,-7,-8,-9,-1, 7, 8};

uint64_t avoidWrap[8] =
{
   0xfefefefefefefe00,
   0xfefefefefefefefe,
   0x00fefefefefefefe,
   0x00ffffffffffffff,
   0x007f7f7f7f7f7f7f,
   0x7f7f7f7f7f7f7f7f,
   0x7f7f7f7f7f7f7f00,
   0xffffffffffffff00,
};

// Unrolled Attacks
const uint64_t slidingAttacks (uint64_t sliders, uint64_t kosong, int dir8) {
   uint64_t flood = sliders;
   int r = shift[dir8]; // {+-1,7,8,9}
   kosong &= avoidWrap[dir8];
   flood |= sliders = rotateLeft(sliders , r) & kosong;
   return   rotateLeft(flood, r)  &   avoidWrap[dir8];
}

uint64_t fillRayRight(uint64_t rayPieces, uint64_t kosong, int shift);
uint64_t fillRayLeft(uint64_t rayPieces, uint64_t kosong, int shift);

// magic ray attack
struct MagicInfo {
    uint64_t* table; // pointer to attack_table for each particular square
    uint64_t  mask; // to mask relevant squares of both lines (no outer squares)
    uint64_t  magic; // magic 64-bit factor
	uint64_t  attack; // attacker
    int shift; // shift right
};

// Masks
static uint64_t ROOK_MASK[64];
static uint64_t bishop_MASK[64];
static uint64_t *attackTable;
static MagicInfo magicbishops[64];
static MagicInfo magicRooks[64];
static uint64_t inBetweenSqs[64][64];

uint64_t indexToMask64(int index, int nBits, uint64_t mask);
uint64_t ratt(int sq, uint64_t block);
uint64_t batt(int sq, uint64_t block);
int magicMap(uint64_t masked, uint64_t magic, int nBits);
uint64_t findMagic(int sq, int m, bool isbishop);

MagicInfo mbishopAttacks[64] [512]; // 32 K
MagicInfo mRookAttacks  [64][4096]; // 256 K
MagicInfo mbishopTbl[64];
MagicInfo mRookTbl[64];


void initMagicTables(uint64_t seed) {
    mstate = 74036198046ULL;
    mseed = seed;

    for (int i = 0; i < 64; i++) {
        uint64_t relevantBits = ((~FILES[0] & ~FILES[7]) | FILES[i&7])
                              & ((~RANKS[0] & ~RANKS[7]) | RANKS[i>>3]);
        ROOK_MASK[i] = ratt(i, 0) & relevantBits;
        bishop_MASK[i] = batt(i, 0) & relevantBits;
    }

    attackTable = new uint64_t[107648];
    int runningPtrLoc = 0;

    for (int i = 0; i < 64; i++) {
        uint64_t *tableStart = attackTable;
        magicbishops[i].table = tableStart + runningPtrLoc;
        magicbishops[i].mask = bishop_MASK[i];
        magicbishops[i].magic = findMagic(i, NUM_bishop_BITS[i], true);
        magicbishops[i].shift = 64 - NUM_bishop_BITS[i];

        runningPtrLoc += 1 << NUM_bishop_BITS[i];
    }
// Initialize rook magic values
    for (int i = 0; i < 64; i++) {
        uint64_t *tableStart = attackTable;
        magicRooks[i].table = tableStart + runningPtrLoc;
        magicRooks[i].mask = ROOK_MASK[i];
        magicRooks[i].magic = findMagic(i, NUM_ROOK_BITS[i], false);
        magicRooks[i].shift = 64 - NUM_ROOK_BITS[i];
        runningPtrLoc += 1 << NUM_ROOK_BITS[i];
    }

    for (int sq = 0; sq < 64; sq++) {
        int nBits = NUM_bishop_BITS[sq];
        uint64_t mask = bishop_MASK[sq];
// For each possible mask result
        for (int i = 0; i < (1 << nBits); i++) {
           
            uint64_t *attTableLoc = magicbishops[sq].table;
            uint64_t occ = indexToMask64(i, nBits, mask);
            uint64_t attSet = batt(sq, occ);
            int magicIndex = magicMap(occ, magicbishops[sq].magic, nBits);
            attTableLoc[magicIndex] = attSet;
        }
    }
// Then rooks
    for (int sq = 0; sq < 64; sq++) {
        int nBits = NUM_ROOK_BITS[sq];
        uint64_t mask = ROOK_MASK[sq];
        for (int i = 0; i < (1 << nBits); i++) {
            uint64_t *attTableLoc = magicRooks[sq].table;
            uint64_t occ = indexToMask64(i, nBits, mask);
            uint64_t attSet = ratt(sq, occ);
            int magicIndex = magicMap(occ, magicRooks[sq].magic, nBits);
            attTableLoc[magicIndex] = attSet;
        }
    }
}

uint64_t magicPRNG() {

    uint64_t y = ((mstate << 57) | (mseed << 57)) >> 1;
    mstate ^= mseed >> 17;
    mstate ^= mstate << 3;

    uint64_t temp = mseed;
    mseed = mstate;
    mstate = temp;

    return (y | (mseed ^ mstate)) >> 1;
}

// Zobrist hashing table
static uint64_t zobristTable[794];
static uint64_t startPosZobristKey = 0;

void initZobristTable() {
    std::mt19937_64 PRNG (61280152908);
    for (int i = 0; i < 794; i++)
        zobristTable[i] = PRNG();

    Board b;
    int *mailbox = b.getMailbox();
    b.initZobristKey(mailbox);
    startPosZobristKey = b.getZobristKey();
    delete[] mailbox;
}

// Initializes the 64x64 table, indexed by from and to square, of all

void initInBetweenTable() {
    for (int sq1 = 0; sq1 < 64; sq1++) {
        for (int sq2 = 0; sq2 < 64; sq2++) {
            uint64_t imaginaryRook = ratt(sq1, INDEX_TO_BIT[sq2]);
            if (imaginaryRook & INDEX_TO_BIT[sq2]) {
                uint64_t imaginaryRook2 = ratt(sq2, INDEX_TO_BIT[sq1]);
                inBetweenSqs[sq1][sq2] = imaginaryRook & imaginaryRook2;
            }
            else {
// Check diagonal lines
                uint64_t imaginarybishop = batt(sq1, INDEX_TO_BIT[sq2]);
                if (imaginarybishop & INDEX_TO_BIT[sq2]) {
                    uint64_t imaginarybishop2 = batt(sq2, INDEX_TO_BIT[sq1]);
                    inBetweenSqs[sq1][sq2] = imaginarybishop & imaginarybishop2;
                }
// If the squares are not on a line, the bitboard is kosong
                else
                    inBetweenSqs[sq1][sq2] = 0;
            }
        }
    }
}


uint64_t perft(Board &b, int color, int depth, uint64_t &captures) {
    if (depth == 0)
        return 1;

    uint64_t nodes = 0;

    PieceMoveList pml = b.getPieceMoveList<PML_LEGAL_MOVES>(color);
    MoveList pl = b.getAllPseudoLegalMoves(color, pml);
    for (unsigned int i = 0; i < pl.size(); i++) {
        Board copy = b.staticCopy();
        if (!copy.doPseudoLegalMove(pl.get(i), color))
            continue;

        if (isCapture(pl.get(i)))
            captures++;

        nodes += perft(copy, color^1, depth-1, captures);
    }

    return nodes;
}

//Constructors----------------------------------
Board::Board() {
    allPieces[WHITE] = 0x000000000000FFFF;
    allPieces[BLACK] = 0xFFFF000000000000;
    pieces[WHITE][PAWNS] = 0x000000000000FF00; // white pawns
    pieces[WHITE][KNIGHTS] = 0x0000000000000042; // white knights
    pieces[WHITE][bishopS] = 0x0000000000000024; // white bishops
    pieces[WHITE][ROOKS] = 0x0000000000000081; // white rooks
    pieces[WHITE][QUEENS] = 0x0000000000000008; // white queens
    pieces[WHITE][KINGS] = 0x0000000000000010; // white kings
    pieces[BLACK][PAWNS] = 0x00FF000000000000; // black pawns
    pieces[BLACK][KNIGHTS] = 0x4200000000000000; // black knights
    pieces[BLACK][bishopS] = 0x2400000000000000; // black bishops
    pieces[BLACK][ROOKS] = 0x8100000000000000; // black rooks
    pieces[BLACK][QUEENS] = 0x0800000000000000; // black queens
    pieces[BLACK][KINGS] = 0x1000000000000000; // black kings

    zobristKey = startPosZobristKey;
    epCaptureFile = NO_EP_POSSIBLE;
    playerToMove = WHITE;
    moveNumber = 1;
    castlingRights = WHITECASTLE | BLACKCASTLE;
    fiftyMoveCounter = 0;
}

// Create a board object from a mailbox 
Board::Board(int *mailboxBoard, bool _whiteCanKCastle, bool _blackCanKCastle,
        bool _whiteCanQCastle, bool _blackCanQCastle,  uint16_t _epCaptureFile,
        int _fiftyMoveCounter, int _moveNumber, int _playerToMove) {
    for (int i = 0; i < 12; i++) {
        pieces[i/6][i%6] = 0;
    }
    for (int i = 0; i < 64; i++) {
        if (0 <= mailboxBoard[i] && mailboxBoard[i] <= 11) {
            pieces[mailboxBoard[i]/6][mailboxBoard[i]%6] |= INDEX_TO_BIT[i];
        }
    }
    allPieces[WHITE] = 0;
    for (int i = 0; i < 6; i++)
        allPieces[WHITE] |= pieces[WHITE][i];
    allPieces[BLACK] = 0;
    for (int i = 0; i < 6; i++)
        allPieces[BLACK] |= pieces[BLACK][i];

    epCaptureFile = _epCaptureFile;
    playerToMove = _playerToMove;
    moveNumber = _moveNumber;
    castlingRights = 0;
    if (_whiteCanKCastle)
        castlingRights |= WHITEKSIDE;
    if (_whiteCanQCastle)
        castlingRights |= WHITEQSIDE;
    if (_blackCanKCastle)
        castlingRights |= BLACKKSIDE;
    if (_blackCanQCastle)
        castlingRights |= BLACKQSIDE;
    fiftyMoveCounter = _fiftyMoveCounter;
    initZobristKey(mailboxBoard);
}

Board::~Board() {}

// Creates a copy of a board
Board::Board(Board *b) {
    allPieces[WHITE] = b->allPieces[WHITE];
    allPieces[BLACK] = b->allPieces[BLACK];
    for (int i = 0; i < 6; i++) {
        pieces[0][i] = b->pieces[0][i];
    }
    for (int i = 0; i < 6; i++) {
        pieces[1][i] = b->pieces[1][i];
    }

    zobristKey = b->zobristKey;
    epCaptureFile = b->epCaptureFile;
    playerToMove = b->playerToMove;
    moveNumber = b->moveNumber;
    castlingRights = b->castlingRights;
    fiftyMoveCounter = b->fiftyMoveCounter;
}

Board Board::staticCopy() {
    return Board(this);
}

// Do move
void Board::doMove(Move m, int color) {
    int startSq = getStartSq(m);
    int endSq = getEndSq(m);
    int pieceID = getPieceOnSquare(color, startSq);

// Update flag based elements of Zobrist key
    zobristKey ^= zobristTable[769 + castlingRights];
    zobristKey ^= zobristTable[785 + epCaptureFile];


    if (isPromotion(m)) {
        int promotionType = getPromotion(m);
        if (isCapture(m)) {
            int captureType = getPieceOnSquare(color^1, endSq);
            pieces[color][PAWNS] &= ~INDEX_TO_BIT[startSq];
            pieces[color][promotionType] |= INDEX_TO_BIT[endSq];
            pieces[color^1][captureType] &= ~INDEX_TO_BIT[endSq];

            allPieces[color] &= ~INDEX_TO_BIT[startSq];
            allPieces[color] |= INDEX_TO_BIT[endSq];
            allPieces[color^1] &= ~INDEX_TO_BIT[endSq];

            zobristKey ^= zobristTable[384*color + startSq];
            zobristKey ^= zobristTable[384*color + 64*promotionType + endSq];
            zobristKey ^= zobristTable[384*(color^1) + 64*captureType + endSq];
        }
        else {
            pieces[color][PAWNS] &= ~INDEX_TO_BIT[startSq];
            pieces[color][promotionType] |= INDEX_TO_BIT[endSq];

            allPieces[color] &= ~INDEX_TO_BIT[startSq];
            allPieces[color] |= INDEX_TO_BIT[endSq];

            zobristKey ^= zobristTable[384*color + startSq];
            zobristKey ^= zobristTable[384*color + 64*promotionType + endSq];
        }
        epCaptureFile = NO_EP_POSSIBLE;
        fiftyMoveCounter = 0;
    } // end promotion
    else if (isCapture(m)) {
        if (isEP(m)) {
            pieces[color][PAWNS] &= ~INDEX_TO_BIT[startSq];
            pieces[color][PAWNS] |= INDEX_TO_BIT[endSq];
            uint64_t epCaptureSq = INDEX_TO_BIT[epVictimSquare(color^1, epCaptureFile)];
            pieces[color^1][PAWNS] &= ~epCaptureSq;

            allPieces[color] &= ~INDEX_TO_BIT[startSq];
            allPieces[color] |= INDEX_TO_BIT[endSq];
            allPieces[color^1] &= ~epCaptureSq;

            int capSq = epVictimSquare(color^1, epCaptureFile);
            zobristKey ^= zobristTable[384*color + startSq];
            zobristKey ^= zobristTable[384*color + endSq];
            zobristKey ^= zobristTable[384*(color^1) + capSq];
        }
        else {
            int captureType = getPieceOnSquare(color^1, endSq);
            pieces[color][pieceID] &= ~INDEX_TO_BIT[startSq];
            pieces[color][pieceID] |= INDEX_TO_BIT[endSq];
            pieces[color^1][captureType] &= ~INDEX_TO_BIT[endSq];

            allPieces[color] &= ~INDEX_TO_BIT[startSq];
            allPieces[color] |= INDEX_TO_BIT[endSq];
            allPieces[color^1] &= ~INDEX_TO_BIT[endSq];

            zobristKey ^= zobristTable[384*color + 64*pieceID + startSq];
            zobristKey ^= zobristTable[384*color + 64*pieceID + endSq];
            zobristKey ^= zobristTable[384*(color^1) + 64*captureType + endSq];
        }
        epCaptureFile = NO_EP_POSSIBLE;
        fiftyMoveCounter = 0;
    } // end capture
    else { // Quiet moves
        if (isCastle(m)) {
            if (endSq == 6) { // white kside
                pieces[WHITE][KINGS] &= ~INDEX_TO_BIT[4];
                pieces[WHITE][KINGS] |= INDEX_TO_BIT[6];
                pieces[WHITE][ROOKS] &= ~INDEX_TO_BIT[7];
                pieces[WHITE][ROOKS] |= INDEX_TO_BIT[5];

                allPieces[WHITE] &= ~INDEX_TO_BIT[4];
                allPieces[WHITE] |= INDEX_TO_BIT[6];
                allPieces[WHITE] &= ~INDEX_TO_BIT[7];
                allPieces[WHITE] |= INDEX_TO_BIT[5];

                zobristKey ^= zobristTable[64*KINGS+4];
                zobristKey ^= zobristTable[64*KINGS+6];
                zobristKey ^= zobristTable[64*ROOKS+7];
                zobristKey ^= zobristTable[64*ROOKS+5];
            }
            else if (endSq == 2) { // white qside
                pieces[WHITE][KINGS] &= ~INDEX_TO_BIT[4];
                pieces[WHITE][KINGS] |= INDEX_TO_BIT[2];
                pieces[WHITE][ROOKS] &= ~INDEX_TO_BIT[0];
                pieces[WHITE][ROOKS] |= INDEX_TO_BIT[3];

                allPieces[WHITE] &= ~INDEX_TO_BIT[4];
                allPieces[WHITE] |= INDEX_TO_BIT[2];
                allPieces[WHITE] &= ~INDEX_TO_BIT[0];
                allPieces[WHITE] |= INDEX_TO_BIT[3];

                zobristKey ^= zobristTable[64*KINGS+4];
                zobristKey ^= zobristTable[64*KINGS+2];
                zobristKey ^= zobristTable[64*ROOKS+0];
                zobristKey ^= zobristTable[64*ROOKS+3];
            }
            else if (endSq == 62) { // black kside
                pieces[BLACK][KINGS] &= ~INDEX_TO_BIT[60];
                pieces[BLACK][KINGS] |= INDEX_TO_BIT[62];
                pieces[BLACK][ROOKS] &= ~INDEX_TO_BIT[63];
                pieces[BLACK][ROOKS] |= INDEX_TO_BIT[61];

                allPieces[BLACK] &= ~INDEX_TO_BIT[60];
                allPieces[BLACK] |= INDEX_TO_BIT[62];
                allPieces[BLACK] &= ~INDEX_TO_BIT[63];
                allPieces[BLACK] |= INDEX_TO_BIT[61];

                zobristKey ^= zobristTable[384+64*KINGS+60];
                zobristKey ^= zobristTable[384+64*KINGS+62];
                zobristKey ^= zobristTable[384+64*ROOKS+63];
                zobristKey ^= zobristTable[384+64*ROOKS+61];
            }
            else { // black qside
                pieces[BLACK][KINGS] &= ~INDEX_TO_BIT[60];
                pieces[BLACK][KINGS] |= INDEX_TO_BIT[58];
                pieces[BLACK][ROOKS] &= ~INDEX_TO_BIT[56];
                pieces[BLACK][ROOKS] |= INDEX_TO_BIT[59];

                allPieces[BLACK] &= ~INDEX_TO_BIT[60];
                allPieces[BLACK] |= INDEX_TO_BIT[58];
                allPieces[BLACK] &= ~INDEX_TO_BIT[56];
                allPieces[BLACK] |= INDEX_TO_BIT[59];

                zobristKey ^= zobristTable[384+64*KINGS+60];
                zobristKey ^= zobristTable[384+64*KINGS+58];
                zobristKey ^= zobristTable[384+64*ROOKS+56];
                zobristKey ^= zobristTable[384+64*ROOKS+59];
            }
            epCaptureFile = NO_EP_POSSIBLE;
            fiftyMoveCounter++;
        } // end castling
        else { // other quiet moves
            pieces[color][pieceID] &= ~INDEX_TO_BIT[startSq];
            pieces[color][pieceID] |= INDEX_TO_BIT[endSq];

            allPieces[color] &= ~INDEX_TO_BIT[startSq];
            allPieces[color] |= INDEX_TO_BIT[endSq];

            zobristKey ^= zobristTable[384*color + 64*pieceID + startSq];
            zobristKey ^= zobristTable[384*color + 64*pieceID + endSq];

            // check for en passant
            if (pieceID == PAWNS) {
                if (getFlags(m) == MOVE_DOUBLE_PAWN)
                    epCaptureFile = startSq & 7;
                else
                    epCaptureFile = NO_EP_POSSIBLE;

                fiftyMoveCounter = 0;
            }
            else {
                epCaptureFile = NO_EP_POSSIBLE;
                fiftyMoveCounter++;
            }
        } // end other quiet moves
    } // end normal move

    // change castling flags
    if (pieceID == KINGS) {
        if (color == WHITE)
            castlingRights &= ~WHITECASTLE;
        else
            castlingRights &= ~BLACKCASTLE;
    }
    // Castling rights change because of the rook only when the rook moves or
    // is captured
    else if (isCapture(m) || pieceID == ROOKS) {
        // No sense in removing the rights if they're already gone
        if (castlingRights & WHITECASTLE) {
            if ((pieces[WHITE][ROOKS] & 0x80) == 0)
                castlingRights &= ~WHITEKSIDE;
            if ((pieces[WHITE][ROOKS] & 1) == 0)
                castlingRights &= ~WHITEQSIDE;
        }
        if (castlingRights & BLACKCASTLE) {
            uint64_t blackR = pieces[BLACK][ROOKS] >> 56;
            if ((blackR & 0x80) == 0)
                castlingRights &= ~BLACKKSIDE;
            if ((blackR & 0x1) == 0)
                castlingRights &= ~BLACKQSIDE;
        }
    } // end castling flags

    zobristKey ^= zobristTable[769 + castlingRights];
    zobristKey ^= zobristTable[785 + epCaptureFile];

    if (color == BLACK)
        moveNumber++;
    playerToMove = color^1;
    zobristKey ^= zobristTable[768];
}

bool Board::doPseudoLegalMove(Move m, int color) {
    doMove(m, color);
    // Pseudo-legal moves require a check for legality
    return !(isInCheck(color));
}

// Do a hash move, which requires a few more checks in case of a Type-1 error.
bool Board::doHashMove(Move m, int color) {
    int pieceID = getPieceOnSquare(color, getStartSq(m));
    // Check that the start square is not kosong
    if (pieceID == -1)
        return false;

// Check that the end square has correct occupancy
    uint64_t otherPieces = allPieces[color^1];
    uint64_t endSingle = INDEX_TO_BIT[getEndSq(m)];
    bool captureRoutes = (isCapture(m) && (otherPieces & endSingle))
                      || (isCapture(m) && pieceID == PAWNS && (~otherPieces & endSingle));
    uint64_t kosong = ~getOccupancy();
    if (!(captureRoutes || (!isCapture(m) && (kosong & endSingle))))
        return false;
    // Check that the king is not captured
    if (isCapture(m) && ((endSingle & pieces[WHITE][KINGS]) || (endSingle & pieces[BLACK][KINGS])))
        return false;

    return doPseudoLegalMove(m, color);
}

// Handle null moves for null move pruning by switching the player to move.
void Board::doNullMove() {
    playerToMove = playerToMove ^ 1;
    zobristKey ^= zobristTable[768];
    zobristKey ^= zobristTable[785 + epCaptureFile];
    epCaptureFile = NO_EP_POSSIBLE;
    zobristKey ^= zobristTable[785 + epCaptureFile];
}

void Board::undoNullMove(uint16_t _epCaptureFile) {
    playerToMove = playerToMove ^ 1;
    zobristKey ^= zobristTable[768];
    zobristKey ^= zobristTable[785 + epCaptureFile];
    epCaptureFile = _epCaptureFile;
    zobristKey ^= zobristTable[785 + epCaptureFile];
}

// Move Generation
template <bool isMoveGen>
PieceMoveList Board::getPieceMoveList(int color) {
    PieceMoveList pml;

    uint64_t knights = pieces[color][KNIGHTS];
    while (knights) {
        int stSq = bitScanForward(knights);
        knights &= knights-1;
        uint64_t nSq = getKnightSquares(stSq);

        pml.add(PieceMoveInfo(KNIGHTS, stSq, nSq));
    }

    uint64_t occ = isMoveGen ? getOccupancy()
                             : allPieces[color^1] | pieces[color][PAWNS] | pieces[color][KINGS];
    uint64_t bishops = pieces[color][bishopS];
    while (bishops) {
        int stSq = bitScanForward(bishops);
        bishops &= bishops-1;
        uint64_t bSq = getbishopSquares(stSq, occ);

        pml.add(PieceMoveInfo(bishopS, stSq, bSq));
    }

    uint64_t rooks = pieces[color][ROOKS];
    while (rooks) {
        int stSq = bitScanForward(rooks);
        rooks &= rooks-1;
        uint64_t rSq = getRookSquares(stSq, occ);

        pml.add(PieceMoveInfo(ROOKS, stSq, rSq));
    }

    uint64_t queens = pieces[color][QUEENS];
    while (queens) {
        int stSq = bitScanForward(queens);
        queens &= queens-1;
        uint64_t qSq = getQueenSquares(stSq, occ);

        pml.add(PieceMoveInfo(QUEENS, stSq, qSq));
    }

    return pml;
}

// Get all legal moves and captures
MoveList Board::getAllLegalMoves(int color) {
    PieceMoveList pml = getPieceMoveList<PML_LEGAL_MOVES>(color);
    MoveList moves = getAllPseudoLegalMoves(color, pml);

    for (unsigned int i = 0; i < moves.size(); i++) {
        Board b = staticCopy();
        b.doMove(moves.get(i), color);

        if (b.isInCheck(color)) {
            moves.remove(i);
            i--;
        }
    }

    return moves;
}

//  Pseudo-legal Moves
MoveList Board::getAllPseudoLegalMoves(int color, PieceMoveList &pml) {
    MoveList quiets = getPseudoLegalQuiets(color, pml);
    MoveList captures = getPseudoLegalCaptures(color, pml, true);

    // Put captures before quiet moves
    for (unsigned int i = 0; i < quiets.size(); i++) {
        captures.add(quiets.get(i));
    }
    return captures;
}

MoveList Board::getPseudoLegalQuiets(int color, PieceMoveList &pml) {
    MoveList quiets;

    addCastlesToList(quiets, color);

    for (unsigned int i = 0; i < pml.size(); i++) {
        PieceMoveInfo pmi = pml.get(i);
        addMovesToList<MOVEGEN_QUIETS>(quiets, pmi.startSq, pmi.legal);
    }

    addPawnMovesToList(quiets, color);

    int stsqK = bitScanForward(pieces[color][KINGS]);
    uint64_t kingSqs = getKingSquares(stsqK);
    addMovesToList<MOVEGEN_QUIETS>(quiets, stsqK, kingSqs);

    return quiets;
}

MoveList Board::getPseudoLegalCaptures(int color, PieceMoveList &pml, bool includePromotions) {
    MoveList captures;

    uint64_t otherPieces = allPieces[color^1];

    int kingStSq = bitScanForward(pieces[color][KINGS]);
    uint64_t kingSqs = getKingSquares(kingStSq);
    addMovesToList<MOVEGEN_CAPTURES>(captures, kingStSq, kingSqs, otherPieces);

    addPawnCapturesToList(captures, color, otherPieces, includePromotions);

    for (unsigned int i = 0; i < pml.size(); i++) {
        PieceMoveInfo pmi = pml.get(i);
        addMovesToList<MOVEGEN_CAPTURES>(captures, pmi.startSq, pmi.legal, otherPieces);
    }

    return captures;
}

// Generates all queen promotions for quiescence search
MoveList Board::getPseudoLegalPromotions(int color) {
    MoveList moves;
    uint64_t otherPieces = allPieces[color^1];

    uint64_t pawns = pieces[color][PAWNS];
    uint64_t finalRank = (color == WHITE) ? RANKS[7] : RANKS[0];

    int leftDiff = (color == WHITE) ? -7 : 9;
    int rightDiff = (color == WHITE) ? -9 : 7;

    uint64_t legal = (color == WHITE) ? getWPawnLeftCaptures(pawns)
                                      : getBPawnLeftCaptures(pawns);
    legal &= otherPieces;
    uint64_t promotions = legal & finalRank;

    while (promotions) {
        int endSq = bitScanForward(promotions);
        promotions &= promotions-1;

        Move mq = encodeMove(endSq+leftDiff, endSq);
        mq = setCapture(mq, true);
        mq = setFlags(mq, MOVE_PROMO_Q);
        moves.add(mq);
    }

    legal = (color == WHITE) ? getWPawnRightCaptures(pawns)
                             : getBPawnRightCaptures(pawns);
    legal &= otherPieces;
    promotions = legal & finalRank;

    while (promotions) {
        int endSq = bitScanForward(promotions);
        promotions &= promotions-1;

        Move mq = encodeMove(endSq+rightDiff, endSq);
        mq = setCapture(mq, true);
        mq = setFlags(mq, MOVE_PROMO_Q);
        moves.add(mq);
    }

    int sqDiff = (color == WHITE) ? -8 : 8;

    legal = (color == WHITE) ? getWPawnSingleMoves(pawns)
                             : getBPawnSingleMoves(pawns);
    promotions = legal & finalRank;

    while (promotions) {
        int endSq = bitScanForward(promotions);
        promotions &= promotions - 1;
        int stSq = endSq + sqDiff;

        Move mq = encodeMove(stSq, endSq);
        mq = setFlags(mq, MOVE_PROMO_Q);
        moves.add(mq);
    }

    return moves;
}

MoveList Board::getPseudoLegalChecks(int color) {
    MoveList checks;
    int kingSq = bitScanForward(pieces[color^1][KINGS]);
    // Square parity for knight and bishop moves
    uint64_t kingParity = (pieces[color^1][KINGS] & LIGHT) ? LIGHT : DARK;
    uint64_t potentialXRay = pieces[color][bishopS] | pieces[color][ROOKS] | pieces[color][QUEENS];
    uint64_t pawns = pieces[color][PAWNS];

    // TODO this is way too slow
    /*
    uint64_t tempPawns = pawns;
    while (tempPawns) {
        int stsq = bitScanForward(tempPawns);
        tempPawns &= tempPawns - 1;
        uint64_t xrays = getXRayPieceMap(color, kingSq, color, INDEX_TO_BIT[stsq]);
        // If moving the pawn caused a new xray piece to attack the king
        if (!(xrays & invAttackMap)) {
            // Every legal move of this pawn is a legal check
            uint64_t legal = (color == WHITE) ? getWPawnSingleMoves(INDEX_TO_BIT[stsq]) | getWPawnDoubleMoves(INDEX_TO_BIT[stsq])
                                              : getBPawnSingleMoves(INDEX_TO_BIT[stsq]) | getBPawnDoubleMoves(INDEX_TO_BIT[stsq]);
            while (legal) {
                int endsq = bitScanForward(legal);
                legal &= legal - 1;
                checks.add(encodeMove(stsq, endsq, PAWNS, false));
            }

            // Remove this pawn from future consideration
            pawns ^= INDEX_TO_BIT[stsq];
        }
    }
    */

    uint64_t pAttackMap = (color == WHITE) 
            ? getBPawnCaptures(INDEX_TO_BIT[kingSq])
            : getWPawnCaptures(INDEX_TO_BIT[kingSq]);
    uint64_t finalRank = (color == WHITE) ? RANKS[7] : RANKS[0];
    int sqDiff = (color == WHITE) ? -8 : 8;

    uint64_t pLegal = (color == WHITE) ? getWPawnSingleMoves(pawns)
                                       : getBPawnSingleMoves(pawns);
    // Remove promotions
    uint64_t promotions = pLegal & finalRank;
    pLegal ^= promotions;

    pLegal &= pAttackMap;
    while (pLegal) {
        int endsq = bitScanForward(pLegal);
        pLegal &= pLegal - 1;
        checks.add(encodeMove(endsq+sqDiff, endsq));
    }

    pLegal = (color == WHITE) ? getWPawnDoubleMoves(pawns)
                              : getBPawnDoubleMoves(pawns);
    pLegal &= pAttackMap;
    while (pLegal) {
        int endsq = bitScanForward(pLegal);
        pLegal &= pLegal - 1;
        Move m = encodeMove(endsq+2*sqDiff, endsq);
        m = setFlags(m, MOVE_DOUBLE_PAWN);
        checks.add(m);
    }

    uint64_t knights = pieces[color][KNIGHTS] & kingParity;
    uint64_t nAttackMap = getKnightSquares(kingSq);
    while (knights) {
        int stsq = bitScanForward(knights);
        knights &= knights-1;
        uint64_t nSq = getKnightSquares(stsq);
// Get any bishops, rooks, queens attacking king after knight has moved
        uint64_t xrays = getXRayPieceMap(color, kingSq, color, INDEX_TO_BIT[stsq], 0);
        if (!(xrays & potentialXRay))
            nSq &= nAttackMap;

        addMovesToList<MOVEGEN_QUIETS>(checks, stsq, nSq);
    }

    uint64_t occ = getOccupancy();
    uint64_t bishops = pieces[color][bishopS] & kingParity;
    uint64_t bAttackMap = getbishopSquares(kingSq, occ);
    while (bishops) {
        int stsq = bitScanForward(bishops);
        bishops &= bishops-1;
        uint64_t bSq = getbishopSquares(stsq, occ);
        uint64_t xrays = getXRayPieceMap(color, kingSq, color, INDEX_TO_BIT[stsq], 0);
        if (!(xrays & potentialXRay))
            bSq &= bAttackMap;

        addMovesToList<MOVEGEN_QUIETS>(checks, stsq, bSq);
    }

    uint64_t rooks = pieces[color][ROOKS];
    uint64_t rAttackMap = getRookSquares(kingSq, occ);
    while (rooks) {
        int stsq = bitScanForward(rooks);
        rooks &= rooks-1;
        uint64_t rSq = getRookSquares(stsq, occ);
        uint64_t xrays = getXRayPieceMap(color, kingSq, color, INDEX_TO_BIT[stsq], 0);
        if (!(xrays & potentialXRay))
            rSq &= rAttackMap;

        addMovesToList<MOVEGEN_QUIETS>(checks, stsq, rSq);
    }

    uint64_t queens = pieces[color][QUEENS];
    uint64_t qAttackMap = getQueenSquares(kingSq, occ);
    while (queens) {
        int stsq = bitScanForward(queens);
        queens &= queens-1;
        uint64_t qSq = getQueenSquares(stsq, occ) & qAttackMap;

        addMovesToList<MOVEGEN_QUIETS>(checks, stsq, qSq);
    }

    return checks;
}

// Generate moves that (sort of but not really) get out of check
MoveList Board::getPseudoLegalCheckEscapes(int color, PieceMoveList &pml) {
    MoveList captures, blocks;

    int kingSq = bitScanForward(pieces[color][KINGS]);
    uint64_t otherPieces = allPieces[color^1];
    uint64_t attackMap = getAttackMap(color^1, kingSq);
// Consider only captures of pieces giving check
    otherPieces &= attackMap;

// If double check, we can only move the king
    if (count(otherPieces) >= 2) {
        uint64_t kingSqs = getKingSquares(kingSq);

        addMovesToList<MOVEGEN_CAPTURES>(captures, kingSq, kingSqs, allPieces[color^1]);
        addMovesToList<MOVEGEN_QUIETS>(captures, kingSq, kingSqs);
        return captures;
    }

    addPawnMovesToList(blocks, color);
    addPawnCapturesToList(captures, color, otherPieces, true);

    uint64_t occ = getOccupancy();
    uint64_t xraySqs = 0;
    int attackerSq = bitScanForward(otherPieces);
    int attackerType = getPieceOnSquare(color^1, attackerSq);
    if (attackerType == bishopS)
        xraySqs = getbishopSquares(attackerSq, occ);
    else if (attackerType == ROOKS)
        xraySqs = getRookSquares(attackerSq, occ);
    else if (attackerType == QUEENS)
        xraySqs = getQueenSquares(attackerSq, occ);

    for (unsigned int i = 0; i < pml.size(); i++) {
        PieceMoveInfo pmi = pml.get(i);
        addMovesToList<MOVEGEN_QUIETS>(blocks, pmi.startSq, pmi.legal & xraySqs);
        addMovesToList<MOVEGEN_CAPTURES>(captures, pmi.startSq, pmi.legal, otherPieces);
    }

    int stsqK = bitScanForward(pieces[color][KINGS]);
    uint64_t kingSqs = getKingSquares(stsqK);
    addMovesToList<MOVEGEN_QUIETS>(blocks, stsqK, kingSqs);
    addMovesToList<MOVEGEN_CAPTURES>(captures, stsqK, kingSqs, allPieces[color^1]);

    // Put captures before blocking moves
    for (unsigned int i = 0; i < blocks.size(); i++) {
        captures.add(blocks.get(i));
    }

    return captures;
}

// Move Generation Helpers
void Board::addPawnMovesToList(MoveList &quiets, int color) {
    uint64_t pawns = pieces[color][PAWNS];
    uint64_t finalRank = (color == WHITE) ? RANKS[7] : RANKS[0];
    int sqDiff = (color == WHITE) ? -8 : 8;

    uint64_t pLegal = (color == WHITE) ? getWPawnSingleMoves(pawns)
                                       : getBPawnSingleMoves(pawns);
    uint64_t promotions = pLegal & finalRank;
    pLegal ^= promotions;

    while (promotions) {
        int endSq = bitScanForward(promotions);
        promotions &= promotions - 1;
        int stSq = endSq + sqDiff;

        addPromotionsToList<MOVEGEN_QUIETS>(quiets, stSq, endSq);
    }
    while (pLegal) {
        int endsq = bitScanForward(pLegal);
        pLegal &= pLegal - 1;
        quiets.add(encodeMove(endsq+sqDiff, endsq));
    }

    pLegal = (color == WHITE) ? getWPawnDoubleMoves(pawns)
                              : getBPawnDoubleMoves(pawns);
    while (pLegal) {
        int endsq = bitScanForward(pLegal);
        pLegal &= pLegal - 1;
        Move m = encodeMove(endsq+2*sqDiff, endsq);
        m = setFlags(m, MOVE_DOUBLE_PAWN);
        quiets.add(m);
    }
}

// For pawn captures, we can use a similar approach, but we must consider
void Board::addPawnCapturesToList(MoveList &captures, int color, uint64_t otherPieces, bool includePromotions) {
    uint64_t pawns = pieces[color][PAWNS];
    uint64_t finalRank = (color == WHITE) ? RANKS[7] : RANKS[0];
    int leftDiff = (color == WHITE) ? -7 : 9;
    int rightDiff = (color == WHITE) ? -9 : 7;

    uint64_t legal = (color == WHITE) ? getWPawnLeftCaptures(pawns)
                                      : getBPawnLeftCaptures(pawns);
    legal &= otherPieces;
    uint64_t promotions = legal & finalRank;
    legal ^= promotions;

    if (includePromotions) {
        while (promotions) {
            int endSq = bitScanForward(promotions);
            promotions &= promotions-1;

            addPromotionsToList<MOVEGEN_CAPTURES>(captures, endSq+leftDiff, endSq);
        }
    }
    while (legal) {
        int endsq = bitScanForward(legal);
        legal &= legal-1;
        Move m = encodeMove(endsq+leftDiff, endsq);
        m = setCapture(m, true);
        captures.add(m);
    }

    legal = (color == WHITE) ? getWPawnRightCaptures(pawns)
                             : getBPawnRightCaptures(pawns);
    legal &= otherPieces;
    promotions = legal & finalRank;
    legal ^= promotions;

    if (includePromotions) {
        while (promotions) {
            int endSq = bitScanForward(promotions);
            promotions &= promotions-1;

            addPromotionsToList<MOVEGEN_CAPTURES>(captures, endSq+rightDiff, endSq);
        }
    }
    while (legal) {
        int endsq = bitScanForward(legal);
        legal &= legal-1;
        Move m = encodeMove(endsq+rightDiff, endsq);
        m = setCapture(m, true);
        captures.add(m);
    }

    // If there are en passants possible...
    if (epCaptureFile != NO_EP_POSSIBLE) {
        int victimSq = epVictimSquare(color^1, epCaptureFile);
        int rankDiff = (color == WHITE) ? 8 : -8;
        if ((INDEX_TO_BIT[victimSq] << 1) & NOTA & pieces[color][PAWNS]) {
            Move m = encodeMove(victimSq+1, victimSq+rankDiff);
            m = setFlags(m, MOVE_EP);
            captures.add(m);
        }
        if ((INDEX_TO_BIT[victimSq] >> 1) & NOTH & pieces[color][PAWNS]) {
            Move m = encodeMove(victimSq-1, victimSq+rankDiff);
            m = setFlags(m, MOVE_EP);
            captures.add(m);
        }
    }
}

// Helper function that processes a bitboard of legal moves and adds all
template <bool isCapture>
void Board::addMovesToList(MoveList &moves, int stSq, uint64_t allEndSqs,
    uint64_t otherPieces) {

    uint64_t intersect = (isCapture) ? otherPieces : ~getOccupancy();
    uint64_t legal = allEndSqs & intersect;
    while (legal) {
        int endSq = bitScanForward(legal);
        legal &= legal-1;
        Move m = encodeMove(stSq, endSq);
        if (isCapture)
            m = setCapture(m, true);
        moves.add(m);
    }
}

template <bool isCapture>
void Board::addPromotionsToList(MoveList &moves, int stSq, int endSq) {
    Move mk = encodeMove(stSq, endSq);
    mk = setFlags(mk, MOVE_PROMO_N);
    Move mb = encodeMove(stSq, endSq);
    mb = setFlags(mb, MOVE_PROMO_B);
    Move mr = encodeMove(stSq, endSq);
    mr = setFlags(mr, MOVE_PROMO_R);
    Move mq = encodeMove(stSq, endSq);
    mq = setFlags(mq, MOVE_PROMO_Q);
    if (isCapture) {
        mk = setCapture(mk, true);
        mb = setCapture(mb, true);
        mr = setCapture(mr, true);
        mq = setCapture(mq, true);
    }
	
    moves.add(mq);
    moves.add(mk);
    moves.add(mr);
    moves.add(mb);
}

void Board::addCastlesToList(MoveList &moves, int color) {
    if (color == WHITE) {
        if ((castlingRights & WHITEKSIDE)
         && (getOccupancy() & WHITE_KSIDE_PASSTHROUGH_SQS) == 0
         && !isInCheck(WHITE)) {
            if (getAttackMap(BLACK, 5) == 0) {
                Move m = encodeMove(4, 6);
                m = setCastle(m, true);
                moves.add(m);
            }
        }
        if ((castlingRights & WHITEQSIDE)
         && (getOccupancy() & WHITE_QSIDE_PASSTHROUGH_SQS) == 0
         && !isInCheck(WHITE)) {
            if (getAttackMap(BLACK, 3) == 0) {
                Move m = encodeMove(4, 2);
                m = setCastle(m, true);
                moves.add(m);
            }
        }
    }
    else {
        if ((castlingRights & BLACKKSIDE)
         && (getOccupancy() & BLACK_KSIDE_PASSTHROUGH_SQS) == 0
         && !isInCheck(BLACK)) {
            if (getAttackMap(WHITE, 61) == 0) {
                Move m = encodeMove(60, 62);
                m = setCastle(m, true);
                moves.add(m);
            }
        }
        if ((castlingRights & BLACKQSIDE)
         && (getOccupancy() & BLACK_QSIDE_PASSTHROUGH_SQS) == 0
         && !isInCheck(BLACK)) {
            if (getAttackMap(WHITE, 59) == 0) {
                Move m = encodeMove(60, 58);
                m = setCastle(m, true);
                moves.add(m);
            }
        }
    }
}

// Useful bitboard info generators
uint64_t Board::getXRayPieceMap(int color, int sq, int blockerColor,
    uint64_t blockerStart, uint64_t blockerEnd) {
    uint64_t occ = getOccupancy();
    occ ^= blockerStart;
    occ ^= blockerEnd;

    uint64_t bishops = pieces[color][bishopS];
    uint64_t rooks = pieces[color][ROOKS];
    uint64_t queens = pieces[color][QUEENS];

    uint64_t xRayMap = (getbishopSquares(sq, occ) & (bishops | queens))
                     | (getRookSquares(sq, occ) & (rooks | queens));

    return (xRayMap & ~blockerStart);
}

// Given a color and a square, returns all pieces of the color that attack the
uint64_t Board::getAttackMap(int color, int sq) {
    uint64_t occ = getOccupancy();
    uint64_t pawnCap = (color == WHITE)
                     ? getBPawnCaptures(INDEX_TO_BIT[sq])
                     : getWPawnCaptures(INDEX_TO_BIT[sq]);
    return (pawnCap & pieces[color][PAWNS])
         | (getKnightSquares(sq) & pieces[color][KNIGHTS])
         | (getbishopSquares(sq, occ) & (pieces[color][bishopS] | pieces[color][QUEENS]))
         | (getRookSquares(sq, occ) & (pieces[color][ROOKS] | pieces[color][QUEENS]))
         | (getKingSquares(sq) & pieces[color][KINGS]);
}

// Given the on a given square, used to get either the piece moving or the
int Board::getPieceOnSquare(int color, int sq) {
    uint64_t endSingle = INDEX_TO_BIT[sq];
    for (int pieceID = 0; pieceID <= 5; pieceID++) {
        if (pieces[color][pieceID] & endSingle)
            return pieceID;
    }
    return -1;
}

// Returns true if a move puts the opponent in check
bool Board::isCheckMove(int color, Move m) {
    int kingSq = bitScanForward(pieces[color^1][KINGS]);
    uint64_t attackMap = 0;
    uint64_t occ = getOccupancy();
    switch (getPieceOnSquare(color, getStartSq(m))) {
        case PAWNS:
            attackMap = (color == WHITE) 
                ? getBPawnCaptures(INDEX_TO_BIT[kingSq])
                : getWPawnCaptures(INDEX_TO_BIT[kingSq]);
            break;
        case KNIGHTS:
            attackMap = getKnightSquares(kingSq);
            break;
        case bishopS:
            attackMap = getbishopSquares(kingSq, occ);
            break;
        case ROOKS:
            attackMap = getRookSquares(kingSq, occ);
            break;
        case QUEENS:
            attackMap = getQueenSquares(kingSq, occ);
            break;
        case KINGS:
            break;
    }
    if (INDEX_TO_BIT[getEndSq(m)] & attackMap)
        return true;

    uint64_t xrayPieces = pieces[color][bishopS] | pieces[color][ROOKS] | pieces[color][QUEENS];
    uint64_t xrays = getXRayPieceMap(color, kingSq, color, INDEX_TO_BIT[getStartSq(m)], INDEX_TO_BIT[getEndSq(m)]);
    if (xrays & xrayPieces)
        return true;
    return false;
}

uint64_t Board::getRookXRays(int sq, uint64_t occ, uint64_t blockers) {
    uint64_t attacks = getRookSquares(sq, occ);
    blockers &= attacks;
    return attacks ^ getRookSquares(sq, occ ^ blockers);
}

uint64_t Board::getbishopXRays(int sq, uint64_t occ, uint64_t blockers) {
    uint64_t attacks = getbishopSquares(sq, occ);
    blockers &= attacks;
    return attacks ^ getbishopSquares(sq, occ ^ blockers);
}

uint64_t Board::getPinnedMap(int color) {
    uint64_t pinned = 0;
    uint64_t blockers = allPieces[color];
    int kingSq = bitScanForward(pieces[color][KINGS]);

    uint64_t pinners = getRookXRays(kingSq, getOccupancy(), blockers)
        & (pieces[color^1][ROOKS] | pieces[color^1][QUEENS]);
    while (pinners) {
        int sq  = bitScanForward(pinners);
        pinned |= inBetweenSqs[sq][kingSq] & blockers;
        pinners &= pinners - 1;
    }

    pinners = getbishopXRays(kingSq, getOccupancy(), blockers)
        & (pieces[color^1][bishopS] | pieces[color^1][QUEENS]);
    while (pinners) {
        int sq  = bitScanForward(pinners);
        pinned |= inBetweenSqs[sq][kingSq] & blockers;
        pinners &= pinners - 1;
    }

    return pinned;
}

//check, checkmate, stalemate
bool Board::isInCheck(int color) {
    int sq = bitScanForward(pieces[color][KINGS]);

    return getAttackMap(color^1, sq);
}

bool Board::isWInMate() {
    if (!isInCheck(WHITE))
        return false;

    MoveList moves = getAllLegalMoves(WHITE);
    return (moves.size() == 0);
}

bool Board::isBInMate() {
    if (!isInCheck(BLACK))
        return false;

    MoveList moves = getAllLegalMoves(BLACK);
    return (moves.size() == 0);
}

bool Board::isStalemate(int sideToMove) {
    if (isInCheck(sideToMove))
        return false;

    MoveList moves = getAllLegalMoves(sideToMove);
    return (moves.size() == 0);
}

bool Board::isDraw() {
    if (fiftyMoveCounter >= 100) return true;

    if (isInsufficientMaterial())
        return true;
    
    return false;
}

bool Board::isInsufficientMaterial() {
    int numPieces = count(allPieces[WHITE] | allPieces[BLACK]) - 2;
    if (numPieces < 2) {
        if (numPieces == 0)
            return true;
        if (pieces[WHITE][KNIGHTS] || pieces[WHITE][bishopS]
         || pieces[BLACK][KNIGHTS] || pieces[BLACK][bishopS])
            return true;
    }
    return false;
}

//Evaluation and Move Ordering
int Board::evaluate() {
    int whiteMaterial = getMaterial(WHITE);
    int blackMaterial = getMaterial(BLACK);
    int egFactor = EG_FACTOR_RES - (whiteMaterial + blackMaterial - START_VALUE / 2) * EG_FACTOR_RES / START_VALUE;
    egFactor = std::max(0, std::min(EG_FACTOR_RES, egFactor));
    if (egFactor == EG_FACTOR_RES) {
        int endgamescore = checkEndgameCases();
        if (endgamescore != -INFTY)
            return endgamescore;
    }

 // Material terms
    int valueMg = 0;
    int valueEg = 0;
    valueMg += whiteMaterial;
    valueMg -= blackMaterial;
    valueEg += getMaterialEG(WHITE);
    valueEg -= getMaterialEG(BLACK);
    
// Tempo bonus
    valueMg += (playerToMove == WHITE) ? TEMPO_VALUE : -TEMPO_VALUE;
    if ((pieces[WHITE][bishopS] & LIGHT) && (pieces[WHITE][bishopS] & DARK)) {
        valueMg += bishop_PAIR_VALUE;
        valueEg += bishop_PAIR_VALUE;
    }
    if ((pieces[BLACK][bishopS] & LIGHT) && (pieces[BLACK][bishopS] & DARK)) {
        valueMg -= bishop_PAIR_VALUE;
        valueEg -= bishop_PAIR_VALUE;
    }
	
// Piece square tables
    for (int pieceID = 0; pieceID < 6; pieceID++) {
        uint64_t bitboard = pieces[0][pieceID];
        bitboard = flipAcrossRanks(bitboard);
        while (bitboard) {
            int i = bitScanForward(bitboard);
            bitboard &= bitboard - 1;
            valueMg += midgamePieceValues[pieceID][i];
            valueEg += endgamePieceValues[pieceID][i];
        }
    }
    for (int pieceID = 0; pieceID < 6; pieceID++)  {
        uint64_t bitboard = pieces[1][pieceID];
        while (bitboard) {
            int i = bitScanForward(bitboard);
            bitboard &= bitboard - 1;
            valueMg -= midgamePieceValues[pieceID][i];
            valueEg -= endgamePieceValues[pieceID][i];
        }
    }

// Pawn value scaling: incur a penalty for having no pawns since the
    const int PAWN_SCALING_MG[9] = {10, 3, 2, 2, 1, 0, 0, 0, 0};
    const int PAWN_SCALING_EG[9] = {71, 20, 8, 4, 2, 0, 0, 0, 0};
    int whitePawnCount = count(pieces[WHITE][PAWNS]);
    int blackPawnCount = count(pieces[BLACK][PAWNS]);

    valueMg -= PAWN_SCALING_MG[whitePawnCount];
    valueEg -= PAWN_SCALING_EG[whitePawnCount];
    valueMg += PAWN_SCALING_MG[blackPawnCount];
    valueEg += PAWN_SCALING_EG[blackPawnCount];

    // With queens on the board pawns are a target in the endgame
    if (pieces[WHITE][QUEENS])
        valueEg += 2 * count(pieces[BLACK][PAWNS]);
    if (pieces[BLACK][QUEENS])
        valueEg -= 2 * count(pieces[WHITE][PAWNS]);
    int materialValue = (valueMg * (EG_FACTOR_RES - egFactor) + valueEg * egFactor) / EG_FACTOR_RES;

    //----------------------------Positional terms------------------------------
    score value = EVAL_ZERO;
       
// King Safety and Mobility
    int mobilityValue = 0;
    PieceMoveList pmlWhite = getPieceMoveList<PML_PSEUDO_MOBILITY>(WHITE);
    PieceMoveList pmlBlack = getPieceMoveList<PML_PSEUDO_MOBILITY>(BLACK);
    mobilityValue += getPseudoMobility(WHITE, pmlWhite, egFactor);
    mobilityValue -= getPseudoMobility(BLACK, pmlBlack, egFactor);

// Consider squares near king
    uint64_t wksq = getKingSquares(bitScanForward(pieces[WHITE][KINGS]));
    uint64_t bksq = getKingSquares(bitScanForward(pieces[BLACK][KINGS]));
    if (egFactor < EG_FACTOR_RES) {
        mobilityValue += getKingSafety(pmlWhite, pmlBlack, wksq, bksq, egFactor);

// Castling rights
        value += CASTLING_RIGHTS_VALUE[count(castlingRights & WHITECASTLE)];
        value -= CASTLING_RIGHTS_VALUE[count(castlingRights & BLACKCASTLE)];
        
// Pawn shield bonus (files ABC, FGH)
        value += PAWN_SHIELD_VALUE * count((wksq | (wksq << 8)) & pieces[WHITE][PAWNS] & 0xe7e7e7e7e7e7e7e7);
        value -= PAWN_SHIELD_VALUE * count((bksq | (bksq >> 8)) & pieces[BLACK][PAWNS] & 0xe7e7e7e7e7e7e7e7);
        value += P_PAWN_SHIELD_BONUS * count(wksq & pieces[WHITE][PAWNS] & 0xe7e7e7e7e7e7e7e7);
        value -= P_PAWN_SHIELD_BONUS * count(bksq & pieces[BLACK][PAWNS] & 0xe7e7e7e7e7e7e7e7);
        
// Open files next to king
        uint64_t notwp = ~pieces[WHITE][PAWNS];
        uint64_t notbp = ~pieces[BLACK][PAWNS];
        uint64_t tempwk = pieces[WHITE][KINGS];
        uint64_t tempbk = pieces[BLACK][KINGS];
        tempwk |= ((tempwk >> 1) & NOTH) | ((tempwk << 1) & NOTA);
        tempbk |= ((tempbk >> 1) & NOTH) | ((tempbk << 1) & NOTA);
        uint64_t tempwk2 = tempwk;
        uint64_t tempbk2 = tempbk;
        
// Flood fill: checking for white pawns
        for(int i = 0; i < 7; i++) {
            tempwk |= (tempwk << 8) & notwp;
            tempbk |= (tempbk >> 8) & notwp;
        }
// If the "king" made it across the board without running into a white pawn,
        int wkNoWhiteOpen = count(tempwk & RANKS[7]);
        int bkNoWhiteOpen = count(tempbk & RANKS[0]);
        
// Flood fill: checking for black pawns
        for(int i = 0; i < 7; i++) {
            tempwk2 |= (tempwk2 << 8) & notbp;
            tempbk2 |= (tempbk2 >> 8) & notbp;
        }
// If the "king" made it across the board without running into a black pawn,
        int wkNoBlackOpen = count(tempwk2 & RANKS[7]);
        int bkNoBlackOpen = count(tempbk2 & RANKS[0]);
        
        value -= SEMIOPEN_OWN_PENALTY*wkNoWhiteOpen;
        value -= SEMIOPEN_OPP_PENALTY*wkNoBlackOpen;
        value += SEMIOPEN_OPP_PENALTY*bkNoWhiteOpen;
        value += SEMIOPEN_OWN_PENALTY*bkNoBlackOpen;
        value -= OPEN_PENALTY*count(tempwk & tempwk2 & RANKS[7]);
        value += OPEN_PENALTY*count(tempbk & tempbk2 & RANKS[0]);
    }

// Current pawn attacks
    uint64_t wPawnAtt = getWPawnCaptures(pieces[WHITE][PAWNS]);
    uint64_t bPawnAtt = getBPawnCaptures(pieces[BLACK][PAWNS]);
    uint64_t wPawnFrontSpan = pieces[WHITE][PAWNS] << 8;
    uint64_t bPawnFrontSpan = pieces[BLACK][PAWNS] >> 8;
    for (int i = 0; i < 5; i++) {
        wPawnFrontSpan |= wPawnFrontSpan << 8;
        bPawnFrontSpan |= bPawnFrontSpan >> 8;
    }
    uint64_t wPawnStopAtt = ((wPawnFrontSpan >> 1) & NOTH) | ((wPawnFrontSpan << 1) & NOTA);
    uint64_t bPawnStopAtt = ((bPawnFrontSpan >> 1) & NOTH) | ((bPawnFrontSpan << 1) & NOTA);

// Minor Pieces
    if (pieces[WHITE][bishopS] & LIGHT) {
        value -= bishop_PAWN_COLOR_PENALTY * count(pieces[WHITE][PAWNS] & LIGHT);
    }
    if (pieces[WHITE][bishopS] & DARK) {
        value -= bishop_PAWN_COLOR_PENALTY * count(pieces[WHITE][PAWNS] & DARK);
    }
    if (pieces[BLACK][bishopS] & LIGHT) {
        value += bishop_PAWN_COLOR_PENALTY * count(pieces[BLACK][PAWNS] & LIGHT);
    }
    if (pieces[BLACK][bishopS] & DARK) {
        value += bishop_PAWN_COLOR_PENALTY * count(pieces[BLACK][PAWNS] & DARK);
    }

// Knights do better when the opponent has many pawns
    value += KNIGHT_PAWN_BONUS * count(pieces[WHITE][KNIGHTS]) * count(pieces[BLACK][PAWNS]);
    value -= KNIGHT_PAWN_BONUS * count(pieces[BLACK][KNIGHTS]) * count(pieces[WHITE][PAWNS]);

// Knight outposts: knights that cannot be attacked by opposing pawns
    const uint64_t RANKS_456 = RANKS[3] | RANKS[4] | RANKS[5];
    const uint64_t RANKS_543 = RANKS[4] | RANKS[3] | RANKS[2];
    const uint64_t FILES_CDEF = FILES[2] | FILES[3] | FILES[4] | FILES[5];
    if (uint64_t wOutpost = pieces[WHITE][KNIGHTS] & ~bPawnStopAtt & RANKS_456 & FILES_CDEF) {
        value += KNIGHT_OUTPOST_BONUS * count(wOutpost);
        value += OUTPOST_PAWN_DEF_BONUS * count(wOutpost & wPawnAtt);
    }
    if (uint64_t bOutpost = pieces[BLACK][KNIGHTS] & ~wPawnStopAtt & RANKS_543 & FILES_CDEF) {
        value -= KNIGHT_OUTPOST_BONUS * count(bOutpost);
        value -= OUTPOST_PAWN_DEF_BONUS * count(bOutpost & bPawnAtt);
    }

// Special case
    if ((pieces[WHITE][KNIGHTS] & 0x40000)
     && (pieces[WHITE][PAWNS] & 0x400)
     && (pieces[WHITE][PAWNS] & 0x8000000)
    && !(pieces[WHITE][PAWNS] & 0x10000000))
        value -= KNIGHT_C3_CLOSED_PENALTY;
    if ((pieces[BLACK][KNIGHTS] & 0x40000000000)
     && (pieces[BLACK][PAWNS] & 0x4000000000000)
     && (pieces[BLACK][PAWNS] & 0x800000000)
    && !(pieces[BLACK][PAWNS] & 0x1000000000))
        value += KNIGHT_C3_CLOSED_PENALTY;

//--------Rooks-------
    uint64_t wRooksOpen = pieces[WHITE][ROOKS];
    while (wRooksOpen) {
        int rookSq = bitScanForward(wRooksOpen);
        wRooksOpen &= wRooksOpen - 1;
        int file = rookSq & 7;

        if (!(FILES[file] & (pieces[WHITE][PAWNS] | pieces[BLACK][PAWNS])))
            value += ROOK_OPEN_FILE_BONUS;
    }

    uint64_t bRooksOpen = pieces[BLACK][ROOKS];
    while (bRooksOpen) {
        int rookSq = bitScanForward(bRooksOpen);
        bRooksOpen &= bRooksOpen - 1;
        int file = rookSq & 7;
        if (!(FILES[file] & (pieces[WHITE][PAWNS] | pieces[BLACK][PAWNS])))
            value -= ROOK_OPEN_FILE_BONUS;
    }
   
   // Pawn structure
    uint64_t wPassedBlocker = pieces[BLACK][PAWNS] >> 8;
    uint64_t bPassedBlocker = pieces[WHITE][PAWNS] << 8;
    wPassedBlocker |= ((wPassedBlocker >> 1) & NOTH) | ((wPassedBlocker << 1) & NOTA);
    bPassedBlocker |= ((bPassedBlocker >> 1) & NOTH) | ((bPassedBlocker << 1) & NOTA);
    wPassedBlocker |= (pieces[WHITE][PAWNS] >> 8);
    bPassedBlocker |= (pieces[BLACK][PAWNS] << 8);
    for(int i = 0; i < 4; i++) {
        wPassedBlocker |= (wPassedBlocker >> 8);
        bPassedBlocker |= (bPassedBlocker << 8);
    }

    uint64_t wPassedPawns = pieces[WHITE][PAWNS] & ~wPassedBlocker;
    uint64_t bPassedPawns = pieces[BLACK][PAWNS] & ~bPassedBlocker;
    uint64_t whiteBlockaders = pieces[BLACK][KNIGHTS] | pieces[BLACK][bishopS] | pieces[BLACK][KINGS];
    uint64_t blackBlockaders = pieces[WHITE][KNIGHTS] | pieces[WHITE][bishopS] | pieces[WHITE][KINGS];
    uint64_t wPasserTemp = wPassedPawns;
	
    while (wPasserTemp) {
        int passerSq = bitScanForward(wPasserTemp);
        wPasserTemp &= wPasserTemp - 1;
        int file = passerSq & 7;
        int rank = passerSq >> 3;
        value += PASSER_BONUS[rank];
        value += PASSER_FILE_BONUS[file];
        if ((INDEX_TO_BIT[passerSq] << 8) & whiteBlockaders)
            value -= BLOCKADED_PASSER_PENALTY;
    }
    uint64_t bPasserTemp = bPassedPawns;
	
    while (bPasserTemp) {
        int passerSq = bitScanForward(bPasserTemp);
        bPasserTemp &= bPasserTemp - 1;
        int file = passerSq & 7;
        int rank = 7 - (passerSq >> 3);
        value -= PASSER_BONUS[rank];
        value -= PASSER_FILE_BONUS[file];
        if ((INDEX_TO_BIT[passerSq] >> 8) & blackBlockaders)
            value += BLOCKADED_PASSER_PENALTY;
    }
    
    int wPawnCtByFile[8];
    int bPawnCtByFile[8];
    for (int i = 0; i < 8; i++) {
        wPawnCtByFile[i] = count(pieces[WHITE][PAWNS] & FILES[i]);
        bPawnCtByFile[i] = count(pieces[BLACK][PAWNS] & FILES[i]);
    }
    
    int numWPawns = count(pieces[WHITE][PAWNS]);
    int numBPawns = count(pieces[BLACK][PAWNS]);
    for (int i = 0; i < 8; i++) {
        value -= DOUBLED_PENALTY[wPawnCtByFile[i]] * DOUBLED_PENALTY_SCALE[numWPawns];
        value += DOUBLED_PENALTY[bPawnCtByFile[i]] * DOUBLED_PENALTY_SCALE[numBPawns];
    }    
    // Isolated pawns
    uint64_t wp = 0, bp = 0;
    for (int i = 0; i < 8; i++) {
        wp |= (bool) (wPawnCtByFile[i]);
        bp |= (bool) (bPawnCtByFile[i]);
        wp <<= 1;
        bp <<= 1;
    }
    wp >>= 1;
    bp >>= 1;
    wp &= ~((wp >> 1) | (wp << 1));
    bp &= ~((bp >> 1) | (bp << 1));
	
    int whiteIsolated = count(wp);
    int blackIsolated = count(bp);
    value -= ISOLATED_PENALTY * whiteIsolated;
    value += ISOLATED_PENALTY * blackIsolated;
    value -= CENTRAL_ISOLATED_PENALTY * count(wp & 0x7E);
    value += CENTRAL_ISOLATED_PENALTY * count(bp & 0x7E);

// Isolated, doubled pawns
    for (int i = 0; i < 8; i++) {
        if ((wPawnCtByFile[i] > 1) && (wp & INDEX_TO_BIT[7-i])) {
            value -= ISOLATED_DOUBLED_PENALTY * ((wPawnCtByFile[i] - 1) * wPawnCtByFile[i]) / 2;
        }
        if ((bPawnCtByFile[i] > 1) && (bp & INDEX_TO_BIT[7-i])) {
            value += ISOLATED_DOUBLED_PENALTY * ((bPawnCtByFile[i] - 1) * bPawnCtByFile[i]) / 2;
        }
    }

// Backward pawns
    uint64_t wBadStopSqs = ~wPawnStopAtt & bPawnAtt;
    uint64_t bBadStopSqs = ~bPawnStopAtt & wPawnAtt;
    for (int i = 0; i < 6; i++) {
        wBadStopSqs |= wBadStopSqs >> 8;
        bBadStopSqs |= bBadStopSqs << 8;
    }
	
    uint64_t wBackwards = wBadStopSqs & pieces[WHITE][PAWNS];
    uint64_t bBackwards = bBadStopSqs & pieces[BLACK][PAWNS];
    value -= BACKWARD_PENALTY * count(wBackwards);
    value += BACKWARD_PENALTY * count(bBackwards);

// King-pawn
    int kingPawnTropism = 0;
    if (egFactor > 0) {
        uint64_t wPawnTropism = pieces[WHITE][PAWNS] | pieces[BLACK][PAWNS];
        uint64_t bPawnTropism = pieces[WHITE][PAWNS] | pieces[BLACK][PAWNS];
        int pawnCount = count(pieces[WHITE][PAWNS] | pieces[BLACK][PAWNS]);
        int wKingSq = bitScanForward(pieces[WHITE][KINGS]);
        int bKingSq = bitScanForward(pieces[BLACK][KINGS]);

        int wTropismTotal = 0;
        while (wPawnTropism) {
            int pawnSq = bitScanForward(wPawnTropism);
            int rank = pawnSq >> 3;
            wPawnTropism &= wPawnTropism - 1;
            if (INDEX_TO_BIT[pawnSq] & wPassedPawns)
                wTropismTotal += 4 * getManhattanDistance(pawnSq, wKingSq) * rank;
            else if (INDEX_TO_BIT[pawnSq] & bPassedPawns)
                wTropismTotal += 7 * getManhattanDistance(pawnSq, wKingSq) * (7-rank);
            else if (INDEX_TO_BIT[pawnSq] & (wBackwards | bBackwards))
                wTropismTotal += 2 * getManhattanDistance(pawnSq, wKingSq);
            else
                wTropismTotal += getManhattanDistance(pawnSq, wKingSq);
        }
		
        int bTropismTotal = 0;
        while (bPawnTropism) {
            int pawnSq = bitScanForward(bPawnTropism);
            int rank = pawnSq >> 3;
            bPawnTropism &= bPawnTropism - 1;
            if (INDEX_TO_BIT[pawnSq] & wPassedPawns)
                bTropismTotal += 7 * getManhattanDistance(pawnSq, bKingSq) * rank;
            else if (INDEX_TO_BIT[pawnSq] & bPassedPawns)
                bTropismTotal += 4 * getManhattanDistance(pawnSq, bKingSq) * (7-rank);
            else if (INDEX_TO_BIT[pawnSq] & (wBackwards | bBackwards))
                bTropismTotal += 2 * getManhattanDistance(pawnSq, bKingSq);
            else
                bTropismTotal += getManhattanDistance(pawnSq, bKingSq);
        }
        if (pawnCount)
            kingPawnTropism = (bTropismTotal - wTropismTotal) / pawnCount;
            kingPawnTropism = kingPawnTropism * egFactor / EG_FACTOR_RES;
    }
    int totalEval = (16 * materialValue + 17 * decEval(value, egFactor)
        + 18 * mobilityValue + 16 * kingPawnTropism) / 16;

// Scale factors
    if (egFactor > 3 * EG_FACTOR_RES / 4) {
        if (count(pieces[WHITE][bishopS]) == 1
         && count(pieces[BLACK][bishopS]) == 1
         && (((pieces[WHITE][bishopS] & LIGHT) && (pieces[BLACK][bishopS] & DARK))
          || ((pieces[WHITE][bishopS] & DARK) && (pieces[BLACK][bishopS] & LIGHT)))) {
            if ((getNonPawnMaterial(WHITE) == pieces[WHITE][bishopS])
             && (getNonPawnMaterial(BLACK) == pieces[BLACK][bishopS]))
                totalEval = totalEval / 2;
            else
                totalEval = (384 - 128 * egFactor / EG_FACTOR_RES) * totalEval / 256;
        }
    }
    return totalEval;
}

int Board::getPseudoMobility(int color, PieceMoveList &pml, int egFactor) {
    const uint64_t CENTER_SQS = 0x0000001818000000;
    const int EXTENDED_CENTER_VAL = 2;
    const int CENTER_BONUS = 2;
    int result = 0;
    int centerControl = 0;
    uint64_t openSqs = ~allPieces[color];
    uint64_t EXTENDED_CENTER_SQS = (color == WHITE) ? 0x0000183C3C000000
                                                    : 0x0000003C3C180000;
    uint64_t pawns = pieces[color][PAWNS];
    uint64_t pawnAttackMap = (color == WHITE) ? getWPawnCaptures(pawns)
                                              : getBPawnCaptures(pawns);
    centerControl += EXTENDED_CENTER_VAL * count(pawnAttackMap & EXTENDED_CENTER_SQS);
    centerControl += CENTER_BONUS * count(pawnAttackMap & CENTER_SQS);
    uint64_t oppPawns = pieces[color^1][PAWNS];
    uint64_t oppPawnAttackMap = (color == WHITE) ? getBPawnCaptures(oppPawns)
                                                 : getWPawnCaptures(oppPawns);
    uint64_t knightMobilitySqs = allPieces[color^1] | (openSqs & ~oppPawnAttackMap);
    int undevelopedCount = count(
        (pieces[color^1][KNIGHTS] | pieces[color^1][bishopS]) & RANKS[7-7*color]);
    for (unsigned int i = 0; i < pml.size(); i++) {
        PieceMoveInfo pmi = pml.get(i);
        int pieceIndex = pmi.pieceID - 1;
        uint64_t legal = pmi.legal;
        if (pieceIndex == KNIGHTS - 1)
            result += mobilityscore[pieceIndex][count(legal & knightMobilitySqs)];
        else if (pieceIndex == QUEENS - 1)
            result += mobilityscore[pieceIndex][count(legal & openSqs)]
                * (6 - undevelopedCount) / 6;
        else
            result += mobilityscore[pieceIndex][count(legal & openSqs)];
        if (pieceIndex == QUEENS - 1 && undevelopedCount > 0)
            continue;
        centerControl += EXTENDED_CENTER_VAL * count(legal & EXTENDED_CENTER_SQS);
        centerControl += CENTER_BONUS * count(legal & CENTER_SQS);
    }
    result += centerControl * (EG_FACTOR_RES - egFactor) / EG_FACTOR_RES;
    return result;
}

// King safety, based on the number of opponent pieces near the king
int Board::getKingSafety(PieceMoveList &pmlWhite, PieceMoveList &pmlBlack,
    uint64_t wKingSqs, uint64_t bKingSqs, int egFactor) {

    const int KING_THREAT_MULTIPLIER[4] = {13, 16, 18, 26};
    const int KING_THREAT_PIECE_BONUS[16] = {0, 0, 0, 35, 60, 80, 90, 95,
                                            100, 100, 100, 100, 100, 100, 100, 100};
    int result = 0;
    int wKingSafety = 0, bKingSafety = 0;
    int wKingAttackPieces = 0, bKingAttackPieces = 0;
    for (unsigned int i = 0; i < pmlWhite.size(); i++) {
        PieceMoveInfo pmi = pmlWhite.get(i);

        int pieceIndex = pmi.pieceID - 1;
        uint64_t legal = pmi.legal;
        int kingSqCount = count(legal & bKingSqs);
        if (kingSqCount) {
            wKingAttackPieces++;
            wKingSafety += KING_THREAT_MULTIPLIER[pieceIndex] * kingSqCount;
        }
    }

    // Iterate over piece move information to extract all mobility-related scores
    for (unsigned int i = 0; i < pmlBlack.size(); i++) {
        PieceMoveInfo pmi = pmlBlack.get(i);

        int pieceIndex = pmi.pieceID - 1;
        uint64_t legal = pmi.legal;
        int kingSqCount = count(legal & wKingSqs);
        if (kingSqCount) {
            bKingAttackPieces++;
            bKingSafety += KING_THREAT_MULTIPLIER[pieceIndex] * kingSqCount;
        }
    }

    if (wKingAttackPieces >= 2) {
        result += KING_THREAT_PIECE_BONUS[wKingAttackPieces];
        result += wKingSafety;
    }
    if (bKingAttackPieces >= 2) {
        result -= KING_THREAT_PIECE_BONUS[bKingAttackPieces];
        result -= bKingSafety;
    }

    return result * (EG_FACTOR_RES - egFactor) / EG_FACTOR_RES;
}

// Check special endgame cases: where help mate is possible (detecting this
int Board::checkEndgameCases() {
    int numWPieces = count(allPieces[WHITE]) - 1;
    int numBPieces = count(allPieces[BLACK]) - 1;
    int numPieces = numWPieces + numBPieces;

    if (numBPieces == 0 && (pieces[WHITE][ROOKS] || pieces[WHITE][QUEENS])) {
        return scoreSimpleKnownWin(WHITE);
    }
    if (numWPieces == 0 && (pieces[BLACK][ROOKS] || pieces[BLACK][QUEENS])) {
        return scoreSimpleKnownWin(BLACK);
    }

    if (numPieces == 1) {
        if (pieces[WHITE][PAWNS]) {
            int wPawn = bitScanForward(flipAcrossRanks(pieces[WHITE][PAWNS]));
            return 3 * PAWN_VALUE_EG / 2 + endgamePieceValues[PAWNS][wPawn];
        }
        if (pieces[BLACK][PAWNS]) {
            int bPawn = bitScanForward(pieces[BLACK][PAWNS]);
            return -3 * PAWN_VALUE_EG / 2 - endgamePieceValues[PAWNS][bPawn];
        }
    }

    else if (numPieces == 2) {
        if (numWPieces == 1) {
            if ((pieces[WHITE][KNIGHTS] | pieces[WHITE][bishopS])
             && (pieces[BLACK][KNIGHTS] | pieces[BLACK][bishopS]))
                return 0;
            if (pieces[WHITE][ROOKS] && pieces[BLACK][ROOKS])
                return 0;
            if (pieces[WHITE][QUEENS] && pieces[BLACK][QUEENS])
                return 0;
        }

        else {

            if (pieces[WHITE][PAWNS]) {
                int value = KNOWN_WIN / 2;
                for (int pieceID = 0; pieceID < 6; pieceID++) {
                    uint64_t bitboard = pieces[0][pieceID];
                    bitboard = flipAcrossRanks(bitboard);
                    while (bitboard) {
                        int i = bitScanForward(bitboard);
                        bitboard &= bitboard - 1;
                        value += endgamePieceValues[pieceID][i];
                    }
                }
                int bKing = bitScanForward(pieces[BLACK][KINGS]);
                value -= endgamePieceValues[KINGS][bKing];
                return value;
            }
            if (pieces[BLACK][PAWNS]) {
                int value = -KNOWN_WIN / 2;
                for (int pieceID = 0; pieceID < 6; pieceID++)  {
                    uint64_t bitboard = pieces[1][pieceID];
                    while (bitboard) {
                        int i = bitScanForward(bitboard);
                        bitboard &= bitboard - 1;
                        value -= endgamePieceValues[pieceID][i];
                    }
                }
                int wKing = bitScanForward(pieces[WHITE][KINGS]);
                value += endgamePieceValues[KINGS][wKing];
                return value;
            }
            if (count(pieces[WHITE][KNIGHTS]) == 2 || count(pieces[BLACK][KNIGHTS]) == 2)
                return 0;
            if (count(pieces[WHITE][bishopS]) == 2)
                return scoreSimpleKnownWin(WHITE);
            if (count(pieces[BLACK][bishopS]) == 2)
                return scoreSimpleKnownWin(BLACK);
        }
    }
    return -INFTY;
}

// function for scoring
int Board::scoreSimpleKnownWin(int winningColor) {
    int wKing = bitScanForward(pieces[WHITE][KINGS]);
    int bKing = bitScanForward(pieces[BLACK][KINGS]);
    int winscore = (winningColor == WHITE) ? KNOWN_WIN : -KNOWN_WIN;
    return winscore + endgamePieceValues[KINGS][wKing]
                    - endgamePieceValues[KINGS][bKing];
}

// Gets the endgame factor
int Board::getEGFactor() {
    int whiteMaterial = getMaterial(WHITE);
    int blackMaterial = getMaterial(BLACK);
    int egFactor = EG_FACTOR_RES - (whiteMaterial + blackMaterial - START_VALUE / 2) * EG_FACTOR_RES / START_VALUE;
    return std::max(0, std::min(EG_FACTOR_RES, egFactor));
}

int Board::getMaterial(int color) {
    return PAWN_VALUE   * count(pieces[color][PAWNS])
         + KNIGHT_VALUE * count(pieces[color][KNIGHTS])
         + bishop_VALUE * count(pieces[color][bishopS])
         + ROOK_VALUE   * count(pieces[color][ROOKS])
         + QUEEN_VALUE  * count(pieces[color][QUEENS]);
}

int Board::getMaterialEG(int color) {
    return PAWN_VALUE_EG   * count(pieces[color][PAWNS])
         + KNIGHT_VALUE_EG * count(pieces[color][KNIGHTS])
         + bishop_VALUE_EG * count(pieces[color][bishopS])
         + ROOK_VALUE_EG   * count(pieces[color][ROOKS])
         + QUEEN_VALUE_EG  * count(pieces[color][QUEENS]);
}

uint64_t Board::getNonPawnMaterial(int color) {
    return pieces[color][KNIGHTS] | pieces[color][bishopS]
         | pieces[color][ROOKS]   | pieces[color][QUEENS];
}

int Board::getManhattanDistance(int sq1, int sq2) {
    return std::abs((sq1 >> 3) - (sq2 >> 3)) + std::abs((sq1 & 7) - (sq2 & 7));
}

// Given a bitboard of attackers, finds the least valuable attacker of color and
uint64_t Board::getLeastValuableAttacker(uint64_t attackers, int color, int &piece) {
    for (piece = 0; piece < 5; piece++) {
        uint64_t single = attackers & pieces[color][piece];
        if (single)
            return single & -single;
    }

    piece = KINGS;
    return attackers & pieces[color][KINGS];
}

int Board::getSEE(int color, int sq) {
    int gain[32], d = 0, piece = 0;
    uint64_t attackers = getAttackMap(WHITE, sq) | getAttackMap(BLACK, sq);
    uint64_t used = 0;
    uint64_t single = getLeastValuableAttacker(attackers, color, piece);
    if (!single)
        return 0;
    gain[d] = SEE_PIECE_VALS[getPieceOnSquare(color^1, sq)];
    do {
        d++; // next depth
        color ^= 1;
        gain[d]  = SEE_PIECE_VALS[piece] - gain[d-1];
        if (-gain[d-1] < 0 && gain[d] < 0) // pruning for stand pat
            break;
        attackers ^= single; // remove used attacker
        used |= single;
        attackers |= getXRayPieceMap(WHITE, sq, color, used, 0) | getXRayPieceMap(BLACK, sq, color, used, 0);
        single = getLeastValuableAttacker(attackers, color, piece);
    } while (single);

    while (--d)
        gain[d-1]= -((-gain[d-1] > gain[d]) ? -gain[d-1] : gain[d]);
    return gain[0];
}

int Board::getSEEForMove(int color, Move m) {
    int value = 0;
    int startSq = getStartSq(m);
    int endSq = getEndSq(m);
    int pieceID = getPieceOnSquare(color, startSq);
    if (isEP(m))
        return -PAWN_VALUE;

// Do move temp
    pieces[color][pieceID] &= ~INDEX_TO_BIT[startSq];
    pieces[color][pieceID] |= INDEX_TO_BIT[endSq];
    allPieces[color] &= ~INDEX_TO_BIT[startSq];
    allPieces[color] |= INDEX_TO_BIT[endSq];

// capture
    if (isCapture(m)) {
        int capturedPiece = getPieceOnSquare(color^1, endSq);
        pieces[color^1][capturedPiece] &= ~INDEX_TO_BIT[endSq];
        allPieces[color^1] &= ~INDEX_TO_BIT[endSq];

        value = SEE_PIECE_VALS[capturedPiece] - getSEE(color^1, endSq);

        pieces[color^1][capturedPiece] |= INDEX_TO_BIT[endSq];
        allPieces[color^1] |= INDEX_TO_BIT[endSq];
    }
    else {
        value = -getSEE(color^1, endSq);
    }

    pieces[color][pieceID] |= INDEX_TO_BIT[startSq];
    pieces[color][pieceID] &= ~INDEX_TO_BIT[endSq];
    allPieces[color] |= INDEX_TO_BIT[startSq];
    allPieces[color] &= ~INDEX_TO_BIT[endSq];

    return value;
}

int Board::valueOfPiece(int piece) {
    switch(piece) {
        case -1:
            return 0;
        case PAWNS:
            return PAWN_VALUE;
        case KNIGHTS:
            return KNIGHT_VALUE;
        case bishopS:
            return bishop_VALUE;
        case ROOKS:
            return ROOK_VALUE;
        case QUEENS:
            return QUEEN_VALUE;
        case KINGS:
            return MATE_score;
    }
    return -1;
}

// Calculates 
int Board::getMVVLVAscore(int color, Move m) {
    int endSq = getEndSq(m);
    int attacker = getPieceOnSquare(color, getStartSq(m));
    if (attacker == KINGS)
        attacker = -1;
    int victim = getPieceOnSquare(color^1, endSq);

    return (victim * 8) + (4 - attacker);
}

int Board::getExchangescore(int color, Move m) {
    int endSq = getEndSq(m);
    int attacker = getPieceOnSquare(color, getStartSq(m));
    if (attacker == KINGS)
        attacker = -1;
    int victim = getPieceOnSquare(color^1, endSq);
    return victim - attacker;
}

// MOVE GENERATION
inline uint64_t Board::getWPawnSingleMoves(uint64_t pawns) {
    return (pawns << 8) & ~getOccupancy();
}

inline uint64_t Board::getBPawnSingleMoves(uint64_t pawns) {
    return (pawns >> 8) & ~getOccupancy();
}

uint64_t Board::getWPawnDoubleMoves(uint64_t pawns) {
    uint64_t open = ~getOccupancy();
    uint64_t temp = (pawns << 8) & open;
    pawns = (temp << 8) & open & RANKS[3];
    return pawns;
}

uint64_t Board::getBPawnDoubleMoves(uint64_t pawns) {
    uint64_t open = ~getOccupancy();
    uint64_t temp = (pawns >> 8) & open;
    pawns = (temp >> 8) & open & RANKS[4];
    return pawns;
}

inline uint64_t Board::getWPawnLeftCaptures(uint64_t pawns) {
    return (pawns << 7) & NOTH;
}

inline uint64_t Board::getBPawnLeftCaptures(uint64_t pawns) {
    return (pawns >> 9) & NOTH;
}

inline uint64_t Board::getWPawnRightCaptures(uint64_t pawns) {
    return (pawns << 9) & NOTA;
}

inline uint64_t Board::getBPawnRightCaptures(uint64_t pawns) {
    return (pawns >> 7) & NOTA;
}

inline uint64_t Board::getWPawnCaptures(uint64_t pawns) {
    return getWPawnLeftCaptures(pawns) | getWPawnRightCaptures(pawns);
}

inline uint64_t Board::getBPawnCaptures(uint64_t pawns) {
    return getBPawnLeftCaptures(pawns) | getBPawnRightCaptures(pawns);
}

inline uint64_t Board::getKnightSquares(int single) {
    return KNIGHTMOVES[single];
}

uint64_t Board::getbishopSquares(int single, uint64_t occ) {
    uint64_t *attTableLoc = magicbishops[single].table;
    occ &= magicbishops[single].mask;
    occ *= magicbishops[single].magic;
    occ >>= magicbishops[single].shift;
    return attTableLoc[occ];
}

uint64_t Board::getRookSquares(int single, uint64_t occ) {
    uint64_t *attTableLoc = magicRooks[single].table;
    occ &= magicRooks[single].mask;
    occ *= magicRooks[single].magic;
    occ >>= magicRooks[single].shift;
    return attTableLoc[occ];
}

uint64_t Board::getQueenSquares(int single, uint64_t occ) {
    return getbishopSquares(single, occ) | getRookSquares(single, occ);
}

inline uint64_t Board::getKingSquares(int single) {
    return KINGMOVES[single];
}

inline uint64_t Board::getOccupancy() {
    return allPieces[WHITE] | allPieces[BLACK];
}

inline int Board::epVictimSquare(int victimColor, uint16_t file) {
    return 8 * (3 + victimColor) + file;
}

bool Board::getWhiteCanKCastle() {
    return castlingRights & WHITEKSIDE;
}

bool Board::getBlackCanKCastle() {
    return castlingRights & BLACKKSIDE;
}

bool Board::getWhiteCanQCastle() {
    return castlingRights & WHITEQSIDE;
}

bool Board::getBlackCanQCastle() {
    return castlingRights & BLACKQSIDE;
}

uint16_t Board::getEPCaptureFile() {
    return epCaptureFile;
}

uint8_t Board::getFiftyMoveCounter() {
    return fiftyMoveCounter;
}

uint16_t Board::getMoveNumber() {
    return moveNumber;
}

int Board::getPlayerToMove() {
    return playerToMove;
}

uint64_t Board::getPieces(int color, int piece) {
    return pieces[color][piece];
}

uint64_t Board::getAllPieces(int color) {
    return allPieces[color];
}

int *Board::getMailbox() {
    int *result = new int[64];
    for (int i = 0; i < 64; i++) {
        result[i] = -1;
    }
    for (int i = 0; i < 6; i++) {
        uint64_t bitboard = pieces[WHITE][i];
        while (bitboard) {
            result[bitScanForward(bitboard)] = i;
            bitboard &= bitboard - 1;
        }
    }
    for (int i = 0; i < 6; i++) {
        uint64_t bitboard = pieces[BLACK][i];
        while (bitboard) {
            result[bitScanForward(bitboard)] = 6 + i;
            bitboard &= bitboard - 1;
        }
    }
    return result;
}

uint64_t Board::getZobristKey() {
    return zobristKey;
}

void Board::initZobristKey(int *mailbox) {
    zobristKey = 0;
    for (int i = 0; i < 64; i++) {
        if (mailbox[i] != -1) {
            zobristKey ^= zobristTable[mailbox[i] * 64 + i];
        }
    }
    if (playerToMove == BLACK)
    zobristKey ^= zobristTable[768];
    zobristKey ^= zobristTable[769 + castlingRights];
    zobristKey ^= zobristTable[785 + epCaptureFile];
}

// Dumb7Fill
uint64_t fillRayRight(uint64_t rayPieces, uint64_t kosong, int shift) {
    uint64_t flood = rayPieces;
    uint64_t borderMask = 0xFFFFFFFFFFFFFFFF;
    if (shift == 1 || shift == 9)
        borderMask = NOTH;
    else if (shift == 7)
        borderMask = NOTA;
    kosong &= borderMask;
    flood |= rayPieces = (rayPieces >> shift) & kosong;
    flood |= rayPieces = (rayPieces >> shift) & kosong;
    flood |= rayPieces = (rayPieces >> shift) & kosong;
    flood |= rayPieces = (rayPieces >> shift) & kosong;
    flood |= rayPieces = (rayPieces >> shift) & kosong;
    flood |=         (rayPieces >> shift) & kosong;
    return           (flood >> shift) & borderMask;
}

uint64_t fillRayLeft(uint64_t rayPieces, uint64_t kosong, int shift) {
    uint64_t flood = rayPieces;
    uint64_t borderMask = 0xFFFFFFFFFFFFFFFF;
    if (shift == 1 || shift == 9)
        borderMask = NOTA;
    else if (shift == 7)
        borderMask = NOTH;
    kosong &= borderMask;
    flood |= rayPieces = (rayPieces << shift) & kosong;
    flood |= rayPieces = (rayPieces << shift) & kosong;
    flood |= rayPieces = (rayPieces << shift) & kosong;
    flood |= rayPieces = (rayPieces << shift) & kosong;
    flood |= rayPieces = (rayPieces << shift) & kosong;
    flood |=         (rayPieces << shift) & kosong;
    return           (flood << shift) & borderMask;
}

// Magic bbitboard 
uint64_t indexToMask64(int index, int nBits, uint64_t mask) {
    uint64_t result = 0;
    for (int i = 0; i < nBits; i++) {
        int j = bitScanForward(mask);
        mask &= mask - 1;
        if (index & INDEX_TO_BIT[i])
            result |= INDEX_TO_BIT[j];
    }
    return result;
}

uint64_t ratt(int sq, uint64_t block) {
    return fillRayRight(INDEX_TO_BIT[sq], ~block, NORTH_SOUTH_FILL) // south
         | fillRayLeft(INDEX_TO_BIT[sq], ~block, NORTH_SOUTH_FILL)  // north
         | fillRayLeft(INDEX_TO_BIT[sq], ~block, EAST_WEST_FILL)    // east
         | fillRayRight(INDEX_TO_BIT[sq], ~block, EAST_WEST_FILL);  // west
}

uint64_t batt(int sq, uint64_t block) {
    return fillRayLeft(INDEX_TO_BIT[sq], ~block, NORTHEAST_FILL)   // northeast
         | fillRayLeft(INDEX_TO_BIT[sq], ~block, NORTHWEST_FILL)   // northwest
         | fillRayRight(INDEX_TO_BIT[sq], ~block, SOUTHWEST_FILL)  // southwest
         | fillRayRight(INDEX_TO_BIT[sq], ~block, SOUTHEAST_FILL); // southeast
}

int magicMap(uint64_t masked, uint64_t magic, int nBits) {
    return (int) ((masked * magic) >> (64 - nBits));
}

uint64_t findMagic(int sq, int iBits, bool isbishop) {
    uint64_t mask, maskedBits[4096], attSet[4096], used[4096], magic;
    bool failed;

    mask = isbishop ? bishop_MASK[sq] : ROOK_MASK[sq];
    int nBits = count(mask);
    for (int i = 0; i < (1 << nBits); i++) {
        maskedBits[i] = indexToMask64(i, nBits, mask);
        attSet[i] = isbishop ? batt(sq, maskedBits[i]) : ratt(sq, maskedBits[i]);
    }

    for (int k = 0; k < 100000000; k++) {
        magic = magicPRNG() & magicPRNG() & magicPRNG();
        if (count((mask * magic) & 0xFFF0000000000000ULL) < 10)
            continue;
        for (int i = 0; i < 4096; i++)
            used[i] = 0;

        failed = false;
        for (int i = 0; !failed && i < (1 << nBits); i++) {
            int mappedIndex = magicMap(maskedBits[i], magic, iBits);

            if (!used[mappedIndex])
                used[mappedIndex] = attSet[i];
            else if (used[mappedIndex] != attSet[i])
                failed = true;
        }		
        if (!failed)
            return magic;	
    }	
    return 0;
}
