/*******************************************************************************                                                                    
Joko Tole UCI Chess Engine Copyright(C) 2016 Bosdot (Javanese)
influenzed by chessprogramming.wikispaces.com and many others 
chess engine with open source code.
--------------------------------------------------------------------------------					 					  					  
 ******************************************************************************/

#ifndef GENERAL_H
#define GENERAL_H

#include <cstdint>
#include <chrono>
#include <string>

#define USE_INLINE_ASM true

#if defined(_WIN64) && defined(_MSC_VER) // No Makefile used
#  include <intrin.h> // MSVC popcnt and bsfq instrinsics
#  define x_64BIT
#endif

#if defined(USE_POPCNT) && defined(__INTEL_COMPILER) && defined(_MSC_VER)
#  include <nmmintrin.h> // Intel header for _mm_popcnt_u64() intrinsic
#endif

#ifdef USE_POPCNT
const bool HasPopCnt = true;
#else
const bool HasPopCnt = false;
#endif

#ifdef IS_64BIT
const bool x64Bit = true;
#else
const bool x64Bit = false;
#endif


const int WHITE = 0;
const int BLACK = 1;
const int PAWNS = 0;
const int KNIGHTS = 1;
const int bishopS = 2;
const int ROOKS = 3;
const int QUEENS = 4;
const int KINGS = 5;

// Material constants
const int PAWN_VALUE = 100, PAWN_VALUE_EG = 103;
const int KNIGHT_VALUE = 330, KNIGHT_VALUE_EG = 412;
const int bishop_VALUE = 337, bishop_VALUE_EG = 439;
const int ROOK_VALUE = 510, ROOK_VALUE_EG = 672;
const int QUEEN_VALUE = 880, QUEEN_VALUE_EG = 1352;
const int MATE_score = 32761;
const int INFTY = 33761;

const int START_VALUE = 8 * PAWN_VALUE + 2 * KNIGHT_VALUE + 2 * bishop_VALUE + 2 * ROOK_VALUE + QUEEN_VALUE;
const int EG_FACTOR_RES = 1000;

// Other values
const int MAX_DEPTH = 100;
const int MAX_MOVES = 256;

// Stuff for timing
typedef std::chrono::high_resolution_clock ChessClock;
typedef std::chrono::high_resolution_clock::time_point ChessTime;

double getTimeElapsed(ChessTime startTime);

// Bitboard methods
int bitScanForward(uint64_t bb);
int bitScanReverse(uint64_t bb);
int count(uint64_t bb);
uint64_t flipAcrossRanks(uint64_t bb);

typedef uint16_t Move;

enum Phase {
  PHASE_ENDGAME,
  PHASE_MIDGAME = 128,
  MG = 0, EG = 1
};

const Move NULL_MOVE = 65;
const uint16_t MOVE_DOUBLE_PAWN = 0x1;
const uint16_t MOVE_EP = 0x5;
const uint16_t MOVE_PROMO_N = 0x8;
const uint16_t MOVE_PROMO_B = 0x9;
const uint16_t MOVE_PROMO_R = 0xA;
const uint16_t MOVE_PROMO_Q = 0xB;
const int PROMO[16] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 1, 2, 3, 4};

inline Move encodeMove(int startSq, int endSq) {
    return (endSq << 6) | startSq;
}

inline Move setCapture(Move m, bool isCapture) {
    return m | (isCapture << 14);
}

inline Move setCastle(Move m, bool isCastle) {
    return m | (isCastle << 13);
}

inline Move setFlags(Move m, uint16_t f) {
    return m | (f << 12);
}

inline int getStartSq(Move m) {
    return (int) (m & 0x3F);
}

inline int getEndSq(Move m) {
    return (int) ((m >> 6) & 0x3F);
}

inline int getPromotion(Move m) {
    return PROMO[m >> 12];
}

inline bool isPromotion(Move m) {
    return (bool) (m >> 15);
}

inline bool isCapture(Move m) {
    return (bool) ((m >> 14) & 1);
}

inline bool isCastle(Move m) {
    return ((m >> 13) == 1);
}

inline bool isEP(Move m) {
    return ((m >> 12) == MOVE_EP);
}

inline uint16_t getFlags(Move m) {
    return (m >> 12);
}

std::string moveToString(Move m);

template <class T> class SearchArrayList {
public:
    T arrayList[MAX_MOVES];
    unsigned int length;

    SearchArrayList() {
        length = 0;
    }
    ~SearchArrayList() {}

    unsigned int size() { return length; }

    void add(T o) {
        arrayList[length] = o;
        length++;
    }

    T get(int i) { return arrayList[i]; }

    void set(int i, T o) { arrayList[i] = o; }

    T remove(int i) {
        T deleted = arrayList[i];
        for(unsigned int j = i; j < length-1; j++) {
            arrayList[j] = arrayList[j+1];
        }
        length--;
        return deleted;
    }

    void swap(int i, int j) {
        T temp = arrayList[i];
        arrayList[i] = arrayList[j];
        arrayList[j] = temp;
    }

    void clear() {
        length = 0;
    }
};

typedef SearchArrayList<Move> MoveList;
typedef SearchArrayList<int> scoreList;

#endif
