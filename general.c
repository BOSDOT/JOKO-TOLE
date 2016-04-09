/*******************************************************************************                                                                    
 Joko Tole 1.0.0 beta Chess Engine Copyright (C) 2016 BOSDOT
 a lot inspired by chessprogramming.wikispaces.com and many 
 others engine stockfish, Gull, glaurung, fruit, etc
 as long as no code is not fail at build engine can be used 
--------------------------------------------------------------------------------					 					  					  
 ******************************************************************************/
#include "general.h"

// Used for bitscan
const int index64[64] = {
0,  47,  1, 56, 48, 27,  2, 60,
57, 49, 41, 37, 28, 16,  3, 61,
54, 58, 35, 52, 50, 42, 21, 44,
38, 32, 29, 23, 17, 11,  4, 62,
46, 55, 26, 59, 40, 36, 15, 53,
34, 51, 20, 43, 31, 22, 10, 45,
25, 39, 14, 33, 19, 30,  9, 24,
13, 18,  8, 12,  7,  6,  5, 63 };

const int lsb_64_table[64] =
{
   63, 30,  3, 32, 59, 14, 11, 33,
   60, 24, 50,  9, 55, 19, 21, 34,
   61, 29,  2, 53, 51, 23, 41, 18,
   56, 28,  1, 43, 46, 27,  0, 35,
   62, 31, 58,  4,  5, 49, 54,  6,
   15, 52, 12, 40,  7, 42, 45, 16,
   25, 57, 48, 13, 10, 39,  8, 44,
   20, 47, 38, 22, 17, 37, 36, 26
};

int bitScanForward(uint64_t bb) {
	#if USE_INLINE_ASM
	asm ("bsfq %1, %0" : "=r" (bb) : "r" (bb));
	return (int) bb;
	#else
	unsigned int folded;
	assert (bb != 0);
	return popCount( (bb & -bb) - 1 );
	#endif
}

int bitScanReverse(uint64_t bb) {
    #if USE_INLINE_ASM
        asm ("bsrq %1, %0" : "=r" (bb) : "r" (bb));
        return (int) bb;
    #else
        uint64_t debruijn64 = 0x03f79d71b4cb0a89;
        bb |= bb >> 1; 
        bb |= bb >> 2;
        bb |= bb >> 4;
        bb |= bb >> 8;
        bb |= bb >> 16;
        bb |= bb >> 32;
        return index64[(int) ((bb * debruijn64) >> 58)];
    #endif
}

int count(uint64_t bb) {
    #if USE_INLINE_ASM
        asm ("popcntq %1, %0" : "=r" (bb) : "r" (bb));
        return (int) bb;
    #else
        bb = bb - ((bb >> 1) & 0x5555555555555555);
        bb = (bb & 0x3333333333333333) + ((bb >> 2) & 0x3333333333333333);
        bb = (((bb + (bb >> 4)) & 0x0F0F0F0F0F0F0F0F) * 0x0101010101010101) >> 56;
        return (int) bb;
    #endif
}

uint64_t flipAcrossRanks(uint64_t bb) {
    #if USE_INLINE_ASM
        asm ("bswapq %0" : "=r" (bb) : "0" (bb));
        return bb;
    #else
        bb = ((bb >>  8) & 0x00FF00FF00FF00FF) | ((bb & 0x00FF00FF00FF00FF) <<  8);
        bb = ((bb >> 16) & 0x0000FFFF0000FFFF) | ((bb & 0x0000FFFF0000FFFF) << 16);
        bb =  (bb >> 32) | (bb << 32);
        return bb;
    #endif
}

unsigned char_BitScanForward64(unsigned long * Index,  unsigned __int64 Mask);
unsigned char _BitScanReverse64(unsigned long * Index,  unsigned __int64 Mask);
unsigned __int64 __lzcnt64(unsigned __int64 value); // AMD K10 only see CPUID

double getTimeElapsed(ChessTime startTime) {
    auto endTime = ChessClock::now();
    std::chrono::duration<double> timeSpan =
        std::chrono::duration_cast<std::chrono::duration<double>>(
        endTime-startTime);
    return timeSpan.count();
}

std::string moveToString(Move m) {
    char startFile = 'a' + (getStartSq(m) & 7);
    char startRank = '1' + (getStartSq(m) >> 3);
	char endFile = 'a' + (getEndSq(m) & 7);
	char endRank = '1' + (getEndSq(m) >> 3);
    std::string moveStr = {startFile, startRank, endFile, endRank};
    if (getPromotion(m)) moveStr += " nbrq"[getPromotion(m)];
    return moveStr;
}

