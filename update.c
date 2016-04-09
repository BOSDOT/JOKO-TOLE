/*******************************************************************************                                                                    
						     Joko Tole Chess Engine
						   Copyright (C) 2016 BOSDOT						
             a lot inspired by chessprogramming.wikispaces.com and 
		    many others engine stockfish, Gull, glaurung, fruit, etc    
		   as long as no code is not fail at build engine can be used 
--------------------------------------------------------------------------------					 					  					  
 ******************************************************************************/
 
uint_64 attacks_table[...]; // ~840K Byte all rook and bishop attacks
 
struct SMagic {
   uint_64* ptr;  // pointer to the attack-table for one particular square
   uint_64 mask;  // to mask of both lines
   uint_64 magic; // magic 64-bit factor
   int shift; // shift right 
};
 
SMagic mBishopTbl[64];
SMagic mRookTbl[64];
 
uint_64 bishopAttacks(uint_64 occ, enumSquare sq) {
   uint_64* aptr = mBishopTbl[sq].ptr;
   occ      &= mBishopTbl[sq].mask;
   occ      *= mBishopTbl[sq].magic;
   occ     >>= mBishopTbl[sq].shift;
   return aptr[occ];
}
 
uint_64 rookAttacks(uint_64 occ, enumSquare sq) {
   uint_64* aptr = mRookTbl[sq].ptr;
   occ      &= mRookTbl[sq].mask;
   occ      *= mRookTbl[sq].magic;
   occ     >>= mRookTbl[sq].shift;
   return aptr[occ];
}