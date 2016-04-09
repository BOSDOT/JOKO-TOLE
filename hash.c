/*******************************************************************************                                                                    
						     Joko Tole Chess Engine
						   Copyright (C) 2016 BOSDOT						
             a lot inspired by chessprogramming.wikispaces.com and 
		    many others engine stockfish, Gull, glaurung, fruit, etc    
		   as long as no code is not fail at build engine can be used 
--------------------------------------------------------------------------------					 					  					  
 ******************************************************************************/

#include "hash.h"

uint64_t packHashData(int depth, Move m, int score, uint8_t nodeType, uint8_t age) {
    uint64_t data = 0;
    data |= (uint8_t) depth;
    data <<= 8;
    data |= age;
    data <<= 8;
    data |= nodeType;
    data <<= 16;
    data |= m;
    data <<= 16;
    data |= (uint16_t) score;

    return data;
}

Hash::Hash(uint64_t MB) {

    uint64_t arrSize = MB << 20;
    arrSize /= sizeof(HashNode);
    table = new HashNode[arrSize];
    size = arrSize;
    keys = 0;
}

Hash::~Hash() {
    delete[] table;
}

// Assumes key has
void Hash::add(Board &b, uint64_t data, int depth, uint8_t age) {
    uint64_t h = b.getZobristKey();
    uint64_t index = h % size;
    HashNode *node = &(table[index]);
    if (node->slot1.zobristKey == 0) {
        keys++;
        node->slot1.setEntry(b, data);
        return;
    }
    else if (node->slot2.zobristKey == 0) {
        keys++;
        node->slot2.setEntry(b, data);
        return;
    }
    else { 
        if ((node->slot1.zobristKey ^ node->slot1.data) == b.getZobristKey()) {
            if (getHashAge(node->slot1.data) != age)
                keys++;
            node->slot1.setEntry(b, data);
        }
        else if ((node->slot2.zobristKey ^ node->slot2.data) == b.getZobristKey()) {
            if (getHashAge(node->slot2.data) != age)
                keys++;
            node->slot2.setEntry(b, data);
        }
// Replace an entry from a previous search space
        else {
            HashEntry *toReplace = NULL;
            int score1 = 128*((int) (age - getHashAge(node->slot1.data)))
                + depth - getHashDepth(node->slot1.data);
            int score2 = 128*((int) (age - getHashAge(node->slot2.data)))
                + depth - getHashDepth(node->slot2.data);
            if (score1 >= score2)
                toReplace = &(node->slot1);
            else
                toReplace = &(node->slot2);
			
// The node must higher depth if from the same search space.
            if (score1 < 0 && score2 < 0)
                toReplace = NULL;

            if (toReplace != NULL) {
                if (getHashAge(toReplace->data) != age)
                    keys++;
                toReplace->setEntry(b, data);
            }
        }
    }
}

// Get the hash entry, if any, associated with a board b.
uint64_t Hash::get(Board &b) {
    uint64_t h = b.getZobristKey();
    uint64_t index = h % size;
    HashNode *node = &(table[index]);

    if((node->slot1.zobristKey ^ node->slot1.data) == b.getZobristKey())
        return node->slot1.data;
    else if ((node->slot2.zobristKey ^ node->slot2.data) == b.getZobristKey())
        return node->slot2.data;

    return 0;
}

uint64_t Hash::getSize() {
    return (2 * size);
}

void Hash::setSize(uint64_t MB) {
    delete[] table;
    
    uint64_t arrSize = MB << 20;
    arrSize /= sizeof(HashNode);
    size = arrSize;

    table = new HashNode[size];
    keys = 0;
}

void Hash::clear() {
    delete[] table;
    table = new HashNode[size];
    keys = 0;
}
