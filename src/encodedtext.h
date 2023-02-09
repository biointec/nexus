/******************************************************************************
 *  Columba: Approximate Pattern Matching using Search Schemes                *
 *  Copyright (C) 2020-2022 - Luca Renders <luca.renders@ugent.be> and        *
 *                            Jan Fostier <jan.fostier@ugent.be>              *
 *                                                                            *
 *  This program is free software: you can redistribute it and/or modify      *
 *  it under the terms of the GNU Affero General Public License as            *
 *  published by the Free Software Foundation, either version 3 of the        *
 *  License, or (at your option) any later version.                           *
 *                                                                            *
 *  This program is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU Affero General Public License for more details.                       *
 *                                                                            *
 * You should have received a copy of the GNU Affero General Public License   *
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.     *
 ******************************************************************************/

#ifndef ENCODEDTEXT_H
#define ENCODEDTEXT_H

#include "alphabet.h"

#include <array>
#include <bitset>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdint.h>
#include <vector>

/**
 * ENCODEDTEXT CLASS
 *
 * Encodes a string where each character takes up B = ceil(log2(S)) bits
 */
template <size_t S> // S is the size of the alphabet (including '$')
class EncodedText {

  private:
    // the number of bits per character
    const static uint64_t B;
    // bitmask where the first B bits are 1, and other bits are 0
    const static uint64_t bitmask;
    // 64-bit word where the bit at index i indicates whether a symbol starting
    // at index i overflows into the next word
    const static uint64_t hasOverflowBits;
    // masks for overflow (either all 0's or all 1's)
    const static std::array<uint64_t, 2> overflowMasks;

    // the encoded text, where each symbol takes up B bits
    std::vector<uint64_t> encodedText;
    // Size (in symbols) of the text
    size_t tSize;

    /**
     * Helper function to get a 64-bit containing either all 1's or all 0's
     * @param index the index in the word
     * @returns all 1's if a symbol starting at index overflows into next word,
     * all 0's if the symbol starting at index does not overflow into the next
     * word
     */
    static uint64_t hasOverflow(uint64_t index) {
        assert(index < 64);
        // A) get the value at the correct bit
        uint64_t maskIndex =
            ((1ull << (63 - index)) & hasOverflowBits) >> (63 - index);
        // B) return the mask
        return overflowMasks[maskIndex];
    }

    /**
     * Encode the letter with charIndex at place index in the text into the
     * bitvector
     * @param charIndex the index of the character in the alphabet
     * @param index the index of this character in the text
     */
    void encodeLetter(const uint64_t charIndex, const uint64_t index) {

        uint64_t bits = charIndex;

        // A) find word index of first bit
        uint64_t w = (index * B) / 64;

        // B) find bit index of first bit
        uint64_t b = (index * B) % 64;

        // C) split bits in bits for w1 and bits for w2

        // find the mask for the bits that are in the second word
        uint64_t mask = hasOverflow(b) & ((-1ull) ^ (-1ull << (B - (64 - b))));

        // create the bits to be set in the second word
        uint64_t bits2 = bits & mask;

        // indicates the number of bits in the next word
        uint64_t maskSize = __builtin_popcountll(mask);
        // shift the bits of 1st word to the right (shift out bits for 2nd word)
        bits >>= maskSize;

        // D) add bits to word 1
        uint64_t startLocation = 64 - B + maskSize;
        uint64_t shift = startLocation - b;
        encodedText[w] |= (bits << shift);

        // E) add bits to word 2 (in case of overflow)
        uint64_t shift2 = 64 - maskSize;
        encodedText[w + 1] |= (bits2 << shift2);
    }

  public:
    /**
     * Default constructor
     */
    EncodedText() {
        tSize = 0;
    }

    /**
     * Constructor: encodes the text given the alphabet
     * @param sigma the alphabet to use
     * @param text the text to encode
     */
    EncodedText(const Alphabet<S>& sigma, const std::string& text)
        : tSize(text.size()) {
        encodedText.resize(((text.size() * B) / 64) + 1);

        for (uint64_t i = 0; i < text.size(); i++) {
            encodeLetter(sigma.c2i(text[i]), i);
        }
    }

    /**
     * Constructor for a text filled with first character of alphabet (= symbol
     * with only 0 bits)
     * @param size the size of the text
     */
    EncodedText(const uint64_t size) : tSize(size) {
        encodedText.resize(((size * B) / 64) + 1);
    }

    // ----------------------------------------------------------------------------
    // DECODING: convert symbolindex to character
    // ----------------------------------------------------------------------------

    /**
     * Gets the letter at index in the text
     * @param sigma the alphabet
     * @param index the index in the text
     * @returns a char which contains the character at index index in the
     * original text according to the passed alphabet
     */
    char decodeLetter(const Alphabet<S>& sigma, const uint64_t index) const {
        uint64_t cIdex = getEncodedLetter(index % tSize);
        return sigma.i2c(cIdex);
    }

    /**
     * Decode the entire text
     * @param sigma the alphabet
     * @returns the decoded text
     */
    std::string decodeText(const Alphabet<S>& sigma) const {
        std::string text;
        text.resize(tSize);

        for (uint64_t i = 0; i < tSize; i++) {
            text[i] = decodeLetter(sigma, i);
        }

        return text;
    }

    /**
     * Decode the entire text
     * @param sigma the alphabet
     * @returns the original text
     */
    std::string decodeSubstring(const Alphabet<S>& sigma, length_t begin,
                                length_t end) const {
        std::string text;
        text.resize(end - begin);

        for (uint64_t i = begin; i < end; i++) {
            text[i - begin] = decodeLetter(sigma, i);
        }

        return text;
    }

    // ----------------------------------------------------------------------------
    // ACCESS OPERATIONS
    // ----------------------------------------------------------------------------

    /**
     * @returns the size (in characters) of the original text
     */
    size_t size() const {
        return tSize;
    }

    /**
     * Operator overloading, gets the characterindex (according to the used
     * alphabet) of the character that was present at index index in the
     * original text
     * @param index the index to find the character of in the original text
     */
    uint64_t operator[](const uint64_t index) const {
        // A) find word index of first bit of this symbol
        uint64_t w = (index * B) / 64;

        // B) find bit index of first bit of this symbol
        uint64_t b = (index * B) % 64;

        // C) get the bits of the symbol that are in the first word
        uint64_t bits = (encodedText[w] & (bitmask >> b)) << b;

        // D) get bitmask for bits that flow over into next word
        uint64_t mask = hasOverflow(b) & ((-1ull) ^ (-1ull >> (B - (64 - b))));

        // E) get bits in next word
        uint64_t bitsNext = encodedText[w + 1] & mask;

        return (bits >> (64 - B)) +
               (bitsNext >> (64 - __builtin_popcountll(mask)));
    }

    /**
     * @brief Set the letter with charIndex at place index in the text into the
     * bitvector
     *
     * @param index the index of this character in the text
     * @param charIndex the index of the character in the alphabet
     */
    void set(const uint64_t index, const uint64_t charIndex) {
        encodeLetter(charIndex, index);
    }

    /**
     * Gets the characterindex (according to the used
     * alphabet) of the character that was present at index index in the
     * original text
     * @param index the index to find the character of in the original text
     */
    uint64_t getEncodedLetter(const uint64_t index) const {
        return operator[](index);
    }

    // ----------------------------------------------------------------------------
    // I/O operations
    // ----------------------------------------------------------------------------

    /**
     * Write encoded text to disk
     * @param filename File name
     */
    void write(const std::string& filename) {
        std::ofstream ofs(filename);
        if (!ofs)
            throw std::runtime_error("Cannot open file: " + filename);

        ofs.write((char*)&tSize, sizeof(size_t));
        size_t vectorSize = encodedText.size();
        ofs.write((char*)&vectorSize, sizeof(size_t));
        ofs.write((char*)encodedText.data(),
                  encodedText.size() * sizeof(uint64_t));
    }

    /**
     * Load encoded text from disk
     * @param filename File name
     */
    bool read(const std::string& filename) {
        std::ifstream ifs(filename);
        if (!ifs)
            return false;

        ifs.read((char*)&tSize, sizeof(size_t));
        size_t vectorSize;
        ifs.read((char*)&vectorSize, sizeof(size_t));
        encodedText.resize(vectorSize);
        ifs.read((char*)&encodedText[0], encodedText.size() * sizeof(uint64_t));

        return true;
    }

    /**
     * @brief Clear the text (free up memory)
     *
     */
    void clear() {
        encodedText.clear();
        encodedText.resize(0);
        encodedText.shrink_to_fit();
        tSize = 0;
    }

    /**
     * @brief Resize the text
     *
     * @param size new size
     */
    void resize(const size_t size) {
        tSize = size;
        encodedText.resize(((size * B) / 64) + 1);
    }

    /**
     * @brief Check whether the text is empty
     *
     * @return true if empty
     * @return false otherwise
     */
    bool empty() {
        return tSize == 0;
    }
};

template <size_t S>
const std::array<uint64_t, 2> EncodedText<S>::overflowMasks = {0ull, -1ull};
template <size_t S> const uint64_t EncodedText<S>::B = ceil(log2(S));
template <size_t S>
const uint64_t EncodedText<S>::bitmask = (-1ull) ^ (-1ull >> B);
template <size_t S>
const uint64_t EncodedText<S>::hasOverflowBits =
    ~(-1ull << (B - 1)); // the bit at index i indicates whether a symbol
                         // starting at index i overflows into the next word

#endif