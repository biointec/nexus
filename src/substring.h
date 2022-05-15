/******************************************************************************
 *  Nexus: Pan-genome compacted de Bruijn graphs with support for approximate *
 *      pattern matching using search schemes                                 *
 *                                                                            *
 *  Copyright (C) 2022 - Lore Depuydt <lore.depuydt@ugent.be>,                *
 *                       Luca Renders <luca.renders@ugent.be> and             *
 *                       Jan Fostier <jan.fostier@ugent.be>                   *
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

#ifndef SUBSTRING_H
#define SUBSTRING_H

#include "range.h"

#include <iostream> // used for printing
// ============================================================================
// ENUMS
// ============================================================================

/**
 * An enum for the direction of the search
 */
enum Direction { FORWARD, BACKWARD };

// ============================================================================
// CLASS SUBSTRING
// ============================================================================

class Substring {
  private:
    // pointer to the string this is a substring of
    const std::string* text;
    // the startIndex of this substring in the text
    length_t startIndex;
    // the endIndex of this substring in the text (non-inclusive)
    length_t endIndex;
    // The direction of this substring
    Direction d;
    // Length of the reference text
    length_t textLength;

  public:
    /**
     * @brief Construct a new Substring object, the start and end index default
     * to 0 and the size of the text
     *
     * @param t the text to point to
     * @param dir the direction (defaults to FORWARD)
     */
    Substring(const std::string& t, Direction dir = FORWARD)
        : startIndex(0), endIndex(t.size()), d(dir), textLength(t.size()) {
        text = &t;
    }

    /**
     * @brief Construct a new Substring object, the direction defaults to
     * FORWARD
     *
     * @param t the text to point to
     * @param start the start index of this substring in t
     * @param end the end index of this substring in t (non-inclusive)
     * @param dir the direction (defaults to FORWARD)
     */
    Substring(const std::string* t, length_t start, length_t end,
              Direction dir = FORWARD)
        : text(t), startIndex(start), endIndex(end), d(dir),
          textLength(t->size()) {
    }

    /**
     * @brief Construct a substring of the text another substring points to
     *
     * @param s pointer to the other substring
     * @param start the start index of this new substring in the original text
     * @param end the end index of this new substring
     */
    Substring(const Substring* s, length_t start, length_t end)
        : text(s->text), startIndex(start), endIndex(end), d(s->d),
          textLength(s->getTextLength()) {
    }

    /**
     * @brief Construct a new substring of the text another substring points to
     *
     * @param s the other substring
     * @param start the start index of this new substring in the original text
     * @param end the end index of this new substring
     */
    Substring(const Substring* s, length_t start, length_t end, Direction dir)
        : text(s->text), startIndex(start), endIndex(end), d(dir),
          textLength(s->getTextLength()) {
    }
    /**
     * @brief Construct a new substring of the text another substring points to
     *
     * @param s the other substring
     * @param start the start index of this new substring in the original text
     * @param end the end index of this new substring
     */
    Substring(const Substring& s, length_t start, length_t end)
        : text(s.text), startIndex(start), endIndex(end), d(s.d),
          textLength(s.getTextLength()) {
    }

    /**
     * Constructs a substring of the text another substring points and sets the
     * @param s, the other substring
     * @param start, the start index of this new substring in the original text
     * @param end, the end index of this new substring
     */
    Substring(const Substring& s, length_t start, length_t end, Direction dir)
        : text(s.text), startIndex(start), endIndex(end), d(dir),
          textLength(s.getTextLength()) {
    }

    /**
     * @brief Set the direction of this substring
     *
     * @param nd the new direction
     */
    void setDirection(Direction nd) {
        d = nd;
    }

    /**
     * @brief Get a substring of a substring, skipping the first skip characters
     * (relative to the direction)
     *
     * @param skip the number of characters to skip
     * @return const Substring - the new substring
     */
    const Substring getSubPiece(length_t skip) const {
        if (d == FORWARD) {
            return Substring(this, startIndex + skip, endIndex);
        } else {
            return Substring(this, startIndex, endIndex - skip);
        }
    }

    /**
     * @brief Operator overloading. Get the character at index i of this
     * substring.
     *
     * @param i the index to get the character from
     * @return char - the character at index i
     */
    char operator[](length_t i) const {
        return (d == FORWARD) ? text->at((startIndex + i) % textLength)
                              : text->at((endIndex - 1 - i) % textLength);
    }

    /**
     * @brief Get the size of this substring
     *
     * @return length_t - the size of this substring
     */
    length_t size() const {
        if (empty()) {
            return 0;
        }
        return endIndex - startIndex;
    }

    /**
     * @brief Get the length of this substring (equals the size)
     *
     * @return length_t - the length of this substring
     */
    length_t length() const {
        return size();
    }

    /**
     * Check if this substring is empty
     * @returns a bool that indicates whether the substring was empty
     */
    bool empty() const {
        return endIndex <= startIndex;
    }

    /**
     * @brief Get the end of this substring
     *
     * @return length_t - the end Index of this substring (non-inclusive)
     */
    length_t end() const {
        return endIndex;
    }

    /**
     * @brief Get the begin of this substring
     *
     * @return length_t - the begin index of the substring
     */
    length_t begin() const {
        return startIndex;
    }

    /**
     * @brief Operator overloading. This object is made equal to another
     * substring.
     *
     * @param other The other substring to which this must be equal
     * @return Substring& - this
     */
    Substring& operator=(const Substring& other) {
        this->text = other.text;
        this->startIndex = other.begin();
        this->endIndex = other.end();
        this->d = other.d;

        return *this;
    }

    /**
     * @brief Convert the Substring to a c++ string
     *
     * @return std::string - the Substring as a c++ string
     */
    std::string tostring() const {
        if (empty()) {
            return "";
        }
        if (endIndex > textLength) {
            std::cout
                << "Warning: substring is longer than the original string "
                   "itself. Only part of the substring is returned."
                << std::endl;
        }
        return text->substr(startIndex, endIndex - startIndex);
    }

    /**
     * @brief Set the end index of this substring
     *
     * @param newEnd the new end index (non-inclusive)
     */
    void setEnd(length_t newEnd) {
        endIndex = newEnd;
    }

    /**
     * @brief Set the begin index of this substring
     *
     * @param n the new begin index
     */
    void setBegin(length_t newBegin) {
        startIndex = newBegin;
    }

    /**
     * @brief Increment the end index
     *
     */
    void incrementEnd() {
        endIndex++;
    }

    /**
     * @brief Decrement the begin index
     *
     */
    void decrementBegin() {
        startIndex--;
    }

    /**
     * @brief Get the length of the reference text
     *
     * @return length_t - the length of the reference text
     */
    length_t getTextLength() const {
        return textLength;
    }
};

#endif