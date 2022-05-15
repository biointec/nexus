/*
 * radix.h
 *
 *  Created on: Oct 21, 2011
 *      Author: marius
 *  Added to github.uconn.edu on Jan 22, 2013
 */

#ifndef RADIX_H_
#define RADIX_H_

#define _NDEBUG
#include <cassert>
#include <vector>

#ifndef _NDEBUG
#define _DEBUG
#endif

#include <math.h>
#define LITTLE_ENDIAN_FLAG

#include "RadixLSDCache.h"
#include "utils.h"

template <class unum> class Radix {
  private:
    typedef unsigned long long word;
    static const int bitsPerWord = 8 * sizeof(word);
#define DIV_BY_WORDSIZE(p) ((p) >> 6)
#define MOD_WORDSIZE(p) ((p)&63)
    static const int maxChar = 256;

    static const bool DoubleNumWord = (sizeof(unum) * 2 <= sizeof(word));
    static const int unumBits = 8 * sizeof(unum);
    const uchar* originalInput;
    const unum length;

    unum* sa;

    vector<uchar> bucketFlag;

    vector<unum> bucket;

    vector<unum> indexBuff1;
    vector<unum> indexBuff2;

    RadixLSDCache<unum, word, unum> sorter;

    const int bucketPiggyBackBits;
    const unum bucketPiggyBackLimit;
    const unum bucketPiggyBackMask;
    const unum kmerLength;

    // Packs the original input into a buffer of words.
    // 'charCode' maps each character of the original input to an index between
    // 0 and alphabet size.
    void packInput(const uchar* originalInput, const unum length,
                   const uint* charCode, const int bitsPerChar,
                   word* packedInput) {
        int remainBits = bitsPerWord;
        word l = 0;
        unum nWords = 0;
        for (unum i = 0; i < length; ++i) {
            word code = charCode[originalInput[i]];
            remainBits -= bitsPerChar;
            if (remainBits >= 0) {
                l = (l << bitsPerChar) | code;
            } else {
                int bitsInNextWord = -remainBits;
                int bitsAtEnd = bitsPerChar - bitsInNextWord;
                l = (l << bitsAtEnd) | (code >> bitsInNextWord);
                packedInput[nWords++] = l;

                // normally we would only need code's final bitsInNextWord bits,
                // but the left shifting will get rid of the extra bits later
                // anyway
                //                l = GET_LAST(code, bitsInNextWord);
                l = code;

                remainBits = bitsPerWord - bitsInNextWord;
            }
        }
        if (remainBits < bitsPerWord)
            packedInput[nWords++] = l << remainBits;
        packedInput[nWords] = 0;
    }

    template <class W> void shiftLeft(uint* code, int bits, W* shiftedCode) {
        for (int i = 0; i < maxChar; ++i)
            shiftedCode[i] = ((W)code[i]) << bits;
    }

    // Sorts all suffixes by their first bits, thus grouping the suffixes into
    // buckets. Returns a vector (call it bStart) of bucket starting positions
    // in the suffix array (i.e. bStart[0] = 0;
    //       bStart[1] = size of first bucket;
    //       bStart[2] = bStart[1] + size of second bucket;
    //  etc.)
    vector<unum> sortByFirstBits(uint* charCode, int bitsPerChar,
                                 int bitsPerPass) {
        word shiftedCode[maxChar];
        shiftLeft(charCode, bitsPerPass - bitsPerChar, shiftedCode);

        unum n = 1 << bitsPerPass;
        vector<unum> buffer(n + 2);

        // Accumulate characters from the original input into a word,
        // and count the frequency of each word
        {
            // bucket size are stored from position 2 of the vector
            unum* bSize = &buffer[2];
            memset(bSize, 0, n * sizeof(unum));
            word w = 0;
            for (unum i = length; i--;) {
                int c = originalInput[i];
                w >>= bitsPerChar;
                w |= shiftedCode[c];
                bSize[w]++;
            }
        }

        // bucket sizes are stored from position 2 so we put a zero in front of
        // them (at position 1) and do a prefix sum to obtain bucket starting
        // positions
        unum* bStart = &buffer[1];
        bStart[0] = 0;
        prefixSum(bStart, bStart, n);

        // Rearrange suffixes into the corresponding buckets given by their
        // first bits.
        {
            word w = 0;
            for (unum i = length; i--;) {
                int c = originalInput[i];
                w >>= bitsPerChar;
                w |= shiftedCode[c];
                sa[bStart[w]++] = i;
            }
        }

        // Bucket starts have been modified during redistribution of suffixes.
        // However, bucketStart[i] now has the value for bucketStart[i+1].
        // To fix this we just shift the whole array by putting a zero in front
        // of it!
        bStart = bStart - 1;
        bStart[0] = 0;

        // bStart is now at the beginning of the buffer so we just return the
        // buffer
        return buffer;
    }

    // Returns an entire word of the packed input starting at the given bit.
    word getWordAtBit(const word* input, unum p) {
        unum w = DIV_BY_WORDSIZE(p);
        int b = MOD_WORDSIZE(p);
        if (b)
            return (input[w] << b) | (input[w + 1] >> (bitsPerWord - b));
        else
            return input[w];
    }

    // For each suffix, copies the next word after skipping bitsToSkip bits.
    void copyWords(unum* suffix, unum len, const word* input, word* wordBuffer,
                   const int bitsPerChar, const int bitsToSkip) {
        for (; len--;) {
            unum bit = (*suffix++) * bitsPerChar + bitsToSkip;
            *wordBuffer++ = getWordAtBit(input, bit);
        }
    }

    // Subdivides a bucket wherever consecutive keys are different.
    // Uses the given flag to mark the start of the sub-buckets.
    void divideBucket(word* key, unum length, unum bStart, uchar flag) {
        uchar* bf = &bucketFlag[bStart];
        bf[0] = flag;
        for (unum i = 1; i < length; ++i)
            if (key[i] != key[i - 1])
                bf[i] = flag;
    }

    // Sorts a set of buckets by their next word after the first bitsToSkip.
    void sortBuckets(const int bitsToSkip, const unum nBuckets,
                     const unum* bStart, const int bitsPerChar,
                     const word* input, word* buffer,
                     RadixLSDCache<unum, word, unum>& sorter) {
        for (unum i = 0; i < nBuckets; ++i) {
            unum start = bStart[i];
            unum bLen = bStart[i + 1] - start;
            if (bLen > 1) {
                unum* lo = sa + start;
                copyWords(lo, bLen, input, buffer, bitsPerChar, bitsToSkip);
                sorter.sort(bLen, lo, buffer);
                divideBucket(buffer, bLen, start, 1);
            }
        }
    }

    // Sorts a set of buckets by their next word after the first bitsToSkip.
    // The start of a bucket is marked by a non-zero flag at that position in
    // the bucketFlag.
    void sortFlagBuckets(const int bitsToSkip, const uchar* bucketFlag,
                         unum start, unum length, const int bitsPerChar,
                         const word* input, word* buffer,
                         RadixLSDCache<unum, word, unum>& sorter) {
        unum bStart = start;
        unum bLen = 1;
        for (unum i = bStart + 1, e = bStart + length; i <= e; ++i) {
            if (bucketFlag[i]) {
                if (bLen > 1) {
                    unum* lo = sa + bStart;

                    copyWords(lo, bLen, input, buffer, bitsPerChar, bitsToSkip);
                    sorter.sort(bLen, lo, buffer);
                    divideBucket(buffer, bLen, bStart, 1);
                }
                bStart = i;
                bLen = 0;
            }
            bLen++;
        }
    }

    word getBucketWord(unum s1, unum D) {
        if (DoubleNumWord) {
            return (((word)bucket[s1]) << unumBits) | bucket[s1 + D];
        } else {
            return bucket[s1];
        }
    }

    void copyBucketNumbers(unum* sa, unum D, unum D2, word* buf, unum bLen) {
        for (; bLen--;)
            *buf++ = getBucketWord((*sa++) + D, D2);
    }

    void assignSubBucketNumbersAndLen(unum bStart, unum bLen) {
        unum currentBStart = bStart;
        for (unum e = bStart + bLen, i = bStart + 1; i <= e; ++i) {
            if (bucketFlag[i]) {
                unum prevBucketLen = i - currentBStart;
                unum buck = ((currentBStart + 1) << bucketPiggyBackBits);
                if (prevBucketLen < bucketPiggyBackLimit)
                    buck |= prevBucketLen;

                for (; currentBStart < i;)
                    bucket[sa[currentBStart++]] = buck;
            }
        }
    }

    void detectPeriods(unum* lo, unum bLen, unum D, unum& periodLength,
                       unum& maxNbPeriods) {
        periodLength = 0;
        maxNbPeriods = 0;
        unum currentPeriodSize = 1;
        for (unum i = 1; i < bLen; ++i) {
            unum diff = lo[i - 1] - lo[i];
            if (diff <= D) { // period
                if (periodLength && diff != periodLength) {
                    periodLength = 0;
                    maxNbPeriods = 0;
                    return; // no good
                }
                periodLength = diff;
                currentPeriodSize++;
            } else {
                if (currentPeriodSize > maxNbPeriods)
                    maxNbPeriods = currentPeriodSize;
                currentPeriodSize = 1;
            }
        }
        if (currentPeriodSize > maxNbPeriods)
            maxNbPeriods = currentPeriodSize;
    }

    bool detectAndTreatPeriods(unum bStart, unum bLen, unum D, unum D2,
                               uchar bucketStartFlag) {
        unum period;
        unum maxNbPeriods;
        detectPeriods(sa + bStart, bLen, D, period, maxNbPeriods);

        if (maxNbPeriods * period > 2 * D)
            return treatPeriodLeaders(bStart, bLen, D, D2, bucketStartFlag,
                                      period);
        else
            return false;
    }

    bool treatPeriodLeaders(unum bStart, unum bLen, unum D, unum D2,
                            uchar bucketStartFlag, unum periodLength) {

        vector<word> bucketBuffVec;
        bucketBuffVec.reserve(bLen);
        word* bufferP = &bucketBuffVec[0];

        indexBuff1.reserve(bLen);
        unum* suffixBuff3P = &indexBuff1[0];

        unum* lo = sa + bStart;
        unum leaders[2] = {0, bLen - 1};
        for (unum i = 0; i < bLen;) {
            unum l = lo[i];
            word leaderEnd = getBucketWord(l + D - periodLength, D2);
            word rest = getBucketWord(l + D, D2);

            if (leaderEnd == rest)
                return false; // can happen with the copy from next optimization

            if (leaderEnd < rest) {
                bufferP[leaders[1]] = rest;
                suffixBuff3P[leaders[1]] = i;
                --leaders[1];
            } else {
                bufferP[leaders[0]] = rest;
                suffixBuff3P[leaders[0]] = i;
                ++leaders[0];
            }

            for (++i; i < bLen && lo[i - 1] - lo[i] == periodLength; ++i)
                ;
        }

        indexBuff2.reserve(bLen);
        unum* dest = &indexBuff2[0];
        unum k = 0;
        if (leaders[0] > 1) {
            sorter.sort(leaders[0], suffixBuff3P, bufferP);
        }

        if (leaders[0])
            k += populate<false>(leaders[0], periodLength, bStart, bLen,
                                 bucketStartFlag, suffixBuff3P, dest, bufferP);

        unum buffStart = leaders[1] + 1;
        leaders[1] = bLen - buffStart;
        if (leaders[1]) {
            word* bf = bufferP + buffStart;
            unum* sb = suffixBuff3P + buffStart;
            reverseArray(bf, leaders[1]);
            reverseArray(sb, leaders[1]);
            if (leaders[1] > 1) {
                sorter.sort(leaders[1], sb, bf);
            }
            k += populate<true>(leaders[1], periodLength, bStart, bLen,
                                bucketStartFlag, sb, dest, bf);
        }

        memcpy(sa + bStart, dest, bLen * sizeof(unum));
        assignSubBucketNumbersAndLen(bStart, bLen);
        return true;
    }

    template <bool isTypeS>
    int populate(unum leaders, unum periodLength, unum bStart, unum bLen,
                 unum bucketStartFlag, unum* suffixBuff, unum* dest,
                 word* buffer) {

        unum* lo = sa + bStart;
        unum end = isTypeS ? bLen - leaders : 0;
        while (leaders > 0) {
            unum start = end;
            unum s = 0;
            for (unum j = 0; j < leaders; ++j) {
                if (j == 0 || buffer[j] != buffer[j - 1]) {
                    bucketFlag[bStart + end] = bucketStartFlag;
                }

                unum i = suffixBuff[j];
                unum suf = lo[i];
                dest[end] = suf;
                ++end;

                if (i + 1 < bLen && suf - lo[i + 1] == periodLength) {
                    suffixBuff[s] = i + 1;
                    buffer[s] = buffer[j];
                    s++;
                }
            }
            leaders = s;

            if (isTypeS)
                end = start - leaders;
        }

        return isTypeS ? bLen - end : end;
    }

    void handleBucketRepeatsOnly(unum bStart, unum bLen, unum D, unum D2,
                                 uchar bucketStartFlag) {
        detectAndTreatPeriods(bStart, bLen, D, D2, bucketStartFlag);
    }

    void handleBucketNoRepeats(unum bStart, unum bLen, unum D, unum D2,
                               uchar bucketStartFlag) {
        vector<word> bucketBuffVec;
        bucketBuffVec.reserve(bLen);
        word* buffer = &bucketBuffVec[0];
        unum* lo = sa + bStart;
        copyBucketNumbers(lo, D, D2, buffer, bLen);
        sorter.sort(bLen, lo, buffer);
        divideBucket(buffer, bLen, bStart, bucketStartFlag);
        assignSubBucketNumbersAndLen(bStart, bLen);
    }

    void handleBucketWithRepeats(unum bStart, unum bLen, unum D, unum D2,
                                 uchar bucketStartFlag) {
        if (!detectAndTreatPeriods(bStart, bLen, D, D2, bucketStartFlag))
            handleBucketNoRepeats(bStart, bLen, D, D2, bucketStartFlag);
    }

    bool adjacentSuffixes(unum p1, unum p2, unum n) {
        unum* a = sa + p1;
        unum* b = sa + p2;
        for (; n--;)
            if (*a++ != *b++ - 1)
                return false;
        return true;
    }

    void copyAdjacentSuffs(unum d, unum s, unum n) {
        unum* a = sa + d;
        unum* b = sa + s;
        for (; n--;)
            *a++ = *b++ - 1;
    }

    void copyFlags(unum d, unum s, unum n, uchar f) {
        uchar* a = &bucketFlag[d];
        uchar* b = &bucketFlag[s];
        for (; n--; ++a)
            if (*b++)
                *a = f;
    }

    unum depthForFlag(uchar flag, unum expectedSorted) {
        if (DoubleNumWord) {
            return ((flag << 1) - 1) * expectedSorted;
        } else {
            return flag * expectedSorted;
        }
    }

    void finalTouches(unum expectedSorted) {
        assignSubBucketNumbersAndLen(0, length);

        const int allowedRepeats = 128;
        do {
            bool done = true;
            unum prevBLen = 0;
            unum prevBStart = 0;
            unum nextBStart = 0;
            unum bStart = 0;
            unum bLen = 0;
            bool copyFromNext = false;
            unum copyDepth = expectedSorted;

            for (unum i = length; i--;) {

                bool prevShouldCopy = false;
                {
                    unum buck = bucket[i];
                    prevBLen = buck & bucketPiggyBackMask;
                    prevBStart = (buck >> bucketPiggyBackBits) - 1;
                    if (prevBLen != 1) {
                        if (prevBStart != bStart) {
                            if (!prevBLen)
                                prevBLen = getBucketLength(prevBStart);
                            prevShouldCopy =
                                prevBLen == bLen &&
                                adjacentSuffixes(prevBStart, bStart, bLen);
                        }
                    }
                }

                if (bLen > 1) {
                    int flag = bucketFlag[bStart];
                    if (flag <= allowedRepeats) {
                        if (copyFromNext) {
                            copyDepth += 1;
                            bool treatPeriods =
                                copyDepth >= sa[bStart] - sa[bStart + 1];

                            copyAdjacentSuffs(bStart, nextBStart, bLen);
                            copyFlags(bStart, nextBStart, bLen, flag + 1);
                            assignSubBucketNumbersAndLen(bStart, bLen);

                            if (treatPeriods) {
                                unum subLen = 1;
                                for (unum i = bStart + bLen; i-- > bStart;) {
                                    if (bucketFlag[i]) {
                                        if (subLen > 1) {
                                            handleBucketRepeatsOnly(
                                                i, subLen, copyDepth,
                                                expectedSorted, flag + 1);
                                        }
                                        subLen = 1;
                                    } else
                                        ++subLen;
                                }
                            }
                        } else {
                            unum D = depthForFlag(flag, expectedSorted);
                            handleBucketWithRepeats(bStart, bLen, D,
                                                    expectedSorted, flag + 1);
                            copyDepth = depthForFlag(flag + 1, expectedSorted);
                        }
                    } else
                        done = false;
                }

                nextBStart = bStart;
                if (bStart == prevBStart) {
                    unum buck = bucket[i];
                    prevBLen = buck & bucketPiggyBackMask;
                    prevBStart = (buck >> bucketPiggyBackBits) - 1;
                    if (!prevBLen)
                        prevBLen = getBucketLength(prevBStart);
                }
                bLen = prevBLen;
                bStart = prevBStart;
                copyFromNext = prevShouldCopy;
            }

            if (done)
                break;

            for (unum i = 0; i < length; ++i) {
                uchar f = bucketFlag[i];
                if (f > 1)
                    bucketFlag[i] = 1;
            }
            if (DoubleNumWord) {
                expectedSorted *= 3;
            } else {
                expectedSorted *= 2;
            }

        } while (1);
    }

    // Sorts the suffixes by their first bits (notice that we use bits instead
    // of input characters). This is done in two or three passes. The first pass
    // bucket sorts by bitsPerFirstPass bits. The second and third pass further
    // subdivides the buckets obtained in the first pass, by the next word of
    // each suffix.
    void inputBasedSort(uint* charCode, int bitsPerChar, int bitsPerFirstPass,
                        bool doThirdPass) {
        vector<unum> bStartVec =
            sortByFirstBits(charCode, bitsPerChar, bitsPerFirstPass);

        bucketFlag = vector<uchar>(length + 1);
        bucketFlag[0] = bucketFlag[length] = 1;

        unum zeros = fixEndingZeros(charCode);
        bStartVec[0] +=
            zeros; // we can now ignore the suffixes that are all zero

        vector<word> inputVec(4 + (length * bitsPerChar / bitsPerWord));
        word* input = &inputVec[0];

        packInput(originalInput, length, charCode, bitsPerChar, input);

        const unum n = 1 << bitsPerFirstPass;

        unum maxBLen = 1;
        for (unum i = 1; i < n; ++i) {
            bucketFlag[bStartVec[i]] = 1;

            unum bLen = bStartVec[i] - bStartVec[i - 1];
            if (bLen > maxBLen) {
                maxBLen = bLen;
            }
        }
        vector<word> bufferVec(maxBLen);
        word* buffer = &bufferVec[0];

        RadixLSDCache<unum, word, unum> sorter;

        sortBuckets(bitsPerFirstPass, n, &bStartVec[0], bitsPerChar, input,
                    buffer, sorter);

        if (doThirdPass) {
            sortFlagBuckets(bitsPerFirstPass + bitsPerWord, &bucketFlag[0], 0,
                            length, bitsPerChar, input, buffer, sorter);
        }
    }

  public:
    Radix(uchar* input, unum n, unum kmerLength = 0)
        : originalInput(input), length(n),
          bucketPiggyBackBits(unumBits - bitsFor(length)),
          bucketPiggyBackLimit(((unum)1) << bucketPiggyBackBits),
          bucketPiggyBackMask(bucketPiggyBackLimit - 1),
          kmerLength(kmerLength) {
    }

    unum* build() {

        vector<uint> charCodeV(maxChar);
        uint* charCode = &charCodeV[0];

        int bitsPerChar;
        // int alphaSize =
        indexAlphabet(originalInput, length, charCode, bitsPerChar);
        const int bitsPerFirstPass = 16; //(18 / bitsPerChar) * bitsPerChar;
        //        cout << "Alphabet size " << alphaSize << endl;
        //        << "(" << bitsPerChar << " bits per char) "<<""
        //                "bitsPerFirstPass 18/bitsPerChar*BitsPerChar " << ((18
        //                / bitsPerChar) * bitsPerChar) <<
        //                "; 16/bitsPerChar*BitsPerChar " <<  ((16 /
        //                bitsPerChar) * bitsPerChar) << endl;

        // If we do kmer sorting, an extra radix pass may be enough.
        // This avoids an extra array for bucket numbers, as we need
        // when doing full suffix sorting.
        bool doExtraRadixPass = false;
        bool doFullSuffixSorting = true;
        if (kmerLength > 0) {
            doFullSuffixSorting = false;
            unum requiredBits = kmerLength * bitsPerChar;
            if (requiredBits > bitsPerFirstPass + bitsPerWord) {
                if (requiredBits < bitsPerFirstPass + 2 * bitsPerWord) {
                    doExtraRadixPass = true;
                } else {
                    // extra pass not enough, default to creating suffix array
                    doFullSuffixSorting = true;
                }
            }
        }

        sa = new unum[length];

        inputBasedSort(charCode, bitsPerChar, bitsPerFirstPass,
                       doExtraRadixPass);

        if (doFullSuffixSorting) {
            //            cout << "Doing full suffix sorting" << endl;
            int minBitsSortedInSecondPass =
                doExtraRadixPass ? 2 * bitsPerWord : bitsPerWord;
            int expectedSorted =
                (bitsPerFirstPass + minBitsSortedInSecondPass) / bitsPerChar;

            unum sentinels =
                100; // add sentinels to avoid if's in getBucketWord
            bucket = vector<unum>(length + sentinels);

            finalTouches(expectedSorted);
        }

        return sa;
    }

  private:
    // Places into singleton buckets all suffixes that contain only zeros.
    // Naturally, these suffixes will be at the beginning of the suffix array.
    // This avoids the need for a special character ($) at the end of the
    // string.
    unum fixEndingZeros(uint* charCode) {
        unum p = 0;
        for (unum i = length - 1;
             (p < length) && charCode[originalInput[i]] == 0; --i)
            bucketFlag[++p] = 1;
        return p;
    }

    unum getBucketLength(unum bStart) {
        unum j = bStart + bucketPiggyBackLimit;
        while (!bucketFlag[j])
            ++j;
        return j - bStart;
    }

#ifdef _DEBUG
    void printBucketStats() {
        int n = 32;
        int bb[n];
        for (int i = 0; i < n; ++i)
            bb[i] = 0;

        int b = 0;
        int mb = 0;
        int singletons = 0;
        for (index i = 0; i <= length; ++i) {
            if (bucketFlag[i]) {
                if (b > mb) {
                    mb = b;
                }
                //              /*
                for (int k = 0x10000, j = 0; k > 0; k >>= 1, ++j)
                    if (b > k) {
                        bb[j]++;
                    }
                //               */
                if (b == 1) {
                    singletons++;
                }

                b = 1;
            } else {
                ++b;
            }
        }
        printf("LARGEST bucket %d\n", mb);
        //      /*
        for (int i = 0x10000, j = 0; i > 0; i >>= 1, ++j)
            if (bb[j] > 0)
                printf("buckets larger than %d : %d\n", i, bb[j]);
        //       */
        printf("singletons %d\n", singletons);
    }

    index largestBuck() {
        index b = 0;
        index mb = 0;
        for (index i = 0; i <= length; ++i)
            if (bucketFlag[i]) {
                if (b > mb)
                    mb = b;
                b = 1;
            } else
                ++b;
        return mb;
    }

    void printSuffix(int s, int len) {
        printf(" s=%d len=%d [", s, len);
        printCharData(originalInput + s, len);
        printf("]\n");
    }

    void printSuffixTranslated(int s, int len) {
        printf("Translated s=%d len=%d: ", s, len);
        for (int i = 0; i < len; ++i)
            printf("%3d ", charCode[originalInput[s + i]]);
        printf("\n");
    }

    void printLargestBucket(index expectedSorted) {
        index lb = largestBuck();
        printf("Largest %d\n", lb);
        index bLen = 0;
        int bStart = 0;
        if (expectedSorted > 1000) {
            expectedSorted = 1000;
        }
        for (index i = 0; i <= length; ++i) {
            if (bucketFlag[i]) {
                if (bLen > 1 && lb == bLen) {
                    index* lo = sa + bStart;
                    index s = *lo;
                    printf("Bucket of size %d starts with", bLen);
                    printSuffix(s, expectedSorted);
                    break;
                }
                bLen = 1;
                bStart = i;
            } else {
                ++bLen;
            }
        }
    }

    bool checkIncreasingSuffixes(int bStart, int bLen, index expectedSorted) {
        for (int i = bStart + 1; i < bStart + bLen; ++i) {
            int j;
            for (j = 0; j < expectedSorted; ++j) {
                if (sa[i] + j >= length || sa[i - 1] + j >= length) {
                    break;
                } else {
                    if (originalInput[sa[i] + j] !=
                        originalInput[sa[i - 1] + j])
                        break;
                }
            }
            if (j < expectedSorted &&
                originalInput[sa[i] + j] < originalInput[sa[i - 1] + j]
                //              && originalInput[mySA[i]] !=
                //              originalInput[mySA[i] + 1]
            ) {
                printf("OHO suffix %d before suffix %d but they differ at %d\n",
                       sa[i - 1], sa[i], j);
                printf("OHO suffix %d: ", sa[i - 1]);
                printSuffix(sa[i - 1], expectedSorted);
                printSuffixTranslated(sa[i - 1], expectedSorted);
                printf("OHO suffix %d: ", sa[i]);
                printSuffix(sa[i], expectedSorted);
                printSuffixTranslated(sa[i], expectedSorted);

                return false;
                break;
            }
        }
        //
        return true;
    }

    void checkSortedBucket(index bStart, index bLen, index expectedSorted) {
        index* lo = sa + bStart;
        for (int i = 1; i < bLen; ++i)
            if (memcmp(originalInput + lo[i - 1], originalInput + lo[i],
                       expectedSorted) != 0) {
                printf("These should not be in the same bucket\n");
                printSuffix(lo[i - 1], expectedSorted);
                printSuffix(lo[i], expectedSorted);
            }
    }

#endif
};

#endif /* RADIX_H_ */