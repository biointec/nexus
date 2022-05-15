/*
 * MergeSorter.h
 *
 *  Created on: Jul 19, 2012
 *      Author: marius
 */

#ifndef MERGESORTER_H_
#define MERGESORTER_H_
#include "utils.h"
#include <vector>
using namespace std;

template <class DataT, class KeyT, class S> class MergeSorter {
  private:
    vector<KeyT> destK;
    vector<DataT> destD;

  public:
    MergeSorter(S maxLength) {
        destD.reserve(maxLength);
        destK.reserve(maxLength);
    }

    ~MergeSorter() {
    }

    void sort(S length, DataT* data, KeyT* key) {
        // make sure we finish back in the original array to avoid an extra copy
        S insSortThr = (logceil(length) & 1) ? 32 : 16;

        S s = 0;
        for (; s + insSortThr < length; s += insSortThr)
            insertSort(data + s, insSortThr, key + s);
        if (length - s > 1)
            insertSort(data + s, length - s, key + s);

        iterativeMergeSort(length, data, key, insSortThr);
    }

    void iterativeMergeSort(S length, DataT* data, KeyT* key,
                            S initialTileSize) {
        // iterative merge sort
        DataT* tD = &destD[0];
        KeyT* tK = &destK[0];
        for (S bSize = initialTileSize; bSize < length; bSize <<= 1) {
            int n = 0;
            for (S s = 0; s < length; s += (bSize << 1)) {
                S i = s;
                S m = s + bSize;
                if (m > length)
                    m = length;
                S j = m;
                S e = m + bSize;
                if (e > length)
                    e = length;
                for (; i < m && j < e; ++n) {
                    if (key[i] <= key[j]) {
                        tK[n] = key[i];
                        tD[n] = data[i];
                        ++i;
                    } else {
                        tK[n] = key[j];
                        tD[n] = data[j];
                        ++j;
                    }
                }
                for (; i < m; ++n, ++i) {
                    tK[n] = key[i];
                    tD[n] = data[i];
                }
                for (; j < e; ++n, ++j) {
                    tK[n] = key[j];
                    tD[n] = data[j];
                }
            }

            DataT* tmp = data;
            data = tD;
            tD = tmp;

            KeyT* tmpk = key;
            key = tK;
            tK = tmpk;
        }
    }
};

#endif /* MERGESORTER_H_ */