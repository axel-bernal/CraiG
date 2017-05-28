/****************************************************************************
 * FeatureVector.h - part of the lless namespace, a general purpose
 *                  linear semi-markov structure prediction library
 *
 *   Copyright (C) 2002-2007  Axel E. Bernal (abernal@seas.upenn.edu)
 *
 *   This program is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU General Public License as
 *   published by the Free Software Foundation; either version 2 of the
 *   License, or (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful, but
 *   WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 *   02111-1307, USA.
 *
 *   The GNU General Public License is contained in the file COPYING.
 *
 *
 ****************************************************************************/

#ifndef _VECTOR_UTILS_H_
#define _VECTOR_UTILS_H_
#include "Utils.h"
#include <map>
#include <list>
#include <vector>
#include <google/dense_hash_map>
#include <google/sparse_hash_map>

#define DENSE_HASH google::dense_hash_map
#define SPARSE_HASH google::sparse_hash_map

namespace lless {

    template <class TClass1, class TClass2>  class TypedVector;

    class VectorUtils {
    public:
        template <class TClass1, class TClass2, class TClass3>
            static void add(TypedVector<TClass1,TClass3> & fv1,
                            const TypedVector<TClass2,TClass3> & fv2) {

            typename TClass2::const_iterator cit = fv2.begin();

            if(fv2.defaultValue()) {
                typename TClass1::iterator it = fv1.begin();

                while(it != fv1.end() && cit != fv2.end()) {
                    ULONG fv1p = fv1.position(it);
                    ULONG fv2p = fv2.position(cit);

                    if(fv1p < fv2p) {
                        fv1[fv1p] += fv2.defaultValue();
                        it++;
                    }
                    else {
                        fv1[fv2p] += fv2.value(cit);
                        cit++;
                        if(fv1p == fv2p)
                            it++;
                    }
                }

                for( ; it != fv1.end(); it++)
                    fv1[fv1.position(it)] += fv2.defaultValue();

                if(fv1.maxSize() > fv1.size())
                    fv1.setDefaultValue(fv1.defaultValue() + fv2.defaultValue());

            }

            for( ; cit != fv2.end(); cit++)
                fv1[fv2.position(cit)] += fv2.value(cit);

        }

        template <class TClass1, class TClass2, class TClass3>
            static void substract(TypedVector<TClass1,TClass3> & fv1,
                                  const TypedVector<TClass2,TClass3> & fv2) {

            typename TClass2::const_iterator cit = fv2.begin();

            if(fv2.defaultValue()) {
                typename TClass1::iterator it = fv1.begin();

                while(it != fv1.end() && cit != fv2.end()) {
                    ULONG fv1p = fv1.position(it);
                    ULONG fv2p = fv2.position(cit);

                    if(fv1p < fv2p) {
                        fv1[fv1p] -= fv2.defaultValue();
                        it++;
                    }
                    else {
                        fv1[fv2p] -= fv2.value(cit);
                        cit++;
                        if(fv1p == fv2p)
                            it++;
                    }
                }

                for( ; it != fv1.end(); it++)
                    fv1[fv1.position(it)] -= fv2.defaultValue();

                if(fv1.maxSize() > fv1.size())
                    fv1.setDefaultValue(fv1.defaultValue() - fv2.defaultValue());

            }

            for( ; cit != fv2.end(); cit++)
                fv1[fv2.position(cit)] -= fv2.value(cit);

        }

        template <class TClass1, class TClass2, class TClass3 >
            static void substractInverse(TypedVector<TClass1,TClass3> & fv1,
                                         const TypedVector<TClass2,TClass3> & fv2) {

            typename TClass2::const_iterator cit = fv2.begin();

            if(fv2.defaultValue()) {
                typename TClass1::iterator it = fv1.begin();

                while(it != fv1.end() && cit != fv2.end()) {
                    ULONG fv1p = fv1.position(it);
                    ULONG fv2p = fv2.position(cit);

                    if(fv1p < fv2p) {
                        fv1[fv1p]  = 1/(1/fv1[fv1p] - fv2.defaultValue());
                        it++;
                    }
                    else {
                        fv1[fv2p]  = 1/(1/fv1[fv2p] - fv2.value(cit));
                        cit++;
                        if(fv1p == fv2p)
                            it++;
                    }
                }

                for( ; it != fv1.end(); it++) {
                    ULONG i = fv1.position(it);
                    fv1[i] = 1/(1/fv1[i] - fv2.defaultValue());
                }

                if(fv1.maxSize() > fv1.size())
                    fv1.setDefaultValue(1/(1/fv1.defaultValue() - fv2.defaultValue()));

            }

            for( ; cit != fv2.end(); cit++) {
                ULONG i = fv2.position(cit);
                fv1[i] = 1/(1/fv1[i] - fv2.value(cit));
            }
        }

        template <class TClass1, class TClass2, class TClass3>
            static void assign(TypedVector<TClass1,TClass3> & fv1,
                               const TypedVector<TClass2,TClass3> & fv2) {

            typename TClass2::const_iterator cit = fv2.begin();
            fv1.reset();

            for( ; cit != fv2.end(); cit++)
                fv1[fv2.position(cit)] = fv2.value(cit);

            if(fv1.maxSize() > fv1.size())
                fv1.setDefaultValue(fv2.defaultValue());

        }

        template <class TClass1, class TClass2, class TClass3>
            static double dot_product(const TypedVector<TClass1,TClass3> & fv1,
                                      const TypedVector<TClass2,TClass3> & fv2) {

            if(fv1.defaultValue() || fv2.defaultValue())
                throw EXCEPTION(BAD_USAGE, "can't * NonZeroSparse TypedVector");

            double result = 0;
            //      cerr << "sizes " << fv.size() << " " << vals.size() << endl;
            if(fv2.size() < fv1.size()) {
                typename TClass2::const_iterator cit = fv2.begin();

                for( ; cit != fv2.end(); cit++)
                    result += fv2.value(cit)*fv1[fv2.position(cit)];
            }
            else {
                typename TClass1::const_iterator cit = fv1.begin();

                for( ; cit != fv1.end(); cit++)
                    result += fv1.value(cit)*fv2[fv1.position(cit)];
                //	  cerr << cit->first << " " << cit->second << " " << result << endl;
            }

            return result;
        }

        template <class TClass1, class TClass2, class TClass3>
            static void average(int n, TypedVector<TClass1,TClass3> & fv1,
                                const TypedVector<TClass2,TClass3> & fv2) {

            if(fv2.defaultValue())
                throw EXCEPTION(BAD_USAGE, "can't average NonZeroSparse TypedVector");

            fv1.reset();
            typename TClass2::const_iterator cit = fv2.begin();

            for( ; cit != fv2.end(); cit++)
                fv1[fv2.position(cit)] += fv2.value(cit)/n;
        }

    };

}

#endif
