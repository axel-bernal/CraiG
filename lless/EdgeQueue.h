/****************************************************************************
 * EdgeQueue.h - part of the lless namespace, a general purpose
 *               linear semi-markov structure prediction library
 *
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
 ****************************************************************************/

#ifndef _EDGE_QUEUE_H_
#define _EDGE_QUEUE_H_

#include "FilterEngine.h"
#include "FSM.h"
#include "TypeDefs.h"

#define QITERATOR int

namespace lless {

    /**
     * EdgeQueue represents a queue which contain pointers to EdgeInst objects
     * and it is internally represented and implemented with a circular buffer.
     */
    class EdgeQueue {
    private:
        EdgeInst **_buffer;
        QITERATOR _head;
        QITERATOR _tail;
        int _size;
        int _tailPos;
        int _headPos;

    public:
        /**
         * Default destructor
         * @param size the size of the circular buffer which is the internal
         * representation of the queue.
         */
        EdgeQueue(int size) {
            _size = size;
            _buffer = new EdgeInst * [_size];
            _head = 0;
            _tail = 0;
            _tailPos = 0;
            _headPos = 0;
        }

        inline int size() {
            return (_head >= _tail) ? _head - _tail + 1 : _size - (_tail - _head) + 1;
        }

        inline bool empty() {
            return _head == _tail;
        }

        inline bool full() {
            return (_tail + 1) % _size == _head;
        }

        inline void setTailPos(int pos) {
            _tailPos = pos;
        }

        inline void setHeadPos(int pos) {
            _headPos = pos;
        }

        inline int tailPos() {
            return _tailPos;
        }

        inline int headPos() {
            return _headPos;
        }

        inline void makeEmpty() {
            _head = _tail;
        }

        inline EdgeInst *peek() {
            if(empty()) return NULL;
            return _buffer[_head];
        }

        inline EdgeInst *peekTail() {
            if(empty()) return NULL;
            return _buffer[(_tail - 1 >= 0) ? _tail - 1: _size - 1];
        }

        inline EdgeInst *dequeue() {
            EdgeInst *sig = peek();
            if(sig)
                _head = (_head + 1) % _size;
            return sig;
        }

        inline void queue(EdgeInst *sig) {
            //    cerr << "before " << _head << " " << _tail << " " << _size << endl;
            if(full()) {
                expand();
            }
            //    cerr << "after " << _head << " " << _tail << " " << _size << endl;
            _buffer[_tail] = sig;
            _tail = (_tail + 1) % _size;
        }

        inline QITERATOR head() {
            return _head;
        }

        inline void restoreHead(QITERATOR head) {
            _head = head;
        }

        inline void expand() {
            int i = 0;
            EdgeInst **_tmp = new EdgeInst * [_size*2];
            // copying over existing elements
            int limit = Utils::max(_tail, _size);
            for(i = _head; i < limit; i++)
                _tmp[i - _head] = _buffer[i];
            if(_tail < _head) // it wrapped around
                for(i = 0; i < _tail; i++)
                    _tmp[limit + i - _head] = _buffer[i];
            delete [] _buffer;
            _buffer = _tmp;
            _head = 0;
            _tail = limit + _tail - _head;
            _size = _size*2;
        }

        ~EdgeQueue() {
            delete [] _buffer;
        }
    };

}

#endif
