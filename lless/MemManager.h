/****************************************************************************
* MemManager.h - part of the lless namespace, a general purpose
*                linear semi-markov structure prediction library
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

#ifndef _MY_MEMMANAGER_H_
#define _MY_MEMMANAGER_H_

#define MEMORY_BLOCKSIZE 0x80000

namespace lless {
  
  //Header file for memory management
  //1. Apply for memory in big chunks
  //2. give out memory to other callers in small amounts
  //3. free all the memory at once
  
  //optmized for speed. system free is very expensive 
  
  //apply for system memory in 512KB chunks
  
  typedef struct _Memory_Linklist {
    unsigned char *ptr;
    struct _Memory_Linklist *next;
  } Memory_Linklist;
  
  class MemManager
    {
    private:
	Memory_Linklist *link, *current;
	unsigned long pos;
    public:
	MemManager();
	~MemManager();
	void *malloc(unsigned long size); //small amount allocate memory
    };

}

#endif
  
  
