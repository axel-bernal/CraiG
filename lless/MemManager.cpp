#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "MemManager.h"

/****************************************************************************
 * MemManager.cpp - part of the lless namespace, a general purpose
 *                    linear semi-markov structure prediction library
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

namespace lless {

    //Constructor
    MemManager::MemManager()
    {
        pos=0;
        link=(Memory_Linklist *) ::malloc(sizeof(Memory_Linklist));
        assert(link);
        current=link;
        current->ptr=(unsigned char *) ::malloc(MEMORY_BLOCKSIZE);
        assert(current->ptr);
        current->next=NULL;
        return;
    }

    //Give out small amounts of memory each time
    void *MemManager::malloc(unsigned long size)
    {
        void *ret;

        // Ask for too much, return empty pointer
        if(size>MEMORY_BLOCKSIZE) return NULL;

        if(size<=MEMORY_BLOCKSIZE-pos) //Sufficient space in the current block
	{
            ret=(void *)(current->ptr+pos);
            pos+=size;
	}
        else //Apply for a new block
	{
            current->next=(Memory_Linklist *) ::malloc(sizeof(Memory_Linklist));
            assert(current->next);
            current=current->next;
            current->ptr=(unsigned char *) ::malloc(MEMORY_BLOCKSIZE);
            assert(current->ptr);
            current->next=NULL;
            ret=(void *)current->ptr;
            pos=size;
	}
        return ret;
    }

    //Free everything together
    MemManager::~MemManager()
    {
        Memory_Linklist *tcurrent, *temp;

        tcurrent=link;
        while(tcurrent)
	{
            temp=tcurrent;
            tcurrent=tcurrent->next;
            free(temp->ptr);
            free(temp);
	}
        return;
    }

}
