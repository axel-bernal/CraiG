/****************************************************************************
* ArgparseUtils.h - part of the craig namespace, a genomics library
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

#ifndef _ARGPARSE_UTILS_H_
#define _ARGPARSE_UTILS_H_
#include "FeatureVector.h"
#include "Utils.h"
#include <numeric>

using namespace lless;

#define CRAIG_VERSION "1.1"
#define AUTHOR_EMAIL "abernal@seas.upenn.edu"
#define NOTICE "Copyright (C) 2003-2007 Axel E. Bernal (abernal@seas.upenn.edu)"
#define PRINT_VERSION(os, pname, desc) (os) << (pname) << " v. " << CRAIG_VERSION << ". " << (desc) << ".\nCopyright (C) 2007 Axel E. Bernal (abernal@seas.upenn.edu)\n"
#define PRINT_DISCLAIMER(os, pname) (os) << (pname) << " comes with NO WARRANTY,\n" << "to the extent permitted by law.\n" << "You may redistribute copies of " << (pname) << " under the\n" << "terms of the GNU General Public License version 2.\n" << "For more information about these matters,\n" << "see the files named COPYING.\n\n";

namespace craig {

  class ArgParseUtils {
       
   public:
     ArgParseUtils() {

     }
     
     static void displayOption(ofstream &fd, string option, string description, int num_colsopt=24, int num_colspage=80) {
       fd << "  " << option;
       int lensofar = option.length() + 2;
       if(lensofar + 1 > num_colsopt) {
	 fd << endl;
	 lensofar = 0;
       }

       int i = lensofar;
       for( ; i < num_colsopt; i++)
	 fd << " ";

       vector<string> tokens;
       Utils::stringSplit<string>(tokens, description, " ");
       for(int j = 0; j < tokens.size(); j++) {
	 if(i + tokens[j].length() > num_colspage) {
	   fd << endl;
	   for(i = 0; i < num_colsopt; i++)
	     fd << " ";
	 }

	 fd << tokens[j];
	 if(j == tokens.size() - 1)  fd << "\n";
	 else  fd << " ";
	 i += tokens[j].length() + 1;
       }
     }
  };
}
#endif
