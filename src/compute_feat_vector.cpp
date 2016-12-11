/*****************************************************************************
 * This is the main of the compute_feat_vector program. It computes the 
 * feature vector associated to both training and validation tag annotation
 * present int he configuration file.
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


#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>
#include "Utils.h"
#include "ContextIMM.h"
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include "Sequence.h"
#include "SequenceUtils.h"
#include "ResourceEngine.h"
#include "Configuration.h"
#include "FeatureEngine.h"
#include "StructureCore.h"
#include "GeneUtils.h"
#include "Lattice.h"
#include "GeneEvaluator.h"
#include "GeneTagPrinter.h"
#include "FSM.h"
#include "TagUtils.h"
#include "InpFile.h"

#define NUM_PHASES 3

using namespace craig;

void printHelp(const char * pname, ::ofstream &fd) {
  fd << "CRAIG v. " << CRAIG_VERSION <<  "tool to compute the feature vector associated to both\ntraining and validation tag annotation present in the\nconfiguration file.\n";
  fd << "Written by Axel E. Bernal (" << AUTHOR_EMAIL << ")\n\n";
  fd << "  usage :" << pname << " [options] CONF_FILE > FEATVECTOR_FILE\n\n";
  ArgParseUtils::displayOption(fd, "CONF_FILE", "Name of the file containing setup information needed to perform learning. See README for details");
  ArgParseUtils::displayOption(fd, "PARAMS_FILE", "Name of the file where model parameters will be stored");
  fd << "optional arguments:\n";
  ArgParseUtils::displayOption(fd, "--version", "Print version name and license information");
  ArgParseUtils::displayOption(fd, "-h --help", "Shows this help message and exit");
  ArgParseUtils::displayOption(fd, "-v --verbose", "Turns on output of debugging information");
  ArgParseUtils::displayOption(fd, "--resources=RES_FILE", "Resources file RES_FILE will be used instead of default one in model directory");
  ArgParseUtils::displayOption(fd, "--filters=FILT_FILE", "Filters file FILT_FILE, will be used instead of default one in model directory");
  ArgParseUtils::displayOption(fd, "--features=FEAT_FILE", "Features file FEAT_FILE, will be used instead of default one in model directory");
  ArgParseUtils::displayOption(fd, "--prefix-evidence=PREFIX_EVIDENCE", "prefix for fetching files used as input resources");
  ArgParseUtils::displayOption(fd, "--strand=STRAND", "One of forward|backward|both. Specifies the strand for decoding/predicting the input sequence [both]");
  ArgParseUtils::displayOption(fd, "FEATVECTOR_FILE", "Feature vectors for training and validation sets.");
  fd << "Report bugs to <" << AUTHOR_EMAIL << ">\n";
}


bool verbose = false;

int main(int argc, char *argv[]) {
  ::string modelPath, confFile, paramFile, prefix_evidence("");
  ::string resourcesFile = "", filtersFile = "", featuresFile = "";
  float minPercCov = 99.9;
  double learnRate = 1;
  int numIter = 100, origStrand = BOTH_STRANDS;
  int kBest = 1, i;
  std::string topology("partial"), strand("both");
  TTrainMethod tm = MIRA;
  bool lossAugmented = false;
  std::string logFile = "";

  try {
    if(!getenv("CRAIG_HOME"))
      throw EXCEPTION(NOT_ANNOTATED,
                          "CRAIG_HOME must be initialized!. See README for details");
    
    modelPath = std::string(getenv("CRAIG_HOME")) + "/models/";

    for(i = 1; i < argc; i++) {
      if(!strncmp(argv[i], "--strand=", 9)) {
        strand = std::string(argv[i] + 9);
      }
      else if(!strncmp(argv[i], "--resources=", 12)) {
	resourcesFile = ::string(argv[i] + 12);
      }
      else if(!strncmp(argv[i], "--filters=", 10)) {
	filtersFile = ::string(argv[i] + 10);
      }
      else if(!strncmp(argv[i], "--features=", 11)) {
	featuresFile = ::string(argv[i] + 11);
      }
      else if(!strncmp(argv[i], "--verbose", 9) 
              || !strncmp(argv[i], "-v", 2)) {  
        verbose = true;  
      }
      else if(!strncmp(argv[i], "--prefix-evidence=", 18)) {
	prefix_evidence = string(argv[i] + 18);
      }
      else if(!strncmp(argv[i], "--help", 6) 
              || !strncmp(argv[i], "-h", 2)) {
        printHelp("compute_feat_vector", (std::ofstream &)cout); 
        exit(0);  
      }
      else if(!strncmp(argv[i], "--version", 9)) { 
        PRINT_VERSION(cerr, "compute_feat_vector", "tool for computing feature vector from tag annotations");
        PRINT_DISCLAIMER(cerr, "compute_feat_vector"); 
        exit(0);  
      }
      else if(!strncmp(argv[i], "--algorithm=", 12)) {
        if(!strncmp(argv[i] + 12, "PERCEPTRON", 5))
          tm = PERCEPTRON;
        else if(!strncmp(argv[i] + 12, "PEGASOS", 7))
          tm = PEGASOS;
        else if(strncmp(argv[i] + 12, "MIRA", 4))
          throw EXCEPTION(BAD_USAGE,
                              std::string("unrecognized c.l. argument") + string(argv[i]));
        
      }
      else
        break;

    }

    if(strand.compare("forward") == 0)
      origStrand = STRAND_FWD;
    else if(strand.compare("backward") == 0)
      origStrand = STRAND_COMP;
    else if(strand.compare("both") != 0)
      throw EXCEPTION(BAD_USAGE, std::string("unrecognized c.l. argument") + strand);

    if(argc - i < 1)
      throw EXCEPTION(BAD_USAGE, "insufficient arguments");
    
    confFile = std::string(argv[i++]);
    
    ::ifstream confStream(confFile.c_str());
    
    if(!confStream.is_open())
      throw EXCEPTION(FILE_UNAVAILABLE, confFile);

    Configuration config(confStream);
    modelPath = config.name().find('/', 0) == std::string::npos ?
      modelPath + config.name() + "." :
      config.name() + ".";

    if(!resourcesFile.length())
      resourcesFile = modelPath + "resources";
    if(!filtersFile.length())
      filtersFile = modelPath + "filters";
    if(!featuresFile.length())
      featuresFile = modelPath + "features";

    std::ifstream reStream(resourcesFile.c_str());
    std::ifstream flStream(filtersFile.c_str());
    std::ifstream ftStream(featuresFile.c_str());
    
    if(!reStream.is_open() || !flStream.is_open() || !ftStream.is_open())
      throw EXCEPTION(FILE_UNAVAILABLE, modelPath + "*");

    map<std::string, std::string> resource_subs;
    resource_subs["PREFIX_EVIDENCE"] = prefix_evidence;
    resource_subs["~"] = modelPath;
    ResourceEngine re(reStream, &resource_subs);

    Sigma *sigma = (Sigma *)re.getResource(std::string("dna-alpha"));
    FSM & fsm = *(FSM *)re.getResource(topology);
    fsm.setParseStrand((TStrand)origStrand);
    FilterEngine fe(re, flStream);
    FeatureEngine fte(fsm, fe, ftStream);

    /* 
     * loading annotSeqs and gene annotations.
     * It's assumed at training time that each Sequence object contains only
     * one gene.
     * One needs to preprocess the input data to fulfil this requirement. 
     * The perl utility splitGenes.pl is provided along the C++ sources in
     * the craig distribution for this purpose.
     */
    
    InpFile trainFile("default", config.trainingSequences(), FASTA, sigma, true);
    list<Sequence *> & trainSequences = (list<Sequence *> &)trainFile.sequences();

    InpFile validFile("default", config.validationSequences(), FASTA, sigma, true);
    list<Sequence *> & validSequences = (list<Sequence *> &)validFile.sequences();
    
    TagUtils::loadSeqTags(config.trainingSet(), trainSequences, fsm);
    TagUtils::loadSeqTags(config.validationSet(), validSequences, fsm);
    
    GeneEvaluator evaluator(fsm);
    GeneTagPrinter printer;

    // Initializing special variables
    TypedFilter<EdgeInst> **signals = new TypedFilter<EdgeInst> * [NUM_EDGE_INSTS];
    signals[START] = (TypedFilter<EdgeInst> *)fe.getFilter("Start-Signal");
    signals[STOP] = (TypedFilter<EdgeInst> *)fe.getFilter("Stop-Signal");
    signals[DONOR] = (TypedFilter<EdgeInst> *)fe.getFilter("Donor-Signal");
    signals[ACCEPTOR] = (TypedFilter<EdgeInst> *)fe.getFilter("Acceptor-Signal");

    TypedFilter<UCHAR> *contexts = (TypedFilter<UCHAR> *)fe.findFilter("GC-Content");
    FT_HiddenSeq *vftHidden = (FT_HiddenSeq *)fte.findFeature("Hidden-Stop-Sequence");

    //Creating Core object
    Lattice lattice(fe, fsm, NUM_PHASES, evaluator, MAX_NUM_EXONS, 
                    INTERGENIC, signals, contexts, vftHidden);
    // DONT FORGET
    GlobalVector params(fte.getFeatures()); // final regularized parameters

    StructureCore core(fsm, lattice, re, fe, fte, 
                       evaluator, printer, AVG_ALL,
                       COMB_EUCLID, (TStrand)origStrand,
                       1, false);


    list<Sequence *>::iterator cit;
    cerr << "Computing parameter vector for organism " << config.name() << endl;

    fte.setParamVector(&params);
    cout << "Training Set\n";
    for(cit = trainSequences.begin(); cit != trainSequences.end(); cit++) {
      Sequence &c = *(*cit);    
      //      cerr << c.id() << endl;
      core.computeGlobalVector(c, c.getTags());
      //      params.print(cout);
    }

    params.store((ofstream &)cout);
    params = 0.0;

    //    cout << "Evaluation Set\n";
    for(cit = validSequences.begin(); cit != validSequences.end(); cit++) {
      Sequence &c = *(*cit);    
      core.computeGlobalVector(c, c.getTags());
    }

    //    params.print(cout);
    //    cerr << "Done!\n";

    reStream.close();
    flStream.close();
    ftStream.close();
    delete [] signals;
    //DONT FORGET
    //    pStream.close();

  } catch(exception *e) { 
    perror(e->what());

    GenLibExcept *gle = (GenLibExcept *)e;
    if(gle->error() == BAD_USAGE)
      cerr << "use --help for more information\n";
    delete e;
    exit(1);
    
  }
}
