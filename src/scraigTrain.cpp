/******************************************************************************
 * This is training part of CRAIG. It trains parameters for a linear 
 * structure model for predicting genes. The method uses MIRA online training
 * algorithm but PERCEPTRON is also available. The  input is the training and 
 * testing sets of gene coordinates with respective fasta sequences. The 
 * output is the gene model
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
#include "SelfLattice.h"
#include "Lattice.h"
#include "GeneEvaluator.h"
#include "GeneTagPrinter.h"
#include "FSM.h"
#include "TagUtils.h"
#include "InpFile.h"
#include "ArgParseUtils.h"

#define NUM_PHASES 3

using namespace craig;

bool verbose = false;

    

void printHelp(const char * pname, ::ofstream &fd) {
  fd << "CRAIG v. " << CRAIG_VERSION << " tool for ";
#ifdef HAVE_MPI 
  fd << "parallel discriminative training of gene taggers using OpenMPI (TM). All file system input is to be expected to be directories to be accessed easily from every created process\n";
#else
  fd << "discriminative training of gene taggers\n";
#endif
  fd << "Written by Axel E. Bernal (" << AUTHOR_EMAIL << ")\n\n";
  fd << "Usage :" << pname << " [Options] CONF_FILE PARAMS_FILE > DEBUG_INFO\n\n";
  fd << "positional arguments:\n";
  ArgParseUtils::displayOption(fd, "CONF_FILE", "Name of the file containing setup information needed to perform learning. See README for details");
  ArgParseUtils::displayOption(fd, "PARAMS_FILE", "Name of the file where model parameters will be stored");
  fd << "optional arguments:\n";
  ArgParseUtils::displayOption(fd, "--version", "Print version name and license information");
  ArgParseUtils::displayOption(fd, "-h --help", "Shows this help message and exit");
  ArgParseUtils::displayOption(fd, "-v --verbose", "Turns on output of debugging information");
  ArgParseUtils::displayOption(fd, "--strand=STRAND", "One of forward|backward|both. Specifies the strand for decoding/predicting the input sequence [both]");
  ArgParseUtils::displayOption(fd, "-utr --model-utrs", "Model will be trained using UTR regions If not specified, make sure tags used as training do not contain UTR-labeled regions");
  ArgParseUtils::displayOption(fd, "--best=K", "Trains with the top K-best segmentations [1]");
  ArgParseUtils::displayOption(fd, "--prefix-evidence=PREFIX_EVIDENCE", "prefix for fetching files used as input resources");

  fd << "loss-function arguments:\n";
  ArgParseUtils::displayOption(fd, "-xl --max-loss", "Turns on max-loss decoding during training");
  ArgParseUtils::displayOption(fd, "--edge-loss-function=LOSS", "Specifies the type of edge loss function. One of NONE, SEGMENT, EDGE, SOFT_EDGE. [NONE]");
  ArgParseUtils::displayOption(fd, "--loss-function=LOSS", "Specifies the type of loss function. One of HAMMING, NONE, CORR_COEF, ZERO_ONE, or F_SCORE. [HAMMING]");
  ArgParseUtils::displayOption(fd, "--edge-loss-factor=FACTOR", "Specifies the amoung of loss caused by onet edge or segment disagreement. [10 for edges, 20 for segments]");
  ArgParseUtils::displayOption(fd, "--add-start-loss=PERCENT", "Amount of additional loss added in percentage over the current loss. This works as a trade-off between Sp/Sn[0.0]");
  ArgParseUtils::displayOption(fd, "--noise-level=NOISE_LEVEL", "Affects the way the loss is computed. Greater than zero if there are false negative instances in the annotation. Max value is 10 [0].");
  fd << "hot restart arguments:\n";
  ArgParseUtils::displayOption(fd, "--restart=LOG_FILE", "Restart training from log-file[\"\"]");
  ArgParseUtils::displayOption(fd, "--added-features=N" ,"Works along the --restart option Useful for training subsets with N new features defined at the end of the feature file which are set with 0-valued parameters on new training data[0]");
  fd << "learning algorithm arguments:\n";
  ArgParseUtils::displayOption(fd, "--learn-rate=LEARN_RATE", "Learning rate for the parameter updating [1.0]"); 
  ArgParseUtils::displayOption(fd, "--algorithm=ALGORITHM", "Training algorithm to be used. One of MIRA, PERCEPTRON, PEGASOS, CWL or ARROW [MIRA]");
  ArgParseUtils::displayOption(fd, "--avg-method=AVG_METHOD", "One of none|last|all Average parameters using all iterations or only the last one[all]");
  ArgParseUtils::displayOption(fd, "--phi=PHI", "Value to control the margin norm for CWL [0.1]");
  ArgParseUtils::displayOption(fd, "--r=R", "Value to control the margin norm for ARROW [1.0]");
  fd << "parameter update arguments:\n";
  ArgParseUtils::displayOption(fd, "--limit=LIMIT", "Number of update passes to perform [100]");
  ArgParseUtils::displayOption(fd, "--oracle-upd=ORACLE_UPDATE", "Update with oracle hypothesis instead of training instance in case the reference annotation is unreachable. The oracle is selected among the top k-best segmentations using either the TOP segmentation, ALL of them, all those BETTER than the oracle or NONE of them [NONE]");
  ArgParseUtils::displayOption(fd, "--oracle-wait=ORACLE_WAIT", "Waits <number> iterations and then start using oracle-upd for the updates[undef]");
  ArgParseUtils::displayOption(fd, "--comb-method=COMB_METHOD", "One of <euclid|k-l>. Specifies the selected combining method for feature vectors provenient from different iterations. Either one which minimizes the euclidean distance or the k-l divergence between them.[euclid]");
  ArgParseUtils::displayOption(fd, "--multi-label-upd=MULTI_LABEL_UPD", "One of ALL|MIN|EXP|SEP|LONGEST|MINLONGER. Controls the updates for multi-labeled sequences. MIN only updates weights for the min-loss expected labeling in. EXP updates weights for a randomly picked labeling sampled from a exponential distribution over labeling loss differences. SEP considers each labeling as a separate instance LONGER updates against the longest labeling, where for labels different from sync state MINLONGER updates like LONGEST for the first iteration(s), then progressively includes additional longer labelings one at a time, and switchs to MIN updates until all labelings are in [ALL]");
  ArgParseUtils::displayOption(fd, "--min-longer-wait=MIN_LONGER_WAIT", "Number of iterations that must pass before adding the next longest labeling not already included in the set of expected labelings for multi-label learning [1]");
  fd << "parallel algorithm arguments:\n";
  ArgParseUtils::displayOption(fd, "--parallel-comb-method=PAR_COMB_METHOD", "One of <euclid|k-l>. Specifies the selected combining method for feature vectors provenient from different nodes in a parallel training procedure. Either one that minimizes the euclidean distance or the k-l divergence between them.[k-l]");
  ArgParseUtils::displayOption(fd, "--block-size=BLOCK_SIZE", "The number of instances sent as a block to each node in parallel when training in a multi-processor architecture[1]");
  ArgParseUtils::displayOption(fd, "-s --share-models", "Share models instead of updates between nodes when training in a multi-processor architecture. Block size in this case is defaulted to 1");
  ArgParseUtils::displayOption(fd, "-bl --balance-load", "Balances the load of each processor so that synchronization barriers do not become processing bottlenecks");
  fd << "input preprocesing arguments:\n";
  ArgParseUtils::displayOption(fd, "-u --add-unreachable", "Learn model including unreachable genes");
  ArgParseUtils::displayOption(fd, "--sample=PERCENT", "Randomly samples a percent of the total number of available training instances and uses it as a new, smaller training set");
  ArgParseUtils::displayOption(fd, "DEBUG_INFO", "Debug information, useful to know when to stop training");
  fd << "Report bugs to <" << AUTHOR_EMAIL << ">\n";


}

int main(int argc, char *argv[]) {
  ::string modelPath, confFile, paramFile, prefix_evidence("");
  double learnRate = 1, phi = 0.1, r= 1.0, addStartLoss = 0.0;
  bool balanceLoad = false;
  int addedFeatures = 0;
  int maxIterations = 100,
    origStrand = BOTH_STRANDS, oracleWait = -1,
    minLongerWait = 1;
  int kBest = 1, i, noiseLevel = 0, blockSize = 1;
  vector<TNodeId3> id34LenSorting;
  id34LenSorting.push_back(EXON); id34LenSorting.push_back(INTRON);
  double sampleSizePerc = 0;
  std::string topology("partial"), strand("both");
  bool addUnreachable = false, modelUTRs = false;
  TTrainMethod tm = MIRA;
  TLossType lf = LF_HAMMING, elf = LF_NONE;
  double edgeLossFactor = 0;
  TCombMethod combMethod = COMB_EUCLID, parCombMethod = COMB_KL;
  TAvgMethod avgMethod = AVG_ALL;
  bool maxLoss = false, shareModels = false;
  TMultiUpd multiUpd = ML_ALL;
  TOracleUpd oracleUpd = OC_NONE;
  std::string logFile = "";

  try {
    if(!getenv("CRAIG_HOME"))
      throw EXCEPTION(NOT_ANNOTATED,
                          "CRAIG_HOME must be initialized!. See README for details");
    
    modelPath = std::string(getenv("CRAIG_HOME")) + "/models/";

    for(i = 1; i < argc; i++) {
      if(!strncmp(argv[i], "--learn-rate=", 13)) {
        sscanf(argv[i] + 13, "%lf", &learnRate);
      }
      else if(!strncmp(argv[i], "--limit=", 8)) {
        sscanf(argv[i] + 8, "%d", &maxIterations);
      }
      else if(!strncmp(argv[i], "--add-unreachable", 17) 
	      || !strncmp(argv[i], "-u", 2)) {  
	addUnreachable = true;
      }
      else if(!strncmp(argv[i], "--oracle-wait=", 14)) {
        sscanf(argv[i] + 14, "%d", &oracleWait);
      }
      else if(!strncmp(argv[i], "--verbose", 9) 
              || !strncmp(argv[i], "-v", 2)) {  
        verbose = true;  
      }
      else if(!strncmp(argv[i], "--max-loss", 10) 
              || !strncmp(argv[i], "-xl", 3)) {
        maxLoss = true;
      }
      else if(!strncmp(argv[i], "--prefix-evidence=", 18)) {
	prefix_evidence = string(argv[i] + 18);
      }
      else if(!strncmp(argv[i], "--multi-label-upd=", 18)) {
        std::string multiUpdString = std::string(argv[i] + 18);
        multiUpd= TypeDefs::stringToTMultiUpd(multiUpdString);
      }
      else if(!strncmp(argv[i], "--min-longer-wait=", 18)) {
        sscanf(argv[i] + 18, "%d", &minLongerWait);
      }
      else if(!strncmp(argv[i], "--noise-level=", 14)) {
        sscanf(argv[i] + 14, "%d", &noiseLevel);
      }
      else if(!strncmp(argv[i], "--block-size=", 13)) {
        sscanf(argv[i] + 13, "%d", &blockSize);
      }
      else if(!strncmp(argv[i], "--added-features=", 17)) {
        sscanf(argv[i] + 17, "%d", &addedFeatures);
      }
      else if(!strncmp(argv[i], "--share-models", 14)
              || !strncmp(argv[i], "-s", 2)) {
        shareModels = true;
      }
      else if(!strncmp(argv[i], "--balance-load", 14)
              || !strncmp(argv[i], "-bl", 3)) {
        balanceLoad = true;
      }
      else if(!strncmp(argv[i], "--comb-method=", 14)) {
        std::string combMethodString = std::string(argv[i] + 14);
        combMethod = TypeDefs::stringToTCombMethod(combMethodString);
      }
      else if(!strncmp(argv[i], "--parallel-comb-method=", 23)) {
        std::string parCombMethodString = std::string(argv[i] + 23);
        parCombMethod = TypeDefs::stringToTCombMethod(parCombMethodString);
      }
      else if(!strncmp(argv[i], "--oracle-upd=", 13)) {
        std::string oracleUpdString = std::string(argv[i] + 13);
        oracleUpd = TypeDefs::stringToTOracleUpd(oracleUpdString);
      }
      else if(!strncmp(argv[i], "--help", 6) 
              || !strncmp(argv[i], "-h", 2)) {
        printHelp("strain_craig", (std::ofstream &)cout); 
        exit(0);  
      }
      else if(!strncmp(argv[i], "--version", 9)) { 
        PRINT_VERSION(cerr, "train_craig", "tool for discriminative training of gene taggers for eukarya");
        PRINT_DISCLAIMER(cerr, "train_craig"); 
        exit(0);  
      }
      else if(!strncmp(argv[i], "--strand=", 9)) {
        strand = std::string(argv[i] + 9);
      }
      else if(!strncmp(argv[i], "--restart=", 10)) {
        logFile = std::string(argv[i] + 10);
      }
      else if(!strncmp(argv[i], "--best=", 7)) {
        sscanf(argv[i] + 7, "%d", &kBest);
      }
      else if(!strncmp(argv[i], "--loss-function=", 16)) {
        std::string lossString = std::string(argv[i] + 16);
        lf = TypeDefs::stringToTLossType(lossString);
      }
      else if(!strncmp(argv[i], "-utr", 4) ||
	      !strncmp(argv[i], "--model-utrs", 12)) {
	modelUTRs = true;
      }
      else if(!strncmp(argv[i], "--edge-loss-function=", 21)) {
        std::string lossString = std::string(argv[i] + 21);
        elf = TypeDefs::stringToTLossType(lossString);
      }
      else if(!strncmp(argv[i], "--edge-loss-factor=", 19)) {
	sscanf(argv[i] + 19, "%lf", &edgeLossFactor);
      }
      else if(!strncmp(argv[i], "--add-start-loss=", 17)) {
        sscanf(argv[i] + 17, "%lf", &addStartLoss);
      }      
      else if(!strncmp(argv[i], "--avg-method=", 13)) {
        std::string avgString = std::string(argv[i] + 13);
        avgMethod = TypeDefs::stringToTAvgMethod(avgString);
      }
      else if(!strncmp(argv[i], "--algorithm=", 12)) {
        std::string tmString = std::string(argv[i] + 12);
        tm = TypeDefs::stringToTTrainMethod(tmString);
      }
      else if(!strncmp(argv[i], "--phi=", 6)) {
        sscanf(argv[i] + 6, "%lf", &phi);
      }
      else if(!strncmp(argv[i], "--sample=", 9)) {
        sscanf(argv[i] + 9, "%lf", &sampleSizePerc);
      }
      else if(!strncmp(argv[i], "--r=", 4)) {
        sscanf(argv[i] + 4, "%lf", &r);
      }
      else
        break;

    }

    if(argc - i < 2)
      throw EXCEPTION(BAD_USAGE, "insufficient arguments");

    // checking sanity
    if(combMethod == COMB_KL && tm != ARROW && tm != CWL)
      throw EXCEPTION(NOT_SUPPORTED, "K-L average only available for ARROW and CWL");

    if((lf != LF_HAMMING || lf != LF_SEGMENT) && maxLoss) 
      throw EXCEPTION(NOT_SUPPORTED, "Cannot combine max-loss decoding and given loss function");

    if(shareModels)
      blockSize = 1;

    if(strand.compare("forward") == 0)
      origStrand = STRAND_FWD;
    else if(strand.compare("backward") == 0)
      origStrand = STRAND_COMP;
    else if(strand.compare("both") != 0)
      throw EXCEPTION(BAD_USAGE, std::string("unrecognized c.l. argument") + strand);
    
    confFile = std::string(argv[i++]);
    paramFile = std::string(argv[i++]);
    
    cerr << "Training model for organism " << confFile << endl;   
    cerr << "Options :";
    for(i = 1; i < argc; i++)
      cerr << " " <<  argv[i];
    cerr << endl;

    cerr << "Process Id " << getpid() << endl;

    ::ifstream confStream(confFile.c_str());
    
    if(!confStream.is_open())
      throw EXCEPTION(FILE_UNAVAILABLE, confFile);

    ::ofstream paramStream(paramFile.c_str());

    if(!paramStream.is_open()) 
      throw EXCEPTION(FILE_UNAVAILABLE, paramFile);
    
    Configuration config(confStream);
    modelPath = config.name().find('/', 0) == std::string::npos ?
      modelPath + config.name() + "." :
      config.name() + ".";

    std::ifstream reStream((modelPath + "resources").c_str());
    std::ifstream flStream((modelPath + "filters").c_str());
    std::ifstream ftStream((modelPath + "features").c_str());
    
    if(!reStream.is_open() || !flStream.is_open() || !ftStream.is_open())
      throw EXCEPTION(FILE_UNAVAILABLE, modelPath + "*");

    map<std::string, std::string> resource_subs;
    resource_subs["~"] = modelPath;
    resource_subs["PREFIX_EVIDENCE"] = prefix_evidence;
    ResourceEngine re(reStream, &resource_subs);

    Sigma *sigma = (Sigma *)re.getResource(std::string("dna-alpha"));
    FSM & fsm = *(FSM *)re.getResource(topology);
    fsm.setParseStrand((TStrand)origStrand);

    if(!modelUTRs) {
      fsm.removeId2Node(ANY_5UTR);
      fsm.removeId2Node(ANY_3UTR);
      fsm.removeId2Node(ANY_UTR_INTRON);
    }

    FilterEngine fe(re, flStream);
    FeatureEngine fte(fsm, fe, ftStream);

    GeneEvaluator evaluator(fsm, lf, elf, noiseLevel, addStartLoss);
    GeneTagPrinter printer;
    
    if(edgeLossFactor)
      evaluator.setEdgeLossFactor(edgeLossFactor);

    // Initializing special variables
    TypedFilter<EdgeInst> **signals = new TypedFilter<EdgeInst> * [NUM_EDGE_INSTS];
    signals[START] = (TypedFilter<EdgeInst> *)fe.getFilter("Start-Signal");
    signals[STOP] = (TypedFilter<EdgeInst> *)fe.getFilter("Stop-Signal");
    signals[DONOR] = (TypedFilter<EdgeInst> *)fe.getFilter("Donor-Signal");
    signals[ACCEPTOR] = (TypedFilter<EdgeInst> *)fe.getFilter("Acceptor-Signal");

    TypedFilter<UCHAR> *contexts = (TypedFilter<UCHAR> *)fe.findFilter("GC-Content");
    //Creating StructureCore object
    SelfLattice lattice(fe, fsm, NUM_PHASES, evaluator, MAX_NUM_EXONS, 
                        INTERGENIC, signals, contexts);

    // DONT FORGET
    GlobalVector params(fte.getFeatures()); // final regularized parameters

#ifdef HAVE_MPI 
    InpDir trainDir("default", config.trainingSequences(), FASTA, sigma);
    
    MPStructureCore core(fsm, lattice, re, fe, fte, 
                         evaluator, printer, avgMethod, combMethod,
                         (TStrand)origStrand, kBest,
                         maxLoss, multiUpd, oracleWait,
			 oracleUpd, addUnreachable);
    

    if(sampleSizePerc)
      core.enableRandomSampling(sampleSizePerc);

    core.setPhi(phi);
    core.setR(r);
    core.setMinLongerWait(minLongerWait);
    core.setId34LenSorting(id34LenSorting);

    // Start up MPI
    int rank, size;
    
    MPI::Init();

    rank = MPI::COMM_WORLD.Get_rank();
    size = MPI::COMM_WORLD.Get_size();
    
    if(size == 1) {
      assert(0);
      throw EXCEPTION(BAD_USAGE, "MPI option requires np > 1 to work well");
    }

    MPI::Group group = MPI::COMM_WORLD.Get_group();
    int *incl = new int [size - 1];
    for(int i = 0; i < size - 1; i++)
      incl[i] = i + 1;

    MPI::Intracomm workerComm = MPI::COMM_WORLD.Create(group.Incl(size - 1, incl));
    
    delete [] incl;
    
    if(rank == 0) {
      /* 
       * loading annotSeqs and gene annotations.
       * It's assumed at training time that each Sequence object contains only
       * one gene.
       * One needs to preprocess the input data to fulfil this requirement. 
       * The perl utility splitGenes.pl is provided along the C++ sources in
       * the craig distribution for this purpose.
       */
      
      InpFile validFile("default", config.validationSequences(), FASTA, sigma, false);
      list<Sequence *> & validSequences = (list<Sequence *> &)validFile.sequences();
      TagUtils::loadSeqTags(config.validationSet(), validSequences, fsm);
      core.setValidationSet(validSequences);

      list<Sequence *> trainSequences;
      vector<string> & seqNames = trainDir.seqNames();

      for(int i = 0; i < seqNames.size(); i++) {
        std::string tagFile(config.trainingSet());
        tagFile += std::string("/") + seqNames[i];
        TagUtils::loadSeqTags(tagFile, trainSequences, fsm, DEFAULT_SET, true);
      }

      core.master(params, 
                  trainSequences,
                  DEFAULT_SET, PRED_SET, 
                  learnRate, 
                  maxIterations, 
                  tm,
                  blockSize,
                  shareModels,
                  balanceLoad,
                  (char *)paramFile.c_str());

      list<Sequence *>::iterator it = trainSequences.begin();
      for( ; it != trainSequences.end(); it++)
        delete (*it);

    }
    else {
      int accumIterations = 0, iteration = 0;
      if(logFile.length()) {
	ostringstream nodestr;
	nodestr << workerComm.Get_rank();
	std::string node = nodestr.str(), pattern = "_NODE_";
	Utils::substitute(logFile, pattern, node, true);
        core.initParamsFromLog(params, (char *)logFile.c_str(),
                               accumIterations,
			       iteration);
      }

      core.worker(trainDir,
		  config.trainingSet(),
		  DEFAULT_SET, PRED_SET, 
		  learnRate, 
		  maxIterations, 
		  workerComm,
		  tm,
		  shareModels,
		  (char *)paramFile.c_str(),
		  iteration,
		  accumIterations);
    }
#else
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

    StructureCore core(fsm, lattice, re, fe, fte, 
                       evaluator, printer, avgMethod, combMethod,
                       (TStrand)origStrand, kBest,
		       maxLoss, multiUpd, oracleWait,
		       oracleUpd, addUnreachable);

    core.setValidationSet(validSequences);
    if(sampleSizePerc)
      core.enableRandomSampling(sampleSizePerc);

    core.setPhi(phi);
    core.setR(r);
    core.setMinLongerWait(minLongerWait);
    core.setId34LenSorting(id34LenSorting);
    
    // Decoding    
    if(logFile.length()) {
      core.reStartTraining(params, 
                           trainSequences, 
                           DEFAULT_SET, PRED_SET, 
                           learnRate, 
                           maxIterations, 
                           addedFeatures,
                           tm, 
                           (char *)logFile.c_str(),
                           (char *)paramFile.c_str());
    }
    else {
      core.startTraining(params, 
                         trainSequences, 
                         DEFAULT_SET, PRED_SET, 
                         learnRate, 
                         maxIterations, 
                         tm, 
                         (char *)paramFile.c_str());
    }
    
#endif
    
#ifdef HAVE_MPI
    if(rank == 0) {
#endif
      
      // Storing parameters
      re.saveResources(paramStream);
      fe.saveFilters(paramStream);
      fte.saveFeatures(paramStream);
      params.store(paramStream);
      paramStream.close();
      cerr << "Done\n";

#ifdef HAVE_MPI
    }
#endif    
    
    // releasing tags
    SeqTags::releaseMem();
    delete [] signals;

    reStream.close();
    flStream.close();
    ftStream.close();
    
#ifdef HAVE_MPI
    cerr << "Process " << rank << " exiting\n";
    MPI::Finalize();
#endif

  } catch(exception *e) { 
    perror(e->what());

    GenLibExcept *gle = (GenLibExcept *)e;
    if(gle->error() == BAD_USAGE)
      cerr << "use --help for more information\n";
    delete e;
    exit(1);
    
  }
}
