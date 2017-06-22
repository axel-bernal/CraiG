/******************************************************************************
 * This is training part of CRAIG+R. It trains parameters for a linear
 * structure model for predicting genes. The method uses MIRA online
 * training algorithm. The  input is the training set of genes and a
 * top-K predictions. For each sequence, the algorithm reranks the
 * predictions, so that the one with the highest quality is ranked
 * higher than any other. The output is the gene model.
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
#include "RerankingCore.h"
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
    fd << "discriminative learning for reranking of gene taggers. All file system input is to be expected to be directories to be accessed easily from every created process\n";
#else
    fd << "discriminative learning for reranking of gene taggers\n";
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
    ArgParseUtils::displayOption(fd, "-utr --model-utrs", "Model will be trained using UTR regions If not specified, make sure tags used as training do not contain UTR-labeled regions");
    ArgParseUtils::displayOption(fd, "--best=K", "Trains with the top K-best segmentations [1]");
    ArgParseUtils::displayOption(fd, "--prefix-evidence=PREFIX_EVIDENCE", "prefix for fetching files used as input resources");

    fd << "loss-function arguments:\n";
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
    ArgParseUtils::displayOption(fd, "--comb-method=COMB_METHOD", "One of <euclid|k-l>. Specifies the selected combining method for feature vectors provenient from different iterations. Either one which minimizes the euclidean distance or the k-l divergence between them.[euclid]");
    ArgParseUtils::displayOption(fd, "DEBUG_INFO", "Debug information, useful to know when to stop training");
    fd << "Report bugs to <" << AUTHOR_EMAIL << ">\n";
}

int main(int argc, char *argv[]) {
    ::string modelPath, confFile, paramFile, prefix_evidence("");
    double learnRate = 1, r = 1.0, phi = 0.1;
    int addedFeatures = 0;
    int maxIterations = 100, origStrand = BOTH_STRANDS;
    int i, noiseLevel = 0;
    std::string topology("partial"), strand("both");
    TAvgMethod avgMethod = AVG_ALL;
    bool modelUTRs = false;
    TCombMethod combMethod = COMB_EUCLID;
    TTrainMethod tm = MIRA;
    TLossType lf = LF_HAMMING, elf = LF_NONE;
    double edgeLossFactor = 0;
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
            else if(!strncmp(argv[i], "--verbose", 9)
                    || !strncmp(argv[i], "-v", 2)) {
                verbose = true;
            }
            else if(!strncmp(argv[i], "--noise-level=", 14)) {
                sscanf(argv[i] + 14, "%d", &noiseLevel);
            }
            else if(!strncmp(argv[i], "--added-features=", 17)) {
                sscanf(argv[i] + 17, "%d", &addedFeatures);
            }
            else if(!strncmp(argv[i], "--prefix-evidence=", 18)) {
                prefix_evidence = string(argv[i] + 18);
            }
            else if(!strncmp(argv[i], "--help", 6)
                    || !strncmp(argv[i], "-h", 2)) {
                printHelp("rtrain_craig", (std::ofstream &)cout);
                exit(0);
            }
            else if(!strncmp(argv[i], "--version", 9)) {
                PRINT_VERSION(cerr, "rtrain_craig", "tool for discriminative reranking of gene taggers for eukarya");
                PRINT_DISCLAIMER(cerr, "rtrain_craig");
                exit(0);
            }
            else if(!strncmp(argv[i], "--restart=", 10)) {
                logFile = std::string(argv[i] + 10);
            }
            else if(!strncmp(argv[i], "--loss-function=", 16)) {
                std::string lossString = std::string(argv[i] + 16);
                lf = TypeDefs::stringToTLossType(lossString);
            }
            else if(!strncmp(argv[i], "--edge-loss-function=", 21)) {
                std::string lossString = std::string(argv[i] + 21);
                elf = TypeDefs::stringToTLossType(lossString);
            }
            else if(!strncmp(argv[i], "--edge-loss-factor=", 19)) {
                sscanf(argv[i] + 19, "%lf", &edgeLossFactor);
            }
            else if(!strncmp(argv[i], "--avg-method=", 13)) {
                std::string avgString = std::string(argv[i] + 13);
                avgMethod = TypeDefs::stringToTAvgMethod(avgString);
            }
            else if(!strncmp(argv[i], "--comb-method=", 14)) {
                std::string combMethodString = std::string(argv[i] + 14);
                combMethod = TypeDefs::stringToTCombMethod(combMethodString);
            }
            else if(!strncmp(argv[i], "--algorithm=", 12)) {
                std::string tmString = std::string(argv[i] + 12);
                tm = TypeDefs::stringToTTrainMethod(tmString);
            }
            else if(!strncmp(argv[i], "--phi=", 6)) {
                sscanf(argv[i] + 6, "%lf", &phi);
            }
            else if(!strncmp(argv[i], "--r=", 4)) {
                sscanf(argv[i] + 4, "%lf", &r);
            }
            else
                break;

        }

        // checking sanity
        if(combMethod == COMB_KL && tm != ARROW && tm != CWL)
            throw EXCEPTION(NOT_SUPPORTED, "K-L average only available for ARROW and CWL");

        if(argc - i < 2)
            throw EXCEPTION(BAD_USAGE, "insufficient arguments");

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
        resource_subs["PREFIX_EVIDENCE"] = prefix_evidence;
        resource_subs["~"] = modelPath;
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

        TagUtils::loadSeqTags(config.trainingSet(), trainSequences, fsm);
        TagUtils::loadSeqTags(config.validationSet(), trainSequences, fsm, TRAIN_SET);

        GeneEvaluator evaluator(fsm, lf, elf, noiseLevel, 0);
        GeneTagPrinter printer;

        if(edgeLossFactor)
            evaluator.setEdgeLossFactor(edgeLossFactor);

        //Creating RerankingCore object
        Lattice lattice(fe, fsm, NUM_PHASES, evaluator,
                        MAX_NUM_EXONS, INTERGENIC);

        GlobalVector params(fte.getFeatures()); // final regularized parameters

        RerankingCore core(fsm, lattice, re, fe, fte, evaluator,
                           printer, avgMethod, combMethod, (TStrand)origStrand,
                           TRAIN_SET);
        core.setPhi(phi);
        core.setR(r);

        cerr << "Reranking model for organism " << confFile << endl;
        cerr << "Process Id " << getpid() << endl;

        if(logFile.length()) {
            core.reStartTraining(params,
                                 trainSequences,
                                 DEFAULT_SET, PRED_SET,
                                 learnRate,
                                 maxIterations,
                                 addedFeatures,
                                 tm,
                                 logFile.c_str(),
                                 paramFile.c_str());
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

        cerr << "Done\n";

        // Storing parameters
        re.saveResources(paramStream);
        fe.saveFilters(paramStream);
        fte.saveFeatures(paramStream);
        params.store(paramStream);
        paramStream.close();
        reStream.close();
        flStream.close();
        ftStream.close();
        // releasing tags
        SeqTags::releaseMem();

    } catch(exception *e) {
        perror(e->what());

        GenLibExcept *gle = (GenLibExcept *)e;
        if(gle->error() == BAD_USAGE)
            cerr << "use --help for more information\n";
        delete e;
        exit(1);

    }
}
