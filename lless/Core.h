/****************************************************************************
 * Core.h - part of the lless namespace, a general purpose
 *          linear semi-markov structure prediction library
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

#ifndef _CORE_ALGORITHMS_H_
#define _CORE_ALGORITHMS_H_
#include "Organism.h"
#include "FeatureEngine.h"
#include "FilterEngine.h"
#include "Lattice.h"
#include "FSM.h"
#include "ContextIMM.h"
#include "Evaluator.h"
#include "TagPrinter.h"
#include "GlobalVector.h"

#define MAX_RESTRICTIONS 1
#define MAX_MULTILABELINGS 1
#define EPS 1e-10
#define MAX_ITER 1000
#define ZERO exp((double)-DBL_MAX_EXP)

namespace lless {

    /**
     * The Core class integrates resources, filters and features and implements
     * all the core subroutines fo  training and predicting models for linear
     * structure prediction. So far there are two training algorithms which
     * have been implemented: k-best perceptron and MIRA - including the
     * hildreth subroutine to solve the quadratic program.
     *
     ***************************************************************************/

    typedef pair<int, double> loss_t;
    class Core {

    protected:
        FSM *_fsm;
        Lattice *_lattice;
        Evaluator *_evaluator;
        FilterEngine *_fe;
        ResourceEngine *_re;
        TagPrinter *_printer;
        vector<double> loss;
        double _phi, _r;
        FeatureEngine *_fte;
        TAvgMethod _avgMethod;
        TCombMethod _combMethod;

        //!  feature vectors for temporary operations
        GlobalVector *incorr_gv;
        vector<GlobalVector *> gv1, gv2;

        //! feature vectors for permanent operations
        GlobalVector *accum_gv,
            *pivot_gv, *sigma_gv, *accsigma_gv;



        int loStrand, upStrand;
        int t; //!< absolute number of iterations during training

    public:
        Core(
            FSM &fsm,
            Lattice &lattice,
            ResourceEngine &re,
            FilterEngine &fe,
            FeatureEngine &fte,
            Evaluator &evaluator,
            TagPrinter &printer,
            TAvgMethod avgMethod,
            TCombMethod combMethod,
            TStrand strand);

        inline FSM *fsm() {
            return _fsm;
        }

        inline FilterEngine *filterEngine() {
            return _fe;
        }

        inline Evaluator *evaluator() {
            return _evaluator;
        }

        inline Lattice *lattice() {
            return _lattice;
        }

        inline void setPhi(double phi) {
            _phi = phi;
        }

        inline void setR(double r) {
            _r = r;
        }

        inline void freetmpFeatVectors() {
            for(int i = 0; i < gv1.size(); i++)
                if(gv1[i]) delete gv1[i];

            gv1.clear();

            for(int i = 0; i < gv2.size(); i++)
                if(gv2[i]) delete gv2[i];

            gv2.clear();

            if(incorr_gv)
                delete incorr_gv;
            incorr_gv = NULL;

            loss.clear();

        }

        double *hildreth(vector<GlobalVector *> &gv, vector<double> &loss, int size);

        void prePSequence(Sequence &);
        void postPSequence(bool deleteSeqResources = true);

        GlobalVector *computeGlobalVector(SeqTags &seqTags,
                                          double learnRate = 1.0);

        GlobalVector *computeGlobalVector(Sequence &annotSeq,
                                          vector<SeqTags> &,
                                          double learnRate = 1.0);

        loss_t computeMinLoss(vector<double> &losses,
                              int numLosses);

        loss_t computeMinLoss(Sequence &c,
                              SeqTags &tags,
                              vector<SeqTags> &vtags,
                              int size = -1);

        vector<loss_t> computeMinLoss(Sequence &c,
                                      vector<SeqTags> &vtagsA,
                                      vector<SeqTags> &vtagsB,
                                      int sizeB = -1);

        vector<double> computeLosses(Sequence &c,
                                     SeqTags &tags,
                                     vector<SeqTags> &vtags,
                                     int size = -1);

        virtual int prepareParamUpdate(Sequence &c,
                                       vector<SeqTags> &expTags,
                                       vector<SeqTags> &predTags,
                                       double learnRate) = 0;

        void updInstancePerceptron(Sequence &,
                                   vector<SeqTags> & expTags,
                                   vector<SeqTags> & predTags,
                                   double learnRate);

        void updInstanceMIRA(Sequence &,
                             vector<SeqTags> & expTags,
                             vector<SeqTags> & predTags,
                             double learnRate);

        void updInstancePEGASOS(Sequence &,
                                vector<SeqTags> & expTags,
                                vector<SeqTags> & predTags,
                                double learnRate);

        void updInstanceCWL(Sequence &c,
                            vector<SeqTags> & expTags,
                            vector<SeqTags> & predTags,
                            double learnRate);

        void updInstanceARROW(Sequence &c,
                              vector<SeqTags> & expTags,
                              vector<SeqTags> & predTags,
                              double learnRate);

        virtual int argmaxAndUpdate(GlobalVector &params, list<Sequence *> &annotSeqs,
                                    vector<bool> &,
                                    TSetType trainSet, TSetType predSet,
                                    double learnRate,
                                    TTrainMethod = MIRA) = 0;

        void initParamsFromLog(GlobalVector &finalParams,
                               const char *logFile,
                               int &iteration,
                               int &accumIterations,
                               int addedFeatures = 0);

        void reStartTraining(GlobalVector & gv, list<Sequence *> &,
                             TSetType, TSetType,
                             double learnRate, int maxIterations,
                             int addedFeatures = 0,
                             TTrainMethod tm = MIRA,
                             const char *logFile = "log_params",
                             const char *paramsFile = "params");
        virtual void startTraining(GlobalVector & gv,
                                   list<Sequence *> &, TSetType, TSetType,
                                   double learnRate, int maxIterations,
                                   TTrainMethod tm = MIRA,
                                   const char *paramsFiles = "params",
                                   int iteration = 0,
                                   int accumIterations = 0) = 0;

        virtual double predictTags(Sequence &c, TSetType,
                                   list<BioFeature> * = NULL) = 0;

        void storeParams(const char *paramsFile, int iteration,
                         GlobalVector *avg);
        void storeLog(const char *paramsFile, int iteration,
                      int accumIterations);
        void storeIteration(const char *paramsFile, int iteration,
                            int accumIterations,
                            GlobalVector *avg,
                            bool saveLog = true);

        virtual ~Core();

    };

}
#endif
