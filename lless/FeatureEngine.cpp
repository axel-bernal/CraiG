#include "FeatureEngine.h"
#include "FT_StateBackground.h"
#include "FT_ScoringSegment.h"
#include "FT_MeanScoringSegment.h"
#include "FT_Edge.h"
#include "FT_EdgeFeature.h"
#include "FT_EdgePhase.h"
#include "FT_EdgeContent.h"
#include "FT_Bin.h"
#include "FT_LinearFitBin.h"
#include "FT_CSplineBin.h"
#include "FT_MultiBin.h"
#include "FT_HarmonicSTerm.h"
#include "FT_HiddenSeq.h"
#include "FT_WWAM.h"
#include "FT_GramWWAM.h"
#include "FT_ConsensusWWAM.h"
#include "FT_PeptideWWAM.h"
#include "FT_PWM.h"
#include "FT_PWMxPWM.h"
#include "FT_CountingSegment.h"
#include "FT_PeriodicCountingSegment.h"
#include "FT_BinnedCountingSegment.h"
#include "FT_BinnedSegment.h"
#include "FT_GreaterThan.h"
#include "FT_SoftMatchWordCounter.h"
#include "FT_LengthSegment.h"
#include "FT_FeatureBag.h"
#include "FT_Composition.h"
#include "FT_Conjunction.h"
#include "FT_ParseScore.h"
#include "FT_EvdAlignment.h"
#include "FT_MinSegmentScore.h"
#include "FT_ChangePtScore.h"

/****************************************************************************
 * FeatureEngine.cpp - part of the lless namespace, a general purpose
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


using namespace lless;

namespace lless {

    /**
     * Registers all Feature objects that can be created dynamically in the
     * Featuyre Header definition file
     */
    void FeatureEngine::registerDefaultFeatures() {
        if(!featureFactory.Register(std::string("FilterWrapper<int>"), Type2Type<FT_FilterWrapper<int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("FilterWrapper<double>"), Type2Type<FT_FilterWrapper<double> >()))
            assert(0);
        if(!featureFactory.Register(std::string("FilterWrapper<UCHAR>"), Type2Type<FT_FilterWrapper<UCHAR> >()))
            assert(0);
        if(!featureFactory.Register(std::string("StateBackground"), Type2Type<FT_StateBackground>()))
            assert(0);
        if(!featureFactory.Register(std::string("FilterXFeature<UCHAR,UCHAR>"), Type2Type<FT_FilterXFeature<UCHAR, double> >()))
            assert(0);
        if(!featureFactory.Register(std::string("FilterXFeature<UCHAR,int>"), Type2Type<FT_FilterXFeature<UCHAR, int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("FilterXFeature<UCHAR,double>"), Type2Type<FT_FilterXFeature<UCHAR, double> >()))
            assert(0);
        if(!featureFactory.Register(std::string("LazyFilterXFeature<UCHAR,int>"), Type2Type<FT_LazyFilterXFeature<UCHAR, int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("LazyFilterXFeature<UCHAR,UCHAR>"), Type2Type<FT_LazyFilterXFeature<UCHAR, UCHAR> >()))
            assert(0);
        if(!featureFactory.Register(std::string("LazyFilterXFeature<UCHAR,USHORT>"), Type2Type<FT_LazyFilterXFeature<UCHAR, USHORT> >()))
            assert(0);
        if(!featureFactory.Register(std::string("LazyFilterXFeature<UCHAR,ULONG>"), Type2Type<FT_LazyFilterXFeature<UCHAR, ULONG> >()))
            assert(0);
        if(!featureFactory.Register(std::string("LazyFilterXFeature<UCHAR,double>"), Type2Type<FT_LazyFilterXFeature<UCHAR, double> >()))
            assert(0);
        if(!featureFactory.Register(std::string("RatioSegment"), Type2Type<FT_RatioSegment>()))
            assert(0);
        if(!featureFactory.Register(std::string("ScoringSegment<int>"), Type2Type< FT_ScoringSegment<int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("ScoringSegment<double>"), Type2Type< FT_ScoringSegment<double> >()))
            assert(0);
        if(!featureFactory.Register(std::string("MeanScoringSegment<double>"), Type2Type< FT_MeanScoringSegment<double> >()))
            assert(0);
        if(!featureFactory.Register(std::string("MeanScoringSignal<double>"), Type2Type< FT_MeanScoringSignal<double> >()))
            assert(0);
        if(!featureFactory.Register(std::string("RampScoringSignal<double>"), Type2Type< FT_RampScoringSignal<double> >()))
            assert(0);
        if(!featureFactory.Register(std::string("StDevScoringSegment<double>"), Type2Type< FT_StDevScoringSegment<double> >()))
            assert(0);
        if(!featureFactory.Register(std::string("WeightedScoringSegment<int>"), Type2Type< FT_WeightedScoringSegment<int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("WeightedScoringSegment<double>"), Type2Type< FT_WeightedScoringSegment<double> >()))
            assert(0);
        if(!featureFactory.Register(std::string("ScoringAtPos<int>"), Type2Type< FT_ScoringAtPos<int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("ScoringAtPos<double>"), Type2Type< FT_ScoringAtPos<double> >()))
            assert(0);
        if(!featureFactory.Register(std::string("CodingDifferential"), Type2Type<FT_CodingDifferential>()))
            assert(0);
        if(!featureFactory.Register(std::string("LengthSegment"), Type2Type<FT_LengthSegment>()))
            assert(0);
        if(!featureFactory.Register(std::string("HarmonicSTerm"), Type2Type<FT_HarmonicSTerm>()))
            assert(0);
        if(!featureFactory.Register(std::string("Bin<UCHAR>"), Type2Type<FT_Bin<UCHAR> >()))
            assert(0);
        if(!featureFactory.Register(std::string("Bin<int>"), Type2Type<FT_Bin<int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("Bin<double>"), Type2Type<FT_Bin<double> >()))
            assert(0);
        if(!featureFactory.Register(std::string("MultiBin<int>"), Type2Type<FT_MultiBin<int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("MultiBin<UCHAR>"), Type2Type<FT_MultiBin<UCHAR> >()))
            assert(0);
        if(!featureFactory.Register(std::string("MultiBin<double>"), Type2Type<FT_MultiBin<double> >()))
            assert(0);
        if(!featureFactory.Register(std::string("RBin<int>"), Type2Type<FT_RBin<int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("RBin<double>"), Type2Type<FT_RBin<double> >()))
            assert(0);
        if(!featureFactory.Register(std::string("CustomBin<int>"), Type2Type<FT_CustomBin<int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("CustomBin<double>"), Type2Type<FT_CustomBin<double> >()))
            assert(0);
        if(!featureFactory.Register(std::string("CustomRBin<int>"), Type2Type<FT_CustomRBin<int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("CustomRBin<double>"), Type2Type<FT_CustomRBin<double> >()))
            assert(0);
        if(!featureFactory.Register(std::string("LinearFitBin<int>"), Type2Type<FT_LinearFitBin<int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("LinearFitBin<UCHAR>"), Type2Type<FT_LinearFitBin<UCHAR> >()))
            assert(0);
        if(!featureFactory.Register(std::string("LinearFitBin<double>"), Type2Type<FT_LinearFitBin<double> >()))
            assert(0);
        if(!featureFactory.Register(std::string("LinearFitCustomBin<int>"), Type2Type<FT_LinearFitCustomBin<int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("LinearFitCustomBin<UCHAR>"), Type2Type<FT_LinearFitCustomBin<UCHAR> >()))
            assert(0);
        if(!featureFactory.Register(std::string("LinearFitCustomBin<double>"), Type2Type<FT_LinearFitCustomBin<double> >()))
            assert(0);
        if(!featureFactory.Register(std::string("CSplineBin<int>"), Type2Type<FT_CSplineBin<int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("CSplineBin<double>"), Type2Type<FT_CSplineBin<double> >()))
            assert(0);
        if(!featureFactory.Register(std::string("CSplineCustomBin<int>"), Type2Type<FT_CSplineCustomBin<int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("CSplineCustomBin<double>"), Type2Type<FT_CSplineCustomBin<double> >()))
            assert(0);
        if(!featureFactory.Register(std::string("CountingAtPos<int>"), Type2Type<FT_CountingAtPos<int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("CountingAtPos<UCHAR>"), Type2Type<FT_CountingAtPos<UCHAR> >()))
            assert(0);
        if(!featureFactory.Register(std::string("CountingAtPos<USHORT>"), Type2Type<FT_CountingAtPos<USHORT> >()))
            assert(0);
        if(!featureFactory.Register(std::string("PeriodicCountingSegment<UCHAR>"), Type2Type<FT_PeriodicCountingSegment<UCHAR> >()))
            assert(0);
        if(!featureFactory.Register(std::string("PeriodicCountingSegment<USHORT>"), Type2Type<FT_PeriodicCountingSegment<USHORT> >()))
            assert(0);
        if(!featureFactory.Register(std::string("PeriodicCountingSegment<int>"), Type2Type<FT_PeriodicCountingSegment<int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("PeriodicCountingSegment<ULONG>"), Type2Type<FT_PeriodicCountingSegment<ULONG> >()))
            assert(0);
        if(!featureFactory.Register(std::string("DenseBinnedSegment<int,int>"), Type2Type<FT_DenseBinnedSegment<int, int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("DenseBinnedSegment<UCHAR,UCHAR>"), Type2Type<FT_DenseBinnedSegment<UCHAR, UCHAR> >()))
            assert(0);
        if(!featureFactory.Register(std::string("DenseBinnedSegment<int,UCHAR>"), Type2Type<FT_DenseBinnedSegment<int,UCHAR> >()))
            assert(0);
        if(!featureFactory.Register(std::string("DenseBinnedSegment<double,int>"), Type2Type<FT_DenseBinnedSegment<double, int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("SparseBinnedSegment<UCHAR,int,UCHAR>"), Type2Type<FT_SparseBinnedSegment<UCHAR, int, UCHAR> >()))
            assert(0);
        if(!featureFactory.Register(std::string("SparseBinnedSegment<USHORT,int,USHORT>"), Type2Type<FT_SparseBinnedSegment<USHORT, int, USHORT> >()))
            assert(0);
        if(!featureFactory.Register(std::string("PhasedSparseBinnedSegment<UCHAR,int,UCHAR>"), Type2Type<FT_PhasedSparseBinnedSegment<UCHAR, int, UCHAR> >()))
            assert(0);
        if(!featureFactory.Register(std::string("PhasedSparseBinnedSegment<USHORT,int,USHORT>"), Type2Type<FT_PhasedSparseBinnedSegment<USHORT, int, USHORT> >()))
            assert(0);
        if(!featureFactory.Register(std::string("BinnedPeriodicCountingSegment<UCHAR,UCHAR>"), Type2Type<FT_BinnedPeriodicCountingSegment<UCHAR, UCHAR> >()))
            assert(0);
        if(!featureFactory.Register(std::string("BinnedPeriodicCountingSegment<UCHAR,UCHAR#>"), Type2Type<FT_BinnedPeriodicCountingSegment<UCHAR, SPARSE_HASH<UCHAR, int> > >()))
            assert(0);
        if(!featureFactory.Register(std::string("BinnedPeriodicCountingSegment<double,int>"), Type2Type<FT_BinnedPeriodicCountingSegment<double, int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("BinnedPeriodicCountingSegment<double,UCHAR>"), Type2Type<FT_BinnedPeriodicCountingSegment<double, UCHAR> >()))
            assert(0);
        if(!featureFactory.Register(std::string("BinnedPeriodicCountingSegmentWUpperLimit<double,UCHAR>"), Type2Type<FT_BinnedPeriodicCountingSegmentWUpperLimit<double, UCHAR> >()))
            assert(0);

        //    if(!featureFactory.Register(std::string("AggregatedPeriodicCountingSegment"), Type2Type<FT_AggregatedPeriodicCountingSegment>()))
        //      assert(0);
        if(!featureFactory.Register(std::string("TrimmedCountingSegment<int>"), Type2Type<FT_TrimmedCountingSegment<int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("TrimmedCountingSegment<UCHAR>"), Type2Type<FT_TrimmedCountingSegment<UCHAR> >()))
            assert(0);
        //    if(!featureFactory.Register(std::string("AggregatedTrimmedCountingSegment"), Type2Type<FT_AggregatedTrimmedCountingSegment>()))
        //      assert(0);
        if(!featureFactory.Register(std::string("PeriodicCountingSegmentWUpperLimit<int>"), Type2Type<FT_PeriodicCountingSegmentWUpperLimit<int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("PeriodicCountingSegmentWUpperLimit<UCHAR>"), Type2Type<FT_PeriodicCountingSegmentWUpperLimit<UCHAR> >()))
            assert(0);
        if(!featureFactory.Register(std::string("HiddenSeq"), Type2Type<FT_HiddenSeq>()))
            assert(0);
        if(!featureFactory.Register(std::string("Edge<int>"), Type2Type<FT_Edge<int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("Edge<UCHAR>"), Type2Type<FT_Edge<UCHAR> >()))
            assert(0);
        if(!featureFactory.Register(std::string("Edge<double>"), Type2Type<FT_Edge<double> >()))
            assert(0);
        if(!featureFactory.Register(std::string("EdgeFeature<int>"), Type2Type<FT_EdgeFeature<int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("EdgeFeature<UCHAR>"), Type2Type<FT_EdgeFeature<UCHAR> >()))
            assert(0);
        if(!featureFactory.Register(std::string("EdgeFeature<double>"), Type2Type<FT_EdgeFeature<double> >()))
            assert(0);
        if(!featureFactory.Register(std::string("EdgePhase"), Type2Type<FT_EdgePhase >()))
            assert(0);
        if(!featureFactory.Register(std::string("NodePhase"), Type2Type<FT_NodePhase >()))
            assert(0);
        if(!featureFactory.Register(std::string("EdgeContent<int>"), Type2Type<FT_EdgeContent<int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("EdgeContent<double>"), Type2Type<FT_EdgeContent<double> >()))
            assert(0);
        if(!featureFactory.Register(std::string("GramWWAM<UCHAR,UCHAR>"), Type2Type<FT_GramWWAM<UCHAR,UCHAR> >()))
            assert(0);
        if(!featureFactory.Register(std::string("GramWWAM<UCHAR,int>"), Type2Type<FT_GramWWAM<UCHAR,int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("BinnedGramWWAM<int,UCHAR,UCHAR>"), Type2Type<FT_BinnedGramWWAM<int,UCHAR,UCHAR> >()))
            assert(0);
        if(!featureFactory.Register(std::string("WWAMUnion<int>"), Type2Type<FT_WWAMUnion<int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("WWAMUnion<UCHAR>"), Type2Type<FT_WWAMUnion<UCHAR> >()))
            assert(0);
        if(!featureFactory.Register(std::string("GramPWM<UCHAR,UCHAR>"), Type2Type<FT_GramPWM<UCHAR,UCHAR> >()))
            assert(0);
        if(!featureFactory.Register(std::string("GramPWM<UCHAR,USHORT>"), Type2Type<FT_GramPWM<UCHAR,USHORT> >()))
            assert(0);
        if(!featureFactory.Register(std::string("GramPWM<UCHAR,int>"), Type2Type<FT_GramPWM<UCHAR,int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("PWMUnion<int>"), Type2Type<FT_PWMUnion<int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("PWMUnion<UCHAR>"), Type2Type<FT_PWMUnion<UCHAR> >()))
            assert(0);
        if(!featureFactory.Register(std::string("ConsensusPWM<int>"), Type2Type<FT_ConsensusPWM<int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("ConsensusPWM<UCHAR>"), Type2Type<FT_ConsensusPWM<UCHAR> >()))
            assert(0);
        if(!featureFactory.Register(std::string("PWMxPWM<int>"), Type2Type<FT_PWMxPWM<int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("PWMxPWM<UCHAR>"), Type2Type<FT_PWMxPWM<UCHAR> >()))
            assert(0);
        if(!featureFactory.Register(std::string("PWMxFeature<UCHAR,double>"), Type2Type<FT_PWMxFeature<UCHAR,double> >()))
            assert(0);
        if(!featureFactory.Register(std::string("PeptideWWAM<UCHAR,int>"), Type2Type<FT_PeptideWWAM<UCHAR,int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("PeptideWWAM<UCHAR,UCHAR>"), Type2Type<FT_PeptideWWAM<UCHAR,UCHAR> >()))
            assert(0);
        if(!featureFactory.Register(std::string("ConsensusWWAM<int>"), Type2Type<FT_ConsensusWWAM<int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("ConsensusWWAM<UCHAR>"), Type2Type<FT_ConsensusWWAM<UCHAR> >()))
            assert(0);
        if(!featureFactory.Register(std::string("FeatureXFeature<UCHAR,int>"), Type2Type<FT_FeatureXFeature<UCHAR, int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("FeatureXFeature<UCHAR,UCHAR>"), Type2Type<FT_FeatureXFeature<UCHAR, UCHAR> >()))
            assert(0);
        if(!featureFactory.Register(std::string("FeatureXFeature<UCHAR,double>"), Type2Type<FT_FeatureXFeature<UCHAR, double> >()))
            assert(0);
        if(!featureFactory.Register(std::string("SoftMatchWordCounter<int>"), Type2Type<FT_SoftMatchWordCounter<int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("SoftMatchWordCounter<double>"), Type2Type<FT_SoftMatchWordCounter<double> >()))
            assert(0);
        if(!featureFactory.Register(std::string("GreaterThan<int>"), Type2Type<FT_GreaterThan<int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("GreaterThan<double>"), Type2Type<FT_GreaterThan<double> >()))
            assert(0);
        if(!featureFactory.Register(std::string("FeatureOFeature<int>"), Type2Type<FT_FeatureOFeature<int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("FeatureOFeature<double>"), Type2Type<FT_FeatureOFeature<double> >()))
            assert(0);
        if(!featureFactory.Register(std::string("FeatureO2Feature<int>"), Type2Type<FT_FeatureO2Feature<int> >()))
            assert(0);
        if(!featureFactory.Register(std::string("FeatureO2Feature<UCHAR>"), Type2Type<FT_FeatureO2Feature<UCHAR> >()))
            assert(0);
        if(!featureFactory.Register(std::string("FeatureO2Feature<double>"), Type2Type<FT_FeatureO2Feature<double> >()))
            assert(0);
        if(!featureFactory.Register(std::string("FilterOFeature<UCHAR,UCHAR>"), Type2Type<FT_FilterOFeature<UCHAR,UCHAR> >()))
            assert(0);
        if(!featureFactory.Register(std::string("FilterOFeature<UCHAR,double>"), Type2Type<FT_FilterOFeature<UCHAR, double> >()))
            assert(0);
        if(!featureFactory.Register(std::string("FilterOFeature<double,double>"), Type2Type<FT_FilterOFeature<double, double> >()))
            assert(0);
        if(!featureFactory.Register(std::string("FeatureBag"), Type2Type<FT_FeatureBag>()))
            assert(0);
        if(!featureFactory.Register(std::string("ParseScore"), Type2Type<FT_ParseScore>()))
            assert(0);
        if(!featureFactory.Register(std::string("AlignmentType"), Type2Type<FT_EvdAlignmentType>()))
            assert(0);
        if(!featureFactory.Register(std::string("AlignmentWeight"), Type2Type<FT_EvdAlignmentWeight>()))
            assert(0);
        if(!featureFactory.Register(std::string("MinSegmentScore"), Type2Type<FT_MinSegmentScore>()))
            assert(0);
        if(!featureFactory.Register(std::string("ChangePtScore"), Type2Type<FT_ChangePtScore>()))
            assert(0);

    }


    void FeatureEngine::preComputeFeatures(int beg, int end, TStrand strand) {
        // Allocating memory for precomputed dot values
        this->fo->allocatePreComputedFeatures(beg, end, strand);
        // preCompute Phi Dot Theta values
        this->fo->preComputeFeatures(beg, end, strand);

        /*
         * Check for features that are bags of features, if so call this
         * function recursively.
         */
        for(unsigned int i = 0; i < this->features.size(); i++) {
            if(!this->features[i]->isFeatureBag())
                continue;

            FT_FeatureBag *f = (FT_FeatureBag *)this->features[i];
            f->doPreComputation(NULL, 0, beg, end, strand);

        }

    }


    void FeatureEngine::deletePreComputedFeatures() {
        /*
         * Check for features that are bags of features, if so call this
         * function recursively.
         */
        for(unsigned int i = 0; i < this->features.size(); i++) {
            if(!this->features[i]->isFeatureBag())
                continue;

            FT_FeatureBag *f = (FT_FeatureBag *)this->features[i];
            f->undoPreComputation();

        }

        /*
         * free PreComputed arrays
         */
        this->fo->deletePreComputedFeatures();

    }

    /**
     * Ties Tag objects to this feature. This means that the feature is part of
     * the model and can eventually learn its parameters and use them for
     * prediction later.
     */
    void FeatureEngine::tieTagsToFeature(Feature *f, std::string tagGroups) {
        //      cerr << "params found for " << f->getName() << " " << tagGroups << endl;
        if(f->paramInd() < 0)
            f->setParamInd(f->ind());

        /*
         * Substitute any occurrance of a predefined node or edge Set
         */
        std::map<std::string, std::string>::iterator it = tagSets.begin();
        for( ; it != tagSets.end(); it++)
            Utils::substitute(tagGroups, it->first, it->second, true);

        int numGEs;
        int pcEntryIndex, subfVal;
        boost::RegEx rExDelim("\\s+");
        boost::RegEx rExCommaDelim("\\s*\\,\\s*");
        boost::RegEx rExTagPCompGroup("^\\s*(\\d+)\\s+(\\S.+\\S)\\s+\\@\\s+(\\d+)\\s*$");
        boost::RegEx rExTagSubFGroup("^\\s*(\\d+)\\s+(\\S.+\\S)\\s+\\|\\s+(\\d+)\\s*$");
        boost::RegEx rExTagGroup("^\\s*(\\d+)\\s+(\\S.+\\S)\\s*$");

        std::vector<std::string> fTagGroups;
        rExCommaDelim.Split(fTagGroups, tagGroups, 0, 100);

        /*
         * Iterate through all the tags associated with this feature
         */
        for(unsigned int i = 0; i < fTagGroups.size(); i++) {
            bool match = rExTagGroup.Match(fTagGroups[i]);
            bool subf_match = rExTagSubFGroup.Match(fTagGroups[i]);
            bool pcomp_match = rExTagPCompGroup.Match(fTagGroups[i]);
            pcEntryIndex  = -1;
            subfVal = -1;
            if(!match && !subf_match && !pcomp_match)
                throw EXCEPTION( PARSE_ERROR, std::string("undefined Tag ") + fTagGroups[i]);

            if(!sscanf(rExTagGroup[1].c_str(), "%d", &numGEs))
                assert(0);

            std::string tags(rExTagGroup[2]);

            if(subf_match) {
                tags = std::string(rExTagSubFGroup[2]);
                if(!sscanf(rExTagSubFGroup[3].c_str(), "%d", &subfVal))
                    assert(0);

            }
            else if(pcomp_match) {
                /*
                 * This feature precomputes its value
                 */
                int entry;
                tags = std::string(rExTagPCompGroup[2]);

                if(!sscanf(rExTagPCompGroup[3].c_str(), "%d", &entry))
                    assert(0);

                std::map<int, TPCEntry>::iterator eit = pcEntries.find(entry);
                if(eit == pcEntries.end()) {
                    assert(0);
                    throw EXCEPTION(PARSE_ERROR,
                                    std::string("undefined precomputed entry ") + rExTagPCompGroup[3]);
                }
                if(pcEntries[entry].arrIndex < 0) {
                    /*
                     * This feature precomoutes its value in a new entry.
                     * Take then the next available slot for precopmutation
                     */
                    pcEntries[entry].arrIndex = usedPCArrIndexes;
                    this->fo->setPreComputedEntry(
                        pcEntries[entry].arrIndex,
                        pcEntries[entry].arrPeriod,
                        pcEntries[entry].collapseFrames,
                        f);
                }

                pcEntryIndex = pcEntries[entry].arrIndex;
            }

            vector<std::string> ftags;
            rExDelim.Split(ftags, tags, 0, 100);
            if((unsigned)numGEs != ftags.size()) {
                assert(0);
                throw EXCEPTION( PARSE_ERROR,
                                 f->getName() + std::string(" has wrong number of tag groups"));
            }

            /*
             * Tie the tag to the feature
             */
            for(int j = 0; j < numGEs; j++) {
                int pcInd = pcEntryIndex;
                bool tied = false;
                Node *node = fsm->findNode(ftags[j]);

                if(node != NULL) { // the tag was a node
                    pcInd = this->fo->organize(f, node, pcEntryIndex, subfVal);
                    tied = true;
                }

                Edge *edge = fsm->findEdge(ftags[j]);
                if(!tied && edge != NULL) { // it was an edge
                    //	  if(edge->hasSignal() && (pcEntryIndex < 0 || !f->signal()))
                    if(edge->hasSignal() && pcEntryIndex >= 0 && !f->signal())
                        throw EXCEPTION( NOT_SUPPORTED,
                                         f->getName() + std::string(" must precompute through a signal filter"));

                    pcInd = this->fo->organize(f, edge, pcEntryIndex, subfVal);
                    tied = true;
                }

                Word *word = fsm->findWord(ftags[j]);
                if(!tied && word != NULL) { // it was a word
                    pcInd = this->fo->organize(f, word, pcEntryIndex, subfVal);
                    tied = true;
                }

                if(!tied) {
                    assert(0);
                    throw EXCEPTION( PARSE_ERROR, std::string("undefined Tag ") + ftags[j]);
                }

                usedPCArrIndexes = (int)Utils::max(pcInd + 1, usedPCArrIndexes);
                if(usedPCArrIndexes >= MAX_PRECOMP_ARRAYS) {
                    assert(0);
                    throw EXCEPTION( OUT_OF_MEMORY, std::string("preComp Array size too short"));
                }
            }
        }
    }


    /**
     * Reads a Feature Header definition from single line
     */
    bool FeatureEngine::readFeatureHeader(int & fInd, std::string & line) {
        boost::RegEx rExCommaDelim("\\s*\\,\\s*");
        boost::RegEx rExDelim("\\s+");
        boost::RegEx rExSparseFeature("^\\s*Sparse(\\S.+)\\s*$");
        boost::RegEx rExPhaseFeature("^\\s*Phase(\\S.+)\\s*$");
        boost::RegEx rExMultiFeature("^\\s*Multi(\\S+)\\((\\S+)\\)(\\s+\\S.+)\\s*$");
        boost::RegEx rExPackedFeature("^\\s*Packed(\\S+)\\((\\S+)\\)(\\s+\\S.+)\\s*$");
        boost::RegEx rExFeatureName("^\\S+\\s+(\\S+)\\s+.+$");
        boost::RegEx rExFeature("^\\s*Feature\\s+(\\S.+)\\s*$");
        boost::RegEx rExFeatAndParams("^\\s*Feature\\s+(\\S.+\\S)\\s+->\\s+(\\S.+)\\s*$");
        boost::RegEx rExPerTagFeature("^\\s*PerTagFeature\\s+(\\S.+\\S)\\s+->\\s+(\\S.+)\\s*$");
        Feature *f = NULL;
        std::vector<std::string> headers;
        std::vector<std::string> descriptions;
        std::string header = line;
        bool multi_match = rExMultiFeature.Match(header);
        bool packed_match = rExPackedFeature.Match(header);

        if(multi_match) {
            std::string anchor = rExMultiFeature[2];
            vector<std::string> *collapsed_fnames = fe->multiFilterEntries(anchor);

            if(!collapsed_fnames)
                collapsed_fnames = multiFeatureEntries(anchor);

            if(!collapsed_fnames) {
                assert(0);
                throw EXCEPTION( PARSE_ERROR, std::string("unknown Object ") + anchor);
            }

            header = rExMultiFeature[1] + rExMultiFeature[3];

            if(!rExFeatureName.Match(header))
                throw EXCEPTION(PARSE_ERROR, header + " is malformed\n");

            multiFeatures[rExFeatureName[1]] = collapsed_fnames;

            for(int j = 0; j < collapsed_fnames->size(); j++) {
                std::ostringstream mfString;
                mfString << j;
                std::string myHeader = header;
                Utils::substitute(myHeader, "$$", mfString.str().c_str(), true);
                headers.push_back(myHeader);
                descriptions.push_back((*collapsed_fnames)[j]);
            }
        }
        else if(packed_match) {
            std::string anchor = rExPackedFeature[2];
            vector<std::string> *collapsed_fnames = fe->multiFilterEntries(anchor);

            if(!collapsed_fnames)
                collapsed_fnames = multiFeatureEntries(anchor);

            if(!collapsed_fnames) {
                assert(0);
                throw EXCEPTION( PARSE_ERROR, std::string("unknown Object ") + anchor);
            }

            header = rExPackedFeature[1] + rExPackedFeature[3];

            if(!rExFeatureName.Match(header))
                throw EXCEPTION(PARSE_ERROR, header + " is malformed\n");

            std::string myHeader, description;

            for(int j = 0; j < collapsed_fnames->size(); j++) {
                std::ostringstream mfString;
                std::string myAnchor = anchor;
                mfString << j;
                Utils::substitute(myAnchor, "$$", mfString.str().c_str(), true);
                myHeader += myAnchor + " ";
                description += (*collapsed_fnames)[j];
            }

            Utils::substitute(header, anchor, myHeader, true);
            headers.push_back(header);
            descriptions.push_back(description);
        }
        else {
            headers.push_back(line);
            descriptions.push_back(" ");
        }

        int i = 0;

        for( ; i < headers.size(); i++) {
            int paramsInd = -1;
            ::vector<std::string> featArgs;
            ::vector<std::string> featTags;

            bool sparseFeat = rExSparseFeature.Match(headers[i]);

            if(sparseFeat)
                headers[i] = rExSparseFeature[1];

            bool phaseFeat = rExPhaseFeature.Match(headers[i]);

            if(phaseFeat)
                headers[i] = rExPhaseFeature[1];

            bool pertag_match = rExPerTagFeature.Match(headers[i]);
            bool match = rExFeature.Match(headers[i]);
            bool param_match = rExFeatAndParams.Match(headers[i]);

            if(!param_match && !match && !pertag_match)
                break;

            std::string args(rExFeature[1]);
            featTags.push_back(rExFeature[1]);

            if(param_match) {
                paramsInd = fInd;
                args = std::string(rExFeatAndParams[1]);
                featTags.clear();
                featTags.push_back(rExFeatAndParams[2]);
            }
            else if(pertag_match) {
                paramsInd = fInd;
                args = std::string(rExPerTagFeature[1]);
                std::string stringTags = rExPerTagFeature[2];
                /*
                 * Substitute any occurrance of a predefined node or edge Set
                 */
                std::map<std::string, std::string>::iterator it = tagSets.begin();
                for( ; it != tagSets.end(); it++)
                    Utils::substitute(stringTags, it->first, it->second, true);

                featTags.clear();
                rExCommaDelim.Split(featTags, stringTags, 0, 100);
            }

            featArgs.clear();
            rExDelim.Split(featArgs, args, 0, 100);

            std::string featName = featArgs[0];

            for(int j = 0 ; j < featTags.size(); j++) {
                int offset = 0;

                if(j > 0) {
                    std::ostringstream added;
                    added << "_" << j;
                    featArgs[0] = featName + added.str();
                }

                f = featureFactory.Create(featArgs[1], fInd, paramsInd, featArgs, offset, fe, this);

                if(f == NULL || (unsigned)offset != featArgs.size()) {
                    assert(0);
                    throw EXCEPTION(BAD_FEATURE_HEADER, headers[i]);
                }

                fInd++;

                if(paramsInd >= 0) {
                    tieTagsToFeature(f, featTags[j]);
                    paramsInd = fInd;
                }

                std::string description = descriptions[i] + "\t" + featTags[j];
                setFeature(featArgs[0], description, f);

                if(phaseFeat)
                    f->makePhaseDependent();

                if(sparseFeat)
                    f->makeSparse();
            }

        }

        if(i && i != headers.size()) {
            assert(0);
            throw EXCEPTION( PARSE_ERROR, std::string("undefined Filter ") + headers[i]);
        }

        return (i == 0 ? false : true);

    }

    /**
     * Reads Feature Header definitions from input file stream
     * featStream
     */
    void FeatureEngine::readHeaderDefinitions(std::ifstream & featStream) {
        assert(featStream);
        Pair<int> pi;

        /*
         * Creating Feature Optimizer object
         */
        fo = new FeatureOptimizer(MAX_PRECOMP_ARRAYS);

        boost::RegEx rExDelim("\\s+");
        boost::RegEx rExComment("^\\s*\\#.*");
        boost::RegEx rExEndofFile("^\\s*\\/\\/\\s*$");
        boost::RegEx rExTagSet("^\\s*\\@define\\s+TagSet\\s+(\\S+)\\s*(\\S.+)$");
        boost::RegEx rExPreCompute("^\\s*PreComputeTo\\s+(\\d+)\\s+(\\d+)\\s+(\\S.+)$");
        boost::RegEx rExFreeze("^\\s*Freeze\\s+(\\S.+)\\s*$");
        boost::RegEx rExFeatureSet("^\\s*\\@define\\s+FeatureSet\\s+(\\S+)\\s+(\\S+)\\s*$");
        boost::RegEx rExEndFeatureSet("^\\s*\\@end\\s+FeatureSet\\s*$");

        int fInd = 0;
        std::string line;

        while(std::getline(featStream, line), !featStream.eof()) {
            //      cerr << "\"" << line << "\"\n";
            featureHeaders += line + std::string("\n");
            if(rExComment.Match(line)) // line was a comment
                continue;

            if(line.length() == 0)
                continue;

            if(rExEndofFile.Match(line))
                break;

            if(rExTagSet.Match(line)) {
                tagSets[rExTagSet[1]] = rExTagSet[2];
                continue;
            }

            if(rExFreeze.Match(line)) {
                if(!rExFreeze[1].compare("all")) {
                    for(int i = 0; i < features.size(); i++)
                        features[i]->freeze();
                }
                else {
                    getFeature(rExFreeze[1])->freeze();
                }
                continue;
            }

            if(rExPreCompute.Match(line)) {
                TPCEntry entry;
                int pceIndex;

                if(!sscanf(rExPreCompute[1].c_str(), "%d", &pceIndex))
                    assert(0);
                if(!sscanf(rExPreCompute[2].c_str(), "%d", &entry.arrPeriod))
                    assert(0);

                entry.collapseFrames = Utils::stringToBoolean(rExPreCompute[3]);
                entry.arrIndex = -1;

                pcEntries[pceIndex] = entry;
                continue;
            }

            if(rExFeatureSet.Match(line)) {
                pair<string, vector<string> > featSet;
                featSet.first = rExFeatureSet[2];
                featureSets[rExFeatureSet[1]] = featSet;

                while(std::getline(featStream, line), !featStream.eof()) {
                    featureHeaders += line + std::string("\n");

                    if(line.length() == 0 || rExComment.Match(line))
                        continue;

                    if(rExEndofFile.Match(line))
                        throw EXCEPTION(PARSE_ERROR, string("after ...") + featureHeaders.substr((int)Utils::max(0, featureHeaders.size()) - 200, 200));

                    if(rExEndFeatureSet.Match(line)) {
                        break;
                    }
                    featureSets[rExFeatureSet[1]].second.push_back(line);
                }

                continue;
            }

            bool isHeader = readFeatureHeader(fInd, line);

            if(!isHeader) {
                ::vector<std::string> featSetVect;
                rExDelim.Split(featSetVect, line, 0, 100);
                map<string, pair<string, vector<string> > >
                    ::iterator it = featureSets.find(featSetVect[0]);

                if(it != featureSets.end()) {
                    pair<string, vector<string> > &fs = it->second;
                    for(int i = 0; i < fs.second.size(); i++) {

                        std::string myline = fs.second[i];
                        Utils::substitute(myline, fs.first, featSetVect[1], true);

                        if(!readFeatureHeader(fInd, myline))
                            throw EXCEPTION(PARSE_ERROR, string("after ...") + featureHeaders.substr((int)Utils::max(0, featureHeaders.size()) - 200, 200));
                    }
                }
                else
                    throw EXCEPTION(PARSE_ERROR, string("after ...") + featureHeaders.substr((int)Utils::max(0, featureHeaders.size()) - 200, 200));
            }
        }
    }

    set<Feature *> FeatureEngine::syncEndFeatures() {
        set<Feature *> features;

        for(int i = 0; i < fsm->numParseNodes(); i++) {
            Node *node = fsm->node((TParseNode)i);
            if(!node->isSyncEnd())
                continue;

            Edge *edge = fsm->edge(node->id(), SYNC_END_STATE);
            vector<Feature *> & edgeFeats = this->fo->features(TR, edge->id());

            for(int j = 0; j < edgeFeats.size(); j++)
                features.insert(edgeFeats[j]);

        }

        return features;
    }

    set<Feature *> FeatureEngine::syncBegFeatures() {
        set<Feature *> features;

        for(int i = 0; i < fsm->numParseNodes(); i++) {
            Node *node = fsm->node((TParseNode)i);
            if(!node->isSyncBeg())
                continue;

            Edge *edge = fsm->edge(SYNC_BEG_STATE, node->id());
            vector<Feature *> & edgeFeats = this->fo->features(TR, edge->id());

            for(int j = 0; j < edgeFeats.size(); j++)
                features.insert(edgeFeats[j]);

        }

        return features;
    }


    FeatureEngine::~FeatureEngine() {
        unsigned int ind;
        deletePreComputedFeatures();

        for(ind = 0; ind < this->features.size(); ind++) {

            if(this->features[ind])
                delete this->features[ind];

            this->features[ind] = NULL;
        }

        features.clear();

        if(fo)
            delete fo;

    }

}
