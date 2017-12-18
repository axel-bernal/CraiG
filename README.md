# CraiG
CRAIG is a suite of tools that use an underlying semi markov CRF model for performing gene model learning and prediction on structured biological sequences. 
CRAIG stands for a CRf-based Ab-Initio Genefinder, however the suite can integrate output from third-party gene predictors, protein alignments, RNA-Seq features and denovo features into the predictive model. The predicted structures are genes, including UTRs, however the libraries are general enough to support any type of structure. The underlying models can be ab initio, de novo or ensemble.

Related Publications

 * Global Discriminative Training for Higher-Accuracy Computational Gene Prediction. Bernal A, Crammer K, Hatzigeorgiou A, Pereira F. PLoS Comput Biol 3(3):e54. 2007

 * Automated Gene-Model Curation using Global Discriminative Learning. Bioinformatics (2012) 28(12): 1571-1578

##  Getting Started
CRAIG's core executables and libraries are written in C++ and are: craigPredict and craigTrain, for predicting and learning models from input structures respectively.

There is a preprocessing script, craigPreprocess.py, that needs to be performed prior to train or predict structures in all cases. This script prepares  model parameters and organises and formats all evidence sources associated with the input sequences to facilitate training and/or prediction. Model parameters and learned gene models will be located in subdirectory CRAIG_HOME/models, if an output directory is not provided (see Subsection 2.2.e below to see when/how to setup this shell variable).

For automated whole-genome improvement of gene annotations, we have provided a processing pipeline  to conveniently preprocess, train and predict gene models given a genome and a set of existing  gene annotations (if any). This pipeline is described in detail in Section 5.
   

## Prerequisites
The following is a list of third party software installed with the main distribution:

### OS/Development level
 * python v2.7.3 or higher
 * gcc version v4.4.2 20091027 (Red Hat 4.4.2-7) (GCC)
 * libtools v2.2 or higher
 * perl v5.1 or higher
 * automake 1.9 or higher
 * autoconf 2.6 or higher
 * doxygen for generating the documentation

### Applications
 * google sparse hash and vector implementations
 * regtools (https://regtools.readthedocs.io/en/latest/) if using RNA-Seq BAMs in the input
 * boost regex libraries 
 * eval software package (only for performance evaluations)
 * gmap 2013-08-19 or later if gsnap is to be used for computing the RNA-Seq alignments

## Installing
 * Setup the environment (the autoconf tools are needed for doing this)
```
tar xvf craig-VERSION.tgz | gunzip
cd craig-VERSION && ./autogen.sh 
```
 * Configure environment
```
./configure --prefix=PREFIX_INSTALLATION [--enable-opt=yes|no] [--enable-mpi=yes|no]
```
    By default configure does not turn on debug information. If you 
      want to turn on this info, run 'configure' with command line 
      option --enable-opt=no
      
      The --enable-opt option turns debug information ON and also
      makes object files at least five times as large and the code 
      three to five times slower, so unless you are thinking in 
      debugging the program, just run configure with the default values.

      The installation script then will copy all the required data files and 
      executables to PREFIX_INSTALLATION, which is the prefix installation 
      directory.
           
      For the mpi version, the option --enable-mpi will install an mpi 
      compliant version of CRAIG. This version should be used when the input
      training data is large enough to take a few days for training (above .5
      Gb will do that). 

      In those cases, splitting the training data in subsets, train the 
      subsets separatedly and merge the parameters computed at the end for
      each subset will run a factor of N times faster, where N is the number
      of subsets. See the craigTrain help to know what options are available
      for spliting the training data and merging the resulting parameters in
      each case. Performance of the mpi vertion will vary but will usually 
      stay competitive when compared to the single processor version.

  * The environment variable CRAIG_HOME needs to be set permanently 
      to the root directory of the installation directory. To do this the
      .bashrc or .bash_profile files located in the $HOME directory need
       to be edited.
      
      To continue the installation this command (in bash only) may be run:
      
      export CRAIG_HOME=$PREFIX_INSTALLATION

      This is needed so that craigTrain and other applications know exactly 
      where to look for model parameters for training and learned gene models
       for predicting. See next section for more information on this issue. 
  * Run command:
      make
      
      This should build all objects files, libraries and executable binaries.

  * Skip this step if you don't have Doxygen installed in your system
      Optionally the following command could be run 
      
      make doc
      
      This command will generate documentation information and it will only 
      work if doxygen has been installed. See section 5 for more details and
      requirements.      

  * For installation, run commands:

      make install; make installcheck
      
      Root access might be needed if the installation tree  permissions
      require it. 
      This step will copy all binary executables in the bin directory and 
      liblless.so, the shared library containing all the lless library 
      rountines, to the lib directory.
      This step will also test the instalation to make sure the built 
      programs have no errors.

  *   The following directories will also need to be added to the PATH variable:
            $CRAIG_HOME/bin
      	    $CRAIG_HOME/perl/bin
      	    $CRAIG_HOME/python/bin

  * The python library numpy needs to be installed as it is not part of
      the python standard library. The comand "yum install numpy" would do
      this if yum is used as install manager

  * Add the directory $CRAIG_HOME/lib to the LD_LIBRARY_PATH environment
      variable. A command like this in the .bashrc or .bash_profile would
      do that:
      
      export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CRAIG_HOME/lib


##  Using the program(s)

### Preprocessing Stage
#### The craigPreprocess.py Script
      There are many types of transcription/translation evidence sources that can be integrated for learning/predicting gene models. This body of evidence needs to be formatted and organised before learning and predicting gene models. The pre-processing script craigPreprocess.py takes care of the formatting and organising in a transparent and user-friendly manner. 

Usage: craigPreprocess.py [OPTIONS] SPECIES ANNOT_FILE FASTA_FILE
      If specified, the option --config-file CONFIG_FILE makes the script to prepare and all the input information necessary to train a new gene model. The input information is then summarized and written in CONFIG_FILE using the following format:

        ValidationSequences FASTA_FILE_FOR_VALIDATION_SET
        ValidationTags TAGS_FOR_VALIDATION_SET
        Sequences FASTA_FILE_FOR_TRAINING_SET
        Tags Parsing TAGS_FOR_TRAINING_SET
        Path PATH_TO_WHERE_CONFIG_FILE_IS_LOCATED
        Name MODEL_NAME_PREFIX
	PrefixFiles EVIDENCE_PREFIX_FOR_TRAINING_SET
	PrefixGenomeFiles EVIDENCE_PREFIX_FOR_TESTING_SET

      EVIDENCE_PREFIX_FOR_TESTING_SET typically prefixes the input evidence for the whole genome. 
      The script craigPreprocess.py supports three different types of gene models -- specified using MODEL/MODEL_NAME_PREFIX above:
      a) Ab initio for models that use only intrinsic features (MODEL = craig, MODEL_NAME_PREFIX = SPECIES), 
      b) Ensemble or Evidence integration models that use alignments to external evidence of transcription/translation that is not RNA-Seq data (MODEL = ecraig, MODEL_NAME_PREFIX = SPECIES-evid)
      c) RNA-seq based models which require at least one evidence source to be of RNA-Seq type (MODEL = ngscraig, MODEL_NAME_PREFIX = SPECIES-rna). 
	
      Before learning a parameter model, the following files must be found in the $CRAIG_HOME/models directory: MODEL_NAME_PREFIX.resources, MODEL_NAME_PREFIX.filters, MODEL_NAME_PREFIX.features, MODEL_NAME_PREFIX.partial.top, MODEL_NAME_PREFIX.complete.top and any other external evidence file with prefix MODEL_NAME_PREFIX that the particular model requires as input. The latter files should be specified in file PRE_CONFIG_FILE, using option --pre-config-file PRE_CONFIG_FILE above.

      If any of the above files MODEL_NAME_PREFIX.* is missing in $CRAIG_HOME/models then the script uses CLOSEST_SPECIES ("generic" by default) instead of SPECIES to look for the same files. CLOSEST_SPECIES is assumed to be closely related to the target. The "generic" species has been designed to work relatively well in most eukaryotic cases, more so if RNA-Seq evidence is provided in the input.

      The file PRE_CONFIG_FILE contains information about all the external evidence sources that are available to the model. Each line in this file refers to one evidence source and has the following fields (tab separated):

      a)  Evidence Type: 'rnaseq' for RNA-Seq evidence sources, 'genepred' for external gene predictions or 'alignment' for blast/protein or EST alignments. The first one corresponds to RNA-Seq library mappings, the second one corresponds to external gene predictions and the third one to alignments to databases of transcripts and proteins. The main difference between the second and third type of evidence is that the third one do not generally provide with translation information but it provides scores that are meaningful, such as percentage of identity or p-value. The final gene model integrates all these evidence sources in an ensemble-type framework for learning.
      b) Evidence SubType: Available options are 'rum' and 'gsnap' for 'rnaseq' types of evidence and 'gtf', 'locs' and 'gff3' for the 'genepred' or 'alignment' type. This field clarifies what the source format is. 
      c) Dataset Id: An identifier for the sequence/annotation dataset referred to by all the sources of evidence.
      d) Sample Id: An identifier for the evidence source. Could be a compound name consisting of the author and stage(hour, day) in case of RNA-Seq or the program and running parameters for gene prediction/alignment.
      e) Full Path to Evidence Source: The full path to the file or directory containing the evidence information. 
      f) Qualifier: This field is either orientation for 'rnaseq' types of evidence (one of N|F|R|FR|RF|FF|RR) or the maximum confidence score for each prediction for 'genepred' types; a score of 0 means CRAIG should not use prediction confidence scores, only annotations.

      An example of how to specify a RNA-Seq evidence source follows:
      "rnaseq  gsnap   me49-9.0  brady_sibley.day0     /gpfs/fs121/h/abernal/GUS/project_home/DJob/gsnap_test/brad_sibley_input/master/mainresult      F"

#### Additional Notes
      There are a few available gene models that were included along the distribution. These models were carefully studied and their feature sets optimized to fit the genome organization in each case. These models are H. sapiens ab initio (human) and H. sapiens ensemble (human-evid), C. elegans ab initio (celegans) and C. elegans ensemble (celegans-evid), A. thaliana ab initio (athaliana) and A. thaliana ensemble (athaliana-evid), T. gondii ab initio(tgondii) and T. gondii RNA-Seq (tgondii-rna), and P. falciparum ab initio (pfalciparum). An automated method to select the closest SPECIES rather than having it as a parameter is feasible using mutual information measures, but it only makes sense when the pool of learned models is larger than what it is right now.

      There are also some restrictions on what types of evidence can be combined. For example, the RNA-Seq based models can freely combine RNA-Seq data with other types of evidence; however, each learnt model can only integrate a single RNA-Seq sample. The reason for this restriction is that the RNA-Seq evidence is stage-specific and different samples will contain contradictory transcriptional evidence about the same genomic locus. 

      An important requirement to learn gene models that integrate RNA-Seq information is that a suitable ab initio model for the target organism must be provided to obtain a good training data set out of the existing input annotations. This ab initio model is assumed to exist in path $CRAIG_HOME/models/SPECIES.params by default, where SPECIES is defined in the previous section. If the model does not exist, the program will check for $CRAIG_HOME/models/CLOSEST_SPECIES.params with CLOSEST_SPECIES is as defined in the previous section. The output from this ab initio model is used together with the set of input gene annotations and the existing RNA-Seq evidence to compute a training set for the final model

### Learning Gene Models
    The configuration file CONFIG_FILE, defined in the previous section, contains all the information needed to start the learning process. FASTA_FILE_FOR_VALIDATION_SET and  FASTA_FILE_FOR_TRAINING_SET are the file names of the input sequence files for training and validation respectively, they should be in fasta format. The files TAGS_FOR_VALIDATION_SET and TAGS_FOR_TRAINING_SET contain the Tag labelings which should be provided for each input sequence. These latter tag files can be computed from *gff3 or *gtf input annotation files by using a sequence of perl scripts provided in the $CRAIG_HOME/perl/bin. This sequence is performed automatically using the craigPreprocess.py script

   For more information related to this point one can refer to craigTrain's API obtained by executing craigTrain -h and/or section 4 in which a pipeline for improving annotation on whole genomes is described in detail.

### Predicting Gene Structures
    Executing craigPredict -h will display a detailed help on how to run the command. The main requirement is to have a trained gene model ready. 

    For predicting genes using ab initio models, having ready an ab initio model parameter file and the input DNA sequences to predict genes is sufficient.
    For predicting genes using models that integrate external evidence, all the transcriptional evidence for an input set of sequences need to be formatted and organised using craigPreprocess.py. The preprocessing step is explained in detail in 3.1. The option --prefix-evidence=PREFIX_EVIDENCE must be set with the appropiate value, typically specified in the configuration file as PrefixGenomeFiles when trying to re-annotate the whole genome (see Section 4).
  

### Range Limits of Important Input Parameters
    Maximum Number of Exons : = 2^8
    Maximum Number of Alternative Splices = 2^8
    Maximum Number of State Phases = 2^2
    Maximum Number of Strands = 2^1
    Maximum Number of States = 2^6
    Maximum Number of Transitions = 2^6


## Using RNA-Seq to Improve Whole-Genome Annotation


   With the advent of RNA-Seq data, there is an urgent need of revising/improving existing gene models most of which were curated, revised or predicted before the availability of RNA-Seq data and the wealth of transcriptional evidence that it provides.

   We have coded an automated pipeline to accomplish this task. The user only needs to provide a contig assembly, an existing annotation and a RNA-Seq dataset -- any other type of evidence could still be specified but would be optional. The output of the pipeline is two sets of gene reannotations: The first set only adds RNA-Seq derived UTR regions to the existing set of gene annotations and the second set corresponds to ngsCRAIG's full denovo gene structure predictions, including UTRs.

   The command to execute is craigReannotate.py. Option -h will provide a detailed information of the command's API. Though this command makes the most sense when RNA-Seq evidence is available (ngscraig), it can also be used to reannotate the genome using any type of model, even ab initio ones. The program craigReannotate.py calls craigPreprocess.py to preprocess the RNA-Seq data and craigTrain to train gene models. After the gene models are obtained, it uses them to predict the reannotations.

   After computing gene models for different RNA-Seq libraries/stages, a merge procedure is necessary to report all available gene predictions in a sound manner. This procedure is called using craigPostprocess.py. This script reports stage-specific alternative splicing events and selects either the best scoring (Default) or the longest UTR call for gene predictions with identical CDS.

## Examples

### Example for reannotating me49 using Hehl day7 RNA-Seq dataset
  * See examples/hehl_day7.preconf for an example of a preconfiguration file for Hehl's day7 dataset
  * See examples/craig_pipeline.sh for all the commands needed to generate annotations for Hehl's day7 dataset using examples/hehl_day7.preconf as preconfiguration file. See comments inside for requirements

### Example for reannotating a genome using Reid's day4 RNA-Seq data in a cluster where distribjob is available
  * Generate task.prop, controller.prop and distribjob directories. Run command buildDJobPropFiles4Craig.pl. The usage for this command is as follows:
For the example at hand the following command would suffice -- provided the RNA-Seq data has been aligned with rum (and not gsnap):

buildDJobPropFiles4Craig.pl --rsq-type rum --species tgondii --rsq-inputdir RNASEQ-DIR --rsq-orientation N --dataset-id tqz_reid --sample-id day4 --prefix-propdir PROP_DIR --prefix-craig-output CRAIGOUTPUT_DIR --fasta-file FASTA_FILE --annot-file GFF3_FILE

Comments: 
  - All File Paths above should be absolute paths, not relative
  - Make sure reads have been mapped and are present in directory RNASEQ-DIR
  - Also make sure that the output directories for --prefix-propdir --prefix-craig-output exist. 
  - The option rsq-orientation can have N,F,R, FR, RF values. The latter two values are for paired reads; R,F stand for the strand orientation, reverse or forward, and N stands for non-stranded samples.
  - Directories CRAIGOUTPUT_DIR.preproc CRAIGOUTPUT_DIR.model CRAIGOUTPUT_DIR.test  will contain the processing data, learned models/predictions and temporary data respectively.
  - If no model files exist for the target, the option --closet-species must be specified

  * Run the distribjob command:
cd PROP_DIR_input && nohup distribjob --propFile controller.prop --numNodes 1 --memoryPerNode 20.0;

Comments: 
  - The memory requirements of 20Gb are only needed for the initialization node only. Please do not use  a smaller size as this could result in insufficient memory problems.
  -  Make sure you are running this in the directory where controller.prop is located. 
  -  Also, make sure to erase the master subdirectory if something seriously wrong happened in the last distribjob run and a run from scratch is needed.

  * Obtaining the output predictions.
The paths for the generated gff3 files with the CRAIG gene models are 
 CRAIGOUTPUT_DIR.model/tgondii-rna.tqz_reid.day4.denovo.chr.gff3 for full CDS and UTR denovo predictions.
  * CRAIGOUTPUT_DIR.model/tgondii-rna.tqz_reid.day4.utr_only.chr.gff3 for UTR only prediction, i.e. only UTR regions are predicted, while the input annotation's CDS is preserved.

  * Running a postprocessing step for Reid's day 3 and day4 libraries. 
   Run post processing command
craigPostprocess.py --use-model-scores --out-dir SUMMARY --list-prefixes CRAIGOUTPUT_DIR_FOR_DAY4.model/tgondii-rna.tqz_reid.day4,CRAIGOUTPUT_DIR_FOR_DAY3.model/tgondii-rna.tqz_reid.day3 CRAIGOUTPUT_DIR_FOR_DAY4.preproc/tgondii-rna.chr.locs CRAIGOUTPUT_DIR_FOR_DAY4..preproc/tgondii-rna.chr.fa

Comments:
  - The path for the final summarized(merged output) w/o alt splicing is SUMMARY/unionized.final.gff3. This file contains the best possible gene model predictions that use all RNA-Seq libraries to predict UTRs and the input gene annotations to predict CDS. This file has no stage-specific alternative splicing however. The last field in the gff3 for mRNA entries contains detail about the experiment/sample id used as support evidence to predict the either UTR.

  - The path SUMMARY/denovo.final.gff3 contains the stage-specific alternative splicing denovo predictions. The last field in the gff3 for mRNA entries contains detail about the experiment/sample id used as support evidence to predict the gene model. Transcripts with UTR variations but otherwise identical CDS will be reported as one transcript only and the last gff3 field will also contain details about the UTR annotations for each sample-id involved.
   The list of prefixes for each RNA-Seq library can be specified in a file if there are many libraries involved.

## Class Documentation

Documentation for the library classes and other programs is provided in the doc directory. There are two files: refman.html.tar.gz and latex/refman.pdf which contain the reference manual documentation for craig and its supporting library llessin html and pdf format.
Doxygen was used to generate documentation from the source, in the style of Javadoc. The documentation is mostly complete, but it needs improvement. We plan to keep working on it in next versions of craig.
  You can execute make doc from the top directory of the installation to generate the doxygen documentation on your own. You need to install doxygen, pdflatex, tar and gzip utilities in order to do this.
 There is a small User's manual with a brief explanation of the main classes of lless and craig and give an overview on the how to build a training environment for a new organism. For more detailed information refer to the files locatd in the doc subdirectory.


## TODO
    The integration of RNA-Seq and ensemble-type evidence sources has not been fully tested. It does not occur in practice too often.
   The human RNA-Seq model has not been fully tested either. The chance never came to be. This is something that is almost ready, but there needs to be some real testing data.
