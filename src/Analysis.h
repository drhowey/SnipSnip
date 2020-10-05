/************************************************************************
 * SnipSnip, version 1.1
 * Copyright 2013
 * Richard Howey
 * Institute of Genetic Medicine, Newcastle University
 *
 * richard.howey@ncl.ac.uk
 * http://www.staff.ncl.ac.uk/richard.howey/
 *
 * This file is part of SnipSnip, the pedigree file processing program for EMIM.
 *
 * SnipSnip is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SnipSnip is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SnipSnip.  If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/


/*! \file Analysis.h
    \brief This file organises the analysis.
    
*/

#ifndef __ANALYSIS
#define __ANALYSIS

#include <string>

using namespace std; // initiates the "std" or "standard" namespace

#include "Data.h"
#include "Model.h"
#include "ModelFitting.h"

//! Returns a string of the run time.
string getTime(const double & t);

//! Class for storing three parameters.
struct threeParameters
{
	double beta0, beta1, beta2;

	threeParameters() : beta0(0), beta1(0), beta2(0) {};
	threeParameters(const double & v1, const double & v2, const double & v3) : beta0(v1), beta1(v2), beta2(v3) {};
};

//! Organises the correect analysis to perform.
class Analysis
{
private:

	//parameters for the options of the analysis to perform
	string filename;
	string outputFileName;
	unsigned int windowBPSize;  //actually half the size of the window, bps to either left or right of tested SNP
	unsigned int windowSNPSize; //actually half the size of the window, number of SNPs to either left or right of tested SNP
	bool partnerResults;
	unsigned int startSNP;
	unsigned int endSNP;
	unsigned int metricType;
	unsigned int modelType;
	bool linear;
	bool lr;	
	SNPWindow * snpWindow;
	Model2SNPs * model2SNPs;
	Model1SNPLogReg * model1SNPLogReg;
	CovariateData * covariateData;
	bool useCovariates;
	unsigned int noCovariates;
	ofstream resultsFile;
	
	//used for fitting models, previous fitted parameters used as starting values
	double prevH0Beta0;
	double prevM2ABeta0, prevM2ABeta1;
	double prevM2PBeta0, prevM2PBeta2;
	double prevM3Beta0, prevM3Beta1, prevM3Beta2;
	list<double> prevCovarH0;
	list<double> prevCovarM2A; 
	list<double> prevCovarM2P;
	list<double> prevCovarM3;

	map<unsigned int, pair<double, double> > initialParasM2;
	map<unsigned int, threeParameters > initialParasM3;
	map<unsigned int, list<double> > initialParasCovarM2;
	map<unsigned int, list<double> > initialParasCovarM3;

	//used for reading in data
	char oneBuffer[1];
	int aBit;

	static const double oneOverSqRoot2;

public:

	Analysis(string & fn, string & ofn, unsigned int & ws, unsigned int & wss, unsigned int & mt, 
		 unsigned int & ss, unsigned int & es, unsigned int & exsnp, bool & lin, double & missingQTValue, bool & l, string & cvfn, string & cvs)
		:  filename(fn), outputFileName(ofn), windowBPSize(ws/2),  windowSNPSize(wss/2), metricType(mt), startSNP(ss), endSNP(es),  linear(lin), lr(l)
	{	
		if(windowSNPSize!=0) windowSNPSize++; //make window SNP size be the number of PPs

		//setup the model to use depending on whether linear regression is used or not
		if(linear) model2SNPs = new TwoSNPLinearRegModel();
		else model2SNPs = new TwoSNPLogRegModel();			

		model1SNPLogReg = new Model1SNPLogReg();
		
		ScorePotentialPartners * scorePotentialPartners;
		if(metricType == 1)	scorePotentialPartners = new PPScoreCorrMulti(model2SNPs);		
		else if(metricType == 2) scorePotentialPartners = new PPScoreCorrDom(model2SNPs);
		else scorePotentialPartners = new PPScoreCorrReces(model2SNPs);
		
		if(linear)
		{
			QuantitiveTraits * quantitiveTraits = new QuantitiveTraits();
			scorePotentialPartners->setupQuantitiveTraits(quantitiveTraits, filename, missingQTValue);		
			model2SNPs->setQuantitiveTraits(quantitiveTraits, missingQTValue);
		}
		else scorePotentialPartners->setupCaseControls(filename);

		if(cvfn != "")
		{			
			covariateData = new CovariateData(cvfn, cvs, filename, missingQTValue);
			useCovariates = true;
			noCovariates = covariateData->covariateDataAllSubjects.begin()->size();
		}
		else
		{
			covariateData = new CovariateData();
			useCovariates = false;
			noCovariates = 0;
		};

		scorePotentialPartners->updateCaseControlWithMissing(covariateData->caseControls, linear);

		if(cvfn != "" || linear) model2SNPs->setCaseControls(scorePotentialPartners->getCaseControls());
		if(cvfn != "" && lr && !linear) model1SNPLogReg->setCaseControls(scorePotentialPartners->getCaseControls());		

		//set covariate data for fitting models
		model2SNPs->setCovariateData(covariateData, useCovariates);
		if(lr && !linear) model1SNPLogReg->setCovariateData(covariateData, useCovariates);

		snpWindow = new SNPWindow(filename, ws, windowSNPSize, startSNP, scorePotentialPartners, exsnp, linear, lr); 
		resultsFile.unsetf(ios::floatfield);            
		resultsFile.precision(10);
	
		setUpInitialParameters();	
	};

	//! Delete analysis things
	~Analysis()
	{		
		delete snpWindow;
		delete model2SNPs;
		delete model1SNPLogReg;
		delete covariateData;
	};

	void runAnalysis();	
	void fitModels(const unsigned int & snpID);
	void fitModelsLR(const unsigned int & snpID); //for stats with anchor + (partner - anchor)
	void fitModelsLinear(const unsigned int & snpID);
	void outputMissingResultsLine(const unsigned int & snpID, const unsigned int & chr, const string & anchorName, const unsigned int & anchorBP);
	void recordSNPnotAnalysed(const unsigned int & snpNo, unsigned int & noSNPsNotAnalysed, const string & snpName, ofstream & notAnalSNPsFile, const double & reason);
	void setUpInitialParameters();
	void updateM2AInitialParameters();
	void updateM2PInitialParameters();
	void updateM3InitialParameters();
	void updateCovarInitialParametersM2A();
	void updateCovarInitialParametersM2P();
	void updateCovarInitialParametersM3();
	void updateCovarParameters(const unsigned int & parametersNo, const unsigned int & modelNo);
	void updatePrevCovarParametersM2A();
	void updatePrevCovarParametersM2P();
	void updatePrevCovarParametersM3();

	unsigned int getNextNoOfMinorAlleles(ifstream & readSNPData, unsigned int & bitCount);
	double getPvalueChiSq1DF(double & chisq);
	double getPvalueFStat(double & fstat, const unsigned int & p, const unsigned int & q);
};


#endif
