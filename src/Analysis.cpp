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


/*! \file Analysis.cpp
    \brief This file contains the methods the various analyse.
    
*/
#include <iostream>
#include <sstream>
#include <math.h>

using namespace std; // initiates the "std" or "standard" namespace

#include "Analysis.h"
#include "Model.h"
#include "ModelFitting.h"
#include "main.h"
#include "cdflib.h"

const double Analysis::oneOverSqRoot2 = 0.7071067811865475244008443621048490392848359376884740365883;

//! Converts an integer to a string
string toString(int & i)
{
	ostringstream aStringStream;
	aStringStream << i;

	return aStringStream.str();
};

//! Returns a string of the run time
string getTime(const double & t)
{
	double time = t;
	int days = 0;
	int hours = 0;
	int minutes = 0;
	int seconds = 0;

	string ans = "";
	days = (int) (time / 86400); time -= days*86400;
	hours = (int) (time / 3600); time -= hours*3600;
	minutes = (int) (time / 60); time -= minutes*60;
	seconds = (int) time;

	if(days == 1) ans += "1 day";
	else if(days > 0) ans += toString(days) + " days";

	if(hours > 0)
	{
		if(days != 0)
		{
			if(minutes == 0 && seconds == 0) ans += " and ";
			else ans += ", ";
		};

		if(hours == 1) ans += "1 hour";
		else ans += toString(hours) + " hours";
	};

	if(minutes > 0)
	{
		if(ans != "")
		{
			if(seconds == 0) ans += " and ";
			else ans += ", ";
		};

		if(minutes == 1) ans += "1 minute";
		else ans += toString(minutes) + " minutes";
	};

	if(seconds > 0)
	{
		if(ans != "")
		{
			ans += " and ";			
		};

		if(seconds == 1) ans += "1 second";
		else ans += toString(seconds) + " seconds";
	};

	if(ans == "") ans = "less than one second";

	return ans;
};

//! Gets the log filename from the results file name
string getSNPNotAnalFileName(string & outFileName)
{
	unsigned int filenameLength = outFileName.length();
	string notAnalFileName;

	//find extension
	unsigned int a = filenameLength - 1;
	while(a > 0)
	{
		if(outFileName.substr(a, 1) == ".") break;
		a--;
	};

	if(a > 0) notAnalFileName = outFileName.substr(0, a) + "-SNPs-not-anal.txt";
	else notAnalFileName = outFileName + "-SNPs-not-anal.txt";

	return notAnalFileName;
};

//! Records SNPs not analysed due to missing data or no suitable partner SNP found.
void Analysis::recordSNPnotAnalysed(const unsigned int & snpNo, unsigned int & noSNPsNotAnalysed, const string & snpName, ofstream & notAnalSNPsFile, const double & reason)
{
	if(noSNPsNotAnalysed == 0)
	{
		string notAnalFileName = getSNPNotAnalFileName(outputFileName);
		notAnalSNPsFile.open(notAnalFileName.c_str());
		notAnalSNPsFile << "POS_IN_FILE SNP REASON\n";				
	};

	noSNPsNotAnalysed++;

	notAnalSNPsFile << snpNo << " " << snpName << " ";
	if(reason == -1) notAnalSNPsFile << "No_data_for_SNP";
	else if(reason == -2) notAnalSNPsFile << "No_case_data_for_SNP";
	else if(reason == -3) notAnalSNPsFile << "No_control_data_for_SNP";
	else if(reason == -4) notAnalSNPsFile << "No_suitable_partner_SNP_found";
	else if(reason == -5) notAnalSNPsFile << "SNP_window_size_of_0";
	notAnalSNPsFile << "\n";

	outputMissingResultsLine(snpNo, snpWindow->getAnchorChromosome(), snpWindow->getAnchorSNPName(), snpWindow->getAnchorSNPBP());
};

//! Runs the chosen analysis.
void Analysis::runAnalysis()
{
	unsigned int noSNPsNotAnalysed = 0;	
	ofstream notAnalSNPsFile;
	
	resultsFile.open(outputFileName.c_str());
	
	string test = "CHISQ";
	if(linear) test = "FSTAT";

	//output header line for results file
	resultsFile << "SNP CHR ID BP PARTNER_ID PARTNER_BP CORRELATION SCORE FIT_STATUS "<<test<<" P";
	if(lr) resultsFile << " FIT_STATUS_LR "<<test<<"_LR P_LR";	
	resultsFile << "\n";

	unsigned int noSNPs = snpWindow->getNoSNPs();

	//loop thro' all the SNPs updating the SNP window as we go along
	unsigned int snpID = startSNP;

	//initial default parameters
	map<unsigned int, double> parameters;
	parameters[1] = 0.05; //beta0
	parameters[2] = 0.05; //beta1	
	parameters[3] = 0; //beta2

	model2SNPs->setNewParameters(parameters);	

	//set a parameter for correlation of last partner SNP	
	model2SNPs->setParameter(4, snpWindow->getLastPartnerScore());	

	//perform analysis for the first window - if a partner SNP was found
	if(snpWindow->getLastPartnerScore() > -1) fitModels(snpID); //all special numbers are less or equal to -1
	else recordSNPnotAnalysed(snpID, noSNPsNotAnalysed, snpWindow->getAnchorSNPName(), notAnalSNPsFile, snpWindow->getLastPartnerScore());
		
	++snpID;

	while(snpID <= noSNPs && (endSNP == 0 || snpID <= endSNP))
	{
		//move to next window and calculate LD between start SNP and other SNPs and set the partner SNP
		//and set the joint genotype counts, if no partner SNP is found skip model fitting
		if(snpWindow->nextWindow())
		{
			//set a parameter for score of last partner SNP			
			model2SNPs->setParameter(4, snpWindow->getLastPartnerScore());

			//perform analysis for the next window
			if(snpWindow->getLastPartnerScore() > -1) fitModels(snpID); 
		}
		else recordSNPnotAnalysed(snpID, noSNPsNotAnalysed, snpWindow->getAnchorSNPName(), notAnalSNPsFile, snpWindow->getLastPartnerScore());

		++snpID;
	};

	resultsFile.close();
	if(noSNPsNotAnalysed != 0) notAnalSNPsFile.close();
	snpWindow->displayWindowStats();
	if(noSNPsNotAnalysed > 0) {out("Number of SNPs not analysed: "); out(noSNPsNotAnalysed); out(" (SNPs in file: "); out(getSNPNotAnalFileName(outputFileName)); out(")\n\n");};
	
};

//! Outputs for any SNP that is not analysed for whatever reason.
void Analysis::outputMissingResultsLine(const unsigned int & snpID, const unsigned int & chr, const string & anchorName, const unsigned int & anchorBP)
{
	resultsFile << snpID << " " << chr << " " << anchorName << " " << anchorBP << " NA NA";
	resultsFile <<" NA NA NA NA NA";
	if(lr) resultsFile << " NA NA NA";
	resultsFile << "\n";
};

unsigned int Analysis::getNextNoOfMinorAlleles(ifstream & readSNPData, unsigned int & bitCount)
{
	int allele1, allele2;
	unsigned int noMinorAlleles = 0;
	int one = '\1';
		
	//read in the next piece of data
	if(bitCount == 9)
	{
		
		readSNPData.read(oneBuffer, 1);
		if(readSNPData.eof())
		{			
			outErr("Error: binary SNP file (.bed) is incomplete!\n");
			exit(0);
		};
		
		aBit = oneBuffer[0];
			
		bitCount = 1;
	};

	allele1 = aBit & one; //read the least significant bit				
	aBit = aBit >> 1; //shift bits to the right
	allele2 = aBit & one; //read the new least significant bit				
	aBit = aBit >> 1; //shift bits to the right for next time

	bitCount += 2;	

	//if genotype is encoded 1/0 then the genotype is missing so do not add it
	if(allele1 == 1 && allele2 == 1)
	{	
		noMinorAlleles = 0;
	}
	else if(allele1 == 0 && allele2 == 1)
	{	
		noMinorAlleles = 1;
	}
	else if(allele1 == 0 && allele2 == 0)
	{	
		noMinorAlleles = 2;
	}
	else
		noMinorAlleles = 3; //denotes missing genotype

	return noMinorAlleles;
};


void Analysis::updateM2AInitialParameters()
{
	if(useCovariates) updateCovarInitialParametersM2A();

	initialParasM2[1] = make_pair(prevM2PBeta0, prevM2PBeta2);
	initialParasM2[2] = make_pair(prevM2ABeta0, prevM2ABeta1);
	initialParasM2[3] = make_pair(prevH0Beta0, 0);	
};

//! Sets up initial parameters for model M2P.
void Analysis::updateM2PInitialParameters()
{	
	if(useCovariates) updateCovarInitialParametersM2P();

	initialParasM2[1] = make_pair(prevM2PBeta0, prevM2PBeta2);
	if(lr)
	{
		initialParasM2[2] = make_pair(prevM2ABeta0, prevM2ABeta1);
		initialParasM2[3] = make_pair(prevH0Beta0, 0);
	}
	else
	{
		initialParasM2[2] = make_pair(prevM2PBeta0, 0);
		initialParasM2[3] = make_pair(0, 0);
	};

};

//! Sets up initial parameters model M3.
void Analysis::updateM3InitialParameters()
{
	if(useCovariates) updateCovarInitialParametersM3();

	if(lr)
	{
		initialParasM3[1] = threeParameters(prevM2PBeta0, prevM2ABeta1, prevM2PBeta2);	
		initialParasM3[2] = threeParameters(prevM3Beta0, prevM3Beta1, prevM3Beta2);
		initialParasM3[3] = threeParameters(prevM2ABeta0, prevM2ABeta1, prevM2PBeta2);
	}
	else
	{		
		initialParasM3[1] = threeParameters(prevM3Beta0, prevM3Beta1, prevM3Beta2);
		initialParasM3[2] = threeParameters(prevM2PBeta0, prevM2PBeta2, prevM2PBeta2);	
		initialParasM3[3] = threeParameters(0, 0, 0);
	};

	for(unsigned int i = 1; i <= 25; ++i) initialParasM3[4+i] = threeParameters(prevM2PBeta0, -0.1*i, prevM2PBeta2);
	for(unsigned int j = 1; j <= 25; ++j) initialParasM3[29+j] = threeParameters(prevM2PBeta0, 0.1*j, prevM2PBeta2);
};

//! Sets up initial parameters model M2A.
void Analysis::updateCovarInitialParametersM2A()
{
	initialParasCovarM2[1] = prevCovarM2P;
	initialParasCovarM2[2] = prevCovarM2A;
};

//! Sets up initial parameters model M2P.
void Analysis::updateCovarInitialParametersM2P()
{
	if(lr)
	{
		initialParasCovarM2[1] = prevCovarM2A;
		initialParasCovarM2[2] = prevCovarM2P;
	}
	else
	{
		initialParasCovarM2[1] = prevCovarM2P;
	};
};

//! Sets up initial parameters model M3.
void Analysis::updateCovarInitialParametersM3()
{
	if(lr)
	{
		initialParasCovarM3[1] = prevCovarM2P;
		initialParasCovarM3[2] = prevCovarM2A;
		initialParasCovarM3[3] = prevCovarM3;
	}
	else
	{
		initialParasCovarM3[1] = prevCovarM2P;
		initialParasCovarM3[2] = prevCovarM3;
	};
	
};

//! Updates covariate parameters.
void Analysis::updateCovarParameters(const unsigned int & parametersNo, const unsigned int & modelNo)
{
	unsigned int uc = 1;
	map<unsigned int, list<double> > initialParasCovar;
	map<unsigned int, list<double> >::const_iterator covars;

	if(modelNo == 2) initialParasCovar = initialParasCovarM2;
	else if(modelNo == 3) initialParasCovar = initialParasCovarM3;
	else if(modelNo == 1) initialParasCovar = initialParasCovarM2; //use model 1 for single SNP log reg

	covars = initialParasCovar.find(parametersNo);
	
	if(covars != initialParasCovar.end())
	{
		for(list<double>::const_iterator ip = covars->second.begin(); ip != covars->second.end(); ++ip, ++uc)
		{
			if(modelNo != 1) model2SNPs->setParameter(4+uc, *ip);
			else model1SNPLogReg->setParameter(4+uc, *ip);
		};
	}
	else
	{
		covars = initialParasCovar.begin();
		//set covariates paras to 0 otherwise
		for(list<double>::const_iterator ip = covars->second.begin(); ip != covars->second.end(); ++ip, ++uc)
		{
			if(modelNo != 1) model2SNPs->setParameter(4+uc, 0);
			else model1SNPLogReg->setParameter(4+uc, *ip);
		};
	};			
};

//! Updates the previous covariate parameters for M2A.
void Analysis::updatePrevCovarParametersM2A()
{
	unsigned int uc = 1; 
	for(list<double>::iterator pc = prevCovarM2A.begin(); pc != prevCovarM2A.end(); ++pc, ++uc) *pc = model2SNPs->getParameter(4+uc); 
};

//! Updates the previous covariate parameters for M2P.
void Analysis::updatePrevCovarParametersM2P()
{
	unsigned int uc = 1; 
	for(list<double>::iterator pc = prevCovarM2P.begin(); pc != prevCovarM2P.end(); ++pc, ++uc) *pc = model2SNPs->getParameter(4+uc); 
};

//! Updates the previous covariate parameters for M3.
void Analysis::updatePrevCovarParametersM3()
{
	unsigned int uc = 1; 
	for(list<double>::iterator pc = prevCovarM3.begin(); pc != prevCovarM3.end(); ++pc, ++uc) *pc = model2SNPs->getParameter(4+uc); 
};

//! Sets up initial parameters.
void Analysis::setUpInitialParameters()
{
	//set initial start values of parameters
	//setup values for fitting model M2A and M2P
	prevH0Beta0 = 0;
	prevM2ABeta0 = 0; prevM2ABeta1 = 0;
	prevM2PBeta0 = 0; prevM2PBeta2 = 0;
	prevM3Beta0 = 0; prevM3Beta1 = 0; prevM3Beta2 = 0;
	
	initialParasM2[5] = make_pair(0, 0);

	//setup starting values for fitting model M3	
	initialParasM3[4] = threeParameters(0, 0, 0);

	//set covariate parameters
	if(useCovariates)
	{
		unsigned int noCovars = covariateData->covariateDataAllSubjects.begin()->size();
		list<double> paraVals;
		for(unsigned int uc = 1; uc <= noCovars; ++uc)
		{
			paraVals.push_back(0);
			prevCovarH0.push_back(0);
			prevCovarM2A.push_back(0);
			prevCovarM2P.push_back(0);
			prevCovarM3.push_back(0);
		};			
		for(unsigned int i = 2; i <= 2; ++i)
		{
			initialParasCovarM2[i] = paraVals;
		};
	};

};

//! Calculates the p-value from a Chi square value with 1 df.
double Analysis::getPvalueChiSq1DF(double & chisq)
{	
	double a = sqrt(chisq)*oneOverSqRoot2;
	int ind = 0;
	return erfc1(&ind, &a);
};

//! Calculates the p-value from a f-statistic with d1 and d2 dfs.
double Analysis::getPvalueFStat(double & fstat, const unsigned int & d1, const unsigned int & d2)
{	
	if(fstat <= 0) return 1.0;
	double x = (double)(d1*fstat)/(double)(d1*fstat + d2);
	double y = 1 - x;
	double a = (double)(d1)/2.0;
	double b = (double)(d2)/2.0;
	double w, w1;
	int ierr;

	bratio(&a, &b, &x, &y, &w, &w1, &ierr);

	return w1;
};

void Analysis::fitModels(const unsigned int & snpID)
{
	if(useCovariates || linear)
	{
		//set SNP data for Anchor and Partner	
		model2SNPs->setAnchorData(snpWindow->getAnchorSNPData());
		model2SNPs->setPartnerData(snpWindow->getPartnerSNPData());
		if(lr) model1SNPLogReg->setAnchorData(snpWindow->getAnchorSNPData());
	};

	if(linear) fitModelsLinear(snpID);
	else fitModelsLR(snpID);
};

//! Fit the models and compare them.
void Analysis::fitModelsLR(const unsigned int & snpID)
{
	FindFit findFit(model2SNPs); 	

	double negLogBeta0 = 0, negLogBeta0_1 = 0, negLogBeta0_2 = 0, negLogBeta0_1_2 = 0;
	map<unsigned int, double> parameters;
	set<unsigned int> parasToFit;

	//used for outputting details of interesting SNPs using -oncc
	double H0Beta0 = -1, M2ABeta0 = -1, M2ABeta1 = -1, M2PBeta0 = -1, M2PBeta2 = -1,
		M3Beta0 = -1, M3Beta1 = -1, M3Beta2 = -1;
	int H0Fit = 1, M2AFit = 0, M2PFit = 0, M3Fit = 0;

	//no SNPs
	parameters[1] =  prevH0Beta0; //beta0
	parameters[2] = 0; //beta1	
	parameters[3] = 0; //beta2	

	//set up covariate parameters, parameter 5 onwards
	if(useCovariates)
	{
		unsigned int noCovars = covariateData->covariateDataAllSubjects.begin()->size();
		for(unsigned int uc = 1; uc <= noCovars; ++uc)
		{
			parameters[4+uc] = 0; 
			parasToFit.insert(4+uc); 
		};
	};

	model2SNPs->setNewParameters(parameters);

	//check data is ok for fitting further models	
	bool anchorDataOK = model2SNPs->checkAnchorSNPData();
	bool partnerDataOK = model2SNPs->checkPartnerSNPData();
	bool anchorPartnerDiffOK = model2SNPs->checkAnchorDiffPartnerSNPData();

	string fitStatusM2A = "Y";
	string fitStatusM2P = "Y";
	string fitStatusM3 = "Y";

	bool fittedOK = false;
	
	//fit standard logistic regression using all available Anchor SNP data
	if(lr)
	{
		model1SNPLogReg->setNewParameters(parameters);
		//fit the model	H0 for beta0 only - using all Anchor data
		FindFit findFitLR(model1SNPLogReg);	
		pair<unsigned int, unsigned int> noCasesControls = snpWindow->getNoCasesControlsAnchor(model1SNPLogReg);

		if(useCovariates)
		{
			parasToFit.insert(1); //beta0

			//set parameters
			model1SNPLogReg->setParameter(1, prevH0Beta0);							

			//set covariate parameters
			//updateCovarParameters(par->first, 2);	

			//fit the model	
			fittedOK = findFitLR.newtonsMethod(negLogBeta0, parasToFit);
			fittedOK = fittedOK && negLogBeta0*0 == 0;

			prevH0Beta0 = model1SNPLogReg->getParameter(1);
		}
		else
		{
			
			pair<double, double> negLogAndBeta0 = model1SNPLogReg->getNegLogBeta0(noCasesControls.first, noCasesControls.second);
			negLogBeta0 = negLogAndBeta0.first;
			prevH0Beta0 = negLogAndBeta0.second;
			fittedOK = true;
		};

		if(noCasesControls.first != 0 && noCasesControls.second != 0)
		{
			if(fittedOK)
			{
				fittedOK = false;				
	
				//fit model M2A
				updateM2AInitialParameters();

				parasToFit.insert(1); //beta0	
				parasToFit.insert(2); //beta1

				//try fitting model M2A with different initial values
				for(map<unsigned int, pair<double, double> >::const_iterator par = initialParasM2.begin(); !fittedOK && par != initialParasM2.end(); ++par)
				{
					//set parameters
					model1SNPLogReg->setParameter(1, par->second.first);			
					model1SNPLogReg->setParameter(2, par->second.second);			

					//set covariate parameters
					if(useCovariates) updateCovarParameters(par->first, 2);	

					//fit the model	
					fittedOK = findFitLR.newtonsMethod(negLogBeta0_1, parasToFit);
					fittedOK = fittedOK && negLogBeta0_1*0 == 0 && negLogBeta0_1 <= (negLogBeta0 + 1e-6);
					M2AFit = par->first;
				};
			};

			if(!fittedOK)
			{
				fitStatusM2A = "N";
				M2AFit = -4;
			};

			//if data is poor then a fit may not be found, and there will be no sign. improvement
			if(fittedOK) 
			{
				prevM2ABeta0 = model1SNPLogReg->getParameter(1);
				prevM2ABeta1 = model1SNPLogReg->getParameter(2);
				if(useCovariates) updatePrevCovarParametersM2A();
			};	

			if(negLogBeta0_1 > negLogBeta0) negLogBeta0_1 = negLogBeta0;

		}
		else
		{
			fittedOK = false;
			fitStatusM2A = "D";
			M2AFit = -1;
		};	
	};

	//fit model M2P
	fittedOK = false;	

	if(partnerDataOK)
	{
		updateM2PInitialParameters();		

		//Remove Anchor SNP from model		
		model2SNPs->setParameter(2, 0);

		parasToFit.insert(1); //beta0
		parasToFit.erase(2); //beta1
		parasToFit.insert(3); //beta2
		
		//try fitting model M2P with different initial values		
		for(map<unsigned int, pair<double, double> >::const_iterator par = initialParasM2.begin(); !fittedOK && par != initialParasM2.end(); ++par)
		{
			//set parameters
			model2SNPs->setParameter(1, par->second.first);			
			model2SNPs->setParameter(3, par->second.second);			

			//set covariate parameters
			if(useCovariates) updateCovarParameters(par->first, 2);			

			//fit the model	
			fittedOK = findFit.newtonsMethod(negLogBeta0_2, parasToFit);
			fittedOK = fittedOK && negLogBeta0_2*0 == 0;

			M2PFit = par->first;

		};
		
		if(!fittedOK)
		{
			fitStatusM2P = "N";
			M2PFit = -4;
		};

		//if data is poor then a fit may not be found, and there will be no sign. improvement
		if(fittedOK) 		
		{
			prevM2PBeta0 = model2SNPs->getParameter(1);
			prevM2PBeta2 = model2SNPs->getParameter(3);
			if(useCovariates) updatePrevCovarParametersM2P();
		};

	}
	else
	{
		//no minor allele data for partner SNP so the best fit model is as for no SNPs
		negLogBeta0_2 = negLogBeta0;

		fitStatusM2P = "D";
		M2PFit = -1;
	};

	//fit model M3
	fittedOK = false;
	
	if(anchorDataOK && partnerDataOK && anchorPartnerDiffOK)
	{
		updateM3InitialParameters();		

		parasToFit.insert(2); //beta1		

		//try fitting model M3 with different initial values
		for(map<unsigned int, threeParameters>::const_iterator par = initialParasM3.begin(); !fittedOK && par != initialParasM3.end(); ++par)
		{
			//set parameters
			model2SNPs->setParameter(1, par->second.beta0);			
			model2SNPs->setParameter(2, par->second.beta1);	
			model2SNPs->setParameter(3, par->second.beta2);

			//set covariate parameters			
			if(useCovariates) updateCovarParameters(par->first, 3);

			//fit the model	
			fittedOK = findFit.newtonsMethod(negLogBeta0_1_2, parasToFit);
			fittedOK = fittedOK && negLogBeta0_1_2*0 == 0 && negLogBeta0_1_2 <= (negLogBeta0_2 + 1e-6);

			M3Fit = par->first;

		};		

		if(!fittedOK)
		{
			fitStatusM3 = "N";
			M3Fit = -4;
		};		
	}
	else
	{
		fitStatusM3 = "D";
		M3Fit = -1;
	};

	//if data is poor then a fit may not be found, and there will be no sign. improvement
	if(fittedOK) 
	{
		prevM3Beta0 = model2SNPs->getParameter(1);
		prevM3Beta1 = model2SNPs->getParameter(2);
		prevM3Beta2 = model2SNPs->getParameter(3);
		if(useCovariates) updatePrevCovarParametersM3();
	};

	if(negLogBeta0_1_2 > negLogBeta0_2) negLogBeta0_1_2 = negLogBeta0_2;

	//Output results
	resultsFile << snpID << " " << snpWindow->getAnchorChromosome() << " " << snpWindow->getAnchorSNPName() << " " << snpWindow->getAnchorSNPBP() << " " <<snpWindow->getPartnerSNPName() << " " << snpWindow->getPartnerSNPBP();
	resultsFile << " " << snpWindow->getLastAnchorPartnerCorr() << " " << model2SNPs->getParameter(4);

	//Output standard SnipSnip test result
	if(fitStatusM2P == "Y" && fitStatusM3 =="Y")
	{
		double chisq = 2*(negLogBeta0_2 - negLogBeta0_1_2);
		resultsFile << " Y " << chisq << " " << getPvalueChiSq1DF(chisq);
	}
	else if(fitStatusM2P == "D" || fitStatusM3 == "D") resultsFile << " D NA NA";
	else resultsFile << " N NA NA";

	//Output standard logistic regression test result
	if(lr)
	{
		if(fitStatusM2A == "Y")
		{
			double chisqLR = 2*(negLogBeta0 - negLogBeta0_1);
			resultsFile << " Y " << chisqLR << " " << getPvalueChiSq1DF(chisqLR);
		}
		else if(fitStatusM2A == "D") resultsFile << " D NA NA";
		else resultsFile << " N NA NA";
	};

	resultsFile << "\n";
};

void Analysis::fitModelsLinear(const unsigned int & snpID)
{
	//get SNP data for Anchor and Partner
	SNPData * anchorSNPData = snpWindow->getAnchorSNPData();
	SNPData * partnerSNPData = snpWindow->getPartnerSNPData();

	double rss0 = 0, rss0_1 = 0, rss0_2 = 0, rss0_1_2 = 0, rss0_1_2_3 = 0;//
	map<unsigned int, double> parameters;
	set<unsigned int> parasToFit;

	//used for outputting details of interesting SNPs using -oncc
	double H0Beta0 = -1, M2ABeta0 = -1, M2ABeta1 = -1, M2PBeta0 = -1, M2PBeta2 = -1,
		M3Beta0 = -1, M3Beta1 = -1, M3Beta2 = -1;
	int H0Fit = 1, M2AFit = 0, M2PFit = 0, M3Fit = 0;

	//no SNPs
	parameters[1] = prevH0Beta0; //beta0
	parameters[2] = 0; //beta1	
	parameters[3] = 0; //beta2

	model2SNPs->setNewParameters(parameters);

	parasToFit.insert(1); //beta0	
		
	prevH0Beta0 = model2SNPs->getParameter(1);		

	//check data is ok for fitting further models	
	string fitStatusM2A = "Y";
	string fitStatusM2P = "Y";
	string fitStatusM3 = "Y";
		
	if(lr)
	{
		//fit model with no SNPs
		model2SNPs->fitModel(rss0, 0, 0, true);

		//fit M2A model
		if(!model2SNPs->fitModel(rss0_1, 1, 0, true)) fitStatusM2A = "D";
	};

	//fit M2P model
	if(!model2SNPs->fitModel(rss0_2, 0, 1)) fitStatusM2P = "D";

	//fit M3 model
	if(!model2SNPs->fitModel(rss0_1_2, 1, 1)) fitStatusM3 = "D";

	//Output results
	resultsFile << snpID << " " << snpWindow->getAnchorChromosome() << " " << snpWindow->getAnchorSNPName() << " " << snpWindow->getAnchorSNPBP() << " " <<snpWindow->getPartnerSNPName() << " " << snpWindow->getPartnerSNPBP();
	resultsFile << " " << snpWindow->getLastAnchorPartnerCorr() << " " << model2SNPs->getParameter(4);

	unsigned int nMinuspf;	
	double fstat;

	//Output standard SnipSnip test result
	//FSTAT = ((RSSr - RSSf)/(pf - pr))/RSSf/(n-pf)
	//RSSr residual sum of sqs reduced model, RSSf full model
	//pr, no of parameters reduced model, pf full model, FSTAT ~ F_{pf-pr, n-pf}
	if(fitStatusM2P == "Y" && fitStatusM3 =="Y")
	{
		if(rss0_1_2 == 0) if(rss0_2 > 0) resultsFile << " Y Inf 0"; else resultsFile << " Y NA 1";
		else
		{
			unsigned int nonMissingSubs = model2SNPs->getNoNonMissingSubjects();
			if(nonMissingSubs > (3 + noCovariates))
			{
				nMinuspf = nonMissingSubs - 3 - noCovariates;
				fstat = (rss0_2 - rss0_1_2)/(rss0_1_2/((double)(nMinuspf)));
				resultsFile << " Y " << fstat << " " << getPvalueFStat(fstat, 1, nMinuspf);	
			}
			else
			{
				 resultsFile << " D NA NA";
			};
		};
	}
	else if(fitStatusM2P == "D" || fitStatusM3 == "D") resultsFile << " D NA NA";
	else resultsFile << " N NA NA";

	//Output standard linear regression test result
	if(lr)
	{
		if(fitStatusM2A == "Y")
		{
			if(rss0_1 == 0) if(rss0 > 0) resultsFile << " Y Inf 0"; else resultsFile << " Y NA 1";
			else
			{
				unsigned int nonMissingSubsLR = model2SNPs->getNoNonMissingSubjectsLR();
				if(nonMissingSubsLR > (2 + noCovariates))
				{
					nMinuspf = nonMissingSubsLR - 2 - noCovariates;
					fstat = (rss0 - rss0_1)/(rss0_1/((double)(nMinuspf)));
					resultsFile << " Y " << fstat << " " << getPvalueFStat(fstat, 1, nMinuspf);
				}
				else
				{
					resultsFile << " D NA NA";
				};
			};
		}
		else if(fitStatusM2A == "D") resultsFile << " D NA NA";
		else resultsFile << " N NA NA";
	};

	resultsFile << "\n";
};
