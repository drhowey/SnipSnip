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


/*! \file main.cpp
    \brief This file reads in the initial input files and options.
    
    This file also outputs usage instructions and program details.
*/

#include <iostream>
#include <ostream>
#include <string>
#include <time.h>

using namespace std; // initiates the "std" or "standard" namespace
 
#include "main.h"
#include "Analysis.h"

bool outputToScreen = true; 
ofstream logFile;

//! Output program title to screen
void header()
{
	out("\nSnipSnip: Imputation without imputation, v1.1\n");
	out("---------------------------------------------\n");
	out("Copyright 2013 Richard Howey, GNU General Public License, v3\n");
	out("Institute of Genetic Medicine, Newcastle University\n\n");
};

//! Output program usage to screen
void usage()
{
		header();
	 	
		out("Usage:\n  ./snipsnip [options] pedigree.bed \n");
		out(" or ./snipsnip -pf parameterfile [pedigree.bed]\n\n");

		out("Options:\n");
		
		out("  -window-size n    -- fix window at n SNPS, n must be even\n");
		out("  -window-size-bp x -- size of window, x, in kB\n");
		out("  -start a          -- start analysis from SNP number a\n");
		out("  -start-end a b    -- start and end analysis from SNP numbers a to b\n");
		out("  -i file.bed       -- input binary pedigree file, file.bed\n");
		out("  -o results.dat    -- output results file, results.dat\n");
		out("  -log results.log  -- log filename, results.log\n");
		out("  -covar covars.dat -- covariate filename, covars.dat\n");
		out("  -covar-number no  -- covariate numbers, no\n");
		out("  -covar-name na    -- covariate names, na\n");
		out("  -linear           -- use linear regression\n");
		out("  -mqtv x           -- missing quantitive trait value for linear regression\n");	
		out("  -dominant         -- use dominant correlation partner metric\n");
		out("  -recessive        -- use recessive correlation partner metric\n");		
		out("  -lr               -- perform standard logistic(linear) regression tests\n");
		out("  -excsnp bp        -- exclude SNP(base pair bp) as partner\n");
		out("  -so               -- suppress output to screen\n\n");

		out("Default Options in Effect:\n");
		out("  -window-size 10\n");		
		out("  -o snipsnipResults.dat\n\n");
};


//! Get an option value either from the parameter file or the command line
void getOptionValue(unsigned int & anUnInt, bool & useParaFile, int & argcount, int & argc, char * argv[], ifstream & readParaFile)
{	
	if(useParaFile)
	{
		if(readParaFile.eof()) return;
		readParaFile >> anUnInt;		
	}
	else
	{
		argcount++; if(argcount >= argc) return;				
		anUnInt = atoi(argv[argcount]);		
	};
};

//! Get an option value either from the parameter file or the command line
void getOptionValue(double & aDouble, bool & useParaFile, int & argcount, int & argc, char * argv[], ifstream & readParaFile)
{
	if(useParaFile)
	{
		if(readParaFile.eof()) return;		
		readParaFile >> aDouble;		
	}
	else
	{
		argcount++; if(argcount >= argc) return;			
		aDouble = atof(argv[argcount]);		
	};
};

//! Get an option value either from the parameter file or the command line
void getOptionValue(string & aString, bool & useParaFile, int & argcount, int & argc, char * argv[], ifstream & readParaFile)
{
	if(useParaFile)
	{
		if(readParaFile.eof()) return;		
		readParaFile >> aString;
	}
	else
	{
		argcount++; if(argcount >= argc) return;		
		aString = argv[argcount];
	};
};

//! Gets the log filename from the results file name
string getDefaultLogFileName(string & outFileName)
{
	unsigned int filenameLength = outFileName.length();
	string logFileName;

	//find extension
	unsigned int a = filenameLength - 1;
	while(a > 0)
	{
		if(outFileName.substr(a, 1) == ".") break;
		a--;
	};

	if(a > 0) logFileName = outFileName.substr(0, a) + ".log";
	else  logFileName = outFileName + ".log";

	return logFileName;
};

//! The start of the program
int main(int argc, char * argv[])
{
	time_t start,end;
	double dif;
	time(&start);

	int argcount = 1;
	string option;
	string filename = "";
	string paraFilename = "";
	string outputFileName = "snipsnipResults.dat";
	string logFilename = "";
	string covariateFilename = "";
	string covariates = "";
	unsigned int windowBPSize = 0;
	unsigned int windowSNPSize = 10;
	unsigned int startSNP = 1;
	unsigned int endSNP = 0;
	unsigned int metricType = 1; //potential partner metric type
		
	bool lr = false;
	bool linear = false;	
	outputToScreen = true;	
	unsigned int excludeSNP = 0;
	double missingQTValue = -9;

	bool useParaFile = false;
	ifstream readParaFile;

	if(argcount < argc) option = argv[argcount];

	//deal with parameter file
	if(option == "-pf")
	{
		argcount++; 
		if(argcount < argc) paraFilename = argv[argcount];

		//open parameter file		
		readParaFile.open(paraFilename.c_str());
		if(!readParaFile.is_open())
		{
			header();
			outErr("Cannot read parameter file: "); outErr(paraFilename); outErr("!\n");
			exit(0);
		};

		argcount++; 
		useParaFile = true;
	};

	//set given options
	while((!useParaFile && argcount < argc && argv[argcount][0] == '-') || (useParaFile && !readParaFile.eof()))
	{

		if(useParaFile)
		{
			//find the start of the next command
			do{
				readParaFile >> option;
				if(option.length() >= 2 && option.substr(0,1) == "-") break;				
			}while(!readParaFile.eof());
		}
		else
		{
			option = argv[argcount];
		};

		if(useParaFile && readParaFile.eof()) break;

		if(option ==  "-window-size-bp" || option ==  "-w") {getOptionValue(windowBPSize, useParaFile, argcount, argc, argv, readParaFile); windowSNPSize = 0;}
		else if(option ==  "-window-size" || option ==  "-ws") getOptionValue(windowSNPSize, useParaFile, argcount, argc, argv, readParaFile);		
		else if(option ==  "-start" || option ==  "-s") getOptionValue(startSNP, useParaFile, argcount, argc, argv, readParaFile);
		else if(option ==  "-i") getOptionValue(filename, useParaFile, argcount, argc, argv, readParaFile);
		else if(option ==  "-covar") getOptionValue(covariateFilename, useParaFile, argcount, argc, argv, readParaFile);
		else if(option ==  "-covar-name" || option ==  "-covar-number") getOptionValue(covariates, useParaFile, argcount, argc, argv, readParaFile);
		else if(option ==  "-log") getOptionValue(logFilename, useParaFile, argcount, argc, argv, readParaFile);	
		else if(option ==  "-linear") linear = true;
		else if(option ==  "-lr") lr = true;		
		else if(option ==  "-mqtv") getOptionValue(missingQTValue, useParaFile, argcount, argc, argv, readParaFile);		
		else if(option ==  "-dominant") metricType = 2;
		else if(option ==  "-recessive") metricType = 3;		
		else if(option ==  "-start-end" || option ==  "-se")
		{	
			getOptionValue(startSNP, useParaFile, argcount, argc, argv, readParaFile);
			getOptionValue(endSNP, useParaFile, argcount, argc, argv, readParaFile);			
		}
		else if(option ==  "-o") getOptionValue(outputFileName, useParaFile, argcount, argc, argv, readParaFile);
		else if(option ==  "-excsnp") getOptionValue(excludeSNP, useParaFile, argcount, argc, argv, readParaFile);
		else if(option == "-so") outputToScreen = false;		
		else
		{
			if(logFilename != "") logFile.open(logFilename.c_str());
			else logFile.open(getDefaultLogFileName(outputFileName).c_str());			

			header();
			if(useParaFile) {outErr("Unrecognised option: "); outErr(option); outErr("\n\n");}
			else {outErr("Unrecognised command line switch: "); outErr(option); outErr("\n\n");};			
    		exit(0);
		};

		if(!useParaFile) argcount++;
	};

	if(argcount < argc) filename = argv[argcount];	
	
	if(logFilename != "") logFile.open(logFilename.c_str());
	else
	{
		logFilename = getDefaultLogFileName(outputFileName);
		logFile.open(logFilename.c_str());	
	};

	if(filename == "")
	{	
		usage();
		if(argc > 1) outErr("\nInput file not set!\n\n");
		exit(0);
	};

	if(filename.length() >= 4 && filename.substr(filename.length()-4) != ".bed")
	{
		header();
		outErr("A binary pedigree file (.bed) is required for SnipSnip!\n\n");
		outErr("Try creating a .bed file from a .ped file using PLINK:\n");
		outErr("plink --file mydata --make-bed\n\n");
		exit(0);
	}
	else if((int)(windowSNPSize/2) != ((windowSNPSize+1)/2))
	{
		header();
		outErr("SNP window size must be even!\n\n");
		outErr("Try "); outErr((windowSNPSize+1)); outErr(" instead of "); outErr((windowSNPSize)); outErr("!\n");
		exit(0);

	}
	else if(startSNP > endSNP && endSNP != 0)
	{
		outErr("The start SNP must not be after the end SNP!\n");
		exit(0);
	};

	//output options to screen	
	header();
	out("Parameters:\n");
	out("Input file: "); out(filename); out("\n");
	out("Output file: "); out(outputFileName); out("\n");
	out("Log file: "); out(logFilename); out("\n");
	if(covariateFilename != "")
	{
			out("Covariate file: "); out(covariateFilename);
			if(covariates != "") {out(" (using covariates: "); out(covariates); out(")");};
			out("\n");
	};
	if(lr)
	{
		if(!linear) out("Also including single-SNP logistic regression results\n");
		else out("Also including single-SNP linear regression results\n");
	};

	if(linear)
	{ 
		out("Using linear regression (with missing QT value "); out(missingQTValue); out(")\n");		
	};

	if(windowSNPSize != 0) {out("Fixed SNP window size: "); out(windowSNPSize); out(" SNPs\n");}
	else { out("SNP window size: "); out(windowBPSize); out(" base pairs (bps)\n");};

	out("Partner metric: ");	
	if(metricType == 1) out("Multiplicative");
	else if(metricType == 2) out("Dominant");
	else if(metricType == 3) out("Recessive");
	
	out("\n");	

	if(startSNP != 1) {out("Start at SNP number: "); out(startSNP); out("\n");};
	if(endSNP != 0) {out("End at SNP number: "); out(endSNP); out("\n");};
	if(excludeSNP > 0) {out("Exclude SNP from being chosen as partner, with BP: "); out(excludeSNP); out("\n");};
	out("\n");
	
	//create analysis option and run analysis
	Analysis anAnalysis(filename, outputFileName, windowBPSize, windowSNPSize, metricType, startSNP, endSNP, excludeSNP, linear, missingQTValue, lr, covariateFilename, covariates);

	anAnalysis.runAnalysis();	

	time(&end);
	dif = difftime(end, start);
	out("Run time: "); out(getTime(dif)); out("\n\n");

	logFile.close();
};

