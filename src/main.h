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


/*! \file main.h
    \brief This file defines a global variable to suppress output to screen. 
    
*/

#ifndef __MAIN
#define __MAIN

extern bool outputToScreen; 
extern ofstream logFile; 

template<typename T>
//! Outputs message to screen and log file
void out(const T & text)
{	
	if(outputToScreen) cout << text;
	logFile << text;
};

template<typename T>
//! Outputs error message to screen and log file
void outErr(const T & text)
{
	cerr << text;
	logFile << text;
};

#endif
