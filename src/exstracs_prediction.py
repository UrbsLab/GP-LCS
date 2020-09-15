"""
Name:        ExSTraCS_Prediction.py
Authors:     Ryan Urbanowicz - Written at Dartmouth College, Hanover, NH, USA
Contact:     ryan.j.urbanowicz@darmouth.edu
Created:     April 25, 2014
Modified:    August 25,2014
Description: Based on a given match set, this module uses a voting scheme to select the phenotype prediction for ExSTraCS.

---------------------------------------------------------------------------------------------------------------------------------------------------------
ExSTraCS V2.0: Extended Supervised Tracking and Classifying System - An advanced LCS designed specifically for complex, noisy classification/data mining tasks,
such as biomedical/bioinformatics/epidemiological problem domains.  This algorithm should be well suited to any supervised learning problem involving
classification, prediction, data mining, and knowledge discovery.  This algorithm would NOT be suited to function approximation, behavioral modeling,
or other multi-step problems.  This LCS algorithm is most closely based on the "UCS" algorithm, an LCS introduced by Ester Bernado-Mansilla and
Josep Garrell-Guiu (2003) which in turn is based heavily on "XCS", an LCS introduced by Stewart Wilson (1995).

Copyright (C) 2014 Ryan Urbanowicz
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABLILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation,
Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
---------------------------------------------------------------------------------------------------------------------------------------------------------
"""

# Import Required Modules-------------------------------
from exstracs_constants import *
import random


# ------------------------------------------------------

class Prediction:
    def __init__(self, population,
                 exploreIter):  # now takes in population ( have to reference the match set to do prediction)  pop.matchSet
        """ Constructs the voting array and determines the prediction decision. """
        self.decision = None
        # self.classCount
        # -------------------------------------------------------
        # DISCRETE PHENOTYPE
        # -------------------------------------------------------
        if cons.env.formatData.discretePhenotype:

            # Temporary voting method:
            num_ones = 0
            num_zeros = 0
            for ref in population.matchSet:
                cl = population.popSet[ref]
                if cl.phenotype == '1':
                    num_ones += 1
                else:
                    num_zeros += 1
            self.decision = '1' if num_ones > num_zeros else '0'

            self.vote = {}

        else:  # ContinuousCode #########################
            if len(population.matchSet) < 1:
                self.decision = None
            else:
                treeCount = 0
                value = 0
                for ref in population.matchSet:
                    cl = population.popSet[ref]
                    value += cl.phenotype * cl.numerosity
                    treeCount += cl.numerosity

                if (cons.popInitGP == 1.0):
                    self.decision = value / treeCount
                    return



    def getDecision(self):
        """ Returns prediction decision. """
        return self.decision

    def getSet(self):
        """ Returns prediction decision. """
        return self.vote
