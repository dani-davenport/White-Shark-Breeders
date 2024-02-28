# White Shark Breeders

This repository contains the code and data for a re-analysisof the publication entitled “Effective number of white shark (Carcharodon carcharias, Linnaeus) breeders is stable over four successive years in the population adjacent to eastern Australia and New Zealand” which was published in Ecology and Evolution (Volume 11, Issue 1, pp 186-198) on 7 October 2020.

# Scripts 
Reanalysis is in file called Final_Report.rmd  
All code scripts reffered to in this report can be found in Scripts/  
  
A brief description of code in "Scripts/":  
To extract r* values from NeEstimator outputs, use script Scripts/Get_RStar.R. To get the file needed, you need to run Ne estimator with the check-box marked for Output Burrows coefficient, then pull the line that says "Weighted Mean of r^2-drift" using this script.   
To make combined estimates, use spreadsheet Scripts/Combining Estimates.xlsx, replace with your own values where indicated  
To estimate age at length in samples, use Scripts/Age_at_Length.R
