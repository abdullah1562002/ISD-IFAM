from FAM import FAM
from iFAM import iFAM
from ISD_iFAM import ISD_iFAM
import pandas as pd
import warnings

# Ignore all future warnings (Related to Syntax Updates)
warnings.simplefilter(action='ignore', category=FutureWarning)
 


# Create an instance of the class
analyzer_FAM = FAM("ISD_Project_Data.xlsx")

# Perform operations
analyzer_FAM.run_analysis_FAM()

# Create an instance of the class
analyzer_iFAM = iFAM("ISD_Project_Data.xlsx")

# Perform operations
analyzer_iFAM.run_analysis_iFAM()

# Create an instance of the class
analyzer_ISD_iFAM = ISD_iFAM("ISD_Project_Data.xlsx")

# Perform operations
analyzer_ISD_iFAM.run_analysis_ISD_iFAM()