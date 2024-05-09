import FAM
import pandas as pd
import warnings

# Ignore all future warnings (Related to Syntax Updates)
warnings.simplefilter(action='ignore', category=FutureWarning)


# Create an instance of the class
analyzer_FAM = FAM.FAM("ISD_Project_Data.xlsx")

# Perform operations
analyzer_FAM.run_analysis_FAM()

# Create an instance of the class
analyzer_iFAM = FAM.iFAM("ISD_Project_Data.xlsx")

# Perform operations
analyzer_iFAM.run_analysis_iFAM()