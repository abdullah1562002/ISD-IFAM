import FAM
import pandas as pd


pd.set_option('future.no_silent_downcasting', True)

# Create an instance of the class
analyzer_FAM = FAM.FAM("ISD_Project_Data.xlsx")

# Perform operations
analyzer_FAM.run_analysis_FAM()

# Create an instance of the class
analyzer_iFAM = FAM.iFAM("ISD_Project_Data.xlsx")

# Perform operations
analyzer_iFAM.run_analysis_iFAM()