import pandas as pd
import pyomo.environ as pe
import pyomo.opt as po
from scipy.optimize import linprog
import numpy as np
import json
import openpyxl
from FAM import FAM
from iFAM import iFAM


class ISD_iFAM(FAM, iFAM):
    
    def __init__(self, excel_file):
        self.excel_file = excel_file
        self.Flight_Schedule = pd.DataFrame()
        self.fleets = pd.DataFrame()
        self.Itenaries = pd.DataFrame()
        self.Stations = []
        self.Nodes = []
        self.ground_links = []
        self.Variables = []
        self.Dummy_variables = []
        self.RON_variables = []
        self.ground_link_variables = [] 
        self.Number_fleets = 0
        self.Number_of_ground_links = 0
        self.Number_Stations = 0
        self.Number_Flights = 0
        self.Number_of_Nodes = 0
        self.Number_of_Variables = 0
        self.Number_of_Coverage_Constraints = 0
        self.Number_of_Resource_Constraints = 0
        self.Number_of_Balance_Constraints = 0
        self.coverage_rhs = []
        self.flight_itenary_incidence_Matrix = pd.DataFrame()
        self.spill_recapture_variables = []
        self.inequality_constraints = pd.DataFrame()
        self.equality_constraints = pd.DataFrame()
        self.inequality_rhs = np.array([])
        self.equality_rhs = np.array([])