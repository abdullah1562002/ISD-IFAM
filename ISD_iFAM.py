import pandas as pd
import pyomo.environ as pe
import pyomo.opt as po
from scipy.optimize import linprog
import numpy as np
import json
import openpyxl
from FAM import FAM
from iFAM import iFAM


class ISD_iFAM(iFAM):
    
    def __init__(self, excel_file):
        self.excel_file = excel_file
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

                                            # ------------------------- Read Change in Demands ------------------------- # 
    def read_demand_change(self):
        self.demand_change = pd.read_excel(self.excel_file, sheet_name="Change_Demand", index_col=None, engine="openpyxl")
        

                                        # ------------------------- Create Optional Itenary List ------------------------- # 
    def create_optional_itenary_list(self):
        self.optional_itenaries = pd.Series()
        for optional_flight in self.Optional_Flights["flight number"]:
            for Itenary in self.Itenaries.index:
                if str(optional_flight) in str(self.Itenaries.loc[Itenary,"Flights"]).split(","):
                    self.optional_itenaries = self.optional_itenaries._append(pd.Series(Itenary+1),ignore_index = True)
                    self.optional_itenaries = self.optional_itenaries.sort_values()
        
        # New_Variables (Zq)
        self.Zq_variables = []
        for value in self.optional_itenaries:
            self.Zq_variables.append(f"Z_{value}")
        
        
                                       # -------------------------     Optimization (Pyomo)   ------------------------- # 
    def optimize_ISD_iFAM_pyomo(self):

        model = pe.ConcreteModel()
 
        model.x = pe.Var(self.Dummy_variables, within = pe.Binary)
        model.RON = pe.Var(self.RON_variables, within = pe.NonNegativeIntegers)
        model.y = pe.Var(self.ground_link_variables, within = pe.NonNegativeIntegers)
        model.t = pe.Var(self.spill_recapture_variables, within = pe.NonNegativeIntegers)
        model.z = pe.Var(self.Zq_variables, within = pe.Binary)


        # Coverage Constraints 
        model.coverage_constraints = pe.ConstraintList()
        for flight in self.Flight_Schedule["flight number"]:
            coverage_sum = sum([model.x[f"X{flight}_{fleet+1}"] for fleet in range(self.Number_fleets)]) 
            if flight in self.Optional_Flights["flight number"].to_list():
                model.coverage_constraints.add(coverage_sum <= 1)
            else:
                model.coverage_constraints.add(coverage_sum == 1)
        
        # Resource Constraints
        model.resource_constraints = pe.ConstraintList()
        for fleet in range(self.Number_fleets):
            resource_sum = sum([model.RON[f"RON_{station}_{fleet+1}"] for station in self.Stations])
            model.resource_constraints.add(resource_sum <= self.fleets["size"][fleet])
        
        # Balance Constraints
        model.balance_constraints = pe.ConstraintList()
        
        new_x = {}
        new_ron = {}
        new_y = {}
        for node in range(self.Number_of_Nodes):
            for fleet in range(self.Number_fleets):
                #Inbounds
                for In in self.Nodes[node]["Inbound"]: 
                    if type(In) is str:
                        if In[0] == "R":
                            new_ron[f"{In}_{fleet+1}"] = 1
                        else: 
                            new_y[f"{In}_{fleet+1}"] = 1
                    else:
                        new_x[f"X{In}_{fleet+1}"] = 1
                #Outbounds
                for Out in self.Nodes[node]["Outbound"]: 
                    if type(Out) is str:
                        if Out[0] == "R":
                            new_ron[f"{Out}_{fleet+1}"] = -1
                        else: 
                            new_y[f"{Out}_{fleet+1}"] = -1
                    else:
                        new_x[f"X{Out}_{fleet+1}"] = -1
                
                
                # Sum all Variables in balance Constraints
                flight_sum = sum(model.x[f"{key}"]*value for key, value in new_x.items())
                Ron_sum = sum(model.RON[f"{key}"]*value for key, value in new_ron.items())
                y_sum = sum(model.y[f"{key}"]*value for key, value in new_y.items())
                balance_sum = sum([flight_sum, Ron_sum, y_sum])

                model.balance_constraints.add(expr= balance_sum == 0)

                new_x.clear()
                new_ron.clear()
                new_y.clear()
        
        # Spill-Recapture Demand Constraints
        model.spill_recapture_Demand_constraints = pe.ConstraintList() 
        counter = 0
        for Itenary in self.I_I_matrix_summation.index:
            spill_recapture_Demand_constraints = sum(model.t[f"{var}"] for var in self.I_I_matrix_summation.loc[Itenary].split(" + "))
            demand_correction = 0
            for q in self.optional_itenaries:
                demand_correction += sum( [-self.demand_change.loc[q-1, counter+1]*(1-model.z[f"Z_{q}"])] )
                
            model.spill_recapture_Demand_constraints.add(expr= sum([spill_recapture_Demand_constraints,demand_correction]) <= self.Itenaries.loc[counter,"Passenger Demand"] )
            counter += 1
  
        # Flight Interaction Constrainst
        
                                                                    # Preform Calculations  
                                                    # (spill - recapture + Seat_f - Zq_terms >= Qf)  multiply by negative
                                                    # (-spill + recapture - Seat_f + Zq_terms <= -Qf)  

        model.flight_interaction_constraints = pe.ConstraintList()
        counter = 0
        
        self.F_I_I_matrix_recapture_summation = self.F_I_I_matrix_recapture_summation.replace(0,'')
        for flight in self.Flight_Schedule.index:
            spill_demand = sum(model.t[f"{var}"]*-1 for var in self.F_I_I_matrix_summation.iloc[flight].split(" + "))
            if self.F_I_I_matrix_recapture_summation.iloc[flight] != '':
                recapture_demand = sum(model.t[f"{var}"]*self.F_I_I_matrix_recapture_summation_probability.iloc[flight] for var in [self.F_I_I_matrix_recapture_summation.iloc[flight]])
            else:
                recapture_demand = 0
            
            seat_f = sum(model.x[f"X{flight+1}_{fleet+1}"]*-self.fleets_capacity[fleet] for fleet in self.fleets.index)
            
            Zq_terms = 0
            for Itenary in self.Itenaries.index:
                for Optional_itenary in self.optional_itenaries:
                    if self.F_I_I_Matrix_symbolic.loc[f"F{flight+1}", f"I{Itenary+1}"] != 0:
                        Zq_terms += sum( [self.demand_change.loc[Optional_itenary-1, Itenary+1]*(1-model.z[f"Z_{Optional_itenary}"])] )

            
            flight_interaction_sum = sum([recapture_demand, spill_demand, seat_f, Zq_terms])
            model.flight_interaction_constraints.add(expr= flight_interaction_sum   <= -self.Unconstrained_Demand_Qf[flight])

        
        # Zq_First Constraint 
        model.Zq_First_constraint  = pe.ConstraintList()
        for optional_itenary in self.optional_itenaries:
            flight_in_itenary = str(self.Itenaries.loc[optional_itenary -1, "Flights"]).split(",")
            for flight in flight_in_itenary:
                X_terms = 0
                for fleet in self.fleets.index:
                    X_terms += -model.x[f"X{flight}_{fleet+1}"]
            
                model.Zq_First_constraint.add(expr=sum( [model.z[f"Z_{optional_itenary}"], X_terms] )    <= 0)

        # Zq_Second Constraint
        model.Zq_Second_constraint  = pe.ConstraintList()
        for optional_itenary in self.optional_itenaries:
            flight_in_itenary = str(self.Itenaries.loc[optional_itenary -1, "Flights"]).split(",")
            X_terms = 0
            for flight in flight_in_itenary:
                for fleet in self.fleets.index:
                    X_terms += model.x[f"X{flight}_{fleet+1}"]
        
            model.Zq_Second_constraint.add(expr=sum( [-model.z[f"Z_{optional_itenary}"], X_terms] )   <= -1+len(flight_in_itenary))

        # Objective_Function
        
        # Minimize Cost ----> C+(S-M)+delta_R 

        #---> Operational Costs (C)
        new_dict = {}
        for flight in self.Flight_Schedule.index:
                for fleet in self.fleets.index:
                    new_dict[f"X{flight+1}_{fleet+1}"] = self.Flight_Schedule.loc[flight, f"e{fleet+1}"] 

        operational_costs = sum(model.x[f"{key}"]*value for key, value in new_dict.items())

        #---> Spilled Costs (S)
        new_dict = {}
        for Itenary in self.Itenaries.index:
            for var in self.I_I_matrix_summation.iloc[Itenary].split(" + "):
                new_dict[var] = self.Itenaries.loc[Itenary, "Fare"]
        
        spilled_costs =sum(model.t[f"{key}"]*value for key, value in new_dict.items())

        #---> Recaptured Costs/Revenue (M)
        new_dict = {}
        for Itenary in self.Itenaries.index:
            for var in self.I_I_matrix_recapture_summation.iloc[Itenary].split(" + "):
                if var != '' and var != 0:
                    new_dict[var] = self.Itenaries.loc[Itenary, "Fare"]*self.I_I_matrix_recapture_summation_probability.iloc[Itenary]
                
        recaptured_costs =sum(model.t[f"{key}"]*-value for key, value in new_dict.items())

        #---> Unconstrained Revenue loss due to flight deletion (delta_R)
        delta_r = 0
        for itenary in self.optional_itenaries:
            term_1 = self.Itenaries.loc[itenary-1, "Passenger Demand"]*self.Itenaries.loc[itenary-1, "Fare"]
            term_2 = sum([self.demand_change.loc[itenary-1, p+1]*self.Itenaries.loc[p, "Fare"] for p in self.Itenaries.index])
            delta_r += sum([(1-model.z[f"Z_{itenary}"])*(term_1-term_2)])
        

        #---> Unconstrained revenue (R) R = fare_p*D_p
        counter = 0
        self.Flight_Schedule["Fare"] = 0
        for flight_Demand in self.Unconstrained_Demand_Qf:
            self.Flight_Schedule.loc[counter, "Fare"] = self.Itenaries.loc[counter, "Fare"]
            counter += 1

        self.Unconstrained_Revenue = 0
        for Itenary in self.Itenaries.index:
            self.Unconstrained_Revenue += self.Itenaries.loc[Itenary, "Fare"] * self.Itenaries.loc[Itenary, "Passenger Demand"]

       
        objective_fun = sum([operational_costs, spilled_costs, recaptured_costs, delta_r])

        model.objective = pe.Objective(sense = pe.minimize, expr = objective_fun)

        with open(f"Outputs/{self.save_folder}/model_{self.save_folder}_output.txt", "w") as output_file:
             model.pprint(output_file)


        solver = po.SolverFactory("gurobi")     
        results = solver.solve(model, tee=False)

        print("##################################################")
        print("##                    ISD-iFAM                  ##")
        print("##################################################")
        print(f"Maxium Profit (iFAM) : {self.Unconstrained_Revenue - pe.value(model.objective)}")
        print(f"Minimum Cost  (iFAM) : {pe.value(model.objective)}")

        for var in model.x:
            print(f"{var} : {pe.value(model.x[f"{var}"])}")

        for var in model.RON:
            print(f"{var} : {pe.value(model.RON[f"{var}"])}")

        for var in model.y:
            print(f"{var} : {pe.value(model.y[f"{var}"])}")

        # for var in model.t:
        #     try: 
        #         print(f"{var} : {pe.value(model.t[f"{var}"])}")
        #     except:
        #         print(f"{var} : 0")

        for var in model.z:
            print(f"{var} : {pe.value(model.z[f"{var}"])}")

                                                 # ------------------------- Retuen Outputs ------------------------- #      
    def run_analysis_ISD_iFAM(self):

        self.save_folder = "ISD_iFAM"
        self.read_flight_schedule()
        self.read_fleets()
        self.read_Itenaries()
        self.read_probability()
        self.identify_stations()
        self.create_nodes()
        self.constraints_matrices()
        self.Calculate_Unconstrained_Demand_Qf()
        self.spilled_recaptured_Demand_Calculation()
        self.constraints_matrices()
        self.Calculate_Unconstrained_Demand_Qf()
        self.read_demand_change()
        self.create_optional_itenary_list()
        self.optimize_ISD_iFAM_pyomo()