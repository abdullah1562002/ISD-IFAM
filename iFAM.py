import pandas as pd
import pyomo.environ as pe
import pyomo.opt as po
from scipy.optimize import linprog
import numpy as np
import json
import openpyxl
from FAM import FAM



class iFAM(FAM):
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
        
        
    def read_probability(self):
        self.recapture_probability = pd.read_excel(self.excel_file, sheet_name="Recapture_Probability", index_col=None, engine="openpyxl")
        self.recapture_probability["Variable_Equivalent"] = self.recapture_probability["Variable"].replace("b","t", regex=True)
        
        
                             # ------------------ Calculate Demand (Change Demand from Itenary to Flight Based) ------------------ #  
    def Calculate_Unconstrained_Demand_Qf(self):
        self.Itenaries_var = [f"I{I_No+1}" for I_No in range(len(self.Itenaries))]
        self.Flights_var = [f"F{F_No+1}" for F_No in range(len(self.Flight_Schedule.index))]
        self.flight_itenary_incidence_Matrix = pd.DataFrame(columns = self.Itenaries_var)
        
        # Create Flight Itenary Incidence Matrix
        new_row = {}
        for flight in self.Flight_Schedule["flight number"]:
            for Itenary in self.Itenaries.index:
                if str(flight) in list(str(self.Itenaries["Flights"][Itenary]).split(",")):
                    new_row[f"I{Itenary + 1}"] = 1
            self.flight_itenary_incidence_Matrix = self.flight_itenary_incidence_Matrix._append(new_row, ignore_index=True).fillna(0)
            new_row.clear()

        self.flight_itenary_incidence_Matrix.index = self.Flights_var
        self.flight_itenary_incidence_Matrix.to_csv(f"Outputs/{self.save_folder}/Flight_Itenary_Incidence_Matrix_{self.save_folder}.csv")

        self.passenger_Demand = np.array(self.Itenaries["Passenger Demand"])
        self.Unconstrained_Demand_Qf = np.matmul(self.flight_itenary_incidence_Matrix.to_numpy() , self.passenger_Demand)
        self.Flight_Schedule["Demand"] = pd.Series(self.Unconstrained_Demand_Qf)
        self.fleets_capacity = np.array(self.fleets["capacity"])
        
        
        # SEAT_f Matrix
        self.SEAT_matrix = pd.DataFrame(columns=self.Variables)
        new_row = {}
                
        for flight in self.Flight_Schedule["flight number"]:
            for fleet in range(self.Number_fleets):
                new_row[f"X{flight}_{fleet+1}"] = self.fleets_capacity[fleet]
            
            
            self.SEAT_matrix = self.SEAT_matrix._append(new_row, ignore_index=True)
            new_row.clear()
        
        self.SEAT_matrix = self.SEAT_matrix.fillna(0)
        self.SEAT_matrix.to_csv(f"Outputs/{self.save_folder}/SEAT_f_matrix{self.save_folder}.csv")

                                             # ------------------ Calculate Spilled and Recaptured Demands ------------------ #  
    def spilled_recaptured_Demand_Calculation(self):
        self.itenary_itenary_spill_recapture_Matrix = pd.DataFrame(columns = self.Itenaries_var)
        self.itenary_itenary_spill_recapture_Matrix[f"I{len(self.Itenaries_var)+1}"] = ""
        
        new_row = {}

        # Create Iteanry - Itenary Incidence Matrix
        for itenary_1 in self.Itenaries.index:
            ref_Itenary = list(self.Itenaries["Itenary"][itenary_1].split("-"))
            for itenary_2 in range(len(self.Itenaries)+1):
                    if itenary_2 < len(self.Itenaries):
                        check_Itenary = list(self.Itenaries["Itenary"][itenary_2].split("-"))
                        if  itenary_1 != itenary_2 and ref_Itenary[0] == check_Itenary[0] and ref_Itenary[-1] == check_Itenary[-1] and  len(ref_Itenary) <= len(check_Itenary):
                            new_row[f"I{itenary_2 + 1}"] = f"t({itenary_1 + 1},{itenary_2 + 1})"
                        else:
                            new_row[f"I{len(self.Itenaries_var)+1}"] = f"t({itenary_1 + 1},{len(self.Itenaries_var)+1})"

                        # Create Spill Recapture Variables t(p,r)
                        self.spill_recapture_variables.append(f"t({itenary_1 + 1},{itenary_2 + 1})")
                    else: 
                        self.spill_recapture_variables.append(f"t({itenary_1 + 1},{itenary_2 + 1})")
            
            self.itenary_itenary_spill_recapture_Matrix = self.itenary_itenary_spill_recapture_Matrix._append(new_row, ignore_index=True).fillna(0)
            new_row.clear()

        self.itenary_itenary_spill_recapture_Matrix.index = self.Itenaries_var
        self.itenary_itenary_spill_recapture_Matrix.to_csv(f"Outputs/{self.save_folder}/Itenary_Itenary_Incidence_Matrix_{self.save_folder}.csv")
        
        # spill recapture variables bounds
        self.spill_recapture_variables_bounds = [(0,None) for i in self.spill_recapture_variables]
        
        ###############################################  Spilled Demand #########################################

        # Take Summation of each Inteary row (symbolic)    
        self.I_I_matrix_summation = pd.Series()
        for Itenary in self.itenary_itenary_spill_recapture_Matrix.index:
            new_list = self.itenary_itenary_spill_recapture_Matrix.loc[Itenary]
            row_sum  = [variable for variable in new_list if variable != 0]
            row_sum = " + ".join(row_sum)

            self.I_I_matrix_summation = self.I_I_matrix_summation._append(pd.Series(row_sum), ignore_index=True)
        
        
        self.itenary_itenary_spill_recapture_Matrix.index = self.Itenaries_var
        self.Variables.extend(self.spill_recapture_variables)
        self.I_I_matrix_summation.index = self.Itenaries_var
        
        
        # Multiply the Summation of t(p,r) to Flight - Itenary - Incidence Matrix
        self.F_I_I_Matrix_symbolic = self.flight_itenary_incidence_Matrix.transpose()

        for Itenary in self.F_I_I_Matrix_symbolic.index:
            self.F_I_I_Matrix_symbolic.loc[Itenary] = self.F_I_I_Matrix_symbolic.loc[Itenary].replace(1, self.I_I_matrix_summation.loc[Itenary])
        
        self.F_I_I_Matrix_symbolic = self.F_I_I_Matrix_symbolic.transpose()
        
        
        # Take Summation of each Flight row (symbolic)    
        self.F_I_I_matrix_summation = pd.Series()
        for Flight in self.F_I_I_Matrix_symbolic.index:
            new_list = self.F_I_I_Matrix_symbolic.loc[Flight]
            row_sum  = [variable for variable in new_list if variable != 0]
            row_sum = " + ".join(row_sum)

            self.F_I_I_matrix_summation = self.F_I_I_matrix_summation._append(pd.Series(row_sum), ignore_index=True)
                
        self.F_I_I_matrix_summation.index = self.Flights_var
        
        # Spilled Demand Matrix 
        self.spilled_Demand_Matrix = pd.DataFrame(columns=self.Variables)
        
        new_row = {}
        for Flight in self.F_I_I_matrix_summation.index:
            elements = self.F_I_I_matrix_summation.loc[Flight].split(" + ")
            for element in elements:
                if element != 0 :
                    new_row[f"{element}"] = 1
                
            
            self.spilled_Demand_Matrix = self.spilled_Demand_Matrix._append(new_row, ignore_index = True).fillna(0)
            new_row.clear()
        
        # check 
        if self.Number_Flights == len(self.spilled_Demand_Matrix):
            self.spilled_Demand_Matrix.to_csv(f"Outputs/{self.save_folder}/spilled_Demand_Matrix_{self.save_folder}.csv")
        else : 
            print("Dimensions Dosent Match")


        ###############################################  Recapture Demand #########################################

        # Take Summation of each Inteary Column (symbolic) 
        self.I_I_matrix_recapture = self.itenary_itenary_spill_recapture_Matrix.iloc[:, :-1]
        self.I_I_matrix_recapture_summation = pd.Series()
        self.I_I_matrix_recapture_summation_probability = pd.Series()

        for Itenary in self.I_I_matrix_recapture.columns:
            new_list = self.itenary_itenary_spill_recapture_Matrix.loc[:, Itenary]
            row_sum  = [variable for variable in new_list if variable != 0]
            row_probability = []

            if len(row_sum) != 0:
                for element in row_sum:
                    for variable in self.recapture_probability.index:
                        if element == self.recapture_probability.loc[variable, "Variable_Equivalent"]:
                            row_probability.append(self.recapture_probability.loc[variable, "Probability"])

            else:
                row_probability.append(0)
            
            row_sum = " + ".join(row_sum)
            
            self.I_I_matrix_recapture_summation = self.I_I_matrix_recapture_summation._append(pd.Series(row_sum), ignore_index=True)
            self.I_I_matrix_recapture_summation_probability = self.I_I_matrix_recapture_summation_probability._append(pd.Series(row_probability), ignore_index =True)

        self.I_I_matrix_recapture_summation.index = self.Itenaries_var
        self.F_I_I_Matrix_recapture_symbolic = self.flight_itenary_incidence_Matrix
        
        
        # Multiply the Summation of t(p,r) to Flight - Itenary - Incidence Matrix
        for Itenary in self.F_I_I_Matrix_recapture_symbolic.columns:
            self.F_I_I_Matrix_recapture_symbolic.loc[:, Itenary] = self.F_I_I_Matrix_recapture_symbolic.loc[:, Itenary].replace(1, self.I_I_matrix_recapture_summation.loc[Itenary])
            self.F_I_I_Matrix_recapture_symbolic = self.F_I_I_Matrix_recapture_symbolic.replace('',0)
        
        
        
        # Take Summation of each Flight row (symbolic) 
        self.F_I_I_matrix_recapture_summation = pd.Series()
        self.F_I_I_matrix_recapture_summation_probability = pd.Series()
        for Flight in self.F_I_I_Matrix_recapture_symbolic.index:
            new_list = self.F_I_I_Matrix_recapture_symbolic.loc[Flight]
            row_sum  = [variable for variable in new_list if variable != 0]
            row_probability = []

            if len(row_sum) != 0:
                for element in row_sum:
                    for variable in self.recapture_probability.index:
                        if element == self.recapture_probability.loc[variable, "Variable_Equivalent"]:
                            row_probability.append(self.recapture_probability.loc[variable, "Probability"])

            else:
                row_probability.append(0)
            
            row_sum = " + ".join(row_sum)

            self.F_I_I_matrix_recapture_summation = self.F_I_I_matrix_recapture_summation._append(pd.Series(row_sum), ignore_index=True)
            self.F_I_I_matrix_recapture_summation_probability = self.F_I_I_matrix_recapture_summation_probability._append(pd.Series(row_probability), ignore_index=True)

        self.F_I_I_matrix_recapture_summation = self.F_I_I_matrix_recapture_summation.replace('',0)     
        self.F_I_I_matrix_recapture_summation.index = self.Flights_var

        
        # Recapture Demand Matrix 
        self.recapture_Demand_Matrix = pd.DataFrame(columns=self.Variables)
        
        new_row = {}
        for Flight in self.F_I_I_matrix_recapture_summation.index:
            elements = str(self.F_I_I_matrix_recapture_summation.loc[Flight]).split(" + ")
            for element in elements:
                if element != "0":
                    for probability in self.recapture_probability.index:
                        if element == self.recapture_probability["Variable_Equivalent"][probability]:
                            new_row[f"{element}"] = self.recapture_probability["Probability"][probability]
                            
                
                else:
                    new_row = pd.DataFrame(np.zeros((1, len(self.Variables))), columns = self.Variables)
                    
            self.recapture_Demand_Matrix = self.recapture_Demand_Matrix._append(new_row, ignore_index = True).fillna(0)
            new_row = {}
        
        # check 
        if self.Number_Flights == len(self.recapture_Demand_Matrix):
            self.recapture_Demand_Matrix.to_csv(f"Outputs/{self.save_folder}/recaptured_Demand_Matrix_{self.save_folder}.csv")
        else : 
            print("Dimensions Dosent Match")
        
                                            # ------------------ Create Flight Interaction Constraints Matrix ------------------ #  
    def Flight_Interaction_constraints(self):
        
                        # Preform Calculations  
        # (spill - recapture + Seat_f >= Qf)  multiply by negative
        # (-spill + recapture - Seat_f <= -Qf)  multiply by negative

        spill = self.spilled_Demand_Matrix.to_numpy()
        recapture = self.recapture_Demand_Matrix.to_numpy()
        Seat_f = self.SEAT_matrix.to_numpy()

        self.Flight_Interaction_constraints_matrix = pd.DataFrame(-spill + recapture -Seat_f, columns=self.Variables)
        self.Flight_Interaction_constraints_matrix.to_csv(f"Outputs/{self.save_folder}/Flight_Interaction_constraints_matrix_{self.save_folder}.csv")
        self.Flight_Interaction_constraints_rhs = np.array(-self.Unconstrained_Demand_Qf)

                                                # ------------------ Create Demand Constraints Matrix ------------------ # 
    def Demand_constraints(self):
        # Demand Constraint Matrix 
        self.Demand_constraint_Matrix = pd.DataFrame(columns=self.Variables)
        
        new_row = {}
        for Itenary in self.I_I_matrix_summation.index:
            elements = self.I_I_matrix_summation.loc[Itenary].split(" + ")
            for element in elements:
                if element != '0' :
                    new_row[f"{element}"] = 1
                
            
            self.Demand_constraint_Matrix = self.Demand_constraint_Matrix._append(new_row, ignore_index = True).fillna(0)
            new_row.clear()
        
        # check 
        if self.Number_Itenaries == len(self.Demand_constraint_Matrix):
            self.Demand_constraint_Matrix.to_csv(f"Outputs/{self.save_folder}/Demand_constraint_Matrix_{self.save_folder}.csv")
        else : 
            print("Dimensions Dosent Match")

        # Demand Constraint RHS
        self.Demand_constraint_rhs = np.array(self.Itenaries["Passenger Demand"])

                                                 # ------------------ Create Objective Function ------------------ # 
    def Objective_Function(self):
        
        # Objective Function ----> Minimizing C+(S-M)
        self.Objective_Function = pd.DataFrame(columns=self.Variables)

        # Operating Cost (C)
        self.Operating_Cost = pd.DataFrame(columns=self.Variables)

        new_row = {}
        for flight in self.Flight_Schedule.index:
            for fleet in self.fleets.index:
                new_row[f"X{flight+1}_{fleet+1}"] = self.Flight_Schedule[f"e{fleet+1}"][flight]

        self.Operating_Cost = self.Operating_Cost._append(new_row, ignore_index = True).fillna(0)

        # Spilled Cost/Revenue (S)
        self.Spilled_Cost = pd.DataFrame(columns=self.Variables)

        new_row = {}
        counter = 0
        for Itenary in self.I_I_matrix_summation:
            elements = Itenary.split(" + ")
            for element in elements:
                if element != '0':
                    new_row[element] = self.Itenaries["Fare"][counter]
            
            counter += 1

        self.Spilled_Cost = self.Spilled_Cost._append(new_row, ignore_index = True).fillna(0)
        # self.Spilled_Cost.to_excel(f"Outputs/{self.save_folder}/Spilled_Cost_{self.save_folder}.xlsx", engine="openpyxl")
        

        # Recaptured Cost/Revenue (M)
        self.Recaptured_Cost = pd.DataFrame(columns=self.Variables)

        new_row = {}
        counter = 0
        for Itenary in self.I_I_matrix_recapture_summation:
            elements = Itenary.split(" + ")
            for element in elements:
                if element != "0":
                    for probability in self.recapture_probability.index:
                        if element == self.recapture_probability["Variable_Equivalent"][probability]:
                            new_row[f"{element}"] = self.recapture_probability["Probability"][probability] * self.Itenaries["Fare"][counter]
        
            counter += 1

        self.Recaptured_Cost = self.Recaptured_Cost._append(new_row, ignore_index = True).fillna(0)
        # self.Recaptured_Cost.to_excel(f"Outputs/{self.save_folder}/Recaptured_Cost_{self.save_folder}.xlsx", engine="openpyxl")

        # Append To Obective Function
        self.Objective_Function = self.Operating_Cost + (self.Spilled_Cost - self.Recaptured_Cost)
        self.Objective_Function.to_csv(f"Outputs/{self.save_folder}/Objective_Function_{self.save_folder}.csv")
        

                                             # ------------------ Create FAM Constraints ------------------ #
    def Create_iFAM_constraints(self):

         # -------------> Append Flight Interaction Constraints To Inequality Matrix & to Inequality RHS
        self.inequality_constraints = pd.concat([self.inequality_constraints, self.Flight_Interaction_constraints_matrix],ignore_index=True)
        self.inequality_rhs = np.concatenate((self.inequality_rhs, self.Flight_Interaction_constraints_rhs))

         # -------------> Append Demand Constraints To Inequality Matrix & to Inequality RHS
        self.inequality_constraints = pd.concat([self.inequality_constraints, self.Demand_constraint_Matrix],ignore_index=True)
        self.inequality_rhs = np.concatenate((self.inequality_rhs, self.Demand_constraint_rhs))

                                             # ------------------ Optimization ------------------ #
    def optimize_iFAM(self):
        
        objective_function = self.Objective_Function.to_numpy()
        inequality_constraints = self.inequality_constraints.to_numpy()
        equality_constraints = self.equality_constraints.to_numpy()
        inequality_rhs = self.inequality_rhs
        equality_rhs = self.equality_rhs
        bounds = []
        bounds.extend(self.coverage_bounds)
        bounds.extend(self.resource_bounds)
        bounds.extend(self.Balnace_bounds)
        bounds.extend(self.spill_recapture_variables_bounds)
        
        # using scipy linprog to optimize
        self.result = linprog(objective_function,
                         A_ub = inequality_constraints,
                         b_ub = inequality_rhs,
                         A_eq = equality_constraints,
                         b_eq = equality_rhs,
                         bounds = bounds, 
                              
        )

        #---> Unconstrained revenue (R) R = fare_p*D_p
        counter = 0
        self.Flight_Schedule["Fare"] = 0
        for flight_Demand in self.Unconstrained_Demand_Qf:
            self.Flight_Schedule.loc[counter, "Fare"] = self.Itenaries.loc[counter, "Fare"]
            counter += 1

        self.Unconstrained_Revenue = 0
        for Itenary in self.Itenaries.index:
            self.Unconstrained_Revenue += self.Itenaries.loc[Itenary, "Fare"] * self.Itenaries.loc[Itenary, "Passenger Demand"]
        


        # Save and Display Solution Output
        self.output = pd.DataFrame(columns=self.Variables)
        self.output.loc[len(self.output)] = self.result.x
        self.output.to_csv(f"Outputs/{self.save_folder}/Optimization Results ({self.save_folder}).csv")
        self.Minimum_Cost = self.result.fun 
        self.Maximum_profit = self.Unconstrained_Revenue - self.Minimum_Cost  

        print("##################################################")
        print("##                     iFAM                     ##")
        print("##################################################")
        print(f"Minimum Cost    (iFAM): {self.Minimum_Cost}")
        print(f"Maximum Profit  (iFAM): {self.Maximum_profit}\n")

    
                                        # -------------------------     Optimization (Pyomo)   ------------------------- # 
    def optimize_iFAM_pyomo(self):

        model = pe.ConcreteModel()
 
        model.x = pe.Var(self.Dummy_variables, within = pe.Binary)
        model.RON = pe.Var(self.RON_variables, within = pe.NonNegativeIntegers)
        model.y = pe.Var(self.ground_link_variables, within = pe.NonNegativeIntegers)
        model.t = pe.Var(self.spill_recapture_variables, within = pe.NonNegativeIntegers)


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
            model.spill_recapture_Demand_constraints.add(expr= spill_recapture_Demand_constraints <= self.Itenaries.loc[counter,"Passenger Demand"] )
            counter += 1

        # Flight Interaction Constrainst
        
                                                                    # Preform Calculations  
                                                    # (spill - recapture + Seat_f >= Qf)  multiply by negative
                                                    # (-spill + recapture - Seat_f <= -Qf)  

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
            
            
            flight_interaction_sum = sum([recapture_demand, spill_demand, seat_f])
            model.flight_interaction_constraints.add(expr= flight_interaction_sum   <= -self.Unconstrained_Demand_Qf[flight])

        # Objective_Function
        
        # Minimize Cost ----> C+(S-M) 

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

        #---> Unconstrained revenue (R) R = fare_p*D_p
        counter = 0
        self.Flight_Schedule["Fare"] = 0
        for flight_Demand in self.Unconstrained_Demand_Qf:
            self.Flight_Schedule.loc[counter, "Fare"] = self.Itenaries.loc[counter, "Fare"]
            counter += 1

        self.Unconstrained_Revenue = 0
        for Itenary in self.Itenaries.index:
            self.Unconstrained_Revenue += self.Itenaries.loc[Itenary, "Fare"] * self.Itenaries.loc[Itenary, "Passenger Demand"]

       
        objective_fun = sum([operational_costs, spilled_costs, recaptured_costs])

        model.objective = pe.Objective(sense = pe.minimize, expr = objective_fun)

        with open(f"Outputs/{self.save_folder}/model_{self.save_folder}_output.txt", "w") as output_file:
             model.pprint(output_file)


        solver = po.SolverFactory("gurobi")     
        results = solver.solve(model, tee=False)

        print("##################################################")
        print("##                     iFAM                     ##")
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

                                             #       Retuen Outputs      #      
    def run_analysis_iFAM(self):

        self.save_folder = "iFAM"
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
        self.Flight_Interaction_constraints()
        self.Demand_constraints()
        self.Objective_Function()
        self.Create_FAM_constraints()
        self.Create_iFAM_constraints()
        self.check_optional_flights()
        self.optimize_iFAM_pyomo()