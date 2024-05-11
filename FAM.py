# FFFFFFFFFFFFFFFFFFFFFF      AAA               MMMMMMMM               MMMMMMMM
# F::::::::::::::::::::F     A:::A              M:::::::M             M:::::::M
# F::::::::::::::::::::F    A:::::A             M::::::::M           M::::::::M
# FF::::::FFFFFFFFF::::F   A:::::::A            M:::::::::M         M:::::::::M
#   F:::::F       FFFFFF  A:::::::::A           M::::::::::M       M::::::::::M
#   F:::::F              A:::::A:::::A          M:::::::::::M     M:::::::::::M
#   F::::::FFFFFFFFFF   A:::::A A:::::A         M:::::::M::::M   M::::M:::::::M
#   F:::::::::::::::F  A:::::A   A:::::A        M::::::M M::::M M::::M M::::::M
#   F:::::::::::::::F A:::::A     A:::::A       M::::::M  M::::M::::M  M::::::M
#   F::::::FFFFFFFFFFA:::::AAAAAAAAA:::::A      M::::::M   M:::::::M   M::::::M
#   F:::::F         A:::::::::::::::::::::A     M::::::M    M:::::M    M::::::M
#   F:::::F        A:::::AAAAAAAAAAAAA:::::A    M::::::M     MMMMM     M::::::M
# FF:::::::FF     A:::::A             A:::::A   M::::::M               M::::::M
# F::::::::FF    A:::::A               A:::::A  M::::::M               M::::::M
# F::::::::FF   A:::::A                 A:::::A M::::::M               M::::::M
# FFFFFFFFFFF  AAAAAAA                   AAAAAAAMMMMMMMM               MMMMMMMM

import pandas as pd
import pyomo.environ as pe
import pyomo.opt as po
from scipy.optimize import linprog
import numpy as np
import json
import openpyxl


class FAM:
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
        self.inequality_constraints = pd.DataFrame()
        self.equality_constraints = pd.DataFrame()
        self.inequality_rhs = np.array([])
        self.equality_rhs = np.array([])

                                                    #  Read Flight Secdule #
    def read_flight_schedule(self):
        self.Flight_Schedule = pd.DataFrame(pd.read_excel(self.excel_file, sheet_name="Flights", index_col=None, engine="openpyxl"))
        self.Number_Flights = len(self.Flight_Schedule)
    
                                                    #  Read Fleets Data   #
    def read_fleets(self):
        self.fleets = pd.DataFrame(pd.read_excel(self.excel_file, sheet_name="Fleets", index_col=None, engine="openpyxl"))
        self.Number_fleets = len(self.fleets)

                                                    #  Read Itenaries Data   #
    def read_Itenaries(self):
        self.Itenaries = pd.DataFrame(pd.read_excel(self.excel_file, sheet_name="Itenaries", index_col=None, engine="openpyxl"))
        self.Number_Itenaries = len(self.Itenaries)
        
        
                                                    #  Identify Stations  #
    def identify_stations(self):
        for flight in self.Flight_Schedule.index:
            if self.Flight_Schedule["from"][flight] not in self.Stations:
                self.Stations.append(self.Flight_Schedule["from"][flight])
        
        self.Number_Stations = len(self.Stations)
                       
                                                    #   Create Null Node  #
    def create_empty_node(self, node_number, station_name):
        new_node = {
            "Node" : node_number,
            "Station" : station_name,
            "Inbound" : [],
            "Outbound" : []
        } 
        return new_node
                                            
                                            #  Update Null Nodes (Interconnection Nodes) #
    def create_nodes(self):
        for station in self.Stations:
            Number_of_Nodes = len(self.Nodes)
            new_node = self.create_empty_node(Number_of_Nodes + 1, station)

                                                    #   Create Station Schedule   #
            station_schedule = self.Flight_Schedule[(self.Flight_Schedule["from"] == station) | (self.Flight_Schedule["to"] == station)]
            station_schedule = station_schedule.reset_index(drop=True)
            station_schedule["Station Time"] = ""
            station_flights = len(station_schedule)

            for flight in station_schedule.index:
                if station == station_schedule["to"][flight]:
                    station_schedule.loc[flight, "Station Time"] = station_schedule.loc[flight, "Arrival"]
                elif station == station_schedule["from"][flight]:
                    station_schedule.loc[flight, "Station Time"] = station_schedule.loc[flight, "Departure"]

            station_schedule = station_schedule.sort_values(['Station Time'])
            station_schedule = station_schedule.reset_index(drop=True)
                                               
                                                #    Open Nodes and Append Flights   #
            for flight in station_schedule.index:
                if station == station_schedule["to"][flight]:
                    if flight == 0:
                        new_node["Inbound"].append(f"RON_{station}")
                        new_node["Inbound"].append(int(station_schedule["flight number"][flight]))
                    elif flight == station_flights - 1:
                        new_node["Inbound"].append(int(station_schedule["flight number"][flight]))
                        new_node["Outbound"].append(f"RON_{station}")
                    else:
                        new_node["Inbound"].append(int(station_schedule["flight number"][flight]))
                elif station == station_schedule["from"][flight]:
                    if flight == 0:
                        new_node["Inbound"].append(f"RON_{station}")
                        if flight != (station_flights - 1) and station == station_schedule["to"][flight + 1]:
                            new_node["Outbound"].append(int(station_schedule["flight number"][flight]))
                            self.Nodes.append(new_node)
                            Number_of_Nodes = len(self.Nodes)
                            new_node = self.create_empty_node(Number_of_Nodes + 1, station)
                        else:
                            new_node["Outbound"].append(int(station_schedule["flight number"][flight]))
                    elif flight == station_flights - 1:
                        if flight != (station_flights - 1) and station == station_schedule["to"][flight + 1]:
                            new_node["Outbound"].append(int(station_schedule["flight number"][flight]))
                            self.Nodes.append(new_node)
                            Number_of_Nodes = len(self.Nodes)
                            new_node = self.create_empty_node(Number_of_Nodes + 1, station)
                        else:
                            new_node["Outbound"].append(int(station_schedule["flight number"][flight]))
                        new_node["Outbound"].append(f"RON_{station}")
                    else:
                        if flight != (station_flights - 1) and station == station_schedule["to"][flight + 1]:
                            new_node["Outbound"].append(int(station_schedule["flight number"][flight]))
                            self.Nodes.append(new_node)
                            Number_of_Nodes = len(self.Nodes)
                            new_node = self.create_empty_node(Number_of_Nodes + 1, station)
                        else:
                            new_node["Outbound"].append(int(station_schedule["flight number"][flight]))

            self.Nodes.append(new_node)
            self.Number_of_Nodes = len(self.Nodes)
        
                                                #   Create Ground Links   #
    
        for Station in self.Stations:
            for node in range(self.Number_of_Nodes - 1):
                if self.Nodes[node]["Station"] == Station and self.Nodes[node + 1]["Station"] == Station:
                        ground_link = f"y({node + 1}-{node + 2})"
                        self.ground_links.append(ground_link)
                        self.Nodes[node]["Outbound"].append(ground_link)
                        self.Nodes[node + 1]["Inbound"].insert(0, ground_link)

            # check ground links = Number_of_Nodes - RONs
            if len(self.ground_links) == self.Number_of_Nodes - self.Number_Stations: 
                self.Number_of_ground_links = len(self.ground_links)

            # File path where you want to save the JSON data
            file_path = "Outputs/Nodes.json"
            # Write the data to the JSON file
            with open(file_path, 'w') as json_file:
                json.dump(self.Nodes, json_file, indent = 4)

            self.Number_of_Variables = (self.Number_Flights + self.Number_Stations + self.Number_of_ground_links)*self.Number_fleets 
            self.Number_of_Coverage_Constraints = self.Number_Flights
            self.Number_of_Resource_Constraints = self.Number_fleets
            self.Number_of_Balance_Constraints = self.Number_of_Nodes*self.Number_fleets

                                            # Construct Eqns (Matrices) #
    def  constraints_matrices(self):
        
        if len(self.Variables) == 0:
            # Coverage Vairables
            for flight in self.Flight_Schedule["flight number"]:
                for fleet in range(self.Number_fleets):
                    self.Variables.append(f"X{flight}_{fleet+1}")
                    self.Dummy_variables.append(f"X{flight}_{fleet+1}")

            # Resources Variables
            for fleet in range(self.Number_fleets):
                for station in self.Stations:
                    self.Variables.append(f"RON_{station}_{fleet+1}")
                    self.RON_variables.append(f"RON_{station}_{fleet+1}")

            # Balance (Interconnection) Variables
            for g_link in self.ground_links:
                for fleet in range(self.Number_fleets):
                    self.Variables.append(f"{g_link}_{fleet+1}")
                    self.ground_link_variables.append(f"{g_link}_{fleet+1}")

        # Coverage Matrix
        self.coverage_matrix = pd.DataFrame(columns=self.Variables)
        new_row = {}
                
        for flight in self.Flight_Schedule["flight number"]:
            for fleet in range(self.Number_fleets):
                new_row[f"X{flight}_{fleet+1}"] = 1
            
            
            self.coverage_matrix = self.coverage_matrix._append(new_row, ignore_index=True)
            new_row.clear()
        
        self.coverage_matrix = self.coverage_matrix.fillna(0)
        self.coverage_matrix.to_csv(f"Outputs/{self.save_folder}/coverage_matrix_{self.save_folder}.csv")
        self.coverage_rhs = np.ones(self.Number_Flights)
        self.coverage_bounds = [(0,1) for i in range(len(self.Dummy_variables))]  

        # Resource Matrix
        self.resource_matrix = pd.DataFrame(columns=self.Variables)
        new_row = {}
        self.resource_bounds = []
        self.resource_rhs = []

        for fleet in range(self.Number_fleets):
            for station in self.Stations:
                new_row[f"RON_{station}_{fleet+1}"] = 1
                # max fleet in ron 
                self.resource_bounds.append((0,self.fleets["size"][fleet]))

            
            self.resource_matrix = self.resource_matrix._append(new_row, ignore_index=True)
            self.resource_rhs.append(self.fleets["size"][fleet])
            new_row.clear()
        
        self.resource_rhs = np.array(self.resource_rhs)
        self.resource_matrix = self.resource_matrix.fillna(0)
        self.resource_matrix.to_csv(f"Outputs/{self.save_folder}/resource_matrix_{self.save_folder}.csv")


        # Balance (Interconnection) Matrix
        self.Balance_matrix = pd.DataFrame(columns=self.Variables)
        new_row = {}

        for node in range(self.Number_of_Nodes):
            for fleet in range(self.Number_fleets):
                
                #Inbounds
                for In in self.Nodes[node]["Inbound"]: 
                    if type(In) is str:
                        new_row[f"{In}_{fleet+1}"] = 1
                    else:
                        new_row[f"X{In}_{fleet+1}"] = 1
                #Outbounds
                for Out in self.Nodes[node]["Outbound"]: 
                    if type(Out) is str:
                        new_row[f"{Out}_{fleet+1}"] = -1
                    else:
                        new_row[f"X{Out}_{fleet+1}"] = -1

            
                self.Balance_matrix = self.Balance_matrix._append(new_row, ignore_index=True).fillna(0)
                new_row.clear()
        
        self.Balance_matrix = self.Balance_matrix.fillna(0)
        self.Balance_matrix.to_csv(f"Outputs/{self.save_folder}/Balance_matrix_{self.save_folder}.csv")
        self.Balance_rhs = np.zeros(len(self.Balance_matrix))
        self.Balnace_bounds = [(0,None) for i in self.ground_link_variables]
        
                                             #  Profit Calaculation  #
    def  profit_calculation(self):
        
        row = {}
        self.profit_function = pd.DataFrame(columns=self.Variables)
        
        # if sheet includes flight based Fares with cost/mile given
        if 'Fare' in self.Flight_Schedule.columns:
            for flight in self.Flight_Schedule.index:
                for fleet in range(self.Number_fleets):
                    row[f"X{flight+1}_{fleet+1}"] =  self.Flight_Schedule["Fare"][flight] * min(self.fleets["capacity"][fleet], self.Flight_Schedule["Demand"][flight]) - self.fleets["cost_per_mile"][fleet] * self.Flight_Schedule["Distance"][flight]

            self.profit_function = self.profit_function._append(row, ignore_index=True).fillna(0)
            self.profit_function.to_csv(f"Outputs/{self.save_folder}/profit_function_{self.save_folder}.csv")
        
        # if sheet includes Itenary Based Fares with cost of flight given
        else :
            for flight in self.Flight_Schedule.index:
                for fleet in range(self.Number_fleets):
                    row[f"X{flight+1}_{fleet+1}"] =  self.Itenaries["Fare"][flight] * min(self.fleets["capacity"][fleet], self.Itenaries["Passenger Demand"][flight]) -  self.Flight_Schedule[f"e{fleet+1}"][flight]

            self.profit_function = self.profit_function._append(row, ignore_index=True).fillna(0)
            self.profit_function.to_csv(f"Outputs/{self.save_folder}/profit_function.csv")

                                             #  Create FAM Constraints  #
    def Create_FAM_constraints(self):
        # -------------> Append Coverage Constraints To Equality Matrix & to Equality RHS
        self.equality_constraints = pd.concat([self.equality_constraints, self.coverage_matrix],ignore_index=True)
        self.equality_rhs = np.concatenate((self.equality_rhs, self.coverage_rhs))

        # -------------> Append Balance Constraints To Equality Matrix & to Equality RHS
        self.equality_constraints = pd.concat([self.equality_constraints, self.Balance_matrix],ignore_index=True)
        self.equality_rhs = np.concatenate((self.equality_rhs, self.Balance_rhs))

         # -------------> Append Resource Constraints To Inequality Matrix & to Inequality RHS
        self.inequality_constraints = pd.concat([self.inequality_constraints, self.resource_matrix],ignore_index=True)
        self.inequality_rhs = np.concatenate((self.inequality_rhs, self.resource_rhs))


                                             #  Check Optional Flights and Update Constraints  #
    def check_optional_flights(self):
        workbook = openpyxl.load_workbook(self.excel_file)
        sheet_name = "Optional_Flights"
        if sheet_name in workbook.sheetnames:
            self.Optional_Flights = pd.read_excel(self.excel_file, sheet_name=sheet_name, index_col=None, engine="openpyxl")
            if len(self.Optional_Flights) != 0:
                self.indices_to_remove = []
                for Flight in self.Optional_Flights["flight number"]:
                    self.indices_to_remove.append(Flight-1)
                    new_inequality_constraint = self.coverage_matrix.iloc[Flight-1]
                    new_inequality_rhs = self.coverage_rhs[Flight-1]
                    self.inequality_constraints = self.inequality_constraints._append(new_inequality_constraint, ignore_index=True)
                    self.inequality_rhs = np.append(self.inequality_rhs, new_inequality_rhs)
                
                self.equality_constraints = self.equality_constraints.drop(self.indices_to_remove)
                self.equality_rhs = np.delete(self.equality_rhs, self.indices_to_remove)
                

            else:
                return 0
        
        else: 
            return 0


                                             #      Optimization     #
    def optimize_FAM(self):
        
        objective_function = self.profit_function.to_numpy()
        inequality_constraints = self.inequality_constraints.to_numpy()
        equality_constraints = self.equality_constraints.to_numpy()
        inequality_rhs = self.inequality_rhs
        equality_rhs = self.equality_rhs
        bounds = []
        bounds.extend(self.coverage_bounds)
        bounds.extend(self.resource_bounds)
        bounds.extend(self.Balnace_bounds)
        
        
        # using scipy linprog to optimize
        self.result = linprog(-objective_function,
                         A_ub = inequality_constraints,
                         b_ub = inequality_rhs,
                         A_eq = equality_constraints,
                         b_eq = equality_rhs,
                         bounds = bounds,
                                                    
        )

        # Save and Display Solution Output
        self.output = pd.DataFrame(columns=self.Variables)
        self.output.loc[len(self.output)] = self.result.x
        self.output.to_csv(f"Outputs/{self.save_folder}/Optimization Results ({self.save_folder}).csv")
        self.Maximum_profit = -self.result.fun 

        print("##################################################")
        print("##                      FAM                     ##")
        print("##################################################")
        print(f"Maximum Profit (FAM): {self.Maximum_profit}\n")
        

    
                                             #      Optimization (Pyomo)     # 
    def optimize_FAM_pyomo(self):

        model = pe.ConcreteModel()
 
        model.x = pe.Var(self.Dummy_variables, within = pe.Binary)
        model.RON = pe.Var(self.RON_variables, within = pe.NonNegativeIntegers)
        model.y = pe.Var(self.ground_link_variables, within = pe.NonNegativeIntegers)

        
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
                
                
                # Constraints
                flight_sum = sum(model.x[f"{key}"]*value for key, value in new_x.items())
                Ron_sum = sum(model.RON[f"{key}"]*value for key, value in new_ron.items())
                y_sum = sum(model.y[f"{key}"]*value for key, value in new_y.items())
                balance_sum = sum([flight_sum, Ron_sum, y_sum])

                model.balance_constraints.add(expr= balance_sum == 0)

                new_x.clear()
                new_ron.clear()
                new_y.clear()

        # Objective_Function
        
        row = {}
        # if sheet includes flight based Fares with cost/mile given
        if 'Fare' in self.Flight_Schedule.columns:
            for flight in self.Flight_Schedule.index:
                for fleet in range(self.Number_fleets):
                    row[f"X{flight+1}_{fleet+1}"] =  self.Flight_Schedule["Fare"][flight] * min(self.fleets["capacity"][fleet], self.Flight_Schedule["Demand"][flight]) - self.fleets["cost_per_mile"][fleet] * self.Flight_Schedule["Distance"][flight]

            self.profit_function = self.profit_function._append(row, ignore_index=True).fillna(0)
            self.profit_function.to_csv(f"Outputs/{self.save_folder}/profit_function_{self.save_folder}.csv")
        
        # if sheet includes Itenary Based Fares with cost of flight given
        else :
            for flight in self.Flight_Schedule.index:
                for fleet in range(self.Number_fleets):
                    row[f"X{flight+1}_{fleet+1}"] =  self.Itenaries["Fare"][flight] * min(self.fleets["capacity"][fleet], self.Itenaries["Passenger Demand"][flight]) -  self.Flight_Schedule[f"e{fleet+1}"][flight]


        objective_fun = sum(model.x[f"{key}"]*value for key, value in row.items())

        model.objective = pe.Objective(sense = pe.maximize, expr = objective_fun)


        solver = po.SolverFactory("gurobi")     
        results = solver.solve(model, tee=False)

        print("##################################################")
        print("##                      FAM                     ##")
        print("##################################################")
        print(f"Maxium Profit (FAM) : {pe.value(model.objective)}")

        for var in self.Dummy_variables:
            print(f"{var} : {pe.value(model.x[f"{var}"])}")

        for var in self.RON_variables:
            print(f"{var} : {pe.value(model.RON[f"{var}"])}")

        for var in self.ground_link_variables:
            print(f"{var} : {pe.value(model.y[f"{var}"])}")
        
        # model.pprint()
    
    
                                             #       Retuen Outputs      #                  
    def run_analysis_FAM(self):

        self.save_folder = "FAM"
        self.read_flight_schedule()
        self.read_fleets()
        self.read_Itenaries()
        self.identify_stations()
        self.create_nodes()
        self.constraints_matrices()
        self.profit_calculation()
        self.Create_FAM_constraints()
        self.check_optional_flights()
        self.optimize_FAM_pyomo()
        
        
