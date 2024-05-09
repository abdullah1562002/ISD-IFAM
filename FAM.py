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
# import pyomo.environ as pe
# import pyomo.opt as po
from scipy.optimize import linprog
import numpy as np
import json


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
        

                                                    # Read Flight Secdule #
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
        self.coverage_matrix.to_excel(f"Outputs/{self.save_folder}/coverage_matrix_{self.save_folder}.xlsx", engine="openpyxl")
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
        
        self.resource_matrix = self.resource_matrix.fillna(0)
        self.resource_matrix.to_excel(f"Outputs/{self.save_folder}/resource_matrix_{self.save_folder}.xlsx", engine="openpyxl")
        
        

        
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
        self.Balance_matrix.to_excel(f"Outputs/{self.save_folder}/Balance_matrix_{self.save_folder}.xlsx", engine="openpyxl")
        self.Balance_rhs = np.zeros(len(self.Balance_matrix))
        self.Balnace_bounds = [(0,None) for i in self.ground_link_variables]
        
        # Equality Matrix
        self.coverage_Balance_matrix = pd.concat([self.coverage_matrix, self.Balance_matrix], ignore_index = True)
        self.coverage_Balance_rhs = np.append(self.coverage_rhs, self.Balance_rhs)
        
        # Constraint_Matrix
        self.constraints_matrix = pd.concat([self.coverage_matrix, self.Balance_matrix, self.resource_matrix], ignore_index = True).fillna(0)
        self.constraints_matrix.to_excel(f"Outputs/{self.save_folder}/constraints_matrix_{self.save_folder}.xlsx", engine="openpyxl")

    
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
            self.profit_function.to_excel(f"Outputs/{self.save_folder}/profit_function_{self.save_folder}.xlsx", engine="openpyxl")
        
        # if sheet includes Itenary Based Fares with cost of flight given
        else :
            for flight in self.Flight_Schedule.index:
                for fleet in range(self.Number_fleets):
                    row[f"X{flight+1}_{fleet+1}"] =  self.Itenaries["Fare"][flight] * min(self.fleets["capacity"][fleet], self.Itenaries["Passenger Demand"][flight]) -  self.Flight_Schedule[f"e{fleet+1}"][flight]

            self.profit_function = self.profit_function._append(row, ignore_index=True).fillna(0)
            self.profit_function.to_excel("Outputs/FAM/profit_function.xlsx", engine="openpyxl")


                                             #      Optimization     #
    def optimize_FAM(self):
        
        objective_function = self.profit_function.to_numpy()
        inequality_constraints = self.resource_matrix.to_numpy()
        equality_constraints = self.coverage_Balance_matrix.to_numpy()
        inequality_rhs = self.resource_rhs
        equality_rhs = self.coverage_Balance_rhs
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

        print(f"Maximum Profit (FAM): {self.Maximum_profit}")
        

    
                                             #      Optimization (Pyomo)     # 
                                             #        Future Upgrade         #
    # def optimize_pyomo(self):

    #     model = pe.ConcreteModel()
 
    #     model.x = pe.Var(self.Dummy_variables, within = pe.Binary)
    #     model.RON = pe.Var(self.RON_variables, within = pe.NonNegativeIntegers)
    #     model.ground_links = pe.Var(self.ground_link_variables, within = pe.NonNegativeIntegers)

    #     # Coverage Constraints 
    #     model.coverage_constraints = pe.ConstraintList()
    #     for flight in self.Flight_Schedule["flight number"]:
    #         coverage_sum = sum([model.x[f"X{flight}_{fleet+1}"] for fleet in range(self.Number_fleets)]) 
    #         model.coverage_constraints.add(coverage_sum == 1)
        
    #     # Resource Constraints
    #     model.resource_constraints = pe.ConstraintList()
    #     for fleet in range(self.Number_fleets):
    #         resource_sum = sum([model.RON[f"RON_{station}_{fleet+1}"] for station in self.Stations])
    #         model.resource_constraints.add(resource_sum == self.fleets["size"][fleet])
        
    #     # Balance Constraints
    #     model.balance_constraints = pe.ConstraintList()
        
    #     new_row = {}
    #     for node in range(self.Number_of_Nodes):
    #         for fleet in range(self.Number_fleets):
                
    #             #Inbounds
    #             for In in self.Nodes[node]["Inbound"]: 
    #                 if type(In) is str:
    #                     new_row[f"{In}_{fleet+1}"] = 1
    #                 else:
    #                     new_row[f"X{In}_{fleet+1}"] = 1
    #             #Outbounds
    #             for Out in self.Nodes[node]["Outbound"]: 
    #                 if type(Out) is str:
    #                     new_row[f"{Out}_{fleet+1}"] = -1
    #                 else:
    #                     new_row[f"X{Out}_{fleet+1}"] = -1

                
    #     model.pprint()
    
    
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
        self.optimize_FAM()

        
    

class iFAM(FAM):
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
        
        
                             #       Calculate Demand (Change Demand from Itenary to Flight Based)     #  
    def Calculate_Unconstrained_Demand_Qf(self):
        self.Itenaries_var = [f"I{I_No+1}" for I_No in range(len(self.Itenaries))]
        self.Flights_var = [f"F{F_No+1}" for F_No in range(len(self.Flight_Schedule.index))]
        self.flight_itenary_incidence_Matrix = pd.DataFrame(columns = self.Itenaries_var)
        
        # Flight Itenary Incidence Matrix
        new_row = {}
        for flight in self.Flight_Schedule["flight number"]:
            for Itenary in self.Itenaries.index:
                if str(flight) in list(str(self.Itenaries["Flights"][Itenary]).split(",")):
                    new_row[f"I{Itenary + 1}"] = 1
            self.flight_itenary_incidence_Matrix = self.flight_itenary_incidence_Matrix._append(new_row, ignore_index=True).fillna(0)
            new_row.clear()

        self.flight_itenary_incidence_Matrix.index = self.Flights_var

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
        self.SEAT_matrix.to_excel("Outputs/iFAM/SEAT_f_matrix.xlsx", engine="openpyxl")

                                             #       Calculate Spilled and Recaptured Demands  #  
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
                        if  itenary_1 != itenary_2 and ref_Itenary[0] == check_Itenary[0] and ref_Itenary[-1] == check_Itenary[-1] and  len(ref_Itenary) <= 2:
                            new_row[f"I{itenary_2 + 1}"] = f"t({itenary_1 + 1},{itenary_2 + 1})"
                        else:
                            new_row[f"I{len(self.Itenaries_var)+1}"] = f"t({itenary_1 + 1},{len(self.Itenaries_var)+1})"

                        # Create Spill Recapture Variables t(p,r)
                        self.spill_recapture_variables.append(f"t({itenary_1 + 1},{itenary_2 + 1})")
                    else: 
                        self.spill_recapture_variables.append(f"t({itenary_1 + 1},{itenary_2 + 1})")
            
            self.itenary_itenary_spill_recapture_Matrix = self.itenary_itenary_spill_recapture_Matrix._append(new_row, ignore_index=True).fillna(0)
            new_row.clear()
        
        # spill recapture variables bounds
        self.spill_recapture_variables_bounds = [(0,None) for i in self.spill_recapture_variables]
        
        ###############################################  Spilled Demand #########################################

        # Take Summation of each Inteary row (symbolic)    
        self.I_I_matrix_summation = pd.Series()
        for Itenary in self.itenary_itenary_spill_recapture_Matrix.index:
            new_list = self.itenary_itenary_spill_recapture_Matrix.iloc[Itenary]
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
            self.spilled_Demand_Matrix.to_excel(f"Outputs/{self.save_folder}/spilled_Demand_Matrix.xlsx", engine="openpyxl")
        else : 
            print("Dimensions Dosent Match")


        ###############################################  Recapture Demand #########################################

        # Take Summation of each Inteary Column (symbolic) 
        self.I_I_matrix_recapture = self.itenary_itenary_spill_recapture_Matrix.iloc[:, :-1]
        self.I_I_matrix_recapture_summation = pd.Series()
        
        for Itenary in self.I_I_matrix_recapture.columns:
            new_list = self.itenary_itenary_spill_recapture_Matrix.loc[:, Itenary]
            row_sum  = [variable for variable in new_list if variable != 0]
            row_sum = " + ".join(row_sum)
            
            self.I_I_matrix_recapture_summation = self.I_I_matrix_recapture_summation._append(pd.Series(row_sum), ignore_index=True)

        self.I_I_matrix_recapture_summation.index = self.Itenaries_var
        self.F_I_I_Matrix_recapture_symbolic = self.flight_itenary_incidence_Matrix
        
        
        # Multiply the Summation of t(p,r) to Flight - Itenary - Incidence Matrix
        for Itenary in self.F_I_I_Matrix_recapture_symbolic.columns:
            self.F_I_I_Matrix_recapture_symbolic.loc[:, Itenary] = self.F_I_I_Matrix_recapture_symbolic.loc[:, Itenary].replace(1, self.I_I_matrix_recapture_summation.loc[Itenary])
            self.F_I_I_Matrix_recapture_symbolic = self.F_I_I_Matrix_recapture_symbolic.replace('',0)
        
        
        
        # Take Summation of each Flight row (symbolic) 
        self.F_I_I_matrix_recapture_summation = pd.Series()
        for Flight in self.F_I_I_Matrix_recapture_symbolic.index:
            new_list = self.F_I_I_Matrix_recapture_symbolic.loc[Flight]
            row_sum  = [variable for variable in new_list if variable != 0]
            row_sum = " + ".join(row_sum)

            self.F_I_I_matrix_recapture_summation = self.F_I_I_matrix_recapture_summation._append(pd.Series(row_sum), ignore_index=True)

        self.F_I_I_matrix_recapture_summation = self.F_I_I_matrix_recapture_summation.replace('',0)      
        self.F_I_I_matrix_recapture_summation.index = self.Flights_var

        # Recapture Demand Matrix 
        self.recapture_Demand_Matrix = pd.DataFrame(columns=self.Variables)
        self.recapture_probability = pd.read_excel(self.excel_file, sheet_name="Recapture_Probability", index_col=None, engine="openpyxl")
        self.recapture_probability["Variable_Equivalent"] = self.recapture_probability["Variable"].replace("b","t", regex=True)
        
        
        
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
            self.recapture_Demand_Matrix.to_excel(f"Outputs/{self.save_folder}/recaptured_Demand_Matrix.xlsx", engine="openpyxl")
        else : 
            print("Dimensions Dosent Match")
        
                                            #   Create Flight Interaction Constraints Matrix  #  
    def Flight_Interaction_constraints(self):
        
                        # Preform Calculations  
        # (spill - recapture + Seat_f >= Qf)  multiply by negative
        # (-spill + recapture - Seat_f <= -Qf)  multiply by negative

        spill = self.spilled_Demand_Matrix.to_numpy()
        recapture = self.recapture_Demand_Matrix.to_numpy()
        Seat_f = self.SEAT_matrix.to_numpy()

        self.Flight_Interaction_constraints_matrix = pd.DataFrame(-spill + recapture -Seat_f, columns=self.Variables)
        self.Flight_Interaction_constraints_matrix.to_excel(f"Outputs/{self.save_folder}/Flight_Interaction_constraints_matrix_{self.save_folder}.xlsx", engine="openpyxl")
        self.Flight_Interaction_constraints_rhs = - self.Unconstrained_Demand_Qf

                                                #   Create Demand Constraints Matrix  # 
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
            self.Demand_constraint_Matrix.to_excel(f"Outputs/{self.save_folder}/Demand_constraint_Matrix_{self.save_folder}.xlsx", engine="openpyxl")
        else : 
            print("Dimensions Dosent Match")

        # Demand Constraint RHS
        self.Demand_constraint_rhs = np.array(self.Itenaries["Passenger Demand"])

                                                 #   Create Objective Function  # 
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
        self.Objective_Function.to_excel(f"Outputs/{self.save_folder}/Objective_Function_{self.save_folder}.xlsx", engine="openpyxl")
        



    
        






    
    
    
    
    
    
    
    def run_analysis_iFAM(self):

        self.save_folder = "iFAM"
        self.read_flight_schedule()
        self.read_fleets()
        self.read_Itenaries()
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
        # self.optimize_iFAM()
    
        
