(* ::Package:: *)

(* ::Title:: *)
(*StoChemSim Package*)


(* ::Text:: *)
(*Installation instructions: Drop this package file, CRNSimulator.m, and the "StoChemSimSequential" directory in the directory opened by SystemOpen@FileNameJoin[{$UserBaseDirectory, "Applications"}]*)


Needs["CRNSimulator`"];
BeginPackage["StoChemSim`", {"CRNSimulator`", "CCompilerDriver`"}];


SimulateDirectSSA::usage =
"result = SimulateDirectSSA[{rxn1, rxn2, ..., conc1, conc2, ...}, Options]
Simulates the given reaction system via Gillespie Direct SSA
Backend is optimized in C++ for computational efficiency
Options include:
	timeEnd (Real), default = Infinity, ending time of simulation (if useIter = False)
	iterEnd (Integer), default = Infinity, ending iteration of simulation (if useIter = True)
	useIter (Boolean), default = False, setting to True uses iterEnd instead of timeEnd
	statesOnly (Boolean), default = False, setting to True avoids simulating reaction times
	finalOnly (Boolean), default = False, only records final state to conserve memory
	outputTS (Boolean), default = True, setting to True outputs result as TimeSeries, setting to False outputs result as List";
	
SimulateBoundedTauLeaping::usage =
"result = SimulateBoundedTauLeaping[{rxn1, rxn2, ..., conc1, conc2, ...}, Options]
Simulates the given reaction system via Soloveichik Bounded Tau Leaping
Backend is optimized in C++ for computational efficiency
Options include:
	timeEnd (Real), default = Infinity, ending time of simulation (if useIter = False)
	iterEnd (Integer), default = Infinity, ending iteration of simulation (if useIter = True)
	useIter (Boolean), default = False, setting to True uses iterEnd instead of timeEnd
	finalOnly (Boolean), default = False, only records final state to conserve memory
	outputTS (Boolean), default = True, setting to True outputs result as TimeSeries, setting to False outputs result as List
	rho (Real), default = 0.25, threshold between 0 and 1 used in calculating epsilon (ignored if epsilon given)
	epsilon (Real), default = 3/(4*numRxns)*(1 - Sqrt[(1 + rho/9)/(1 + rho)]),
		threshold between 0 and 1 used in calculating firing bounds for each reaction
Note: epsilon must be well chosen to guarantee correct functionality of Bounded Tau Leaping
A rho value can be provided instead, and a valid epsilon will be calculated given that all reactions are unimolecular or bimolecular
See https://arxiv.org/pdf/0803.1030.pdf for more details";

PlotLastSimulation::usage =
"PlotLastSimulation[Options]
Plots the last simulation ran
Uses same Options as ListLinePlot";

GetRuntimeInfo::usage=
"Retrieves runtime data collected from previous simulation as Association object
Frontend runtime includes time spent converting rxnsys to NumericArray objects in Wolfram
Interface runtime includes datatype conversion time in Library Link interface
Backend runtime includes the simulation runtime in C++";


Begin["`Private`"];


(*Compile C++ part of package*)
Quiet@LibraryUnload["StoChemSimInterface"];
packageDir = $InputFileName//DirectoryName;
If[packageDir === "", packageDir = NotebookDirectory[]];
pathInterface = FileNameJoin[
	{packageDir, "StoChemSimSequential", "interface.cpp"}, 
	OperatingSystem -> $OperatingSystem];
CreateLibrary[{pathInterface}, "StoChemSimInterface", "Language" -> "C++"];


(*Error given when timeEnd value is invalid*)
timeenderrMsg = "Error: timeEnd (`1`) must be a real number greater than zero"
SimulateDirectSSA::timeenderr = timeenderrMsg
SimulateBoundedTauLeaping::timeenderr = timeenderrMsg

(*Error given when iterEnd value is invalid*)
iterenderrMsg = "Error: iterEnd (`1`) must be an integer greater than zero"
SimulateDirectSSA::iterenderr = iterenderrMsg
SimulateBoundedTauLeaping::iterenderr = iterenderrMsg

(*Error given when epsilon value is invalid*)
epsilonerrMsg = "Error: epsilon (`1`) must be a real number between 0 and 1"
SimulateBoundedTauLeaping::epsilonerr = epsilonerrMsg

(*Error given when no simulations have been run*)
nosimerrMsg = "Error: No simulations have been run"
GetRuntimeInfo::nosimerr = nosimerrMsg
PlotLastSimulation::nosimerr = nosimerrMsg

(*Error given when plotting is attempted and no states were saved*)
finalonlyerrMsg = "Error: Cannot plot when finalOnly = True"
PlotLastSimulation::finalonlyerr = finalonlyerrMsg


(*Creates initial state vector using conc objects*)
GetInitCounts[concs_, spcs_] := (Plus @@ Cases[concs, conc[#, count_] :> count])& /@ spcs

(*Creates reactant count vectors by obtaining reactant coefficients from each reaction*)
GetReactCounts[rxnls_, spcs_] := Module[
	{numRxnls = Count[rxnls, rxnl[__]],
	reactants = Cases[rxnls, rxnl[r_, _, _] :> r],
	reactantsStoich, spcsIndexMapping},
	reactantsStoich = Table[0, {numRxnls},{Length[spcs]}];
	spcsIndexMapping = AssociationThread[spcs, Range[Length[spcs]]];
	MapIndexed[(reactantsStoich[[First[#2], spcsIndexMapping[#1]]]++)&, reactants, {2}];
	reactantsStoich
]

(*Creates product count vectors by obtaining product coefficients from each reaction*)
GetProdCounts[rxnls_, spcs_] := Module[
	{numRxnls = Count[rxnls, rxnl[__]],
	products = Cases[rxnls, rxnl[_, p_, _] :> p],
	productsStoich, spcsIndexMapping},
	productsStoich = Table[0, {numRxnls},{Length[spcs]}];
	spcsIndexMapping = AssociationThread[spcs, Range[Length[spcs]]];
	MapIndexed[(productsStoich[[First[#2], spcsIndexMapping[#1]]]++)&, products, {2}];
	productsStoich
]

(*Creates rate vector by obtaining reaction rate from each reaction*)
GetRates[rxnls_] := Cases[rxnls, rxnl[_, _, k_] :> k]


(*Loads C++ library and loads interface functions with specified argument types*)
library = LibraryLoad["StoChemSimInterface"]

(*Entire backend implementation for Direct SSA, has no return type*)
DirectSSABackend = LibraryFunctionLoad[library, "directSSAInterface",
	{LibraryDataType[NumericArray],
	LibraryDataType[NumericArray],
	LibraryDataType[NumericArray],
	LibraryDataType[NumericArray],
	Real,
	Integer,
	True|False,
	True|False,
	True|False,
	True|False},
	"Void"
]

(*Entire backend implementation for BTL, has no return type*)
BTLBackend = LibraryFunctionLoad[library, "BTLInterface",
	{LibraryDataType[NumericArray],
	LibraryDataType[NumericArray],
	LibraryDataType[NumericArray],
	LibraryDataType[NumericArray],
	Real,
	Integer,
	True|False,
	True|False,
	True|False,
	Real},
	"Void"
]

(*Once a simulation has been run, these functions obtain the simulation results*)
GetStates = LibraryFunctionLoad[library, "getStates", {}, LibraryDataType[NumericArray]]
GetTimes = LibraryFunctionLoad[library, "getTimes", {}, LibraryDataType[NumericArray]]
GetRuntimes = LibraryFunctionLoad[library, "getRuntimes", {}, LibraryDataType[NumericArray]]


(*Specifies options that can be passed into SimulateDirectSSA with default values*)
Options[SimulateDirectSSA] = {
	"timeEnd" -> Infinity,
	"iterEnd" -> Infinity,
	"useIter" -> False,
	"statesOnly" -> False,
	"finalOnly" -> False,
	"outputTS" -> True
}

(*Main function that runs Direct SSA simulation on given rxnsys*)
SimulateDirectSSA[rxnsys_, OptionsPattern[]] := Module[
	(*Define local variables*)
	{initCounts, reactCounts, prodCounts, rates,
	initCountsNA, reactCountsNA, prodCountsNA, ratesNA,
	infTime, infIter, frontendTime, runtimes},
	
	frontendTime = Timing[Module[{},
	(*Initialize global variables that are stored for PlotLastSimulation*)
	concs = Cases[ExpandConcs[rxnsys], conc[_, _]];
	rxnls = Cases[RxnsToRxnls[rxnsys], rxnl[_List, _List, _]];
	spcs = SpeciesInRxnsys[rxnsys];
	timeEnd = OptionValue["timeEnd"];
	iterEnd = OptionValue["iterEnd"];
	useIter = OptionValue["useIter"];
	statesOnly = OptionValue["statesOnly"];
	finalOnly = OptionValue["finalOnly"];
	outputTS = OptionValue["outputTS"];
	BTL = False;
	
	(*Exception handling block*)
	(*Set infTime flag if timeEnd is infinity or invalid*)
	If[timeEnd === Infinity || timeEnd <= 0 || !NumericQ[timeEnd],
		timeEndR = 1000000.0; infTime = True,
		timeEndR = N[timeEnd]; infTime = False
	];
	(*If timeEnd is invalid, give error*)
	If[((timeEnd =!= Infinity && !NumericQ[timeEnd]) || timeEnd <= 0.0) && useIter === False, 
		Message[SimulateDirectSSA::timeenderr, timeEnd]
	];
	(*Set infIter flag if iterEnd is infinity or invalid*)
	If[iterEnd === Infinity || iterEnd <= 0 || !IntegerQ[iterEnd],
		iterEndI = 1000000; infIter = True,
		iterEndI = Round[iterEnd]; infIter = False];
	(*If iterEnd is invalid, give error*)
	If[((iterEnd =!= Infinity && !IntegerQ[iterEnd]) || iterEnd <= 0) && useIter === True,
		Message[SimulateDirectSSA::iterenderr, iterEnd]
	];
	(*Set inf flag based on userIter, infTime, and infIter, which determines if simulation runs to completion or not*)
	If[(infTime === True && useIter === False) || (infIter === True && useIter === True),
		inf = True,
		inf = False
	];
	
	(*Determine simulation parameters from rxnsys and convert to numeric arrays of correct datatypes*)
	initCounts = GetInitCounts[concs, spcs];
	reactCounts = GetReactCounts[rxnls, spcs];
	prodCounts = GetProdCounts[rxnls, spcs];
	rates = GetRates[rxnls];
	initCountsNA = NumericArray[initCounts, "Integer32"];
	reactCountsNA = NumericArray[reactCounts, "Integer32"];
	prodCountsNA = NumericArray[prodCounts, "Integer32"];
	ratesNA = NumericArray[rates, "Real64"];
	];][[1]];
	
	(*Run simulation via C++ library*)
	DirectSSABackend[initCountsNA, reactCountsNA, prodCountsNA, ratesNA, timeEndR, iterEndI, inf, useIter, statesOnly, finalOnly];
	(*Output format depends on outputTS and statesOnly flags*)
	If[outputTS && !finalOnly,
		If[statesOnly,
			simulationResult = TimeSeries[Normal[GetStates[]], {0, Length[Normal[GetStates[]]]-1}],
			simulationResult = TimeSeries[Normal[GetStates[]], {Normal[GetTimes[]]}]
		],
		If[statesOnly,
			simulationResult = Normal[GetStates[]],
			simulationResult = {Normal[GetTimes[]], Normal[GetStates[]]}
		]
	];
	runtimes = Normal[GetRuntimes[]];
	runtimeInfo =  <|
		"frontend conversion" -> frontendTime,
		"interface conversion" -> runtimes[[1]],
		"backend preprocessing" -> runtimes[[2]],
		"backend algorithm" -> runtimes[[3]]
	|>;
	simulationResult
]


(*Specifies options that can be passed into SimulateBoundedTauLeaping with default values*)
Options[SimulateBoundedTauLeaping] = {
	"timeEnd" -> Infinity,
	"iterEnd" -> Infinity,
	"useIter" -> False,
	"finalOnly" -> False,
	"outputTS" -> True,
	"rho" -> Null,
	"epsilon" -> Null
}

(*Main function that runs Bounded Tau Leaping simulation on given rxnsys*)
SimulateBoundedTauLeaping[rxnsys_, OptionsPattern[]] := Module[
	(*Define local variables*)
	{initCounts, reactCounts, prodCounts, rates,
	initCountsNA, reactCountsNA, prodCountsNA, ratesNA,
	infTime, infIter, frontendTime, runtimes},
	
	frontendTime = Timing[Module[{},
	(*Initialize global variables that are stored for PlotLastSimulation*)
	concs = Cases[ExpandConcs[rxnsys], conc[_, _]];
	rxnls = Cases[RxnsToRxnls[rxnsys], rxnl[_List, _List, _]];
	spcs = SpeciesInRxnsys[rxnsys];
	timeEnd = OptionValue["timeEnd"];
	iterEnd = OptionValue["iterEnd"];
	useIter = OptionValue["useIter"];
	statesOnly = False;
	finalOnly = OptionValue["finalOnly"];
	outputTS = OptionValue["outputTS"];
	rho = OptionValue["rho"];
	epsilon = OptionValue["epsilon"];
	BTL = True;
	(*If no rho or epsilon set, use well-defined epsilon (equivalent to p=0.25)*)
	If[rho === Null && epsilon === Null,
		rho = 0.25
	];
	If[epsilon === Null,
		epsilon = 3/(4*Length[rxnls])*(1 - Sqrt[(1 + rho/9)/(1 + rho)])
	];
	
	(*Exception handling block*)
	(*Set infTime flag if timeEnd is infinity or invalid*)
	If[timeEnd === Infinity || timeEnd <= 0 || !NumericQ[timeEnd],
		timeEndR = 1000000.0; infTime = True,
		timeEndR = N[timeEnd]; infTime = False
	];
	(*If timeEnd is invalid, give error*)
	If[((timeEnd =!= Infinity && !NumericQ[timeEnd]) || timeEnd <= 0.0) && useIter === False, 
		Message[SimulateBoundedTauLeaping::timeenderr, timeEnd]
	];
	(*Set infIter flag if iterEnd is infinity or invalid*)
	If[iterEnd === Infinity || iterEnd <= 0 || !IntegerQ[iterEnd],
		iterEndI = 1000000; infIter = True,
		iterEndI = Round[iterEnd]; infIter = False];
	(*If iterEnd is invalid, give error*)
	If[((iterEnd =!= Infinity && !IntegerQ[iterEnd]) || iterEnd <= 0) && useIter === True,
		Message[SimulateBoundedTauLeaping::iterenderr, iterEnd]
	];
	(*Set inf flag based on userIter, infTime, and infIter, which determines if simulation runs to completion or not*)
	If[(infTime === True && useIter === False) || (infIter === True && useIter === True),
		inf = True,
		inf = False
	];
	(*If epsilon is invalid, give error*)
	If[!NumericQ[epsilon] || epsilon < 0 || epsilon > 1,
		Message[SimulateBoundedTauLeaping::epsilonerr, epsilon];
		rho = 0.25;
		epsilon = 3/(4*Length[rxnls])*(1 - Sqrt[(1 + rho/9)/(1 + rho)])
	];
	
	(*Determine simulation parameters from rxnsys and convert to numeric arrays of correct datatypes*)
	initCounts = GetInitCounts[concs, spcs];
	reactCounts = GetReactCounts[rxnls, spcs];
	prodCounts = GetProdCounts[rxnls, spcs];
	rates = GetRates[rxnls];
	initCountsNA = NumericArray[initCounts, "Integer32"];
	reactCountsNA = NumericArray[reactCounts, "Integer32"];
	prodCountsNA = NumericArray[prodCounts, "Integer32"];
	ratesNA = NumericArray[rates, "Real64"];
	];][[1]];
	
	(*Run simulation via C++ library*)
	BTLBackend[initCountsNA, reactCountsNA, prodCountsNA, ratesNA, timeEndR, iterEndI, inf, useIter, finalOnly, epsilon];
	(*Output format depends on outputTS flag*)
	If[outputTS && !finalOnly,
		simulationResult = TimeSeries[Normal[GetStates[]], {Normal[GetTimes[]]}],
		simulationResult = {Normal[GetTimes[]], Normal[GetStates[]]}
	];
	runtimes = Normal[GetRuntimes[]];
	runtimeInfo =  <|
		"frontend conversion" -> frontendTime,
		"interface conversion" -> runtimes[[1]],
		"backend preprocessing" -> runtimes[[2]],
		"backend algorithm" -> runtimes[[3]]
	|>;
	simulationResult
]


(*Retrieves stored runtime results*)
GetRuntimeInfo[] := Module[{},
	If[Head[runtimeInfo] === Symbol,
		Message[GetRuntimeInfo::nosimerr],
		runtimeInfo
	]
]

(*Plots the last simulation ran using the same options as ListLinePlot*)
PlotLastSimulation[opts:OptionsPattern[ListLinePlot]] := Module[
	(*Define local variables*)
	{ts, title, xLabel},
	
	(*If no simulation ran, give error*)
	If[Head[simulationResult] === Symbol,
		Message[PlotLastSimulation::nosimerr],
		(*If no states saved, give error*)
		If[finalOnly === True,
			Message[PlotLastSimulation::finalonlyerr],
			(*Cases for formatting TimeSeries vs. List object*)
			If[outputTS === True,
				ts = simulationResult,
				(*Cases for formatting with times or no times*)
				If[statesOnly === True,
					ts = TimeSeries[simulationResult, {0, Length[simulationResult]-1}],
					ts = TimeSeries[simulationResult[[2]], {simulationResult[[1]]}]
				]
			];
			(*Set appropriate title*)
			If[BTL,
				title = "BTL Simulation",
				title = "Direct SSA Simulation"
			];
			(*Set appropriate x axis label*)
			If[statesOnly,
				xLabel = "Iterations",
				xLabel = "Time [s]"
			];
			(*Plot result with user given options overriding preset options*)
			ListLinePlot[
				ts,
				opts,
				PlotRange -> {0, All},
				PlotLegends -> spcs,
				AxesLabel -> {xLabel, "Molecular Count"},
				PlotLabel -> title,
				MaxPlotPoints -> 500
			]
		]
	]
]


End[]
EndPackage[]
