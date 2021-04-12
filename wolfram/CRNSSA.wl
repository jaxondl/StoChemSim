(* ::Package:: *)

(* ::Title:: *)
(*CRNSSA Package*)


Needs["CRNSimulator`"]
BeginPackage["CRNSSA`", {"CRNSimulator`"}]


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
Note: implementation is not complete
Options include:
	timeEnd (Real), default = Infinity, ending time of simulation (if useIter = False)
	iterEnd (Integer), default = Infinity, ending iteration of simulation (if useIter = True)
	useIter (Boolean), default = False, setting to True uses iterEnd instead of timeEnd
	finalOnly (Boolean), default = False, only records final state to conserve memory
	outputTS (Boolean), default = True, setting to True outputs result as TimeSeries, setting to False outputs result as List
	epsilon (Real), default = 0.0309/numRxns, threshold between 0 and 1 used in calculating firing bounds for each reaction";

PlotLastSimulation::usage =
"PlotLastSimulation[Options]
Plots the last simulation ran
Uses same Options as ListLinePlot"


Begin["`Private`"]


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

(*Error given when plotting is attempted with no simulations run*)
nosimerrMsg = "Error: No simulations have been run"
PlotLastSimulation::nosimerr = nosimerrMsg

(*Error given when plotting is attempted and no states were saved*)
finalonlyerrMsg = "Error: Cannot plot when finalOnly = True"
PlotLastSimulation::finalonlyerr = finalonlyerrMsg


(*Creates initial state vector using conc objects*)
GetInitCounts[concs_, spcs_] := (Plus @@ Cases[concs, conc[#, count_] :> count])& /@ spcs

(*Creates reactant count vectors by obtaining reactant coefficients from each reaction*)
GetReactCounts[rxnls_, spcs_] := Outer[Coefficient[#1, #2]&, Cases[rxnls, rxnl[r_List, _List, _] :> Plus @@ r], spcs]

(*Creates product count vectors by obtaining product coefficients from each reaction*)
GetProdCounts[rxnls_, spcs_] := Outer[Coefficient[#1, #2]&, Cases[rxnls, rxnl[_List, p_List, _] :> Plus @@ p], spcs]

(*Creates rate vector by obtaining reaction rate from each reaction*)
GetRates[rxnls_] := Cases[rxnls, rxnl[_, _, k_] :> k]


(*Loads C++ library and loads interface functions with specified argument types*)
library = LibraryLoad["directSSAinterface"];

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
];

(*Once DirectSSABackend has been run, these functions obtain the simulation results*)
GetStates = LibraryFunctionLoad[library, "getStates", {}, LibraryDataType[NumericArray]];
GetTimes = LibraryFunctionLoad[library, "getTimes", {}, LibraryDataType[NumericArray]];
GetDebugs = LibraryFunctionLoad[library, "getDebugs", {}, LibraryDataType[NumericArray]];


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
	infTime, infIter},
	
	(*Initialize global variables that are stored for PlotLastSimulation*)
	initTime = Timing[Module[{},
	concs = Quiet[Cases[ExpandConcs[rxnsys], conc[_, _]]];
	rxnls = Quiet[Cases[RxnsToRxnls[rxnsys], rxnl[_List, _List, _]]];
	spcs = Quiet[SpeciesInRxnsys[rxnsys]];
	timeEnd = OptionValue["timeEnd"];
	iterEnd = OptionValue["iterEnd"];
	useIter = OptionValue["useIter"];
	statesOnly = OptionValue["statesOnly"];
	finalOnly = OptionValue["finalOnly"];
	outputTS = OptionValue["outputTS"];
	];][[1]];
	
	(*Exception handling block*)
	exceptionTime = Timing[Module[{},
	CheckSyntaxErrors[rxnsys];
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
	];][[1]];
	
	(*Determine simulation parameters from rxnsys and convert to numeric arrays of correct datatypes*)
	naTime = Timing[Module[{},
	initCounts = GetInitCounts[concs, spcs];
	reactCounts = GetReactCounts[rxnls, spcs];
	prodCounts = GetProdCounts[rxnls, spcs];
	rates = GetRates[rxnls];
	initCountsNA = NumericArray[initCounts, "Integer32"];
	reactCountsNA = NumericArray[reactCounts, "Integer64"];
	prodCountsNA = NumericArray[prodCounts, "Integer64"];
	ratesNA = NumericArray[rates, "Real64"];
	];][[1]];
	
	(*Run simulation via C++ library*)
	backendTime = Timing[Module[{},
	DirectSSABackend[initCountsNA, reactCountsNA, prodCountsNA, ratesNA, timeEndR, iterEndI, inf, useIter, statesOnly, finalOnly];
	];][[1]];
	(*Output format depends on outputTS and statesOnly flags*)
	resultTime = Timing[Module[{},
	If[outputTS,
		If[statesOnly,
			simulationResult = TimeSeries[Normal[GetStates[]], {0, Length[Normal[GetStates[]]]-1}],
			simulationResult = TimeSeries[Normal[GetStates[]], {Normal[GetTimes[]]}]
		],
		If[statesOnly,
			simulationResult = Normal[GetStates[]],
			simulationResult = {Normal[GetStates[]], Normal[GetTimes[]]}
		]
	];
	];][[1]];
	{simulationResult, {initTime, exceptionTime, naTime, backendTime, resultTime, Normal[GetDebugs[]]}}
]


(*Specifies options that can be passed into SimulateBoundedTauLeaping with default values*)
Options[SimulateBoundedTauLeaping] = {
	"timeEnd" -> Infinity,
	"iterEnd" -> Infinity,
	"useIter" -> False,
	"finalOnly" -> False,
	"outputTS" -> True,
	"epsilon" -> Null
}

(*Main function that runs Bounded Tau Leaping simulation on given rxnsys*)
SimulateBoundedTauLeaping[rxnsys_, OptionsPattern[]] := Module[
	(*Define local variables*)
	{initCounts, reactCounts, prodCounts, rates,
	initCountsNA, reactCountsNA, prodCountsNA, ratesNA,
	infTime, infIter},
	
	(*Initialize global variables that are stored for PlotLastSimulation*)
	concs = Quiet[Cases[ExpandInits[rxnsys], conc[_, _]]];
	rxnls = Quiet[Cases[RxnsToRxnls[rxnsys], rxnl[_List, _List, _]]];
	spcs = Quiet[SpeciesInRxnsys[rxnsys]];
	timeEnd = OptionValue["timeEnd"];
	iterEnd = OptionValue["iterEnd"];
	useIter = OptionValue["useIter"];
	statesOnly = False;
	finalOnly = OptionValue["finalOnly"];
	outputTS = OptionValue["outputTS"];
	epsilon = OptionValue["epsilon"];
	(*If no epsilon set, use well-defined epsilon (equivalent to p=0.1)*)
	If[epsilon === Null,
		epsilon = 0.0309/Length[rxnls]
	];
	
	(*Exception handling block*)
	CheckSyntaxErrors[rxnsys];
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
		epsilon = 0.0309/Length[rxnls]
	];
	
	(*Determine simulation parameters from rxnsys and convert to numeric arrays of correct datatypes*)
	initCounts = GetInitCounts[concs, spcs];
	reactCounts = GetReactCounts[rxnls, spcs];
	prodCounts = GetProdCounts[rxnls, spcs];
	rates = GetRates[rxnls];
	initCountsNA = NumericArray[initCounts, "Integer32"];
	reactCountsNA = NumericArray[reactCounts, "Integer64"];
	prodCountsNA = NumericArray[prodCounts, "Integer64"];
	ratesNA = NumericArray[rates, "Real64"];
	
(*	(*Run simulation via C++ library*)
	BoundedTauLeapingBackend[initCountsNA, reactCountsNA, prodCountsNA, ratesNA, timeEndR, iterEndI, inf, useIter, finalOnly, epsilon];
	(*Output format depends on outputTS flag*)
	If[outputTS,
		simulationResult = TimeSeries[Normal[GetStates[]], {Normal[GetTimes[]]}],
		simulationResult = {Normal[GetStates[]], Normal[GetTimes[]]}
	]*)
	(*Temporary return while backend unimplemented*)
	epsilon
]


(*Plots the last simulation ran using the same options as ListLinePlot*)
PlotLastSimulation[opts:OptionsPattern[ListLinePlot]] := Module[
	(*Define local variables*)
	{xLabel, ts},
	
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
					ts = TimeSeries[simulationResult[[1]], {simulationResult[[2]]}]
				]
			];
			(*Set appropriate x axis label*)
			If[statesOnly === True,
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
				PlotLabel -> "CRN Simulation",
				MaxPlotPoints -> 500
			]
		]
	]
]


End[]
EndPackage[]
