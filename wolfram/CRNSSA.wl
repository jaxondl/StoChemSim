(* ::Package:: *)

(* ::Title:: *)
(*CRNSSA Package*)


BeginPackage["CRNSSA`"]


rxn::usage =
"rxn[reactants expression, products expression, rate]
A chemical reaction with a reactants expression,
products expression, and rate";

revrxn::usage =
"revrxn[reactants expression, products expression, forward rate, backward rate]
A reversible chemical reaction with a reactants expression,
products expression, forwards rate, and backward rate";

init::usage =
"init[species, quantity] sets the initial quantity of a species
init[species list, quantity] sets the initial quantity of every species in a species list";

GetSpecies::usage =
"GetSpecies[rxnsys] returns the species that appear in rxnsys
GetSpecies[rxnsys, pattern] returns the species in rxnsys that match the specified pattern";

DirectSSA::usage =
"result = DirectSSA[{rxn1, rxn2, ..., init1, init2, ...}, Options]
Simulates the given reaction system via Gillespie Direct SSA
Backend is optimized in C++ for computational efficiency
Options include:
	timeEnd (Real), default = Infinity, ending time of simulation (if useIter = False)
	iterEnd (Integer), default = Infinity, ending iteration of simulation (if useIter = True)
	useIter (Boolean), default = False, setting to True uses iterEnd instead of timeEnd
	statesOnly (Boolean), default = False, setting to True avoids simulating reaction times
	finalOnly (Boolean), default = False, only records final state to conserve memory
	outputTS (Boolean), default = True, setting to True outputs result as TimeSeries, setting to False outputs result as List";
	
BoundedTauLeaping::usage =
"result = BoundedTauLeaping[{rxn1, rxn2, ..., init1, init2, ...}, Options]
Simulates the given reaction system via Soloveichik Bounded Tau Leaping
Backend is optimized in C++ for computational efficiency
Options include:
	timeEnd (Real), default = Infinity, ending time of simulation (if useIter = False)
	iterEnd (Integer), default = Infinity, ending iteration of simulation (if useIter = True)
	useIter (Boolean), default = False, setting to True uses iterEnd instead of timeEnd
	finalOnly (Boolean), default = False, only records final state to conserve memory
	outputTS (Boolean), default = True, setting to True outputs result as TimeSeries, setting to False outputs result as List
	epsilon (Real), default = 0.0309/numRxns, threshold between 0 and 1 using in calculating firing bounds for each reaction";

PlotLastSimulation::usage =
"PlotLastSimulation[Options]
Plots the last simulation ran
Uses same Options as ListLinePlot"


Begin["`Private`"]


(*Handle rxn cases with no products or reactants*)
rxn[Except[1, _Integer], p_, k_] := rxn[1, p, k]
rxn[r_, Except[1, _Integer], k_] := rxn[r, 1, k]
(*Create two rxn objects from revrxn*)
revrxn[r_, p_, kf_, kb_] := Sequence[rxn[r, p, kf], rxn[p, r, kb]]
(*Create multiple init objects when list of species given*)
init[s_List, n_] := Sequence@@(init[#, n]& /@ s)


(*Function to obtain list of objects that don't match expected pattern in rxnsys*)
GetUnkObjs[rxnsys_] := Cases[rxnsys, Except[rxn[_, _, _] | init[_, _]]]
(*Warning given when unknown objects present in rxnsys*)
rxnsyswarnMsg = "Warning: Unknown objects detected in rxnsys. These will be ignored by Wolfram pattern matching: `1`"
GetSpecies::rxnsyswarn = rxnsyswarnMsg
DirectSSA::rxnsyswarn = rxnsyswarnMsg
BoundedTauLeaping::rxnsyswarn = rxnsyswarnMsg
(*Error given when timeEnd value is invalid*)
timeenderrMsg = "Error: timeEnd (`1`) must be a real number greater than zero"
DirectSSA::timeenderr = timeenderrMsg
BoundedTauLeaping::timeenderr = timeenderrMsg
(*Error given when iterEnd value is invalid*)
iterenderrMsg = "Error: iterEnd (`1`) must be an integer greater than zero"
DirectSSA::iterenderr = iterenderrMsg
BoundedTauLeaping::iterenderr = iterenderrMsg
(*Error given when epsilon value is invalid*)
epsilonerrMsg = "Error: epsilon (`1`) must be a real number between 0 and 1"
BoundedTauLeaping::epsilonerr = epsilonerrMsg
(*Error given when plotting is attempted with no simulations run*)
nosimerrMsg = "Error: No simulations have been run"
PlotLastSimulation::nosimerr = nosimerrMsg
(*Error given when plotting is attempted and no states were saved*)
finalonlyerrMsg = "Error: Cannot plot when finalOnly = True"
PlotLastSimulation::finalonlyerr = finalonlyerrMsg


(*Function to get all species types in rxnsys, obtained from both rxn and init objects*)
GetSpecies[rxnsys_] := Module[{unkObjs = GetUnkObjs[rxnsys]},
	If[Length[unkObjs] =!= 0, Message[GetSpecies::rxnsyswarn, unkObjs]];
	Sort[Union[
		Cases[Cases[rxnsys, rxn[r_, p_, _] :> Sequence[r, p]] /. Times | Plus -> Sequence, s_Symbol | s_Symbol[__]],
		Cases[rxnsys, init[x_, _] :> x]
	]]
]
(*Allows usage of GetSpecies with a specific pattern that returned species must match*)
GetSpecies[rxnsys_, pattern_] := Cases[GetSpecies[rxnsys], pattern]


(*Creates initial state vector using init objects*)
GetInitCounts[inits_, spcs_] := (Plus @@ Cases[inits, init[#, count_] :> count])& /@ spcs
(*Creates reactant count vectors by obtaining reactant coefficients from each reaction*)
GetReactCounts[rxns_, spcs_] := Outer[Coefficient[#1, #2]&, Cases[rxns, rxn[r_, _, _] :> r], spcs]
(*Creates product count vectors by obtaining product coefficients from each reaction*)
GetProdCounts[rxns_, spcs_] := Outer[Coefficient[#1, #2]&, Cases[rxns, rxn[_, p_, _] :> p], spcs]
(*Creates rate vector by obtaining reaction rate from each reaction*)
GetRates[rxns_] := Cases[rxns, rxn[_, _, k_] :> k]


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


(*Specifies options that can be passed into DirectSSA with default values*)
Options[DirectSSA] = {
	"timeEnd" -> Infinity,
	"iterEnd" -> Infinity,
	"useIter" -> False,
	"statesOnly" -> False,
	"finalOnly" -> False,
	"outputTS" -> True
}
(*Main funciton that runs Direct SSA simulation on given rxnsys*)
DirectSSA[rxnsys_, OptionsPattern[]] := Module[
	(*Define local variables*)
	{initCounts, reactCounts, prodCounts, rates,
	initCountsNA, reactCountsNA, prodCountsNA, ratesNA,
	infTime, infIter, unkObjs = GetUnkObjs[rxnsys]},
	
	(*Initialize global variables that are stored for PlotLastSimulation*)
	inits = Cases[rxnsys, init[_, _]];
	rxns = Cases[rxnsys, rxn[_, _, _]];
	spcs = Quiet[GetSpecies[rxnsys]];
	timeEnd = OptionValue["timeEnd"];
	iterEnd = OptionValue["iterEnd"];
	useIter = OptionValue["useIter"];
	statesOnly = OptionValue["statesOnly"];
	finalOnly = OptionValue["finalOnly"];
	outputTS = OptionValue["outputTS"];
	
	(*Exception handling block*)
	(*If unknown objects present, give warning*)
	If[Length[unkObjs] =!= 0,
		Message[DirectSSA::rxnsyswarn, unkObjs]
	];
	(*Set infTime flag if timeEnd is infinity or invalid*)
	If[timeEnd === Infinity || timeEnd <= 0 || !NumericQ[timeEnd],
		timeEndR = 1000000.0; infTime = True,
		timeEndR = N[timeEnd]; infTime = False
	];
	(*If timeEnd is invalid, give error*)
	If[((timeEnd =!= Infinity && !NumericQ[timeEnd]) || timeEnd <= 0.0) && useIter === False, 
		Message[DirectSSA::timeenderr, timeEnd]
	];
	(*Set infIter flag if iterEnd is infinity or invalid*)
	If[iterEnd === Infinity || iterEnd <= 0 || !IntegerQ[iterEnd],
		iterEndI = 1000000; infIter = True,
		iterEndI = Round[iterEnd]; infIter = False];
	(*If iterEnd is invalid, give error*)
	If[((iterEnd =!= Infinity && !IntegerQ[iterEnd]) || iterEnd <= 0) && useIter === True,
		Message[DirectSSA::iterenderr, iterEnd]
	];
	(*Set inf flag based on userIter, infTime, and infIter, which determines if simulation runs to completion or not*)
	If[(infTime === True && useIter === False) || (infIter === True && useIter === True),
		inf = True,
		inf = False
	];
	
	(*Determine simulation parameters from rxnsys and convert to numeric arrays of correct datatypes*)
	initCounts = GetInitCounts[inits, spcs];
	reactCounts = GetReactCounts[rxns, spcs];
	prodCounts = GetProdCounts[rxns, spcs];
	rates = GetRates[rxns];
	initCountsNA = NumericArray[initCounts, "Integer32"];
	reactCountsNA = NumericArray[reactCounts, "Integer64"];
	prodCountsNA = NumericArray[prodCounts, "Integer64"];
	ratesNA = NumericArray[rates, "Real64"];
	
	(*Run simulation via C++ library*)
	DirectSSABackend[initCountsNA, reactCountsNA, prodCountsNA, ratesNA, timeEndR, iterEndI, inf, useIter, statesOnly, finalOnly];
	(*Output format depends on outputTS and statesOnly flags*)
	If[outputTS,
		If[statesOnly,
			simulationResult = TimeSeries[Normal[GetStates[]], {0, Length[Normal[GetStates[]]]-1}],
			simulationResult = TimeSeries[Normal[GetStates[]], {Normal[GetTimes[]]}]
		],
		If[statesOnly,
			simulationResult = Normal[GetStates[]],
			simulationResult = {Normal[GetStates[]], Normal[GetTimes[]]}
		]
	]
]


(*Specifies options that can be passed into BoundedTauLeaping with default values*)
Options[BoundedTauLeaping] = {
	"timeEnd" -> Infinity,
	"iterEnd" -> Infinity,
	"useIter" -> False,
	"finalOnly" -> False,
	"outputTS" -> True,
	"epsilon" -> Null
}
(*Main funciton that runs Direct SSA simulation on given rxnsys*)
BoundedTauLeaping[rxnsys_, OptionsPattern[]] := Module[
	(*Define local variables*)
	{initCounts, reactCounts, prodCounts, rates,
	initCountsNA, reactCountsNA, prodCountsNA, ratesNA,
	infTime, infIter, unkObjs = GetUnkObjs[rxnsys]},
	
	(*Initialize global variables that are stored for PlotLastSimulation*)
	inits = Cases[rxnsys, init[_, _]];
	rxns = Cases[rxnsys, rxn[_, _, _]];
	spcs = Quiet[GetSpecies[rxnsys]];
	timeEnd = OptionValue["timeEnd"];
	iterEnd = OptionValue["iterEnd"];
	useIter = OptionValue["useIter"];
	statesOnly = False;
	finalOnly = OptionValue["finalOnly"];
	outputTS = OptionValue["outputTS"];
	epsilon = OptionValue["epsilon"];
	(*If no epsilon set, use well-defined epsilon (equivalent to p=0.1)*)
	If[epsilon === Null,
		epsilon = 0.0309/Length[rxns]
	];
	
	(*Exception handling block*)
	(*If unknown objects present, give warning*)
	If[Length[unkObjs] =!= 0,
		Message[BoundedTauLeaping::rxnsyswarn, unkObjs]
	];
	(*Set infTime flag if timeEnd is infinity or invalid*)
	If[timeEnd === Infinity || timeEnd <= 0 || !NumericQ[timeEnd],
		timeEndR = 1000000.0; infTime = True,
		timeEndR = N[timeEnd]; infTime = False
	];
	(*If timeEnd is invalid, give error*)
	If[((timeEnd =!= Infinity && !NumericQ[timeEnd]) || timeEnd <= 0.0) && useIter === False, 
		Message[BoundedTauLeaping::timeenderr, timeEnd]
	];
	(*Set infIter flag if iterEnd is infinity or invalid*)
	If[iterEnd === Infinity || iterEnd <= 0 || !IntegerQ[iterEnd],
		iterEndI = 1000000; infIter = True,
		iterEndI = Round[iterEnd]; infIter = False];
	(*If iterEnd is invalid, give error*)
	If[((iterEnd =!= Infinity && !IntegerQ[iterEnd]) || iterEnd <= 0) && useIter === True,
		Message[BoundedTauLeaping::iterenderr, iterEnd]
	];
	(*Set inf flag based on userIter, infTime, and infIter, which determines if simulation runs to completion or not*)
	If[(infTime === True && useIter === False) || (infIter === True && useIter === True),
		inf = True,
		inf = False
	];
	(*If epsilon is invalid, give error*)
	If[!NumericQ[epsilon] || epsilon < 0 || epsilon > 1,
		Message[BoundedTauLeaping::epsilonerr, epsilon];
		epsilon = 0.0309/Length[rxns]
	];
	
	(*Determine simulation parameters from rxnsys and convert to numeric arrays of correct datatypes*)
	initCounts = GetInitCounts[inits, spcs];
	reactCounts = GetReactCounts[rxns, spcs];
	prodCounts = GetProdCounts[rxns, spcs];
	rates = GetRates[rxns];
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
