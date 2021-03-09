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

PlotLastSimulation::usage =
"PlotLastSimulation[Options]
Plots the last simulation ran
Uses same Options as ListLinePlot"


Begin["`Private`"]


rxn[Except[1, _Integer], p_, k_] := rxn[1, p, k]
rxn[r_, Except[1, _Integer], k_] := rxn[r, 1, k]
revrxn[r_, p_, kf_, kb_] := Sequence[rxn[r, p, kf], rxn[p, r, kb]]
init[s_List, n_] := Sequence@@(init[#, n]& /@ s)


GetUnkObjs[rxnsys_] := Cases[rxnsys, Except[rxn[_, _, _] | init[_, _]]]
rxnsyswarnMsg = "Warning: Unknown objects detected in rxnsys. These will be ignored by Wolfram pattern matching: `1`"
GetSpecies::rxnsyswarn = rxnsyswarnMsg
DirectSSA::rxnsyswarn = rxnsyswarnMsg
timeenderrMsg = "Error: timeEnd (`1`) must be a real number greater than zero"
DirectSSA::timeenderr = timeenderrMsg
iterenderrMsg = "Error: iterEnd (`1`) must be an integer greater than zero"
DirectSSA::iterenderr = iterenderrMsg
nosimerrMsg = "Error: No simulations have been run"
PlotLastSimulation::nosimerr = nosimerrMsg
finalonlyerrMsg = "Error: Cannot plot when finalOnly = True"
PlotLastSimulation::finalonlyerr = finalonlyerrMsg


GetSpecies[rxnsys_] := Module[{unkObjs = GetUnkObjs[rxnsys]},
	If[Length[unkObjs] =!= 0, Message[GetSpecies::rxnsyswarn, unkObjs]];
	Sort[Union[
		Cases[Cases[rxnsys, rxn[r_, p_, _] :> Sequence[r, p]] /. Times | Plus -> Sequence, s_Symbol | s_Symbol[__]],
		Cases[rxnsys, init[x_, _] :> x]
	]]
]
GetSpecies[rxnsys_, pattern_] := Cases[GetSpecies[rxnsys], pattern]


GetInitCounts[inits_, spcs_] := (Plus @@ Cases[inits, init[#, count_] :> count])& /@ spcs
GetReactCounts[rxns_, spcs_] := Outer[Coefficient[#1, #2]&, Cases[rxns, rxn[r_, _, _] :> r], spcs]
GetProdCounts[rxns_, spcs_] := Outer[Coefficient[#1, #2]&, Cases[rxns, rxn[_, p_, _] :> p], spcs]
GetRates[rxns_] := Cases[rxns, rxn[_, _, k_] :> k]


library = LibraryLoad["directSSAinterface"];
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
GetStates = LibraryFunctionLoad[library, "getStates", {}, LibraryDataType[NumericArray]];
GetTimes = LibraryFunctionLoad[library, "getTimes", {}, LibraryDataType[NumericArray]];


Options[DirectSSA] = {
	"timeEnd" -> Infinity,
	"iterEnd" -> Infinity,
	"useIter" -> False,
	"statesOnly" -> False,
	"finalOnly" -> False,
	"outputTS" -> True
}
DirectSSA[rxnsys_, OptionsPattern[]] := Module[
	{initCounts, reactCounts, prodCounts, rates,
	initCountsNA, reactCountsNA, prodCountsNA, ratesNA,
	infTime, infIter, unkObjs = GetUnkObjs[rxnsys]},
	
	inits = Cases[rxnsys, init[_, _]];
	rxns = Cases[rxnsys, rxn[_, _, _]];
	spcs = Quiet[GetSpecies[rxnsys]];
	timeEnd = OptionValue["timeEnd"];
	iterEnd = OptionValue["iterEnd"];
	useIter = OptionValue["useIter"];
	statesOnly = OptionValue["statesOnly"];
	finalOnly = OptionValue["finalOnly"];
	outputTS = OptionValue["outputTS"];
	
	If[Length[unkObjs] =!= 0,
		Message[DirectSSA::rxnsyswarn, unkObjs]
	];
	If[timeEnd === Infinity || timeEnd <= 0 || !NumericQ[timeEnd],
		timeEndR = 1000000.0; infTime = True,
		timeEndR = N[timeEnd]; infTime = False
	];
	If[((timeEnd =!= Infinity && !NumericQ[timeEnd]) || timeEnd <= 0.0) && useIter === False, 
		Message[DirectSSA::timeenderr, timeEnd]
	];
	If[iterEnd === Infinity || iterEnd <= 0 || !IntegerQ[iterEnd],
		iterEndI = 1000000; infIter = True,
		iterEndI = Round[iterEnd]; infIter = False];
	If[((iterEnd =!= Infinity && !IntegerQ[iterEnd]) || iterEnd <= 0) && useIter === True,
		Message[DirectSSA::iterenderr, iterEnd]
	];
	If[(infTime === True && useIter === False) || (infIter === True && useIter === True),
		inf = True,
		inf = False
	];
	
	initCounts = GetInitCounts[inits, spcs];
	reactCounts = GetReactCounts[rxns, spcs];
	prodCounts = GetProdCounts[rxns, spcs];
	rates = GetRates[rxns];
	initCountsNA = NumericArray[initCounts, "Integer32"];
	reactCountsNA = NumericArray[reactCounts, "Integer64"];
	prodCountsNA = NumericArray[prodCounts, "Integer64"];
	ratesNA = NumericArray[rates, "Real64"];
	
	DirectSSABackend[initCountsNA, reactCountsNA, prodCountsNA, ratesNA, timeEndR, iterEndI, inf, useIter, statesOnly, finalOnly];
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


PlotLastSimulation[opts:OptionsPattern[ListLinePlot]] := Module[
	{xLabel, ts},
	
	If[Head[simulationResult] === Symbol,
		Message[PlotLastSimulation::nosimerr],
		If[finalOnly === True,
			Message[PlotLastSimulation::finalonlyerr],
			If[outputTS === True,
				ts = simulationResult,
				If[statesOnly === True,
					ts = TimeSeries[simulationResult, {0, Length[simulationResult]-1}],
					ts = TimeSeries[simulationResult[[1]], {simulationResult[[2]]}]
				]
			];
			If[statesOnly === True,
				xLabel = "Iterations",
				xLabel = "Time [s]"
			];
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
