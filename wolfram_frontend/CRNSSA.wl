(* ::Package:: *)

(* ::Title:: *)
(*CRN SSA Package*)


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
"{states, times} = DirectSSA[{rxn1, rxn2, ..., init1, init2, ...}, end time]
Simulates the given reaction system via Gillespie Direct SSA
Backend is optimized in C++ for computational efficiency
Note: implementation is not complete";


Begin["`Private`"]


rxn[Except[1, _Integer], p_, k_] := rxn[1, p, k]
rxn[r_, Except[1, _Integer], k_] := rxn[r, 1, k]
revrxn[r_, p_, kf_, kb_] := Sequence[rxn[r, p, kf], rxn[p, r, kb]]
init[s_List, n_] := Sequence@@(init[#, n]& /@ s)


GetUnkObjs[rxnsys_] := Cases[rxnsys, Except[rxn[_, _, _] | init[_, _]]]
rxnsyswarnMsg = "Warning: Unknown objects detected in rxnsys. These will be ignored by Wolfram pattern matching: `1`"
GetSpecies::rxnsyswarn = rxnsyswarnMsg
DirectSSA::rxnsyswarn = rxnsyswarnMsg
tenderrMsg = "Error: tEnd (`1`) must be a number"
DirectSSA::tenderr = tenderrMsg


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


library = LibraryLoad["interface"];
DirectBackend = LibraryFunctionLoad[library, "CRN_SSA",
	{LibraryDataType[NumericArray],
	LibraryDataType[NumericArray],
	LibraryDataType[NumericArray],
	LibraryDataType[NumericArray],
	Real},
	"Void"];
GetStates = LibraryFunctionLoad[library, "getStates", {}, LibraryDataType[NumericArray]];
GetTimes = LibraryFunctionLoad[library, "getTimes", {}, LibraryDataType[NumericArray]];
GetDebug = LibraryFunctionLoad[library, "getDebug", {}, LibraryDataType[NumericArray]];


DirectSSA[rxnsys_] := DirectSSA[rxnsys, Infinity]

DirectSSA[rxnsys_, tEnd_] := Module[
	{inits = Cases[rxnsys, init[_, _]],
	rxns = Cases[rxnsys, rxn[_, _, _]],
	spcs = Quiet[GetSpecies[rxnsys]],
	initCounts, reactCounts, prodCounts, rates,
	initCountsNA, reactCountsNA, prodCountsNA, ratesNA, tEndR,
	unkObjs = GetUnkObjs[rxnsys]},
	
	If[Length[unkObjs] =!= 0, Message[DirectSSA::rxnsyswarn, unkObjs]];
	If[tEnd === Infinity, tEndR = 1000000.0, tEndR = N[tEnd]];
	If[!NumericQ[tEndR], Message[DirectSSA::tenderr, tEndR]];
	
	initCounts = GetInitCounts[inits, spcs];
	reactCounts = GetReactCounts[rxns, spcs];
	prodCounts = GetProdCounts[rxns, spcs];
	rates = GetRates[rxns];
	
	initCountsNA = NumericArray[initCounts, "Integer32"];
	reactCountsNA = NumericArray[reactCounts, "Integer64"];
	prodCountsNA = NumericArray[prodCounts, "Integer64"];
	ratesNA = NumericArray[rates, "Real64"];
	
	DirectBackend[initCountsNA, reactCountsNA, prodCountsNA, ratesNA, tEndR];
	{GetStates[], GetTimes[]}
]


End[]
EndPackage[]
