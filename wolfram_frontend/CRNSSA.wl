(* ::Package:: *)

(* ::Title:: *)
(*CRN SSA Package*)


BeginPackage["CRNSSA`"]


rxn::usage = "Chemical reaction, expressed as rxn[reactants, products, rate]";
revrxn::usage = "Reversible reaction, expressed as revrxn[reactants, products, forward rate, backward rate]";
init::usage = "Initial molecular counts, expressed as init[species or species list, quantity]";
GetSpecies::usage =
"GetSpecies[rxnsys] returns the species in rxnsys
GetSpecies[rxnsys, pattern] returns the species in rxnsys that match the specified pattern"
DirectSSA::usage =
"DirectSSA[{rxn1, rxn2, ..., init1, init2, ...}, end time]
Simulates the given reaction system via Gillespie direct SSA
Backend is optimized in C++ for computational efficiency
Note: implementation not complete";


Begin["`Private`"]


rxn[Except[1, _Integer], p_, k_] := rxn[1, p, k]
rxn[r_, Except[1, _Integer], k_] := rxn[r, 1, k]
revrxn[r_, p_, kf_, kb_] := Sequence[rxn[r, p, kf], rxn[p, r, kb]]
init[s_List, n_] := Sequence@@(init[#, n]& /@ s)


GetSpecies[rxnsys_] := Sort[Union[
	Cases[Cases[rxnsys, rxn[r_, p_, _] :> Sequence[r, p]] /. Times | Plus -> Sequence, s_Symbol | s_Symbol[__]],
	Cases[rxnsys, init[x_, _] :> x]
]]
GetSpecies[rxnsys_, pattern_] := Cases[GetSpecies[rxnsys], pattern]


GetInitCounts[inits_, spcs_] := (Plus @@ Cases[inits, init[#, count_] :> count])& /@ spcs
GetReactCounts[rxns_, spcs_] := Outer[Coefficient[#1, #2]&, Cases[rxns, rxn[r_, _, _] :> r], spcs]
GetProdCounts[rxns_, spcs_] := Outer[Coefficient[#1, #2]&, Cases[rxns, rxn[_, p_, _] :> p], spcs]
GetRates[rxns_] := Cases[rxns, rxn[_, _, k_] :> k]


DirectSSA[rxnsys_, tEnd_] := Module[
	{inits = Cases[rxnsys, init[_, _]],
	rxns = Cases[rxnsys, rxn[_, _, _]],
	spcs = GetSpecies[rxnsys],
	initCounts, reactCounts, prodCounts, rates,
	initCountsNA, reactCountsNA, prodCountsNA, ratesNA},
	
	initCounts = GetInitCounts[inits, spcs];
	reactCounts = GetReactCounts[rxns, spcs];
	prodCounts = GetProdCounts[rxns, spcs];
	rates = GetRates[rxns];
	
	initCountsNA = NumericArray[initCounts, "Integer64"];
	reactCountsNA = NumericArray[reactCounts, "Integer64"];
	prodCountsNA = NumericArray[prodCounts, "Integer64"];
	ratesNA = NumericArray[rates, "Integer64"];
	
	{initCountsNA, reactCountsNA, prodCountsNA, ratesNA, tEnd}
]


End[]
EndPackage[]
