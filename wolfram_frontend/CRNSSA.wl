(* ::Package:: *)

(* ::Title:: *)
(*CRN SSA Package*)


BeginPackage["CRNSSA`"]


rxn::usage = "Chemical reaction, expressed as rxn[reactants, products, rate]";
revrxn::usage = "Reversible reaction, expressed as revrxn[reactants, products, forward rate, backward rate]";
init::usage = "Initial molecular counts, expressed as init[species or species list, quantity]";
DirectSSA::usage =
"Simulates a reaction system via Gillespie direct SSA, expressed as
DirectSSA[{rxn1, rxn2, ..., init1, init2, ...}, end time]
Note: implementation not complete";


Begin["`Private`"]


(* is rxn definition for no reactants necessary, and would we need a similar definition for products? *)
rxn[Except[1, _Integer], p_, k_] := rxn[1, p, k]
revrxn[r_, p_, kf_, kb_] := Sequence[rxn[r, p, kf], rxn[p, r, kb]]
init[s_List, n_] := Sequence@@(init[#, n]& /@ s)


(* is s_Symbol[__] necessary? why would a species be a head with arguments? *)
GetSpecies[rxnsys_] := Union[
 Cases[Cases[rxnsys, rxn[r_, p_, _] :> Sequence[r, p]] /. Times | Plus -> Sequence, s_Symbol | s_Symbol[__]],
 Cases[rxnsys, init[x_, _] :> x]]


DirectSSA[rxnsys_, tend_] := {rxnsys, tend, GetSpecies[rxnsys]}


End[]
EndPackage[]
