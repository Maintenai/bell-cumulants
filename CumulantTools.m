
BeginPackage["CumulantTools`"]

PartialBellY::usage = 
  "PartialBellY[n, k, a] computes the partial Bell polynomial B_{n,k} for sequence a.";
BellY::usage = 
  "BellY[n, a] computes the complete Bell polynomial of order n for sequence a.";
Cumulant::usage = 
  "Cumulant[n, moments] returns the n-th cumulant in terms of moments.";
AllCumulants::usage = 
  "AllCumulants[n, moments] returns a list of cumulants up to order n.";
GeneratingFunctionMoments::usage = 
  "GeneratingFunctionMoments[Z, β, n] returns the first n moments by differentiating log[Z].";

Begin["`Private`"]

ClearAll[PartialBellY, BellY, Cumulant, AllCumulants, GeneratingFunctionMoments];

(* Partial Bell Polynomial (Exponential Type) *)
PartialBellY[0, 0, _] := 1;
PartialBellY[_, 0, _] := 0;
PartialBellY[0, _, _] := 0;
PartialBellY[1, 1, a_] := a[[1]];
PartialBellY[n_, k_, a_] := Sum[
  Binomial[n - 1, j - 1] a[[j]] PartialBellY[n - j, k - 1, a],
  {j, 1, n - k + 1}
];

(* Complete Bell Polynomial *)
BellY[n_, a_] := Sum[
  PartialBellY[n, k, a],
  {k, 1, n}
];

(* Cumulant from moments using Bell polynomial expansion *)
Cumulant[n_, moments_List] := Module[
  {μ = Prepend[moments, 0]},
  (-1)^(n - 1) (n - 1)! Sum[
    BellY[k, Table[μ[[i]]/i!, {i, 1, k}]] StirlingS2[n, k],
    {k, 1, n}
  ]
];

(* List of all cumulants up to order n *)
AllCumulants[n_, moments_List] := Table[
  Cumulant[k, moments],
  {k, 1, n}
];

(* Compute moments from symbolic log Z(β) *)
GeneratingFunctionMoments[Z_, β_, n_] := Module[
  {logZ = Log[Z], μ},
  μ = Table[D[logZ, {β, k}], {k, 1, n}];
  Simplify[μ]
];

End[]
EndPackage[]
