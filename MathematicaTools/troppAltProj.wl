(* ::Package:: *)

(* A rewrite of "troppAltProj.m", available in
   https://github.com/gnikylime/GameofSloanes *)


(* Initial seed *)
Options[InitMat] = {WorkingPrecision -> MachinePrecision};
InitMat[d_, n_, opts : OptionsPattern[]] := 
 Module[{\[Tau], T, done, X, t, x, cor, prec},
  prec = OptionValue[WorkingPrecision];
  \[Tau] = 0.9;
  T = 10000;
  done = 0;
  While[done == 0,
   X = RandomReal[NormalDistribution[], d, WorkingPrecision -> prec] + 
     I RandomReal[NormalDistribution[], d, WorkingPrecision -> prec];
   X = List /@ Normalize[X];
   t = 0;
   While[t < T,
    x = RandomReal[NormalDistribution[], d, WorkingPrecision -> prec] + 
      I RandomReal[NormalDistribution[], d, WorkingPrecision -> prec];
    x = Normalize[x];
    cor = Max@Abs[X\[ConjugateTranspose] . x];
    If[cor < \[Tau],
     X = MapThread[Append, {X, x}];
     t = 0;
     If[Dimensions[X][[2]] == n, t = T];,
     t++;
     ]
    ];
   If[Dimensions[X][[2]] == n, done = 1];
   ];
  X
  ]

(* Alternating projections *)
Options[troppAltProj] = {InitialPrecision -> MachinePrecision, 
   Iterations -> 30000, UpdateInterval -> Infinity};
(* With random initial seed *)
troppAltProj[d_, n_, \[Mu]_, opts : OptionsPattern[]] := Module[{X},
  X = InitMat[d, n, WorkingPrecision -> OptionValue[InitialPrecision]];
  troppAltProj[d, n, \[Mu], X, opts]
  ]
(* With initial seed given as argument Y *)
troppAltProj[d_, n_, \[Mu]_, Y_?MatrixQ, opts : OptionsPattern[]] := 
 Module[
  {X, prec, G, T, t, i, j, dd, V, ind, Vs, first, sumfirst, k, 
   gammapart, D, newD, F, U, lastprec, currentprec},
  X = Y;
  G = X\[ConjugateTranspose] . X;
  T = OptionValue[Iterations];
  If[OptionValue[UpdateInterval] != Infinity,
   lastprec = Precision[X];
   currentprec = Precision[G];
   Print[CurrentDate[], " - Starting loop", " - Precision of G is ", 
    Precision[G], " - Lost ", lastprec - currentprec, 
    " digits of precision"];
   ];
  For[t = 1, t <= T, t++,
   If[Divisible[t, OptionValue[UpdateInterval]],
    lastprec = currentprec;
    currentprec = Precision[G];
    Print[CurrentDate[], " - Iteration ", t, " - Precision of G is ", 
     Precision[G], " - Lost ", lastprec - currentprec, 
     " digits of precision"]
    ];
   For[i = 1, i <= n, i++,
    For[j = 1, j <= n, j++,
     If[i == j, G[[i, j]] = 1,
      If[Abs@G[[i, j]] > \[Mu], 
       G[[i, j]] = G[[i, j]]/Abs[G[[i, j]]]*\[Mu]]
      ]
     ]
    ];
   G = (G + G\[ConjugateTranspose])/2;
   {dd, V} = Chop@Eigensystem[G];
   ind = Ordering[dd, All, Greater];
   Vs = V[[ind]];
   first = dd[[1 ;; d]];
   sumfirst := Plus @@ first;
   If[sumfirst > n,
    While[sumfirst > n,
     k = Plus @@ (If[# > 0, 1, 0] & /@ dd);
     gammapart = (sumfirst - n)/k;
     first = (first - 
         gammapart) (If[# > 0, 1, 0] & /@ (first - gammapart));
     ],
    first += (n - sumfirst)/d;
    ];
   newD = ConstantArray[0, {Length[dd], Length[dd]}];
   newD[[1 ;; d, 1 ;; d]] = DiagonalMatrix[first];
   G = Transpose[Vs] . newD . Vs\[Conjugate];
   ];
  D = DiagonalMatrix@Diagonal[G];
  G = # . G . # &@Inverse[Sqrt[D]];
  G = (G + G\[ConjugateTranspose])/2;
  {U, D, V} = SingularValueDecomposition[G, UpTo[d]];
  (U . Sqrt[D])\[ConjugateTranspose]
  ]
