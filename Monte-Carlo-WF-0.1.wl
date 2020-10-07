(* ::Package:: *)

(* ::Title:: *)
(*Monte Carlo wave function library documentation*)


(* ::Section:: *)
(*Notation and Conventions*)


(* ::Subsection::Closed:: *)
(*Initialization*)


(* ::Input:: *)
(*(*These operators come from a dependence: the QuLib library*)*)


(* ::Input::Initialization:: *)
(*Loads the pauli operators and bosonic operators from que qulib matrices*)
LoadBosonicOperators[nmax = 30];
LoadPauliMatrices[];


(* ::Subsection::Closed:: *)
(*Notation*)


(* ::Input:: *)
(*(*These functions will not be loaded in the .m file*)*)


(* ::Input::Initialization:: *)
(*Aesthetic Way to write tensor products. The function MatrixQ decides wether we should use the function for matrices, kron[], or for vectors, vkron[]*)
a_\[CircleTimes]b_:=If[MatrixQ[a]==True,kron[a,b],kronV[a,b]]

(*This makes an arbitrary number of products define, such as Subscript[|0\[RightAngleBracket], Qubit]\[CircleTimes]Subscript[|0\[RightAngleBracket], Qubit]\[CircleTimes]Subscript[|0\[RightAngleBracket], Qubit]\[CircleTimes]Subscript[|1\[RightAngleBracket], Qubit], etc...*)
a_\[CircleTimes]b_\[CircleTimes]c__:=(a\[CircleTimes]b)\[CircleTimes]c


(* ::Input::Initialization:: *)
(*Fock States and Squeezed States. Usage Example: Fock[2] or Subscript[|2\[RightAngleBracket], Fock]*)

Fock[n_]:=IdentityMatrix[nmax][[n+1]]
Subscript[Ket[n_],Fock]:=Fock[n+1]


(* ::Subsection::Closed:: *)
(*Special Functions*)


(* ::Input::Initialization:: *)
(*Wiener increment for a certain variance Var_ = \[CapitalDelta]t and list size. This is used in diffusive SSEs and SMEs*)
Wiener[Var_, Size_]:=RandomVariate[NormalDistribution[0,Sqrt[Var]],Size]


(* ::Input::Initialization:: *)
(*Innovation.  This term appearns in SMEs*)
\[ScriptCapitalH][\[Rho]_, O_]:=O.\[Rho]+\[Rho].O\[ConjugateTranspose]-Tr[(O+O\[ConjugateTranspose]).\[Rho]]\[Rho]
Subscript[\[ScriptCapitalH], \[Rho]_][O_]:=O.\[Rho]+\[Rho].O\[ConjugateTranspose]-Tr[(O+O\[ConjugateTranspose]).\[Rho]]\[Rho]


(* ::Input::Initialization:: *)
(*Expected value of an operator*)
QAvg[\[Psi]_, O_]:= Conjugate[\[Psi]].O.\[Psi]
QAvg\[Rho][\[Rho]_, O_]:= Tr[\[Rho].O]

(*More compact, but slower, way of writing it with LefAngleBracket and RightAngleBracket;
Usage Example: Subscript[\[LeftAngleBracket]\[Sigma]z\[RightAngleBracket], {0,1}],  which should be -1, Subscript[\[LeftAngleBracket]a\[ConjugateTranspose].a\[RightAngleBracket], Fock[10]] , which should be 9, etc*)
Subscript[\[LeftAngleBracket]O_\[RightAngleBracket], \[Rho]\[Psi]_]:=If[MatrixQ[\[Rho]\[Psi]]==True,Tr[\[Rho]\[Psi].O],Conjugate[\[Rho]\[Psi]].O.\[Rho]\[Psi]]


(* ::Input::Initialization:: *)
(*Commutators and Anti-Commuators. These can be written with DoubleStruckCapitalC and A, respectively*)
\[DoubleStruckCapitalC][A_,B_]:=A.B-B.A
\[DoubleStruckCapitalA][A_,B_]:=A.B+B.A


(* ::Input::Initialization:: *)
(*Dissipator*)
\[DoubleStruckCapitalD][\[Rho]_, L_]:=L.\[Rho].L\[ConjugateTranspose] - 1/2 \[DoubleStruckCapitalA][L\[ConjugateTranspose].L,\[Rho]]


(* ::Input::Initialization:: *)
(*Q - function*)
Q[\[Alpha]_, O_]:=1/\[Pi] CoherentState[\[Alpha]].O.CoherentState[\[Alpha]]


(* ::Input::Initialization:: *)
(*Draws a Bloch Sphere*)
BlochSphere:=Module[{},{
SphericalPlot3D[1,{t,0,Pi},{p,0,2*Pi},PlotStyle->Directive[White,Opacity[0.3],Specularity[White,10]],Mesh->1,MeshStyle->Opacity[0.3]],
Graphics3D[{
{Black,PointSize[0.025],{Point[{0,0,0}]}},
{Thickness-> 0.005,Black,Line[{{-1,0,0},{1,0,0}}],Line[{{0,-1,0},{0,1,0}}],Line[{{0,0,-1},{0,0,1}}]}/.Line->Composition[Arrow,Tube],
{Thickness-> 0.005,Black,Line[{{1,0,0},{-1,0,0}}],Line[{{0,1,0},{0,-1,0}}],Line[{{0,0,1},{0,0,-1}}]}/.Line->Composition[Arrow,Tube]
}]
}
]


(* ::Input::Initialization:: *)
(*Axis Labels for the Bloch Sphere*)
BlochLabel:=Module[{},{AxesLabel-> {Style["(\!\(\*SubscriptBox[\(\[Sigma]\), \(x\)]\))",16,Bold],Style["(\!\(\*SubscriptBox[\(\[Sigma]\), \(y\)]\))",16,Bold],Style["(\!\(\*SubscriptBox[\(\[Sigma]\), \(z\)]\))",16,Bold]},AxesStyle->Thickness[0.01],LabelStyle->Directive[Bold]}]


(* ::Section:: *)
(*Documentation*)


(* ::Subsection::Closed:: *)
(*Photon Counting*)


(* ::Input:: *)
(*(* These two functions are implementations of the photo detection protocol. In this implementation we use two Krauss operators;*)
(*The first function is an implementation for the Bosonic case, the second function is an implementation for the qubit case.*)
(**)
(*> Arguments;*)
(*	- \[Psi]0: the initial wave-function;*)
(*	- H: the Hamiltonian;*)
(*	- nsteps: the number of steps;*)
(*	- \[CapitalDelta]t: the temporal step size;*)
(* *)
(*> Returns;*)
(*	- A list \[Psi] of the type {\[Psi]0, \[Psi]1, \[Psi]2, ..., \[Psi]nsteps} where \[Psi]i represents the wave function at the i-th iteration;*)
(**)*)


(* ::Input::Initialization:: *)
PhotoDetectionBoson[\[Psi]0_,H_, nsteps_,\[CapitalDelta]t_] := Module[{\[Psi]},
K1 = \[DoubleStruckCapitalI] -\[CapitalDelta]t(I H + 1/2 a\[ConjugateTranspose].a);
K2 = Sqrt[\[CapitalDelta]t] a;
\[ScriptCapitalE]1 = K1\[ConjugateTranspose].K1;\[ScriptCapitalE]2 = K2\[ConjugateTranspose].K2;
\[Psi]={\[Psi]0};
rands = RandomReal[{0,1},nsteps];
Do[
p1 =Re[Conjugate[(\[Psi][[-1]])].\[ScriptCapitalE]1.(\[Psi][[-1]])];
p2 = Re[Conjugate[(\[Psi][[-1]])].\[ScriptCapitalE]2.(\[Psi][[-1]])];
If[x<=p1,AppendTo[\[Psi],K1.(\[Psi][[-1]])/Sqrt[p1]],AppendTo[\[Psi],K2.(\[Psi][[-1]])/Sqrt[p2]]]
,{x,rands}];
\[Psi]
]

PhotoDetectionQubit[\[Psi]0_,H_, nsteps_,\[CapitalDelta]t_] := Module[{\[Psi]},
K1 = IdentityMatrix[2] - \[CapitalDelta]t(I H + 1/2 \[Sigma]p.\[Sigma]m);
K2 = Sqrt[\[CapitalDelta]t] \[Sigma]m;
\[ScriptCapitalE]1 = K1\[ConjugateTranspose].K1; \[ScriptCapitalE]2 = K2\[ConjugateTranspose].K2;
\[Psi]={\[Psi]0};
rands = RandomReal[{0,1},nsteps];
Do[
p1 =Re[Conjugate[(\[Psi][[-1]])].\[ScriptCapitalE]1.(\[Psi][[-1]])];
p2 = Re[Conjugate[(\[Psi][[-1]])].\[ScriptCapitalE]2.(\[Psi][[-1]])];
If[x<=p1,AppendTo[\[Psi],K1.(\[Psi][[-1]])/Sqrt[p1]],AppendTo[\[Psi],K2.(\[Psi][[-1]])/Sqrt[p2]]]
,{x,rands}];

(*Returns a list \[Psi] with the WF vector at each time step, like {(Subscript[\[Psi], 0](0)
Subscript[\[Psi], 1](0)

),(Subscript[\[Psi], 0](\[CapitalDelta]t)
Subscript[\[Psi], 1](\[CapitalDelta]t)

),(Subscript[\[Psi], 0](2\[CapitalDelta]t)
Subscript[\[Psi], 1](2\[CapitalDelta]t)

),...}*)
\[Psi]
]


PhotoDetectionBoson::usage="PhotoDetectionBoson[\[Psi]0_, H_, nsteps_, \[CapitalDelta]t_] 
Returns the stochastic evolution given an initial wave function, the Hamiltonian and the integration parameters. 
The implementation is based on the application of the Krauss Operators for the bosonic case";

PhotoDetectionQubit::usage="PhotoDetectionQubit[\[Psi]0_, H_, nsteps_, \[CapitalDelta]t_] 
Returns the stochastic evolution given an initial wave function, the Hamiltonian and the integration parameters. 
The implementation is based on the application of the Krauss Operators for the qubitcase";


(* ::Input:: *)
(*(* The same as the function above, but the implementation is done through the (normalized) SSE itself, for an arbitrary jump operator.*)
(**)
(*> Arguments;*)
(*	- \[Psi]0: the initial wave-function;*)
(*	- H: the Hamiltonian;*)
(*	- c: the jump operator;*)
(*	- nsteps: the number of steps;*)
(*	- \[CapitalDelta]t: the temporal step size;*)
(* *)
(*> Returns;*)
(*	- A list \[Psi] of the type {\[Psi]0, \[Psi]1, \[Psi]2, ..., \[Psi]nsteps} where \[Psi]i represents the wave function at the i-th iteration;*)
(**)*)


(* ::Input::Initialization:: *)
PhotoDetection[\[Psi]0_,H_, c_, nsteps_,\[CapitalDelta]t_] := Module[{\[Psi]},
\[Psi]={\[Psi]0};
rands = RandomReal[{0,1},nsteps];
Do[
p =\[CapitalDelta]t;
If[x<=p,dN=1,dN=0];
d\[Psi]= - \[CapitalDelta]t*(I H + 1/2 (c\[ConjugateTranspose].c-Re[QAvg[\[Psi][[-1]],c\[ConjugateTranspose].c]]*IdentityMatrix[Length[c]] )).\[Psi][[-1]]+dN*(c/Sqrt[Re[QAvg[\[Psi][[-1]],c\[ConjugateTranspose].c]]]-IdentityMatrix[Length[c]]).\[Psi][[-1]];
AppendTo[\[Psi],d\[Psi]+\[Psi][[-1]]];
If[Norm[\[Psi][[-1]]]>1,\[Psi][[-1]]=Normalize[\[Psi][[-1]]]];
,{x,rands}];
\[Psi]
]


(* ::Input::Initialization:: *)
PhotoDetection::usage="PhotoDetection[\[Psi]0_,H_, c_, nsteps_,\[CapitalDelta]t_] 
Returns the stochastic evolution given an initial wave function, the Hamiltonian, the jump operator and the integration parameters. 
The implementation is based on a jump-type Stochastic Schr\[ODoubleDot]dinger Equation";


(* ::Input:: *)
(*(* The same as the functions above, but this time we keep track of the jump record. The implementation is once again done trough Krauss Operators.*)
(**)
(*> Arguments;*)
(*	- \[Psi]0: the initial wave-function;*)
(*	- H: the Hamiltonian;*)
(*	- c: the jump operator;*)
(*	- nsteps: the number of steps;*)
(*	- \[CapitalDelta]t: the temporal step size;*)
(* *)
(*> Returns;*)
(*	- A list \[Psi] of the type {{\[Psi]0, dN0}, {\[Psi]1, dN1}, {\[Psi]2,dN2} , ... } where \[Psi]i represents the wave function at the i-th iteration and dNi is the jump record;*)
(**)*)


(* ::Input::Initialization:: *)
PhotoDetectionJumpRecord[\[Psi]0_,H_, c_, nsteps_,\[CapitalDelta]t_] := Module[{\[Psi]},
K1 = \[DoubleStruckCapitalI] -\[CapitalDelta]t(I H + 1/2 c\[ConjugateTranspose].c);
K2 = Sqrt[\[CapitalDelta]t] c;
\[ScriptCapitalE]1 = K1\[ConjugateTranspose].K1;\[ScriptCapitalE]2 = K2\[ConjugateTranspose].K2;
\[Psi]={{\[Psi]0,0}};
rands = RandomReal[{0,1},nsteps];
Do[
p1 =Re[Conjugate[(\[Psi][[-1, 1]])].\[ScriptCapitalE]1.(\[Psi][[-1, 1]])];
p2 = Re[Conjugate[(\[Psi][[-1, 1]])].\[ScriptCapitalE]2.(\[Psi][[-1, 1]])];
If[x<=p1,AppendTo[\[Psi],{K1.(\[Psi][[-1, 1]])/Sqrt[p1],0}],AppendTo[\[Psi],{K2.(\[Psi][[-1, 1]])/Sqrt[p2],1}]]
,{x,rands}];
\[Psi]
]


(* ::Input::Initialization:: *)
PhotoDetectionJumpRecord::usage="PhotoDetectionJumpRecord[\[Psi]0_, H_, c_, nsteps_, \[CapitalDelta]t_] 
The same as PhotoDetection[...] but it returns elements of the type {\[Psi]t,dN(t)} with the measurement record dN";


(* ::Subsection::Closed:: *)
(*Statistics*)


(* ::Input:: *)
(*(* Plots a customized ListPlot which appears as a bar chart with the corresponding probability pn = |<n|\[Psi]>(|^2) for each fock state given a wave function \[Psi]*)
(**)
(*> Arguments;*)
(*	- \[Psi]: the wave-function;*)
(*	- range: the dimension of the bosonic space;*)
(*	- style: plot style;*)
(*	- extras: any other graphic objects to be added to the ListPlot;*)
(* *)
(*> Returns;*)
(*	- A ListPlot with the appearence of a bar chart;*)
(**)*)


(* ::Input::Initialization:: *)
(*Histogram for the Fock Basos probabilities of an arbitrary wave function*)
PhotonHistogram[\[Psi]_, range_, style_:{},extras_:{}] := Module[{Probabilities},

(*Computes |<n|\[Psi]>|^2*)
Probabilities=Table[ {i-1,Abs[(Fock[i].\[Psi])]^2},{i,2,range}];

ListPlot[
Probabilities,
(*We mark each point with a filled square and shade the region below. The thickness of each bar is inversely proportional to the plot range*)
PlotStyle->{style},PlotMarkers->{"\[FilledSquare]",Scaled[1.8/range]}, PlotRange-> {{0,range},{0,1.02Max[Probabilities[[All,2]]]}}
,Filling->Axis,FillingStyle->Directive[Thickness-> 1/range,Opacity[.25]],
extras
]
]


(* ::Input::Initialization:: *)
PhotonHistogram::usage="PhotonHistogram[\[Psi]_, range_, style_:{},extras_:{}] 
Plots a customized ListPlot which appears as a bar chart with the corresponding probability p = |\[LeftAngleBracket]n|\[Psi]\[RightAngleBracket]\!\(\*SuperscriptBox[\(|\), \(2\)]\) for each fock state given a wave function \[Psi]";


(* ::Input:: *)
(*(* This is similar to the PhotoDetection[...] and similar functions, but instead of the wave function, we record the time interval between consecutive jumps*)
(**)
(*> Arguments;*)
(*	- \[Psi]0: the initial wave-function;*)
(*	- H: the Hamiltonian;*)
(*	- c: the jump operator;*)
(*	- ensemblesize: the number of realizations we perform for this experiment;*)
(*	- nsteps: the number of steps;*)
(*	- \[CapitalDelta]t: the temporal step size;*)
(* *)
(*> Returns;*)
(*	- A list of the type {Subscript[t, 1], Subscript[t, 2], ...} with the time at which each jump occurs. This distribution can be plotted with a histogram*)
(**)*)


(* ::Input::Initialization:: *)
(*Time interval between two consecutive jumps in a histogram*)
(*Note that we use a double blank __ after style. This way we can use as many plot options as we want*)
StochasticJumpTimes[\[Psi]0_,H_,c_, ensemblesize_,steps_,\[CapitalDelta]t_] := Module[{\[Psi], Intervals,Poslist},

Intervals={};

(*Loop for several trajectories. We record the the intervals in the list above*)
Do[{
(*Wave function evolutions*)
\[Psi]=PhotoDetectionJumpRecord[\[Psi]0,H ,c,steps,\[CapitalDelta]t];
 (*The first entry in the vector is the wave function itself, the second entry is the binary stochastic variable N = 0,1*)
Poslist=Flatten@Position[\[Psi][[ All, 2]],1];
(*Appends the time differences into the list. We have to flatten it otherwise the data would come out as {{Subscript[t, 1]},{Subscript[t, 2]},...}*)
AppendTo[Intervals,Flatten@Table[\[CapitalDelta]t(Poslist[[i]]-Poslist[[i-1]]),{i,2,Length@Poslist}]];
},{j,1,ensemblesize}];
Flatten@Intervals
]


(* ::Input::Initialization:: *)
StochasticJumpTimes::usage="StochasticJumpTimes[\[Psi]0_,H_,c_, ensemblesize_,steps_,\[CapitalDelta]t_] 
Returns a list of the type {\!\(\*SubscriptBox[\(t\), \(1\)]\), \!\(\*SubscriptBox[\(t\), \(2\)]\), ...} with the time at which each jump occurs for a photon counting protocol";


(* ::Input:: *)
(*(* Plots a histogram with the number of jumps within a certain time interval. This time interval is determined by \[Psi], whose time interval is determined by \[CapitalDelta]t * nsteps;*)
(**)
(*> Arguments;*)
(*	- \[Psi]:  A list \[Psi] of the type {{\[Psi]0, dN0}, {\[Psi]1, dN1}, {\[Psi]2,dN2} , ... } where \[Psi]i represents the wave function at the i-th iteration and dNi is the jump record, *)
(*given by PhotoDetectionJumpRecord[...];*)
(*	- ensemblesize: the number of realizations we perform for this experiment;*)
(*	- style: plot style;*)
(*	- extras: any other graphic objects to be added to the ListPlot;*)
(* *)
(*> Returns;*)
(*	- Returns a list with {Histogram, Average Numer of Jumps, Total Number of Jumps}*)
(**)*)


(* ::Input::Initialization:: *)
(*We should pay attenction to the fact that the WF evolution should be called inside a Hold[], otherwise we'll use the same trajectory in every loop*)
JumpDistribution[\[Psi]_, ensemblesize_, style_:{}, extra_:{}] := Module[{HistoClicks={},Average,Poslist,AverageSquared,Histo,TotalClicks},

(*Runs the simulation over several trajectories*)
Do[
(*PhotoDetectionJumpRecord returns a function of the type {{Subscript[\[Psi], 0],Subscript[N, 0]},{Subscript[\[Psi], 1],Subscript[N, 1]},..}. The second entry is the binary variable corresponding to No Jump/Jump*)
Poslist=ReleaseHold[\[Psi]][[All, 2]]; (*Here we use Release Hold otherwise the moduele doesn't evaluate the function*)
(*A list with all the instants where a jump occurs*)
PosClicklist=Flatten@Position[Poslist,1];
(*The total number of Jumps*)
NumberofClicks=Length@Flatten@Position[Poslist,1];
(*Stores the number of Clicks in this trajectory*)
AppendTo[HistoClicks,NumberofClicks];
,{j,1,ensemblesize}];
(*Total number of jumps across all the realizations*)
TotalClicks=Total[Table[Count[HistoClicks,k],{k,0,Max[HistoClicks]+1}]];
(*Average number of Jumps in the interval T*)
Average=Total[HistoClicks]/ensemblesize;
(*Average of Squares*)
AverageSquared=Total[HistoClicks^2]/ensemblesize;
(*Plots the Histogram*)
Histo=Histogram[HistoClicks,{"Raw",Max[HistoClicks]+1},style,extra,PlotRange->{{0,Max[HistoClicks]+1},Automatic}];
(*Returns a list with {Histogram, Average Numer of Jumps, Total Number of Jumps (for purposes of normalization, if necessary)}*)
{Histo,Average,AverageSquared,TotalClicks}
]


(* ::Input::Initialization:: *)
JumpDistribution::usage="JumpDistribution[\[Psi]_, ensemblesize_, style_:{}, extra_:{}] 
Plots a histogram with the number of jumps within a certain time interval given a
\[Psi] of the type {{\[Psi]0, dN0}, {\[Psi]1, dN1}, {\[Psi]2,dN2} , ... } where \[Psi]i represents the wave function at the i-th iteration and dNi is the jump record ";


(* ::Subsection::Closed:: *)
(*SSEs*)


(* ::Input:: *)
(*(* These two functions are implementations of the homodyne detection protocol;*)
(*The first function is an implementation for the Bosonic case, the second function is an implementation for the qubit case.*)
(**)
(*> Arguments;*)
(*	- \[Psi]0: the initial wave-function;*)
(*	- H: the Hamiltonian;*)
(*	- nsteps: the number of steps;*)
(*	- \[CapitalDelta]t: the temporal step size;*)
(* *)
(*> Returns;*)
(*	- A list \[Psi] of the type {\[Psi]0, \[Psi]1, \[Psi]2, ..., \[Psi]nsteps} where \[Psi]i represents the wave function at the i-th iteration;*)
(**)*)


(* ::Input::Initialization:: *)
(*Homodyne Measurement for a phase \[CapitalPhi]*)
BosonHomodyneSSE[\[CapitalPhi]_,\[Psi]0_,H_, nsteps_,\[CapitalDelta]t_] := Module[{\[Psi], \[Gamma]=1},
W=RandomVariate[NormalDistribution[0,Sqrt[\[CapitalDelta]t]],nsteps];
\[Psi]={\[Psi]0};
stoevo\[Theta]={};stoevo\[Phi]={};
Do[
(*Stochastic Schrodinger Equation for Homodyne Measurements*)
\[CapitalDelta]\[Psi]=-\[CapitalDelta]t(\[Gamma]/2 a\[ConjugateTranspose].a +I H ).\[Psi][[-1]] +E^(-I \[CapitalPhi]) (\[Gamma] \[CapitalDelta]t QAvg[\[Psi][[-1]],E^(+ I \[CapitalPhi])  a\[ConjugateTranspose]+E^(-I \[CapitalPhi]) a]+Sqrt[\[Gamma]]dW)a.\[Psi][[-1]];

(*Normalizes the WF*)
AppendTo[\[Psi],(\[Psi][[-1]]+\[CapitalDelta]\[Psi])/Norm[\[Psi][[-1]]+\[CapitalDelta]\[Psi]]];
,{dW,W}];
\[Psi] 
]

(*Homodyne Measurement for a phase \[CapitalPhi] for a qubit*)
QubitHomodyneSSE[\[CapitalPhi]_,\[Psi]0_,H_, nsteps_,\[CapitalDelta]t_] := Module[{\[Psi], \[Gamma]=1},
W=RandomVariate[NormalDistribution[0,Sqrt[\[CapitalDelta]t]],nsteps];
\[Psi]={\[Psi]0};
Do[
(*Stochastic Schrodinger Equation for Homodyne Measurements*)
\[CapitalDelta]\[Psi]=-\[CapitalDelta]t(\[Gamma]/2 \[Sigma]p.\[Sigma]m +I H ).\[Psi][[-1]] +E^(-I \[CapitalPhi]) (\[Gamma] \[CapitalDelta]t QAvg[\[Psi][[-1]],E^(+ I \[CapitalPhi])  \[Sigma]p+E^(-I \[CapitalPhi]) \[Sigma]m]+Sqrt[\[Gamma]]dW)\[Sigma]m.\[Psi][[-1]];

(*Normalizes the WF*)
AppendTo[\[Psi],(\[Psi][[-1]]+\[CapitalDelta]\[Psi])/Norm[\[Psi][[-1]]+\[CapitalDelta]\[Psi]]];
,{dW,W}];
\[Psi] 
]


(* ::Input::Initialization:: *)
BosonHomodyneSSE::usage="BosonHomodyneSSE[\[Psi]0_, H_, nsteps_, \[CapitalDelta]t_] 
Returns the stochastic evolution for a homodyne measurement protocol in the fock space given an initial wave function, the Hamiltonian and the integration parameters.";

QubitHomodyneSSE::usage="QubitHomodyneSSE[\[Psi]0_, H_, nsteps_, \[CapitalDelta]t_] 
Returns the stochastic evolution for a homodyne measurement protocol in que qubit space given an initial wave function, the Hamiltonian and the integration parameters.";


(* ::Input:: *)
(*(* Implementation of the homodyne detection protocol for an arbitrary measurement c*)
(**)
(*> Arguments;*)
(*	- \[Psi]0: the initial wave-function;*)
(*	- H: the Hamiltonian;*)
(*	- c: measurement operator;*)
(*	- nsteps: the number of steps;*)
(*	- \[CapitalDelta]t: the temporal step size;*)
(* *)
(*> Returns;*)
(*	- A list \[Psi] of the type {\[Psi]0, \[Psi]1, \[Psi]2, ..., \[Psi]nsteps} where \[Psi]i represents the wave function at the i-th iteration;*)
(**)*)


(* ::Input::Initialization:: *)
(*Homodyne Measurement for a phase \[CapitalPhi]*)
HomodyneSSE[\[CapitalPhi]_,\[Psi]0_,H_,c_, nsteps_,\[CapitalDelta]t_] := Module[{\[Psi], \[Gamma]=1},
W=RandomVariate[NormalDistribution[0,Sqrt[\[CapitalDelta]t]],nsteps];
\[Psi]={\[Psi]0};
stoevo\[Theta]={};stoevo\[Phi]={};
Do[
(*Stochastic Schrodinger Equation for Homodyne Measurements*)
\[CapitalDelta]\[Psi]=-\[CapitalDelta]t(\[Gamma]/2 c\[ConjugateTranspose].c +I H ).\[Psi][[-1]] +E^(-I \[CapitalPhi]) (\[Gamma] \[CapitalDelta]t QAvg[\[Psi][[-1]],E^(+ I \[CapitalPhi])  c\[ConjugateTranspose]+E^(-I \[CapitalPhi]) c]+Sqrt[\[Gamma]]dW)c.\[Psi][[-1]];

(*Normalizes the WF*)
AppendTo[\[Psi],(\[Psi][[-1]]+\[CapitalDelta]\[Psi])/Norm[\[Psi][[-1]]+\[CapitalDelta]\[Psi]]];
,{dW,W}];
\[Psi] 
]



(* ::Input::Initialization:: *)
HomodyneSSE::usage="HomodyneSSE[\[CapitalPhi]_,\[Psi]0_,H_,c_, nsteps_,\[CapitalDelta]t_]
Returns the stochastic evolution for a homodyne measurement protocol for an arbitrary operator c given an initial wave function, the Hamiltonian and the integration parameters.";


(* ::Input:: *)
(*(* Plot of the two angles in the Bloch sphere for a qubit undergoing a homodyne measurement*)
(**)
(*> Arguments;*)
(*	- \[Psi]0: the initial wave-function;*)
(*	- H: the Hamiltonian;*)
(*	- nsteps: the number of steps;*)
(*	- \[CapitalDelta]t: the temporal step size;*)
(* *)
(*> Returns;*)
(*	- A plot of Cos[\[Theta]] and \[Phi]*)
(**)*)


(* ::Input::Initialization:: *)
(*Homodyne Measurement for a phase \[CapitalPhi]*)
QubitHomodyneSSEPlot[\[CapitalPhi]_,\[Psi]0_,H_, nsteps_,\[CapitalDelta]t_] := Module[{\[Psi], \[Gamma]=1},
W=RandomVariate[NormalDistribution[0,Sqrt[\[CapitalDelta]t]],nsteps];
\[Psi]={\[Psi]0};
stoevo\[Theta]={};stoevo\[Phi]={};
Do[
(*Stochastic Schrodinger Equation for Homodyne Measurements*)
\[CapitalDelta]\[Psi]=-\[CapitalDelta]t(\[Gamma]/2 \[Sigma]p.\[Sigma]m +I H ).\[Psi][[-1]] +E^(-I \[CapitalPhi]) (\[Gamma] \[CapitalDelta]t QAvg[\[Psi][[-1]],E^(+ I \[CapitalPhi])  \[Sigma]p+E^(-I \[CapitalPhi]) \[Sigma]m]+Sqrt[\[Gamma]]dW)\[Sigma]m.\[Psi][[-1]];

(*Normalizes the WF*)
AppendTo[\[Psi],(\[Psi][[-1]]+\[CapitalDelta]\[Psi])/Norm[\[Psi][[-1]]+\[CapitalDelta]\[Psi]]];

(*Store \[Theta] = 2 Tan^-1(|(C0/C1)|) and \[Phi] = Arg(C0.C1^*)*)
AppendTo[stoevo\[Theta], Cos[2ArcTan[Abs[\[Psi][[-1,1]]/\[Psi][[-1,2]]]]]];
AppendTo[stoevo\[Phi], Arg[\[Psi][[-1,1]]*Conjugate[\[Psi][[-1,2]]]]];

,{dW,W}];

(*Plot the results*)
GraphicsRow[{
Show[ListLinePlot[stoevo\[Theta],DataRange->{0,nsteps*\[CapitalDelta]t}, PlotLabel-> "Cos \[Theta]"]],
Show[ListLinePlot[stoevo\[Phi],DataRange->{0,nsteps*\[CapitalDelta]t},PlotRange-> {-\[Pi],\[Pi]}, PlotLabel-> "\[Phi]"]]
}]
]


(* ::Input::Initialization:: *)
QubitHomodyneSSEPlot::usage="QubitHomodyneSSEPlot[\[CapitalPhi]_,\[Psi]0_,H_, nsteps_,\[CapitalDelta]t_]
Plot of Cos[\[Theta]] and \[Phi] for a qubit undergoing a homodyne measurement.;


(* ::Subsection::Closed:: *)
(*SSEs Plots*)


(* ::Input:: *)
(*(* Builds a list with stochastic evolution of a certain protocol. Each row represent a time evolution, and reach column represents the emsemble at a time t;*)
(**)
(*> Arguments;*)
(*	- function: the wave-function evolution \[Psi], i.e., the list returned by functions which implement SSEs such as PhotoDetection;*)
(*   - size: the number of trajectories/the ensemble size;*)
(**)
(*> Returns;*)
(*	- A list \[Psi] of the type {{\[Psi]00, \[Psi]01, \[Psi]02, ...}, {\[Psi]10, \[Psi]11, \[Psi]12...}, ..., {\[Psi](size)0, \[Psi](size)1, ....}} where \[Psi]ij represents the wave function at the j-th iteration for the i-th trajectory;*)
(*	In matrix form:*)
(*(\[NoBreak]\[Psi]00	\[Psi]10	.	\[Psi](size)0*)
(*\[Psi]01	\[Psi]11	.	\[Psi](size)1*)
(*.	.	.	\[Psi](size)2*)
(*\[Psi]0n	\[Psi]1n	.	\[Psi](size)3*)
(**)
(*\[NoBreak])*)
(*Note that \[Psi]ij does not represent  a w.f. component but rather the WHOLE wave funtion, if we take into account the k-th wave function component itself we will need three indices! As in \[Psi]ijk.*)
(**)*)


(* ::Input::Initialization:: *)
(* This function must be called as StochasticEnsemble[Hold["Measurement Protocol"], "Ensemble Size"]
The Hold[] is necessary because otherwise the function only evaluates the measurement once and it'll use a single realization. The ReleaseHold in the function
body prevents this, guaranteeing that we'll get k different realizations of the system, building an Ensemble with a different trajectory for each experiment. *)
StochasticEnsemble[function_,size_]:=Module[{},Table[ReleaseHold[function],{k,1,size}]]


(* ::Input::Initialization:: *)
StochasticEnsemble::usage="StochasticEnsemble[function_,size_]
Returns an ensenble of several realizations of SSE/SME realizations.";


(* ::Input:: *)
(*(* Plots the average of several trajectories for measurement scheme;*)
(**)
(*> Arguments;*)
(*	- \[Psi]Ensemble: a list of the type  {{\[Psi]00, \[Psi]01, \[Psi]02, ...}, {\[Psi]10, \[Psi]11, \[Psi]12...}, ...,} with an ensemble fo several trajectories ;*)
(*	- obs: the observable whose expected value will be plotted;*)
(*	- style: any additional options for the ListLinePlot;*)
(* *)
(*> Returns;*)
(*	- A ListLinePlot with the average;*)
(**)*)


(* ::Input::Initialization:: *)
(*Average over an Ensemble*)
StochasticAveragePlot[\[Psi]Ensemble_,obs_,style_] := Module[{Eavg},
Eavg=Table[Mean[Table[Re[QAvg[\[Psi]Ensemble[[All,j]][[i]],obs]],{i,1,Length@\[Psi]Ensemble}]],{j,1,Length@Transpose[\[Psi]Ensemble]}] ;
ListLinePlot[Eavg, style]
]


(* ::Input::Initialization:: *)
StochasticAveragePlot::usage="StochasticAveragePlot[\[Psi]Ensemble_,obs_,style_] 
Plots the expected value of an observable, averaged over several trajectories";


(* ::Input:: *)
(*(* Plots all the trajectories for a measurement scheme;*)
(**)
(*> Arguments;*)
(*	- \[Psi]Ensemble: a list of the type  {{\[Psi]00, \[Psi]01, \[Psi]02, ...}, {\[Psi]10, \[Psi]11, \[Psi]12...}, ...,} with an ensemble fo several trajectories ;*)
(*	- obs: the observable whose expected value will be plotted;*)
(*	- style: any additional options for the ListLinePlot;*)
(* *)
(*> Returns;*)
(*	- A ListLinePlot with all the trajectories;*)
(**)*)


(* ::Input::Initialization:: *)
(*Plot all the evolutions in an Ensemble*)
StochasticPlot[\[Psi]Ensemble_,obs_,style_] := Module[{},
Table[ListLinePlot[Table[Re[QAvg[\[Psi]Ensemble[[All,i]][[j]],obs]],{i,1,Length@Transpose[\[Psi]Ensemble]}] , PlotStyle->{Opacity[0.1,Gray]}, style],{j,1,Length@\[Psi]Ensemble}]
]


(* ::Input::Initialization:: *)
StochasticPlot::usage="StochasticPlot[\[Psi]Ensemble_,obs_,style_] 
Plots the expected value of an observable, for all the trajectories";


(* ::Input:: *)
(*(* Plots a single trajectory for a measurement scheme;*)
(**)
(*> Arguments;*)
(*	- \[Psi]Ensemble: a list of the type  {{\[Psi]00, \[Psi]01, \[Psi]02, ...}, {\[Psi]10, \[Psi]11, \[Psi]12...}, ...,} with an ensemble fo several trajectories ;*)
(*	- obs: the observable whose expected value will be plotted;*)
(*	- style: any additional options for the ListLinePlot;*)
(* *)
(*> Returns;*)
(*	- A ListLinePlot with all the trajectories;*)
(**)*)


(* ::Input::Initialization:: *)
(*Average over an Ensemble*)
StochasticHighlight[\[Psi]Ensemble_,obs_,style_] := Module[{},
ListLinePlot[Table[Re[QAvg[\[Psi]Ensemble[[All,i]][[1]],obs]],{i,1,Length@Transpose[\[Psi]Ensemble]}] ,PlotStyle->{RGBColor[0,0.55,0.31]}, style]
]


(* ::Input::Initialization:: *)
StochasticHighlight::usage="StochasticHighlight[\[Psi]Ensemble_,obs_,style_] 
Plots the expected value of an observable for a single trajectory";


(* ::Input:: *)
(*(* Combines all the plots above*)
(**)
(*> Arguments;*)
(*	- \[Psi]Ensemble: a list of the type  {{\[Psi]00, \[Psi]01, \[Psi]02, ...}, {\[Psi]10, \[Psi]11, \[Psi]12...}, ...,} with an ensemble fo several trajectories ;*)
(*	- obs: the observable whose expected value will be plotted;*)
(*	- style: any additional options for the ListLinePlot;*)
(* *)
(*> Returns;*)
(*	- A ListLinePlot with the average, all the trajectories and a single highlighted trajectory;*)
(**)*)


(* ::Input::Initialization:: *)
(*Plots the three curves above*)
StochasticFullPlot[\[Psi]Ensemble_,obs_, range_] := Module[{},
RandomPlot = StochasticPlot[\[Psi]Ensemble,obs,{DataRange->range}] ;
AveragePlot = StochasticAveragePlot[\[Psi]Ensemble,obs,{DataRange->range,PlotStyle->{Red}}] ;
HighlightedPlot = StochasticHighlight[\[Psi]Ensemble,obs,{DataRange->range,PlotStyle->{Red}}] ;
Show[RandomPlot,HighlightedPlot,AveragePlot]
]


(* ::Input::Initialization:: *)
StochasticFullPlot::usage="StochasticFullPlot[\[Psi]Ensemble_,obs_,style_] 
Plots the expected value of an observable, averaged over several trajectories, for all the trajectories and for a single trajectory";


(* ::Input:: *)
(*(* Plots the expected value of an observable in a photon counting protocol, highliting the instants in which jumps occur*)
(**)
(*> Arguments;*)
(*	- \[Psi]0: the initial wave-function;*)
(*	- H: the Hamiltonian;*)
(*	- c: the jump operator;*)
(*	- steps: the number of steps;*)
(*	- \[CapitalDelta]t: the temporal step size;*)
(*	- obs-: the observable to be plotted*)
(**)
(*> Returns;*)
(*	- A ListLinePlot which highligths the jump instants with dashed lines;*)
(**)*)


(* ::Input::Initialization:: *)
(*Stochastic Plot with Jump Markings*)
StochasticJumpPlot[\[Psi]0_,H_,c_, steps_,\[CapitalDelta]t_,obs_] := Module[{\[Psi]},
\[Psi]=PhotoDetectionJumpRecord[\[Psi]0,H,c , steps,\[CapitalDelta]t];
(*The wave funciton is divided as (\[NoBreak](\[NoBreak]Subscript[\[Psi], 0]
Subscript[\[Psi], 1]

\[NoBreak])
N

\[NoBreak]), the third entry is the binary result of the record*)
JumpRecord=\[Psi][[All,2]];
(*Vector of the point in time at which each jump occurs*)
JumpInstants=Flatten@Position[JumpRecord,1];
(*Plot the jumps. The range of the y axis is between the minimum and the maximum of the observable. Here we use a Map to calculate the minimum and
maximum of the expected value. The map works as f[#]&/@{a,b,c,d,e}={f[a],f[b],...} *)
Clicks=Table[ContourPlot[x==JumpInstants[[i]]\[CapitalDelta]t,{x,0,steps*\[CapitalDelta]t},{y,Min[Re[QAvg[#,obs]]&/@\[Psi][[All,1]]],Max[Re[QAvg[#,obs]]&/@\[Psi][[All,1]]]},ContourStyle-> {Black, Dashed}],{i,1,Length@JumpInstants}];
TrajectoryPlot=StochasticHighlight[{\[Psi][[All,1]]},obs,{DataRange->{0,steps*\[CapitalDelta]t},PlotRange-> {Min[Re[QAvg[#,obs]]&/@\[Psi][[All,1]]],Max[Re[QAvg[#,obs]]&/@\[Psi][[All,1]]]}}];
Show[TrajectoryPlot,If[Length@JumpInstants>0,Clicks,{}]]
]


(* ::Input::Initialization:: *)
StochasticJumpPlot::usage="StochasticJumpPlot[\[Psi]0_,H_,c_, steps_,\[CapitalDelta]t_,obs_] 
Plots the expected value of an observable in a photon counting protocol, highliting the instants in which jumps occur";


(* ::Subsection::Closed:: *)
(*Bloch Sphere Plots*)


(* ::Input:: *)
(*(* This function is used to plot the trajectory of a Qubit on the sphere, given the wave function evolution for a certain time interval*)
(**)
(*> Arguments;*)
(*	- \[Psi]: stochastic evolution of the w.f. in the form {\[Psi]0, \[Psi]1, \[Psi]2, ..., \[Psi]nsteps} where \[Psi]i represents the wave function at the i-th iteration.*)
(*We can get this through PhotoDetection[], QubitHomodyneSSE[], and other similar functions;*)
(*	- extra: any other graphic objects to be drawn, such as parametric plots, arrows and so on;*)
(*	- Style3D: any other style choice, such as frame, color and so on;*)
(* *)
(*> Returns;*)
(*	- A graphich object cointaining the bloch sphere with the trajectory described by \[Psi] and the additional parameters*)
(**)*)


(* ::Input::Initialization:: *)
(*Static Trajectory on the Bloch sphere given the evolution of the wave function*)
BlochTrajectory[\[Psi]_,extra_,Style3D_]:=Module[{},

(*Coordinates in the Bloch Sphere*)
trajX=Table[Re[QAvg[\[Psi][[i]],\[Sigma]x]],{i,1,steps}];
trajY=Table[Re[QAvg[\[Psi][[i]],\[Sigma]y]],{i,1,steps}];
trajZ=Table[Re[QAvg[\[Psi][[i]],\[Sigma]z]],{i,1,steps}]; 

(*The coordinates for the Qubit path on the block sphere*)
Blochcoordinates=Table[{trajX[[t]], trajY[[t]], trajZ[[t]]},{t, 1, steps}];
(*The points above are interpolated*)
InterpolatedCurve=Interpolation[Table[{{t},Blochcoordinates[[t]]},{t,1,steps}]];
(*And then plotted parametrically*)
Trajectory=ParametricPlot3D[InterpolatedCurve[t],{t,1,steps}];

Show[
extra,(*Any extra drawings enter as a parameter in the function*)
Trajectory, (*The trajectory itself*)
BlochSphere, (*The Bloch Sphere*)
BlochLabel, (*Labels in the bloch sphere*)
Style3D
]
]


(* ::Input::Initialization:: *)
BlochTrajectory::usage="BlochTrajectory[\[Psi]_,extra_,Style3D_]
Given a list with the wave function at different time steps, plots the trajectory of a qubit on the Bloch Sphere";


(* ::Input:: *)
(*(* This function is the animated version of the implementation above*)
(**)
(*> Arguments;*)
(*	- \[Psi]: stochastic evolution of the w.f. in the form {\[Psi]0, \[Psi]1, \[Psi]2, ..., \[Psi]nsteps} where \[Psi]i represents the wave function at the i-th iteration.*)
(*We can get this through PhotoDetection[], QubitHomodyneSSE[], and other similar functions;*)
(*	- extra: any other graphic objects to be drawn, such as parametric plots, arrows and so on;*)
(*	- Style3D: any other style choice, such as frame, color and so on;*)
(* *)
(*> Returns;*)
(*	- An animation of the graphich object cointaining the bloch sphere with the trajectory described by \[Psi] and the additional parameters*)
(**)*)


(* ::Input::Initialization:: *)
(*Animates a stochastic trajectory on the Bloch Sphere*)
BlochAnimation[\[Psi]_,extra_,Style3D_]:=Module[{},

(*Coordinates in the Bloch Sphere*)
trajX=Table[Re[QAvg[\[Psi][[i]],\[Sigma]x]],{i,1,steps}];
trajY=Table[Re[QAvg[\[Psi][[i]],\[Sigma]y]],{i,1,steps}];
trajZ=Table[Re[QAvg[\[Psi][[i]],\[Sigma]z]],{i,1,steps}]; 

Animate[

(*The coordinates for the Qubit path on the block sphere*)
Blochcoordinates=Table[{trajX[[t]], trajY[[t]], trajZ[[t]]},{t, 1, steps}];
(*The points above are interpolated*)
InterpolatedCurve=Interpolation[Table[{{t},Blochcoordinates[[t]]},{t,1,steps}]];
(*And then plotted parametrically*)
Trajectory=ParametricPlot3D[InterpolatedCurve[t],{t,1,i}];

Show[
extra, (*Any extra drawings enter as a parameter in the function*)
Trajectory, (*The trajectory itself*)
BlochSphere, (*The Bloch Sphere*)
Show[Graphics3D[{Red,PointSize[0.025],{Point[{trajX[[i]],trajY[[i]],trajZ[[i]]}]}}]], (*Highlights the current point in the animation*)
BlochLabel, (*Labels in the bloch sphere*)
Style3D
],
{i,2,steps,1},AnimationRunning-> False,AnimationRate->50
]
]


(* ::Input::Initialization:: *)
BlochAnimation::usage="BlochAnimation[\[Psi]_,extra_,Style3D_]
Given a list with the wave function at different time steps, plots and animates the trajectory of a qubit on the Bloch Sphere";
