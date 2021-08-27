(* ::Package:: *)

BeginPackage["ela`"]

Formatting::usage="Formatting[abundance_table,sample_metadata,normalization_flag(0 or 1),thresholds]: "
SimpleSA::usage="SimpleSA[cooccurrence_table,convergence_threshold]"
ELA::usage="ELA[h,j]"
SteepestDescent::usage="SteepestDescent[{community_composition,energy},h,j]"
CIntegerDigits::usage="CIntegerDigits[composition_id,max_number_of_species]"
FindTippingPoint::usage="FindTippingPoint[community_composition1,composition2,h,j,maximum]"
Bi::usage="Bi[{community_composition,energy},h,j]"
Energy::usage="Energy[community_composition,h,j]"
GraphObj::usage="GraphObj[output_of_ELA]"
PGraphObj::usage=""
Assort::usage=""
DisconnectivityGraph::usage="DisconnectivityGraph[graph_object,max_number_of_species,plot_label]"
ELPruning::usage="ELAPruning[output_of_ELA,pruning_level]"
FctRescale::usage="FctRescale[environmental_factors_matrix]"
Rescl::usage=""
StochasticApproximation::usage="StochasticApproximation[cooccurrence_table,convergence_threshold]"
OnestepHBS::usage=""
Logprior::usage=""
GradELA::usage="GradELA[{h,j,g},reference_environment,steps,pruning_level,target]"
Renv::usage="Renv[generated_environmental_conditions,scale_coefficients]"
SSDobj::usage="SSDobj[output_of_GradELA]"
SSDiagram::usage="SSDiagram[ssd_object,x-label]"

Unprotect[Formatting];
Unprotect[SimpleSA];
Unprotect[ELA];
Unprotect[SteepestDescent];
Unprotect[CIntegerDigits];
Unprotect[FindTippingPoint];
Unprotect[Bi];
Unprotect[Energy];
Unprotect[GraphObj];
Unprotect[PGraphObj];
Unprotect[Assort];
Unprotect[DisconnectivityGraph];
Unprotect[ELPruning];
Unprotect[FctRescale];
Unprotect[Rescl];
Unprotect[StochasticApproximation];
Unprotect[OnestepHBS];
Unprotect[Logprior];
Unprotect[GradELA];
Unprotect[Renv];
Unprotect[SSDobj];
Unprotect[SSDiagram];

Begin["ela`Private`"]




(*\:30c7\:30fc\:30bf\:6574\:5f62*)
Formatting[abtable_,metadata_,normalizeq_,parameters_]:=Module[{
},
ath=parameters[[1]];(*\:76f8\:5bfe\:983b\:5ea6\:306e\:95be\:5024; \:3053\:306e\:5024\:3088\:308a\:5c11\:306a\:3044\:5834\:54080\:3068\:3059\:308b*)
minocth=parameters[[2]];(*\:751f\:8d77\:7387\:306e\:4e0b\:9650; \:3053\:306e\:5024\:3088\:308a\:751f\:8d77\:7387\:304c\:9ad8\:3044\:7a2e\:306f\:9664\:304f*)
maxocth=parameters[[3]];(*\:751f\:8d77\:7387\:306e\:4e0a\:9650; \:3053\:306e\:5024\:3088\:308a\:751f\:8d77\:7387\:304c\:9ad8\:3044\:7a2e\:306f\:9664\:304f*)
(*\:30e9\:30d9\:30eb\:60c5\:5831\:306e\:53d6\:5f97*)
taxonlabel=Drop[abtable[[1]],1];
absamplelabel=Drop[abtable[[All,1]],1];
mdsamplelabel=Drop[metadata[[All,1]],1];
sharedsamplelabel=Intersection[mdsamplelabel,absamplelabel];
(*\:751f\:8d77\:884c\:5217\:306e\:751f\:6210*)
patab=Flatten[Position[absamplelabel,#]&/@sharedsamplelabel];
patmd=Flatten[Position[mdsamplelabel,#]&/@sharedsamplelabel];
abnum=((Drop[#,1]&/@Drop[abtable,1])/.{"NA"->0})[[patab]];(*\:30b5\:30f3\:30d7\:30eb\:306e\:4e26\:3073\:3092\:7d71\:4e00*)
If[normalizeq==1,abnum=(N[#/(Total[#]/.{0.->1,0->1})])&/@abnum];(*\:30b5\:30f3\:30d7\:30eb\:5185\:3067\:30ce\:30fc\:30de\:30e9\:30a4\:30ba (\:30aa\:30d7\:30b7\:30e7\:30f3)*)
aboc=Map[(If[#>ath,1,0])&,abnum,{2}];(*\:76f8\:5bfe\:983b\:5ea6\:306e\:95be\:5024ath\:30670/1\:5909\:63db*)
occrit=Position[N[Mean[aboc]],_?(#>maxocth||#<minocth&)];(*\:751f\:8d77\:7387\:306e\:7b97\:51fa//\:751f\:8d77\:7387\:306b\:5f93\:3063\:305f\:9664\:53bb\:5bfe\:8c61\:306e\:6c7a\:5b9a*)
(*\:30aa\:30d6\:30b8\:30a7\:30af\:30c8\:751f\:6210*)
ocmatrix=Transpose[Delete[Transpose[aboc],occrit]];(*\:4e0a\:3067\:9078\:629e\:3057\:305f\:5206\:985e\:7fa4\:3092\:9664\:53bb*)
abmatrix=Transpose[Delete[Transpose[abnum],occrit]];(*\:4e0a\:3067\:9078\:629e\:3057\:305f\:5206\:985e\:7fa4\:3092\:9664\:53bb*)
fctmatrix=(Drop[#,1]&/@Drop[metadata,1])[[patmd]];(*\:74b0\:5883\:30d1\:30e9\:30e1\:30fc\:30bf*)
samplelabel=sharedsamplelabel;(*\:30b5\:30f3\:30d7\:30eb\:30e9\:30d9\:30eb*)
taxonlabel=Delete[taxonlabel,occrit];(*\:5206\:985e\:7fa4\:30e9\:30d9\:30eb*)
factorlabel=Drop[metadata[[1]],1];(*\:74b0\:5883\:56e0\:5b50\:30e9\:30d9\:30eb*)
Print[
"Processed "<>ToString[Length[ocmatrix]]<>" samples.
Relative abundance threshold = "<>ToString[ath]<>",
Occurrence threshold (lower) = "<>ToString[minocth]<>",
Occurrence threshold (upper) = "<>ToString[maxocth]<>".
Selected "<>ToString[Length[ocmatrix[[1]]]]<>" out of "<>ToString[Length[abnum[[1]]]]<>" species."
];
{ocmatrix,abmatrix,fctmatrix,samplelabel,taxonlabel,factorlabel}
]

(*\:74b0\:5883\:56e0\:5b50\:3092\:542b\:307e\:306a\:3044\:30d1\:30e9\:30e1\:30fc\:30bf\:63a8\:5b9a*)
SimpleSA[ystr_,qth_]:=Module[{maxit=100000,
nlocation=Length[ystr],nspecies=Length[ystr[[1]]],
ystats=Transpose[ystr].ystr,ysim,lp,
beta,delbeta,delalpha,learningrate0,tt,alpha,logmodel,logmo,
learningrate,momentum,ysimstats,ydif,betagrad,alphagrad,
uod,upds,upd,updmon,upcrit},
ysim=Table[0,{nlocation},{nspecies}];
lp=N[(Total[ystr]+1)/(Length[ystr]+2)];alpha=Log[lp/(1-lp)];
beta=Table[0,{nspecies},{nspecies}];delbeta=0;delalpha=0;
learningrate0=0.1;
tt=0;upd=Infinity;upds={};updmon={};
While[tt<maxit,
tt++;
logmodel=Compile[{{y,_Integer,2}},1/(1+Exp[Transpose[-alpha-Transpose[y.beta]]])];
logmo=logmodel[ysim];
ysim=OnestepHBS[ysim,logmo];
learningrate=learningrate0;
momentum=0.3;
ysimstats=Transpose[ysim].ysim;
ydif=ystats-ysimstats;
betagrad=(ydif + Logprior[beta,0.5])/nlocation*Abs[IdentityMatrix[nspecies]-1];
alphagrad=(Diagonal[ydif]+Logprior[alpha,2])/nlocation;
delbeta=(1-momentum)*betagrad*learningrate+momentum*delbeta;
delalpha=(1-momentum)*alphagrad*learningrate+momentum*delalpha;
beta=beta+delbeta;
alpha=alpha+delalpha;
upd=Max[Abs[Flatten[Join[delalpha,delbeta]]]];
upds=Append[upds,upd];
upds=Take[upds,-Min[Length[upds],1000]];
updmon=Append[updmon,upd];
upcrit=Abs[Mean[(#[[2]]-#[[1]])&/@Partition[upds,2,1]]];
If[tt>100&&upcrit<qth,Break[]]
];If[tt==maxit,
Print["Exceeded maximum iteration length(\!\(\*SuperscriptBox[\(10\), \(5\)]\)). Threshold value may be too low."]];
{alpha,beta,updmon}
]

(*\:30a8\:30cd\:30eb\:30ae\:30fc\:30e9\:30f3\:30c9\:30b9\:30b1\:30fc\:30d7\:89e3\:6790*)
ELA[h_,j_]:=Module[{cEnergy,xinit,xt,n=Length[h],minsets,stablestates,tipps,
tippingenergies,tippingstates,ssenergies},
cEnergy=Compile[ {{x,_Integer,1}},-x.h-x.(x.j)/2];
xinit:={xt=RandomInteger[1,n],cEnergy[xt]};
minsets=ParallelTable[SteepestDescent[xinit,h,j],{20000},DistributedContexts->Automatic];
stablestates=Union[FromDigits[#,2]&/@minsets[[All,1]]];
ssenergies=cEnergy[CIntegerDigits[#,n]]&/@stablestates;
tipps=Partition[ParallelTable[
If[w[[2]]>w[[1]],
FindTippingPoint[stablestates[[w[[1]]]],stablestates[[w[[2]]]],h,j,10000],
{Infinity,Infinity}],{w,Flatten[Table[{x,y},{x,1,Length[stablestates]},{y,1,Length[stablestates]}],1]},
DistributedContexts->Automatic],
Length[stablestates]];
tippingenergies=tipps[[All,All,1]];
tippingstates=tipps[[All,All,2]];
{stablestates,ssenergies,tippingstates,tippingenergies}]

(*\:6700\:6025\:964d\:4e0b\:6cd5*)
SteepestDescent[x_,h_,j_]:=Module[{y,term,ystat,yene,energi,ytj,edif,sign,yenenew,
mine,mp},
y=x;term=0;
While[term!=1,
ystat=y[[1]];yene=y[[2]];
energi=Table[
ytj=ystat[[pp]];(*presence/absence of selected species*)
edif=-h[[pp]]-ystat.j[[pp]];(*energy difference*)
sign=-2*ytj+1;
yenenew=yene+sign*edif;(*energy of candidate state*)
yenenew,{pp,1,Length[ystat]}];
mine=Min[energi];
If[mine<yene,
mp=Position[energi,mine][[1,1]];
y={ReplacePart[ystat,mp->Abs[ystat[[mp]]-1]],mine},
term=1]
];
y]

(*\:6574\:6570\:3092\:7fa4\:96c6\:7d44\:6210(\:4e8c\:5024\:30d9\:30af\:30c8\:30eb)\:306b\:5909\:63db*)
CIntegerDigits[ssid_,n_]:=Module[{id},id=IntegerDigits[ssid,2];
Join[Table[0,{n-Length[id]}],id]]

(*\:5206\:6c34\:5dba\:7279\:5b9a*)
FindTippingPoint[s1_,s2_,h_,j_,tmax_]:=Module[{n=Length[h],ss1,ss2,dif,pd,seq,bde,tt,rpl,
seqnew,samplepath,pathandenergy,bastrack,crossbas,sortedcb,tipping,bdenew,sss1,sss2,sam,pe,bas,
minbde,mintip,cEnergy},
cEnergy=Compile[ {{x,_Integer,1}},-x.h-x.(x.j)/2];
ss1=CIntegerDigits[s1,n];
ss2=CIntegerDigits[s2,n];
sss1=ss1;sss2=ss2;
dif=ss2-ss1;
pd=Flatten[Position[dif,_?(Abs[#]==1&)]];
seq=RandomSample[pd];
If[Length[seq]<=7,
(*If length(seq) is less than 8,
enumerate all paths*)
minbde=Infinity;
Map[(
samplepath=FoldList[ReplacePart[#1,#2->Abs[#1[[#2]]-1]]&,ss1,#];
pathandenergy={cEnergy[#],#}&/@samplepath;
sortedcb=Sort[pathandenergy][[-1]];
bdenew=sortedcb[[1]];
If[bdenew<minbde,
tipping=sortedcb[[2]];
{minbde,mintip}={bdenew,tipping}];
)&,Permutations[seq]],
(*If length(seq) is greater than 7,
run Metropolis-Hastings algorithm to find the tipping point*)
bde=Infinity;
tt=0;minbde=Infinity;
While[tt<tmax,
tt++;
rpl=RandomChoice[Range[1,Length[seq]],2];
seqnew=ReplacePart[seq,{rpl[[1]]->seq[[rpl[[2]]]],rpl[[2]]->seq[[rpl[[1]]]]}];
samplepath=FoldList[ReplacePart[#1,#2->Abs[#1[[#2]]-1]]&,ss1,seqnew];
pathandenergy={cEnergy[#],#}&/@samplepath;
sortedcb=Sort[pathandenergy][[-1]];
bdenew=sortedcb[[1]];
If[bdenew<minbde,
tipping=sortedcb[[2]];
{minbde,mintip}={bdenew,tipping}];
If[RandomReal[]<Min[1,Exp[bde]/Exp[bdenew]],
seq=seqnew;bde=bdenew]
];
];
{minbde,FromDigits[mintip,2]}
]

(*\:6700\:6025\:964d\:4e0b\:6cd5\:306b\:3088\:308b\:30d9\:30a4\:30b7\:30f3\:306e\:7279\:5b9a*)
Bi[xi_,h_,j_]:={FromDigits[#[[1]],2],#[[2]]}&@SteepestDescent[xi,h,j](*<-Bi*)

(*\:30a8\:30cd\:30eb\:30ae\:30fc\:306e\:8a08\:7b97*)
Energy[stats_,h_,j_]:=-stats.h-stats.(stats.j)/2;

(*Disconnectivity graph \:63cf\:753b\:30aa\:30d6\:30b8\:30a7\:30af\:30c8\:751f\:6210*)
GraphObj[eladata_]:=Module[{ssenergies,tippingenergies,rre,energyall,eee,sta,tip,rips,tempmat,
ccc,ese,eps,ori,hee,cee,rcee,xpositions,xx,yy,pp,nodes2xposi,comp,ff,graphinfo,ea},
ssenergies=eladata[[2]];tippingenergies=eladata[[4]];
energyall=(ea=DeleteCases[Union[Join[ssenergies,Flatten[tippingenergies]]],Infinity]);
rre=Union[Join[(#[[2]]->#[[1]])&/@Transpose[{eladata[[1]],eladata[[2]]}],
(#[[2]]->#[[1]])&/@Transpose[{Flatten[eladata[[3]]],Flatten[eladata[[4]]]}]]];
eee=Reverse[Table[
sta=Position[ssenergies,_?(#<=ww&)];
tip=Position[tippingenergies,_?(#<=ww&),{2},Heads->False];
rips=(If[Length[#]==1,Flatten[{#,#}],#])&/@Join[sta,tip];
tempmat=Normal[SparseArray[(#->1)&/@rips]];
ccc=ConnectedComponents[AdjacencyGraph[tempmat+Transpose[tempmat]]];
{ww,DeleteCases[ccc,_?(Length[#]==1&&!MemberQ[sta,#]&)]},
{ww,energyall}]];
ese=#[[-1]]&/@Split[eee,(Sort[Sort/@#1[[2]]]==Sort[Sort/@#2[[2]]]&)];
eps=ese[[All,2]];
ori=eps[[1,1]];
hee=(comp={#}&/@Complement[ori,Flatten[#]];
Join[#,comp])&/@eps;
cee=Drop[FoldList[Assort[#1,#2]&,hee[[1]],hee],1];
rcee=Reverse[cee];
xpositions=Reverse[Drop[FoldList[If[Length[#2]==Length[#1[[1]]],#1,xx=#1;yy=#2;
pp=((Position[yy,#][[1]])&/@Flatten[xx[[1]]])[[All,1]];
{({#[[1]]})&/@yy,(N[Total[#[[All,2]]]/Length[#[[All,2]]]]
)&/@Split[Transpose[{pp,xx[[2]]}],(#1[[1]]==#2[[1]]&)]}
]&,{rcee[[1]],Range[Length[rcee[[1]]]]},rcee][[All,2]],1]];
nodes2xposi=(#[[1]]->#[[2]])&/@Union[Transpose[{Flatten[cee,1],Flatten[xpositions]}]];
ff=Flatten[(Transpose[{#[[2]],Table[#[[1]],{Length[#[[2]]]}],#[[2]]/.nodes2xposi}])&/@ese,1];
graphinfo=(Sort[#,(#1[[2]]<=#2[[2]]&)][[1]])&/@Split[Sort[ff],(Length[#1[[1]]]==Length[Union[Join[#1[[1]],#2[[1]]]]]&)];
(Join[#,{#[[2]]/.rre}])&/@graphinfo
]

(*Disconnectivity graph \:63cf\:753b\:30aa\:30d6\:30b8\:30a7\:30af\:30c8\:751f\:6210: \:4e26\:5217\:7248*)
PGraphObj[eladata_]:=Module[{ssenergies,tippingenergies,rre,energyall,eee,sta,tip,rips,tempmat,
ccc,ese,eps,ori,hee,cee,rcee,xpositions,xx,yy,pp,nodes2xposi,comp,ff,graphinfo,ea},
ssenergies=eladata[[2]];tippingenergies=eladata[[4]];
energyall=(ea=DeleteCases[Union[Join[ssenergies,Flatten[tippingenergies]]],Infinity]);
rre=Union[Join[(#[[2]]->#[[1]])&/@Transpose[{eladata[[1]],eladata[[2]]}],
(#[[2]]->#[[1]])&/@Transpose[{Flatten[eladata[[3]]],Flatten[eladata[[4]]]}]]];
eee=Reverse[ParallelTable[
sta=Position[ssenergies,_?(#<=ww&)];
tip=Position[tippingenergies,_?(#<=ww&),{2},Heads->False];
rips=(If[Length[#]==1,Flatten[{#,#}],#])&/@Join[sta,tip];
tempmat=Normal[SparseArray[(#->1)&/@rips]];
ccc=ConnectedComponents[AdjacencyGraph[tempmat+Transpose[tempmat]]];
{ww,DeleteCases[ccc,_?(Length[#]==1&&!MemberQ[sta,#]&)]},
{ww,energyall},DistributedContexts->Automatic]];
ese=#[[-1]]&/@Split[eee,(Sort[Sort/@#1[[2]]]==Sort[Sort/@#2[[2]]]&)];
eps=ese[[All,2]];
ori=eps[[1,1]];
hee=(comp={#}&/@Complement[ori,Flatten[#]];
Join[#,comp])&/@eps;
cee=Drop[FoldList[Assort[#1,#2]&,hee[[1]],hee],1];
rcee=Reverse[cee];
xpositions=Reverse[Drop[FoldList[If[Length[#2]==Length[#1[[1]]],#1,xx=#1;yy=#2;
pp=((Position[yy,#][[1]])&/@Flatten[xx[[1]]])[[All,1]];
{({#[[1]]})&/@yy,(N[Total[#[[All,2]]]/Length[#[[All,2]]]]
)&/@Split[Transpose[{pp,xx[[2]]}],(#1[[1]]==#2[[1]]&)]}
]&,{rcee[[1]],Range[Length[rcee[[1]]]]},rcee][[All,2]],1]];
nodes2xposi=(#[[1]]->#[[2]])&/@Union[Transpose[{Flatten[cee,1],Flatten[xpositions]}]];
ff=Flatten[(Transpose[{#[[2]],Table[#[[1]],{Length[#[[2]]]}],#[[2]]/.nodes2xposi}])&/@ese,1];
graphinfo=(Sort[#,(#1[[2]]<=#2[[2]]&)][[1]])&/@Split[Sort[ff],
(Length[#1[[1]]]==Length[Union[Join[#1[[1]],#2[[1]]]]]&)];
(Join[#,{#[[2]]/.rre}])&/@graphinfo
]

(*\:968e\:5c64\:7684\:306a\:4e26\:3073\:66ff\:3048*)
Assort[parent_,daughter_]:=Module[{orderkey,aa},
orderkey=(aa=#;{Position[parent,_?(SubsetQ[#,aa]&),1,Heads->False][[1,1]],-Length[aa],aa[[1]]})&/@daughter;
Sort[Transpose[{orderkey,daughter}]][[All,2]]]

(*Disconnectivity graph*)
DisconnectivityGraph[graphobj_,n_,label_]:=Module[{range,jun,jen,hgn,cc,graph,comc,aa,bb,ea},
range={Min[graphobj[[All,2]]],Max[graphobj[[All,2]]]};
jun=Sort[Select[graphobj,(Length[#[[1]]]==1&)],(#1[[3]]<=#2[[3]]&)];
jen=Sort[Select[graphobj,(Length[#[[1]]]>1&)],(#1[[3]]<=#2[[3]]&)];
hgn=Show[{
cc=0;
ListPlot[(cc++;
Labeled[{#[[3]],#[[2]]},
Style[Subscript["C",ToString[#[[4]]]],18,ColorData[97][cc],Bold,FontFamily->"Helvetica"]])&/@jun,PlotStyle->Black],
cc=0;ListPlot[(cc++;
Labeled[{#[[3]],#[[2]]},
Style[Subscript["C",ToString[#[[4]]]],18,Black,Bold,FontFamily->"Helvetica"]])&/@jen,PlotStyle->Black]},PlotRange->All];
graph=Show[If[Length[graphobj]!=1,
Show[ListPlot[(If[Length[ToExpression[#[[1]]]]==1,
Labeled[{#[[3]],#[[2]]},
comc=CIntegerDigits[#[[4]],n];
PieChart[Table[1,{Length[comc]}],ChartStyle->(comc/.{0->Black,1->White}),ImageSize->40]],{#[[3]],#[[2]]}])&/@graphobj,
PlotStyle->{Gray,PointSize[0.001]},GridLines->Automatic],Graphics[(aa=#;
bb=Select[graphobj,(#[[1]]!=aa[[1]]&&SubsetQ[#[[1]],aa[[1]]]&)][[1]];
Line[{aa[[{3,2}]],{aa[[3]],bb[[2]]},bb[[{3,2}]]}]
)&/@Drop[graphobj,-1]],Axes->{False,True},AspectRatio->1,PlotRange->{All,range},
AxesLabel->{None,"Energy"},LabelStyle->Directive[18,Black,FontFamily->"Helvetica"],PlotRangeClipping->False,ImageSize->350],
Show[ListPlot[(If[Length[ToExpression[#[[1]]]]==1,Labeled[{#[[3]],#[[2]]},
comc=CIntegerDigits[(#[[4]]),n];
PieChart[Table[1,{Length[comc]}],ChartStyle->(comc/.{0->Black,1->White}),ImageSize->40]],{#[[3]],#[[2]]}])&/@graphobj,
PlotStyle->{Gray,PointSize[0.001]},GridLines->Automatic],Axes->{False,True},AspectRatio->1,PlotRange->{All,range},
AxesLabel->{None,"Energy"},
LabelStyle->Directive[18,Black,FontFamily->"Helvetica"],PlotRangeClipping->False,ImageSize->350]],hgn,
PlotLabel->label];
graph]

(*\:30a8\:30cd\:30eb\:30ae\:30fc\:30e9\:30f3\:30c9\:30b9\:30b1\:30fc\:30d7\:306e\:679d\:5208\:308a*)
ELPruning[elsummary_,pss_]:=Module[{tss,tssen,tti,ttien,ssrep,pp,te,se,rat,pax,paxmax,mist},
{tss,tssen,tti,ttien}=elsummary;
rat=DeleteCases[Tuples[Range[Length[tss]],2],_?(#[[1]]>=#[[2]]&)];
mist={};pax=(te=Extract[ttien,#];se=tssen[[#]];
{mist=Join[mist,te-se];Min[te-se],Sort[Transpose[{se,#}]][[All,2]]})&/@rat;
paxmax=Max[mist];
ssrep=elsummary[[1]];
While[Min[pax[[All,1]]]<pss*paxmax,
rat=DeleteCases[Tuples[Range[Length[tss]],2],_?(#[[1]]>=#[[2]]&)];
mist={};pax=(te=Extract[ttien,#];se=tssen[[#]];
{mist=Join[mist,te-se];Min[te-se],Sort[Transpose[{se,#}]][[All,2]]})&/@rat;
paxmax=Max[mist];If[Min[pax[[All,1]]]>=pss*paxmax,Break[]];
pp=Sort[pax][[1,2]];
ssrep=ssrep/.(tss[[pp[[2]]]]->tss[[pp[[1]]]]);
tss=Drop[tss,{pp[[2]]}];
tssen=Drop[tssen,{pp[[2]]}];
tti=Drop[#,{pp[[2]]}]&/@Drop[tti,{pp[[2]]}];
ttien=Drop[#,{pp[[2]]}]&/@Drop[ttien,{pp[[2]]}];
];
{{tss,tssen,tti,ttien},
(#[[1]]->#[[2]])&/@Transpose[{elsummary[[1]],ssrep}]}]

(*\:74b0\:5883\:56e0\:5b50\:306e\:30b9\:30b1\:30fc\:30eb\:5909\:63db*)
FctRescale[fm_]:={Transpose[Rescl/@Transpose[fm]],{Min[#],Max[#]}&/@Transpose[fm]}
Rescl[aa_]:=Module[{amin=Min[aa],mm},
N[(mm=aa-amin)/Max[mm]]]

(*\:74b0\:5883\:56e0\:5b50\:3092\:542b\:3080\:5834\:5408\:306e\:30d1\:30e9\:30e1\:30fc\:30bf\:63a8\:5b9a*)
StochasticApproximation[ystr_,env_,qth_]:=Module[{maxit=100000,
nlocation=Length[ystr],nspecies=Length[ystr[[1]]],nenvironment=Length[env[[1]]],
ystats=Transpose[ystr].ystr,yenvstats=Transpose[env].ystr,ysim,lp,learningrate0,learningrate,tt,logmodel,logmo,
momentum,ysimstats,ysimenvstats,ydif,yenvdif,h,g,j,hgrad,jgrad,ggrad,dh,dj,dg,hg,upds,upd,updmon,upcrit},
ysim=Table[0,{nlocation},{nspecies}];g=Table[0,{nenvironment},{nspecies}];
lp=N[(Total[ystr]+1)/(Length[ystr]+2)];h=Log[lp/(1-lp)];
j=Table[0,{nspecies},{nspecies}];dj=0;dh=0;dg=0;
learningrate0=0.1;tt=0;upd=Infinity;upds={};updmon={};
While[tt<maxit,
tt++;
hg=Transpose[Transpose[env.g]+h];
logmodel=Compile[{{y,_Integer,2}},1/(1+Exp[-hg-y.j])];
logmo=logmodel[ysim];
ysim=OnestepHBS[ysim,logmo];
learningrate=learningrate0;
momentum=0.3;
ysimstats=Transpose[ysim].ysim;
ysimenvstats=Transpose[env].ysim;
ydif=ystats-ysimstats;
jgrad=(ydif + Logprior[j,0.5])/nlocation*Abs[IdentityMatrix[nspecies]-1];
hgrad=(Diagonal[ydif]+Logprior[h,2])/nlocation;
yenvdif=yenvstats-ysimenvstats;
ggrad=(yenvdif+Logprior[g,2])/nlocation;
dj=(1-momentum)*jgrad*learningrate+momentum*dj;
dh=(1-momentum)*hgrad*learningrate+momentum*dh;
dg=(1-momentum)*ggrad*learningrate+momentum*dg;
j=j+dj;
h=h+dh;
g=g+dg;
upd=Max[Abs[Flatten[Join[dh,dg,dj]]]];
upds=Append[upds,upd];
upds=Take[upds,-Min[Length[upds],1000]];
updmon=Append[updmon,upd];
upcrit=Abs[Mean[(#[[2]]-#[[1]])&/@Partition[upds,2,1]]];
If[tt>100&&upcrit<qth,Break[]]
];If[tt==maxit,
Print["Exceeded maximum iteration length(\!\(\*SuperscriptBox[\(10\), \(5\)]\)). Threshold value may be too low."]];
{h,j,g,updmon}
]

(*1\:30b9\:30c6\:30c3\:30d7\:306e\:71b1\:6d74\:6cd5*)
OnestepHBS[y_,logmo_]:=Module[{position,ne},MapThread[(position=RandomInteger[{1,Length[y[[1]]]}];
ne=If[RandomReal[]<#2[[position]],1,0];
ReplacePart[#1,position->ne])&,{y,logmo}]]

(*logistic prior*)
Logprior[x_,a_]:=-N[Tanh[(x)/a/2]/a]

(*\:8907\:6570\:306e\:74b0\:5883\:6761\:4ef6\:306b\:6e21\:308b\:30a8\:30cd\:30eb\:30ae\:30fc\:30e9\:30f3\:30c9\:30b9\:30b1\:30fc\:30d7\:5206\:6790*)
GradELA[hjg_,refenv_,steps_,prn_,eid_]:=Module[
{medenv,tgte,emin,emax,de,baseenv,hsum,ela,gela},
de=N[1/(steps-1)];
gela=Transpose[Table[
baseenv=ReplacePart[refenv,eid->ee];
hsum=baseenv.hjg[[3]]+hjg[[1]];
ela=ELPruning[ELA[hsum,hjg[[2]]],prn][[1]];
{baseenv,ela,
PGraphObj[ela]},
{ee,0,1,de}]];
Print["GradELA generated "<>ToString[steps]<>" energy landscapes at reference environment: "<>ToString[refenv]];
gela
]

(*\:74b0\:5883\:56e0\:5b50\:306e\:30b9\:30b1\:30fc\:30eb\:5fa9\:5143*)
Renv[genv_,rsccoeffs_]:=Transpose[(#[[1]]*(#[[2,2]]-#[[2,1]])+#[[2,1]])&/@Transpose[{Transpose[genv],rsccoeffs}]];

(*\:5b89\:5b9a\:72b6\:614b\:30c0\:30a4\:30a2\:30b0\:30e9\:30e0 \:63cf\:753b\:30aa\:30d6\:30b8\:30a7\:30af\:30c8*)
SSDobj[gradela_]:=Module[{esum0,pp,ee,graphobj,stablestates,ssenergies,
tippingstates,tippingenergies,ss,rss,psi,sid,ti,rti,pti,tid,sscomid,ii,
ticomid,esum},
esum0=Map[(
pp=Position[Union/@Transpose[gradela[[1]]],_?(Length[#]!=1&),{1},Heads->False][[1,1]];
ee=#[[1,pp]];
graphobj=#[[3]];
stablestates=#[[2,1]];
ssenergies=#[[2,2]];
tippingstates=#[[2,3]];
tippingenergies=#[[2,4]];
ss=Sort[Select[graphobj,(Length[#[[1]]]==1&)],(#1[[3]]<=#2[[3]]&)];
rss=(psi=Position[ssenergies,#[[2]]][[1]];
sid=Extract[stablestates,psi];
ReplacePart[#,1->sid])&/@ss;
ti=Sort[Select[graphobj,(Length[#[[1]]]>1&)],(#1[[3]]<=#2[[3]]&)];
rti=(pti=Position[tippingenergies,#[[2]]][[1]];
tid=Extract[tippingstates,pti];
ReplacePart[#,1->tid])&/@ti;
{{ee,rss[[All,2]],Flatten[rss[[All,1]]]},{ee,rti[[All,2]],rti[[All,1]]}})&,Transpose[gradela]
];
sscomid=MapIndexed[({#2[[1]],#1[[1,1]]})&,
(Sort/@Split[Sort[Flatten[MapIndexed[(ii=#2;({#,ii[[1]]})&/@#1)&,esum0[[All,1,3]]],1]],(#1[[1]]==#2[[1]]&)])];
ticomid=MapIndexed[({#2[[1]],#1[[1,1]]})&,
(Sort/@Split[Sort[Flatten[MapIndexed[(ii=#2;({#,ii[[1]]})&/@#1)&,esum0[[All,2,3]]],1]],(#1[[1]]==#2[[1]]&)])];
{{sscomid,ticomid},{esum0[[All,1]]/.((#[[2]]->#[[1]])&/@sscomid),
esum0[[All,2]]/.((#[[2]]->#[[1]])&/@ticomid)}}
]

(*\:5b89\:5b9a\:72b6\:614b\:30c0\:30a4\:30a2\:30b0\:30e9\:30e0*)
SSDiagram[ssdobj_,label_]:=Module[{ssdi,xrep,dm,linelem,stp,cc,xrange,yrange,mi,ma},
Show[Table[
ssdi=ssdobj[[2,x]];xrep=(#[[1]]->#[[2]])&/@ssdobj[[1,x]];dm=(#[[2]]-#[[1]])&@ssdi[[All,1]];
xrange=(mi=Min[Flatten[ssdobj[[2,{1,2},All,1]]]];
ma=Max[Flatten[ssdobj[[2,{1,2},All,1]]]];
{mi-0.02*(ma-mi),ma+0.02*(ma-mi)});
yrange=(mi=Min[Flatten[ssdobj[[2,{1,2},All,2]]]];
ma=Max[Flatten[ssdobj[[2,{1,2},All,2]]]];
{mi-0.02*(ma-mi),ma+0.02*(ma-mi)});
linelem=Split[Sort[Flatten[Transpose[{#[[3]],
Table[#[[1]],{Length[#[[2]]]}],#[[2]]}]&/@ssdi,1]],(#1[[1]]==#2[[1]]&&#2[[2]]-#1[[2]]==dm&)];
stp=linelem[[All,-1]];cc=0;
Show[(cc++;
If[Length[#]>1,
Quiet[ListPlot[
If[x==1,Labeled[#[[All,{2,3}]],Style[Subscript[If[x==1,"C","b"],
Style[ToString[#[[1,1]]/.xrep],14]],16,ColorData[1][cc],FontFamily->"Helvetica"],Top],#[[All,{2,3}]]],
PlotStyle->{ColorData[1][cc],If[x==1,Thick,Dashed]},Joined->True]],
ListPlot[If[x==1,Labeled[#[[{2,3}]],Subscript[If[x==1,"C","b"],
Style[ToString[#[[1]]/.xrep],14]],Top]&/@#,#[[All,{2,3}]]],Joined->False,
PlotStyle->ColorData[1][cc]]])&/@linelem,
AspectRatio->1,Frame->True,FrameTicks->True,AxesOrigin->{0,0},
FrameStyle->Directive[20,Black,FontFamily->"Helvetica"],Frame->True,Axes->True,GridLines->Automatic,
FrameLabel->{label,"Energy"},ImageSize->450,PlotRange->{xrange,yrange}],{x,1,2}]]]




End[]

Protect[Formatting];
Protect[SimpleSA];
Protect[ELA];
Protect[SteepestDescent];
Protect[CIntegerDigits];
Protect[FindTippingPoint];
Protect[Bi];
Protect[Energy];
Protect[GraphObj];
Protect[PGraphObj];
Protect[Assort];
Protect[DisconnectivityGraph];
Protect[ELPruning];
Protect[FctRescale];
Protect[Rescl];
Protect[StochasticApproximation];
Protect[OnestepHBS];
Protect[Logprior];
Protect[GradELA];
Protect[Renv];
Protect[SSDobj];
Protect[SSDiagram];




EndPackage[]
