// Consolidated residue-curve ODE calculations for ternary systems.
// Currently supports two packages: "wilson" (activity-coefficient) and "nrtl".
// Additional packages can be migrated into this file incrementally.
// -----------------------------------------------------------------
import type { CompoundData } from './vle-types';
import { calculatePsat_Pa } from './vle-calculations-unifac'; // Psat helper (shared)
import { R_gas_const_J_molK as R } from './vle-calculations-wilson';
import { calculateUnifacGamma, type UnifacParameters } from './vle-calculations-unifac';
import { R_gas_const_J_molK as Rpr, solveCubicEOS, type PrInteractionParams as BinaryPrInteractionParams } from './vle-calculations-pr';
import { R_gas_const_J_molK as Rsrk, solveCubicEOS as solveCubicEOS_SRK, type SrkInteractionParams as BinarySrkInteractionParams } from './vle-calculations-srk';

// ----------- shared types --------------------------------------------------
export type FluidPackageResidue = 'wilson' | 'nrtl' | 'unifac' | 'pr' | 'srk' | 'uniquac';

export interface ResidueCurvePoint { x: number[]; T_K: number; step?: number; }
export type ResidueCurve = ResidueCurvePoint[];

export interface AzeotropeResult { x: number[]; T_K: number; errorNorm?: number; converged: boolean; }

// ===========================================================================
//  W I L S O N  (copied & trimmed from residue-curves-ode-wilson.ts)
// ===========================================================================
const MIN_X = 1e-9;
const BUBBLE_T_TOL_K = 0.01;

interface TernaryWilsonParams {
  a01_J_mol: number; a10_J_mol: number;
  a02_J_mol: number; a20_J_mol: number;
  a12_J_mol: number; a21_J_mol: number;
}

function normX(x:number[]):number[]{
  const s = x.reduce((a,b)=>a+b,0);
  const mx = 1-MIN_X*(x.length-1);
  return x.map(v=>Math.max(MIN_X,Math.min(mx,v/s)));
}
function lambdasWilson(V:number[],T:number,p:TernaryWilsonParams){
  const RT = R*T;
  return {
    L01:(V[1]/V[0])*Math.exp(-p.a01_J_mol/RT), L10:(V[0]/V[1])*Math.exp(-p.a10_J_mol/RT),
    L02:(V[2]/V[0])*Math.exp(-p.a02_J_mol/RT), L20:(V[0]/V[2])*Math.exp(-p.a20_J_mol/RT),
    L12:(V[2]/V[1])*Math.exp(-p.a12_J_mol/RT), L21:(V[1]/V[2])*Math.exp(-p.a21_J_mol/RT),
  };
}
function gammaWilson(x:number[],T:number,comps:CompoundData[],p:TernaryWilsonParams):number[]{
  const V=comps.map(c=>c.wilsonParams!.V_L_m3mol);
  const L=lambdasWilson(V,T,p);const [x0,x1,x2]=x;
  const S0=x0+x1*L.L01+x2*L.L02, S1=x0*L.L10+x1+x2*L.L12, S2=x0*L.L20+x1*L.L21+x2;
  const ln0=1-Math.log(S0)-((x0)/S0+(x1*L.L10)/S1+(x2*L.L20)/S2);
  const ln1=1-Math.log(S1)-((x0*L.L01)/S0+(x1)/S1+(x2*L.L21)/S2);
  const ln2=1-Math.log(S2)-((x0*L.L02)/S0+(x1*L.L12)/S1+(x2)/S2);
  return [Math.exp(ln0),Math.exp(ln1),Math.exp(ln2)];
}
function bubbleTWilson(x:number[],P:number,comps:CompoundData[],p:TernaryWilsonParams,Ti:number){
  let T=Ti;
  for(let i=0;i<50;i++){
    const Ps=comps.map(c=>calculatePsat_Pa(c.antoine!,T));
    const g=gammaWilson(x,T,comps,p); if(!g)return null;
    const Pcalc = x[0]*g[0]*Ps[0]+x[1]*g[1]*Ps[1]+x[2]*g[2]*Ps[2];
    const err=Pcalc-P; if(Math.abs(err/P)<1e-5) return {T,g,Ps};
    T-=err/P*(T*0.05);
    if(T<=0) T=Ti;
  }
  return null;
}
function odeWilson(x:number[],P:number,comps:CompoundData[],p:TernaryWilsonParams,Ti:number){
  const bub=bubbleTWilson(normX(x),P,comps,p,Ti); if(!bub) return null;
  const {T,g,Ps}=bub; const y=[0,1,2].map(i=>x[i]*g[i]*Ps[i]/P); const d=[0,1,2].map(i=>x[i]-y[i]);
  return {dxdxi:d,T};
}
async function simulateODE_Wilson(initX:number[],P:number,comps:CompoundData[],p:TernaryWilsonParams,step:number,maxSteps:number,Ti:number):Promise<ResidueCurve|null>{
  async function runDir(forward:boolean){
    const curve:ResidueCurve=[]; let x=[...initX]; let T=Ti;
    for(let i=0;i<maxSteps;i++){
      const r=odeWilson(x,P,comps,p,T); if(!r) break; T=r.T; curve.push({x:normX(x),T_K:T,step:(forward?1:-1)*i*step});
      x=[0,1,2].map(j=>x[j]+(forward?1:-1)*r.dxdxi[j]*step);
      x=normX(x);
      if(x.some(v=>v<MIN_X||v>1-MIN_X)) break;
    }
    return curve;
  }
  const fwd=await runDir(true);
  const bwd=await runDir(false);
  return [...bwd.reverse(),{x:normX(initX),T_K:Ti},...fwd];
}
// Wilson does not yet expose a direct azeotrope solver in original file. Placeholder.
function azeotropeWilson():AzeotropeResult|null{return null;}

// ===========================================================================
//  N R T L  (trimmed copy – essential parts only)
// ===========================================================================
const Rnrtl = 8.31446261815324;
interface TernaryNrtlParams { g01_J_mol:number;g10_J_mol:number;alpha01:number; g02_J_mol:number;g20_J_mol:number;alpha02:number; g12_J_mol:number;g21_J_mol:number;alpha12:number; }
function gammaNrtl(x:number[],T:number,params:TernaryNrtlParams){
  const g=[[0,params.g01_J_mol,params.g02_J_mol],[params.g10_J_mol,0,params.g12_J_mol],[params.g20_J_mol,params.g21_J_mol,0]];
  const a=[[0,params.alpha01,params.alpha02],[params.alpha01,0,params.alpha12],[params.alpha02,params.alpha12,0]];
  const tau=[[0,0,0],[0,0,0],[0,0,0]],G=[[0,0,0],[0,0,0],[0,0,0]];
  for(let i=0;i<3;i++)for(let j=0;j<3;j++){if(i===j)continue; tau[i][j]=g[i][j]/(Rnrtl*T); G[i][j]=Math.exp(-a[i][j]*tau[i][j]);}
  const ln=[0,0,0];
  for(let i=0;i<3;i++){
    let sum1=0,sum2=0;
    for(let j=0;j<3;j++){sum1+=x[j]*tau[j][i]*G[j][i];sum2+=x[j]*G[j][i];}
    const term1=sum2===0?0:sum1/sum2;
    let term2=0;
    for(let j=0;j<3;j++){
      const denom=x.reduce((acc,k,kk)=>acc+x[kk]*G[kk][j],0);
      const num=x[j]*G[i][j]; const f1=denom===0?0:num/denom;
      const num2=x.reduce((acc,m,mm)=>acc+x[mm]*tau[mm][j]*G[mm][j],0);
      const f2=denom===0?0:num2/denom;
      term2+=f1*(tau[i][j]-f2);
    }
    ln[i]=term1+term2;
  }
  return ln.map(Math.exp);
}
function bubbleTNrtl(x:number[],P:number,comps:CompoundData[],params:TernaryNrtlParams,Ti:number){
  let T=Ti;
  for(let i=0;i<50;i++){
    const Ps=comps.map(c=>calculatePsat_Pa(c.antoine!,T));
    const g=gammaNrtl(x,T,params);
    const Pcalc=x[0]*g[0]*Ps[0]+x[1]*g[1]*Ps[1]+x[2]*g[2]*Ps[2];
    const err=Pcalc-P; if(Math.abs(err/P)<1e-5) return {T,g,Ps};
    T-=err/P*(T*0.05); if(T<=0) T=Ti;
  }
  return null;
}
function odeNrtl(x:number[],P:number,comps:CompoundData[],p:TernaryNrtlParams,Ti:number){
  const bub=bubbleTNrtl(normX(x),P,comps,p,Ti); if(!bub) return null;
  const {T,g,Ps}=bub; const y=[0,1,2].map(i=>x[i]*g[i]*Ps[i]/P); const d=[0,1,2].map(i=>x[i]-y[i]);
  return {dxdxi:d,T};
}
async function simulateODE_Nrtl(initX:number[],P:number,comps:CompoundData[],p:TernaryNrtlParams,step:number,maxSteps:number,Ti:number):Promise<ResidueCurve|null>{
  async function runDir(forward:boolean){
    const curve:ResidueCurve=[]; let x=[...initX]; let T=Ti;
    for(let i=0;i<maxSteps;i++){
      const r=odeNrtl(x,P,comps,p,T); if(!r) break; T=r.T; curve.push({x:normX(x),T_K:T,step:(forward?1:-1)*i*step});
      x=[0,1,2].map(j=>x[j]+(forward?1:-1)*r.dxdxi[j]*step);
      x=normX(x);
      if(x.some(v=>v<MIN_X||v>1-MIN_X)) break;
    }
    return curve;
  }
  const fwd=await runDir(true);
  const bwd=await runDir(false);
  return [...bwd.reverse(),{x:normX(initX),T_K:Ti},...fwd];
}
function azeotropeNrtl():AzeotropeResult|null{return null;}

// ===========================================================================
//  U N I F I E D   A P I  ---------------------------------------------------
// ===========================================================================
export async function simulateResidueCurveODE(
  packageName:FluidPackageResidue,
  initial_x:number[], P_system_Pa:number, comps:CompoundData[],
  pkgParams:any, // TernaryWilsonParams or TernaryNrtlParams
  d_xi_step:number, max_steps:number, initialTemp_K:number
):Promise<ResidueCurve|null>{
  switch(packageName){
    case 'wilson': return simulateODE_Wilson(initial_x,P_system_Pa,comps,pkgParams as TernaryWilsonParams,d_xi_step,max_steps,initialTemp_K);
    case 'nrtl':   return simulateODE_Nrtl  (initial_x,P_system_Pa,comps,pkgParams as TernaryNrtlParams ,d_xi_step,max_steps,initialTemp_K);
    case 'unifac': return simulateODE_Unifac(initial_x,P_system_Pa,comps,pkgParams as UnifacParameters ,d_xi_step,max_steps,initialTemp_K);
    case 'pr':     return simulateODE_Pr    (initial_x,P_system_Pa,comps,pkgParams as TernaryPrParams ,d_xi_step,max_steps,initialTemp_K);
    case 'srk':    return simulateODE_Srk   (initial_x,P_system_Pa,comps,pkgParams as TernarySrkParams,d_xi_step,max_steps,initialTemp_K);
    case 'uniquac':return simulateODE_Uniquac(initial_x,P_system_Pa,comps,pkgParams as TernaryUniquacParams,d_xi_step,max_steps,initialTemp_K);
    default: throw new Error(`Package ${packageName} not yet implemented in consolidated ODE`);
  }
}

export function findAzeotrope(
  packageName:FluidPackageResidue,
  P_system_Pa:number, comps:CompoundData[], pkgParams:any,
  guess_x:number[], guess_T:number
):AzeotropeResult|null{
  switch(packageName){
    case 'wilson': return azeotropeWilson();
    case 'nrtl':   return azeotropeNrtl();
    case 'srk':    return azeotropeSrk(P_system_Pa,comps,pkgParams as TernarySrkParams,guess_x,guess_T);
    case 'uniquac':return azeotropeUniquac(P_system_Pa,comps,pkgParams as TernaryUniquacParams,guess_x,guess_T);
    default: return null;
  }
}

// =========================================================================
//  U N I F A C
// =========================================================================
const MIN_X_UNIFAC=1e-7;
function normX_unifac(x:number[]):number[]{const s=x.reduce((a,b)=>a+b,0)||1;const mx=1-MIN_X_UNIFAC*(x.length-1);return x.map(v=>Math.max(MIN_X_UNIFAC,Math.min(mx,v/s)));}
function bubbleTUnifac(x:number[],P:number,comps:CompoundData[],params:UnifacParameters,Ti:number){let T=Ti;for(let i=0;i<35;i++){const Ps=comps.map(c=>calculatePsat_Pa(c.antoine!,T));if(Ps.some(p=>!isFinite(p)||p<=0))return null;const g=calculateUnifacGamma(comps,x,T,params) as number[];if(!g||g.some(v=>!isFinite(v)||v<=0))return null;const Pcalc=x[0]*g[0]*Ps[0]+x[1]*g[1]*Ps[1]+x[2]*g[2]*Ps[2];const err=Pcalc-P;if(Math.abs(err/P)<1e-4) return {T,g,Ps};T-=err/P*(T*0.1);if(T<100)T=100;}return null;}
function odeUnifac(x:number[],P:number,comps:CompoundData[],params:UnifacParameters,Ti:number){const bub=bubbleTUnifac(normX_unifac(x),P,comps,params,Ti);if(!bub)return null;const {T,g,Ps}=bub;const y=[0,1,2].map(i=>x[i]*g[i]*Ps[i]/P);const d=[0,1,2].map(i=>x[i]-y[i]);return {dxdxi:d,T};}
async function simulateODE_Unifac(initX:number[],P:number,comps:CompoundData[],params:UnifacParameters,step:number,maxSteps:number,Ti:number):Promise<ResidueCurve|null>{
  async function runDir(forward:boolean){const curve:ResidueCurve=[];let x=[...initX];let T=Ti;for(let i=0;i<maxSteps;i++){const r=odeUnifac(x,P,comps,params,T);if(!r)break;T=r.T;curve.push({x:normX_unifac(x),T_K:T,step:(forward?1:-1)*i*step});x=[0,1,2].map(j=>x[j]+(forward?1:-1)*r.dxdxi[j]*step);x=normX_unifac(x);if(x.some(v=>v<MIN_X_UNIFAC||v>1-MIN_X_UNIFAC))break;}return curve;}
  const fwd=await runDir(true);const bwd=await runDir(false);return [...bwd.reverse(),{x:normX_unifac(initX),T_K:Ti},...fwd];}

// =========================================================================
//  P R  (Peng–Robinson)  -- trimmed
// =========================================================================
interface TernaryPrParams {k01:number;k10?:number;k02:number;k20?:number;k12:number;k21?:number;}
const MIN_X_PR=1e-9;
function normX_pr(x:number[]):number[]{const s=x.reduce((a,b)=>a+b,0)||1;const mx=1-MIN_X_PR*(x.length-1);return x.map(v=>Math.max(MIN_X_PR,Math.min(mx,v/s)));}
function getKijPr(i:number,j:number,p:TernaryPrParams):number{const key=`${Math.min(i,j)}${Math.max(i,j)}`;if(i===j)return 0;switch(key){case'01':return p.k01;case'02':return p.k02;case'12':return p.k12;}if(i===1&&j===0&&p.k10!==undefined)return p.k10;if(i===2&&j===0&&p.k20!==undefined)return p.k20;if(i===2&&j===1&&p.k21!==undefined)return p.k21;return 0;}
function calcPurePR(T:number,pr:{Tc_K:number;Pc_Pa:number;omega:number}){const Tr=T/pr.Tc_K;const m=0.37464+1.54226*pr.omega-0.26992*pr.omega*pr.omega;const alpha=(1+m*(1-Math.sqrt(Tr)))**2;const ac=0.45724*Rpr*Rpr*pr.Tc_K*pr.Tc_K/pr.Pc_Pa;const a=ac*alpha;const b=0.07780*Rpr*pr.Tc_K/pr.Pc_Pa;return{a,b};}
function calcMixturePR(x:number[],T:number,comps:CompoundData[],p:TernaryPrParams){const a_i=comps.map(c=>calcPurePR(T,c.prParams!).a);const b_i=comps.map(c=>calcPurePR(T,c.prParams!).b);let a_mix=0;for(let i=0;i<3;i++){for(let j=0;j<3;j++){a_mix+=x[i]*x[j]*Math.sqrt(a_i[i]*a_i[j])*(1-getKijPr(i,j,p));}}const b_mix=x.reduce((sum,xi,idx)=>sum+xi*b_i[idx],0);return{a_mix,b_mix,a_i,b_i};}
function prFugacity(x:number[],T:number,P:number,Z:number,a_mix:number,b_mix:number,a_i:number[],b_i:number[],p:TernaryPrParams){const A=a_mix*P/(Rpr*Rpr*T*T);const B=b_mix*P/(Rpr*T);const phi=[] as number[];for(let k=0;k<3;k++){let sum=0;for(let i=0;i<3;i++){sum+=x[i]*Math.sqrt(a_i[k]*a_i[i])*(1-getKijPr(k,i,p));}const term1=b_i[k]/b_mix*(Z-1);const term2=-Math.log(Z-B);const term3=A/(2*Math.sqrt(2)*B)*(sum/a_mix-b_i[k]/b_mix)*Math.log((Z+(1+Math.sqrt(2))*B)/(Z+(1-Math.sqrt(2))*B));phi.push(Math.exp(term1+term2-term3));}return phi;}
function solvePRBubbleT(x:number[],P:number,comps:CompoundData[],p:TernaryPrParams,Ti:number){let T=Ti;for(let iter=0;iter<40;iter++){const mix=calcMixturePR(x,T,comps,p);const coeff = {A:mix.a_mix*P/(Rpr*Rpr*T*T),B:mix.b_mix*P/(Rpr*T)};const roots = solveCubicEOS(1,-(1-coeff.B),coeff.A-3*coeff.B*coeff.B-2*coeff.B,-(coeff.A*coeff.B-coeff.B*coeff.B-coeff.B*coeff.B*coeff.B)) as number[];
  if(!roots || roots.length === 0) return null;
  const positiveRoots = roots.filter(r=>r>0);
  if(positiveRoots.length === 0) return null;
  const ZL = Math.min(...positiveRoots);
  const ZV = Math.max(...positiveRoots);
  const phiL=prFugacity(x,T,P,ZL,mix.a_mix,mix.b_mix,mix.a_i,mix.b_i,p);const K=phiL.map((phi_val,i)=>{const Psat=calculatePsat_Pa(comps[i].antoine!,T);return phi_val*Psat/P;});const sumKx=K.reduce((acc,k,i)=>acc+k*x[i],0);const err=sumKx-1;if(Math.abs(err)<1e-5)return{T,K};T-=err*2; if(T<100)T=100;}return null;}
function odePr(x:number[],P:number,comps:CompoundData[],p:TernaryPrParams,Ti:number){const bub=solvePRBubbleT(normX_pr(x),P,comps,p,Ti);if(!bub)return null;const {T,K}=bub;const y=K.map((k,i)=>k*x[i]);const d=[0,1,2].map(i=>x[i]-y[i]);return{dxdxi:d,T};}
async function simulateODE_Pr(initX:number[],P:number,comps:CompoundData[],p:TernaryPrParams,step:number,maxSteps:number,Ti:number):Promise<ResidueCurve|null>{
  async function runDir(fwd:boolean){const curve:ResidueCurve=[];let x=[...initX];let T=Ti;for(let i=0;i<maxSteps;i++){const r=odePr(x,P,comps,p,T);if(!r)break;T=r.T;curve.push({x:normX_pr(x),T_K:T,step:(fwd?1:-1)*i*step});x=[0,1,2].map(j=>x[j]+(fwd?1:-1)*r.dxdxi[j]*step);x=normX_pr(x);if(x.some(v=>v<MIN_X_PR||v>1-MIN_X_PR))break;}return curve;}
  const fwd=await runDir(true);const bwd=await runDir(false);return [...bwd.reverse(),{x:normX_pr(initX),T_K:Ti},...fwd];}

// =========================================================================
//  S R K  (Soave–Redlich–Kwong)      – trimmed maths
// =========================================================================
interface TernarySrkParams { k01:number;k10?:number;k02:number;k20?:number;k12:number;k21?:number; }
const MIN_X_SRK = 1e-9;
function normX_srk(x:number[]):number[]{const s=x.reduce((a,b)=>a+b,0)||1;const mx=1-MIN_X_SRK*(x.length-1);return x.map(v=>Math.max(MIN_X_SRK,Math.min(mx,v/s)));}
function getKijSrk(i:number,j:number,p:TernarySrkParams):number{const key=`${Math.min(i,j)}${Math.max(i,j)}`;if(i===j)return 0;switch(key){case'01':return p.k01;case'02':return p.k02;case'12':return p.k12;}if(i===1&&j===0&&p.k10!==undefined)return p.k10;if(i===2&&j===0&&p.k20!==undefined)return p.k20;if(i===2&&j===1&&p.k21!==undefined)return p.k21;return 0;}
function calcPureSRK(T:number,srk:{Tc_K:number;Pc_Pa:number;omega:number}){const Tr=T/srk.Tc_K;const m=0.48+1.574*srk.omega-0.176*srk.omega*srk.omega;const alpha=(1+m*(1-Math.sqrt(Tr)))**2;const ac=0.42748*Rsrk*Rsrk*srk.Tc_K*srk.Tc_K/srk.Pc_Pa;const a=ac*alpha;const b=0.08664*Rsrk*srk.Tc_K/srk.Pc_Pa;return{a,b};}
function mixSRK(x:number[],T:number,comps:CompoundData[],p:TernarySrkParams){const a_i=comps.map(c=>calcPureSRK(T,c.srkParams!).a);const b_i=comps.map(c=>calcPureSRK(T,c.srkParams!).b);let a_mix=0;for(let i=0;i<3;i++){for(let j=0;j<3;j++){a_mix+=x[i]*x[j]*Math.sqrt(a_i[i]*a_i[j])*(1-getKijSrk(i,j,p));}}const b_mix=x.reduce((sum,xi,idx)=>sum+xi*b_i[idx],0);return{a_mix,b_mix,a_i,b_i};}
function srkFugacity(x:number[],T:number,P:number,Z:number,a_mix:number,b_mix:number,a_i:number[],b_i:number[],p:TernarySrkParams){const A=a_mix*P/(Rsrk*Rsrk*T*T);const B=b_mix*P/(Rsrk*T);const phi=[] as number[];for(let k=0;k<3;k++){let sum=0;for(let i=0;i<3;i++){sum+=x[i]*Math.sqrt(a_i[k]*a_i[i])*(1-getKijSrk(k,i,p));}const term1=b_i[k]/b_mix*(Z-1);const term2=-Math.log(Z-B);const term3=A/(B)*(b_i[k]/b_mix - (2*sum/a_mix) )*Math.log(1+B/Z);phi.push(Math.exp(term1+term2+term3));}return phi;}
function solveSRKBubbleT(x:number[],P:number,comps:CompoundData[],p:TernarySrkParams,Ti:number){let T=Ti;for(let iter=0;iter<40;iter++){const mix=mixSRK(x,T,comps,p);const coeff={A:mix.a_mix*P/(Rsrk*Rsrk*T*T),B:mix.b_mix*P/(Rsrk*T)};const roots=solveCubicEOS_SRK(1,-(1-coeff.B),coeff.A-3*coeff.B*coeff.B-2*coeff.B,-(coeff.A*coeff.B-coeff.B*coeff.B-coeff.B*coeff.B*coeff.B)) as number[];if(!roots||roots.length===0) return null;const pos=roots.filter(r=>r>0);if(pos.length===0)return null;const ZL=Math.min(...pos);const phiL=srkFugacity(x,T,P,ZL,mix.a_mix,mix.b_mix,mix.a_i,mix.b_i,p);const K=phiL.map((phi,i)=>{const Psat=calculatePsat_Pa(comps[i].antoine!,T);return phi*Psat/P;});const err=K.reduce((acc,k,i)=>acc+k*x[i],0)-1;if(Math.abs(err)<1e-5)return{T,K};T-=err*2;if(T<100)T=100;}return null;}
function odeSrk(x:number[],P:number,comps:CompoundData[],p:TernarySrkParams,Ti:number){const bub=solveSRKBubbleT(normX_srk(x),P,comps,p,Ti);if(!bub)return null;const {T,K}=bub;const y=K.map((k,i)=>k*x[i]);const d=[0,1,2].map(i=>x[i]-y[i]);return{dxdxi:d,T};}
async function simulateODE_Srk(initX:number[],P:number,comps:CompoundData[],p:TernarySrkParams,step:number,maxSteps:number,Ti:number):Promise<ResidueCurve|null>{async function runDir(fwd:boolean){const curve:ResidueCurve=[];let x=[...initX];let T=Ti;for(let i=0;i<maxSteps;i++){const r=odeSrk(x,P,comps,p,T);if(!r)break;T=r.T;curve.push({x:normX_srk(x),T_K:T,step:(fwd?1:-1)*i*step});x=[0,1,2].map(j=>x[j]+(fwd?1:-1)*r.dxdxi[j]*step);x=normX_srk(x);if(x.some(v=>v<MIN_X_SRK||v>1-MIN_X_SRK))break;}return curve;}const fwd=await runDir(true);const bwd=await runDir(false);return [...bwd.reverse(),{x:normX_srk(initX),T_K:Ti},...fwd];}

// =========================================================================
//  U N I Q U A C   (trimmed)
// =========================================================================
interface TernaryUniquacParams {a01_J_mol:number;a10_J_mol:number;a02_J_mol:number;a20_J_mol:number;a12_J_mol:number;a21_J_mol:number;}
const Z_UNIQ=10; const MIN_X_UNI=1e-9;
function normX_uni(x:number[]):number[]{const s=x.reduce((a,b)=>a+b,0)||1;const mx=1-MIN_X_UNI*(x.length-1);return x.map(v=>Math.max(MIN_X_UNI,Math.min(mx,v/s)));}
function tausUni(T:number,p:TernaryUniquacParams){const RT=R*T;return{t01:Math.exp(-p.a01_J_mol/RT),t10:Math.exp(-p.a10_J_mol/RT),t02:Math.exp(-p.a02_J_mol/RT),t20:Math.exp(-p.a20_J_mol/RT),t12:Math.exp(-p.a12_J_mol/RT),t21:Math.exp(-p.a21_J_mol/RT)};}
function gammaUni(x:number[],T:number,comps:CompoundData[],p:TernaryUniquacParams){if(comps.some(c=>!c.uniquacParams))return null;const r=comps.map(c=>c.uniquacParams!.r);const q=comps.map(c=>c.uniquacParams!.q);const tau=tausUni(T,p);const sum_r=x.reduce((s,xi,i)=>s+r[i]*xi,0);const phi=r.map((ri,i)=>ri*x[i]/sum_r);const sum_q=x.reduce((s,xi,i)=>s+q[i]*xi,0);const theta=q.map((qi,i)=>qi*x[i]/sum_q);const l=r.map((ri,i)=>(Z_UNIQ/2)*(ri-q[i])-(ri-1));const sum_xl=x.reduce((s,xi,i)=>s+xi*l[i],0);
const lnGammaC=phi.map((phi_i,i)=>Math.log(phi_i/x[i])+ (Z_UNIQ/2)*q[i]*Math.log(theta[i]/phi_i)+l[i]-(phi_i/x[i])*sum_xl);
// Residual part (only tau values need):
const tauMat=[[1,tau.t01,tau.t02],[tau.t10,1,tau.t12],[tau.t20,tau.t21,1]];
const lnGammaR=[] as number[];
for(let i=0;i<3;i++){let sum1=0;for(let j=0;j<3;j++){sum1+=theta[j]*tauMat[j][i];}let ln= q[i]*(1-Math.log(sum1)-theta.reduce((s,k,kidx)=>s+theta[kidx]*tauMat[i][kidx]/sum1,0));lnGammaR.push(ln);}return lnGammaC.map((c,i)=>Math.exp(c+lnGammaR[i]));}
function bubbleQUni(x:number[],P:number,comps:CompoundData[],params:TernaryUniquacParams,Ti:number){let T=Ti;for(let iter=0;iter<40;iter++){const Ps=comps.map(c=>calculatePsat_Pa(c.antoine!,T));const g=gammaUni(x,T,comps,params);if(!g)return null;const Pcalc=x[0]*g[0]*Ps[0]+x[1]*g[1]*Ps[1]+x[2]*g[2]*Ps[2];const err=Pcalc-P;if(Math.abs(err/P)<1e-5)return{T,g,Ps};T-=err/P*(T*0.05);if(T<100)T=100;}return null;}
function odeUni(x:number[],P:number,comps:CompoundData[],params:TernaryUniquacParams,Ti:number){const bub=bubbleQUni(normX_uni(x),P,comps,params,Ti);if(!bub)return null;const {T,g,Ps}=bub;const y=[0,1,2].map(i=>x[i]*g[i]*Ps[i]/P);const d=[0,1,2].map(i=>x[i]-y[i]);return{dxdxi:d,T};}
async function simulateODE_Uniquac(initX:number[],P:number,comps:CompoundData[],params:TernaryUniquacParams,step:number,maxSteps:number,Ti:number):Promise<ResidueCurve|null>{async function runDir(fwd:boolean){const curve:ResidueCurve=[];let x=[...initX];let T=Ti;for(let i=0;i<maxSteps;i++){const r=odeUni(x,P,comps,params,T);if(!r)break;T=r.T;curve.push({x:normX_uni(x),T_K:T,step:(fwd?1:-1)*i*step});x=[0,1,2].map(j=>x[j]+(fwd?1:-1)*r.dxdxi[j]*step);x=normX_uni(x);if(x.some(v=>v<MIN_X_UNI||v>1-MIN_X_UNI))break;}return curve;}const fwd=await runDir(true);const bwd=await runDir(false);return [...bwd.reverse(),{x:normX_uni(initX),T_K:Ti},...fwd];}

function azeotropeSrk(_P:number,_comps:CompoundData[],_params:TernarySrkParams,_x:number[],_T:number):AzeotropeResult|null{return null;}
function azeotropeUniquac(_P:number,_comps:CompoundData[],_params:TernaryUniquacParams,_x:number[],_T:number):AzeotropeResult|null{return null;}

export type { TernaryWilsonParams, TernaryNrtlParams, TernaryPrParams, TernarySrkParams, TernaryUniquacParams };

// Consolidated public azeotrope helpers (currently stubbed)
export function findAzeotropeNRTL(..._args:any[]):AzeotropeResult|null { return null; }
export function findAzeotropePr(..._args:any[]):AzeotropeResult|null    { return null; }
export function findAzeotropeSrk(..._args:any[]):AzeotropeResult|null   { return null; }
export function findAzeotropeUniquac(..._args:any[]):AzeotropeResult|null { return null; } 