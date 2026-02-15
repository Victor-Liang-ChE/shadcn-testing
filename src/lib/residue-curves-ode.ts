// Consolidated residue-curve ODE calculations for ternary systems.
// Currently supports two packages: "wilson" (activity-coefficient) and "nrtl".
// Additional packages can be migrated into this file incrementally.
// -----------------------------------------------------------------
import type { CompoundData } from './vle-types';

// Module-level verbose flag: set true during azeotrope solving to enable bubble-T logging
let _verboseBubbleT = false;
import {
  R_gas_const_J_molK as R,
  calculatePsat_Pa,
  calculateUnifacGamma,
  solveCubicEOS,
  type UnifacParameters,
} from './vle-calculations';

// ----------- shared types --------------------------------------------------
export type FluidPackageResidue = 'wilson' | 'nrtl' | 'unifac' | 'pr' | 'srk' | 'uniquac';

export interface ResidueCurvePoint { x: number[]; T_K: number; step?: number; }
export type ResidueCurve = ResidueCurvePoint[];

export interface AzeotropeResult { x: number[]; T_K: number; errorNorm?: number; converged: boolean; }

export type ResidueODEPoint = { d:number[]; T_K:number } | null;

// ===========================================================================
//  W I L S O N  (copied & trimmed from residue-curves-ode-wilson.ts)
// ===========================================================================
const MIN_X = 1e-9;

interface TernaryWilsonParams {
  A01: number; B01: number; A10: number; B10: number;
  A02: number; B02: number; A20: number; B20: number;
  A12: number; B12: number; A21: number; B21: number;
}

function normX(x:number[]):number[]{
  const s = x.reduce((a,b)=>a+b,0);
  const mx = 1-MIN_X*(x.length-1);
  return x.map(v=>Math.max(MIN_X,Math.min(mx,v/s)));
}
function lambdasWilson(V:number[],T:number,p:TernaryWilsonParams){
  // HYSYS Wilson Aij are in cal/mol. Use R_cal=1.9872 and swap convention:
  // L01 uses A10/B10 (reverse) to match HYSYS export format.
  const R_cal = 1.9872; // cal·mol⁻¹·K⁻¹
  return {
    L01:(V[1]/V[0])*Math.exp(-(p.A01/(R_cal*T)+p.B01)), L10:(V[0]/V[1])*Math.exp(-(p.A10/(R_cal*T)+p.B10)),
    L02:(V[2]/V[0])*Math.exp(-(p.A02/(R_cal*T)+p.B02)), L20:(V[0]/V[2])*Math.exp(-(p.A20/(R_cal*T)+p.B20)),
    L12:(V[2]/V[1])*Math.exp(-(p.A12/(R_cal*T)+p.B12)), L21:(V[1]/V[2])*Math.exp(-(p.A21/(R_cal*T)+p.B21)),
  };
}
function gammaWilson(x:number[],T:number,comps:CompoundData[],p:TernaryWilsonParams):number[]{
  const V=comps.map(c=>c.wilsonParams!.V_L_m3mol);
  if (V.some(v => v == null || v <= 0)) {
    console.error(`[gammaWilson] FAILED: Component is missing required Wilson parameter 'V_L_m3mol'.`);
    return []; // Return empty array to signal critical failure
  }
  const L=lambdasWilson(V,T,p);const [x0,x1,x2]=x;
  const S0=x0+x1*L.L01+x2*L.L02, S1=x0*L.L10+x1+x2*L.L12, S2=x0*L.L20+x1*L.L21+x2;
  
  // Check for non-positive or very small denominators before log/division
  if (S0 <= 1e-20 || S1 <= 1e-20 || S2 <= 1e-20) {
    console.error(`[gammaWilson] FAILED: Potential division by zero or log of non-positive. S-values: ${S0}, ${S1}, ${S2}`);
    return [];
  }

  const ln0=1-Math.log(S0)-((x0)/S0+(x1*L.L10)/S1+(x2*L.L20)/S2);
  const ln1=1-Math.log(S1)-((x0*L.L01)/S0+(x1)/S1+(x2*L.L21)/S2);
  const ln2=1-Math.log(S2)-((x0*L.L02)/S0+(x1*L.L12)/S1+(x2)/S2);
  
  const gammas = [Math.exp(ln0),Math.exp(ln1),Math.exp(ln2)];

  if (gammas.some(g => !isFinite(g))) {
    console.error(`[gammaWilson] FAILED: Non-finite gamma value produced. Inputs: x=${JSON.stringify(x)}, T=${T}`);
    return [];
  }
  return gammas;
}
function bubbleTWilson(x:number[],P:number,comps:CompoundData[],p:TernaryWilsonParams,Ti:number){
  let T=Ti;
  let Tprev = T;
  let errprev = 0;
  // console.log(`[bubbleTWilson] START: x=${JSON.stringify(x)}, P=${P}, Ti=${Ti}`);
  for(let i=0;i<80;i++){
    const Ps=comps.map(c=>calculatePsat_Pa(c.antoine!,T));
    const g=gammaWilson(x,T,comps,p); if(g.length === 0) {
        console.error(`[bubbleTWilson] iter ${i}: gammaWilson returned failure`);
        return null;
    }
    const Pcalc = x[0]*g[0]*Ps[0]+x[1]*g[1]*Ps[1]+x[2]*g[2]*Ps[2];
    const err=Pcalc-P;

    if(Math.abs(err/P)<1e-5) {
        if(_verboseBubbleT) console.log(`[bubbleTWilson] CONVERGED: T=${T.toFixed(2)}, gamma=[${g.map(v=>v.toFixed(4))}], Ps=[${Ps.map(v=>v.toFixed(0))}]`);
        return {T,g,Ps};
    }

    let dT: number;
    if (i === 0) {
      dT = -Math.sign(err/P) * 5.0;
    } else {
      const denom = err - errprev;
      dT = denom === 0 ? -Math.sign(err/P) * 2.0 : -err * (T - Tprev) / denom;
    }
    
    Tprev = T;
    errprev = err;
    T += Math.max(-20, Math.min(20, dT));

    if(!isFinite(T) || T <= 0) {
        // console.warn(`[bubbleTWilson] iter ${i}: Temperature became invalid (T=${T}). Resetting to Ti=${Ti}.`);
        T = Ti;
    }
  }
  // console.error(`[bubbleTWilson] FAILED to converge after 50 iterations for x=${JSON.stringify(x)}`);
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

// ===========================================================================
//  N R T L  (trimmed copy – essential parts only)
// ===========================================================================
const Rnrtl = 1.9872; // cal·mol⁻¹·K⁻¹ - HYSYS NRTL table uses cal/mol for Aij in this dataset
interface TernaryNrtlParams {
  A01:number;B01:number;A10:number;B10:number;alpha01:number;
  A02:number;B02:number;A20:number;B20:number;alpha02:number;
  A12:number;B12:number;A21:number;B21:number;alpha12:number;
}
function gammaNrtl(x:number[],T:number,params:TernaryNrtlParams){
  // g[i][j] = A{i}{j} – the swap from HYSYS DB conventions is already applied
  // in the fetch function (fetchNrtlParameters), so params are in model convention.
  const g=[[0,params.A01,params.A02],[params.A10,0,params.A12],[params.A20,params.A21,0]];
  const bMat=[[0,params.B01,params.B02],[params.B10,0,params.B12],[params.B20,params.B21,0]];
  const a=[[0,params.alpha01,params.alpha02],[params.alpha01,0,params.alpha12],[params.alpha02,params.alpha12,0]];
  const tau=[[0,0,0],[0,0,0],[0,0,0]],G=[[1,0,0],[0,1,0],[0,0,1]];
  for(let i=0;i<3;i++)for(let j=0;j<3;j++){if(i===j)continue; tau[i][j]=g[i][j]/(Rnrtl*T)+bMat[i][j]; G[i][j]=Math.exp(-a[i][j]*tau[i][j]);}
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
  let Tprev = T;
  let errprev = 0;
  for(let i=0;i<80;i++){
    const Ps=comps.map(c=>calculatePsat_Pa(c.antoine!,T));
    const g=gammaNrtl(x,T,params);
    const Pcalc=x[0]*g[0]*Ps[0]+x[1]*g[1]*Ps[1]+x[2]*g[2]*Ps[2];
    const err=Pcalc-P; if(Math.abs(err/P)<1e-5) {
      if(_verboseBubbleT) console.log(`[bubbleTNrtl] CONVERGED: T=${T.toFixed(2)}, gamma=[${g.map(v=>v.toFixed(4))}], Ps=[${Ps.map(v=>v.toFixed(0))}]`);
      return {T,g,Ps};
    }
    
    let dT: number;
    if (i === 0) {
      dT = -Math.sign(err/P) * 5.0;
    } else {
      const denom = err - errprev;
      dT = denom === 0 ? -Math.sign(err/P) * 2.0 : -err * (T - Tprev) / denom;
    }
    
    Tprev = T;
    errprev = err;
    T += Math.max(-20, Math.min(20, dT));

    if(!isFinite(T) || T <= 0) T = Ti;
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

// =========================================================================
//  U N I F A C
// =========================================================================
const MIN_X_UNIFAC=1e-7;
function normX_unifac(x:number[]):number[]{const s=x.reduce((a,b)=>a+b,0)||1;const mx=1-MIN_X_UNIFAC*(x.length-1);return x.map(v=>Math.max(MIN_X_UNIFAC,Math.min(mx,v/s)));}
function bubbleTUnifac(x:number[],P:number,comps:CompoundData[],params:UnifacParameters,Ti:number){
  let T=Ti;
  let Tprev = T;
  let errprev = 0;
  for(let i=0;i<80;i++){
    const Ps=comps.map(c=>calculatePsat_Pa(c.antoine!,T));
    if(Ps.some(p=>!isFinite(p)||p<=0))return null;
    const g=calculateUnifacGamma(comps,x,T,params) as number[];
    if(!g||g.some(v=>!isFinite(v)||v<=0))return null;
    const Pcalc=x[0]*g[0]*Ps[0]+x[1]*g[1]*Ps[1]+x[2]*g[2]*Ps[2];
    const err=Pcalc-P;if(Math.abs(err/P)<1e-4) {
      if(_verboseBubbleT) console.log(`[bubbleTUnifac] CONVERGED: T=${T.toFixed(2)}, gamma=[${g.map(v=>v.toFixed(4))}], Ps=[${Ps.map(v=>v.toFixed(0))}]`);
      return {T,g,Ps};
    }

    let dT: number;
    if (i === 0) {
      dT = -Math.sign(err/P) * 5.0;
    } else {
      const denom = err - errprev;
      dT = denom === 0 ? -Math.sign(err/P) * 2.0 : -err * (T - Tprev) / denom;
    }
    
    Tprev = T;
    errprev = err;
    T += Math.max(-20, Math.min(20, dT));

    if(!isFinite(T) || T < 100) T = Math.max(100, Ti);
  }
  return null;
}
function odeUnifac(x:number[],P:number,comps:CompoundData[],params:UnifacParameters,Ti:number){const bub=bubbleTUnifac(normX_unifac(x),P,comps,params,Ti);if(!bub)return null;const {T,g,Ps}=bub;const y=[0,1,2].map(i=>x[i]*g[i]*Ps[i]/P);const d=[0,1,2].map(i=>x[i]-y[i]);return{dxdxi:d,T};}
async function simulateODE_Unifac(initX:number[],P:number,comps:CompoundData[],params:UnifacParameters,step:number,maxSteps:number,Ti:number):Promise<ResidueCurve|null>{
  async function runDir(forward:boolean){
    const curve:ResidueCurve=[];
    let x=[...initX];
    let T=Ti;
    for(let i=0;i<maxSteps;i++){
      const r=odeUnifac(x,P,comps,params,T);
      if(!r)break;
      T=r.T;
      curve.push({x:normX_unifac(x),T_K:T,step:(forward?1:-1)*i*step});
      x=[0,1,2].map(j=>x[j]+(forward?1:-1)*r.dxdxi[j]*step);
      x=normX_unifac(x);
      if(x.some(v=>v<MIN_X_UNIFAC||v>1-MIN_X_UNIFAC))break;
    }
    return curve;
  }
  const fwd=await runDir(true);
  const bwd=await runDir(false);
  return [...bwd.reverse(),{x:normX_unifac(initX),T_K:Ti},...fwd];
}

// =========================================================================
//  P R  (Peng–Robinson)  -- trimmed
// =========================================================================
interface TernaryPrParams {k01:number;k10?:number;k02:number;k20?:number;k12:number;k21?:number;}
const MIN_X_PR=1e-9;
function normX_pr(x:number[]):number[]{const s=x.reduce((a,b)=>a+b,0)||1;const mx=1-MIN_X_PR*(x.length-1);return x.map(v=>Math.max(MIN_X_PR,Math.min(mx,v/s)));}
function getKijPr(i:number,j:number,p:TernaryPrParams):number{
    // Handle asymmetric kij values first (explicitly provided opposites)
    if (i === 1 && j === 0 && p.k10 !== undefined) return p.k10;
    if (i === 2 && j === 0 && p.k20 !== undefined) return p.k20;
    if (i === 2 && j === 1 && p.k21 !== undefined) return p.k21;

    const key = `${Math.min(i,j)}${Math.max(i,j)}`;
    if (i === j) return 0;

    switch(key){
        case '01': return p.k01;
        case '02': return p.k02;
        case '12': return p.k12;
    }
    return 0;
}
function calcPurePR(T:number,pr:{Tc_K:number;Pc_Pa:number;omega:number}){const Tr=T/pr.Tc_K;const m=0.37464+1.54226*pr.omega-0.26992*pr.omega*pr.omega;const alpha=(1+m*(1-Math.sqrt(Tr)))**2;const ac=0.45724*R*R*pr.Tc_K*pr.Tc_K/pr.Pc_Pa;const a=ac*alpha;const b=0.07780*R*pr.Tc_K/pr.Pc_Pa;return{a,b};}
function calcMixturePR(x:number[],T:number,comps:CompoundData[],p:TernaryPrParams){const a_i=comps.map(c=>calcPurePR(T,c.prParams!).a);const b_i=comps.map(c=>calcPurePR(T,c.prParams!).b);let a_mix=0;for(let i=0;i<3;i++){for(let j=0;j<3;j++){a_mix+=x[i]*x[j]*Math.sqrt(a_i[i]*a_i[j])*(1-getKijPr(i,j,p));}}const b_mix=x.reduce((sum,xi,idx)=>sum+xi*b_i[idx],0);return{a_mix,b_mix,a_i,b_i};}
function prFugacity(x:number[],T:number,P:number,Z:number,a_mix:number,b_mix:number,a_i:number[],b_i:number[],p:TernaryPrParams){const A=a_mix*P/(R*R*T*T);const B=b_mix*P/(R*T);const phi=[] as number[];for(let k=0;k<3;k++){let sum=0;for(let i=0;i<3;i++){sum+=x[i]*Math.sqrt(a_i[k]*a_i[i])*(1-getKijPr(k,i,p));}const term1=b_i[k]/b_mix*(Z-1);const term2=-Math.log(Z-B);const term3=A/(2*Math.sqrt(2)*B)*((2*sum/a_mix)-b_i[k]/b_mix)*Math.log((Z+(1+Math.sqrt(2))*B)/(Z+(1-Math.sqrt(2))*B));phi.push(Math.exp(term1+term2-term3));}return phi;}
function solvePRBubbleT(x:number[],P:number,comps:CompoundData[],p:TernaryPrParams,Ti:number){
  // Abort if operating essentially super-critical for any component
  if (comps.some(c => P > 0.95 * c.prParams!.Pc_Pa)) return null;

  let T = Ti;
  let Tprev = T;
  let errprev = 0;
  // Initialize K from Raoult's law
  let K = comps.map(c => { const Psat = calculatePsat_Pa(c.antoine!, T); return Psat / P; });
  if(_verboseBubbleT) console.log(`[solvePRBubbleT] START: x=[${x.map(v=>v.toFixed(4))}], P=${P.toFixed(0)}, Ti=${Ti.toFixed(2)}`);

  for(let iter=0;iter<80;iter++){
    // Liquid phase: EOS with composition x → smallest valid Z
    const mixL = calcMixturePR(x, T, comps, p);
    const AL = mixL.a_mix*P/(R*R*T*T), BL = mixL.b_mix*P/(R*T);
    const rootsL = solveCubicEOS(1, -(1-BL), AL-3*BL*BL-2*BL, -(AL*BL-BL*BL-BL*BL*BL)) as number[];
    if(!rootsL || rootsL.length === 0) { console.error(`[solvePRBubbleT] iter ${iter}: No liquid EOS roots`); return null; }
    const validL = rootsL.filter(r => r > BL);
    if(validL.length === 0) { console.error(`[solvePRBubbleT] iter ${iter}: No valid liquid roots (Z>B). roots=${JSON.stringify(rootsL)}, B=${BL}`); return null; }
    const ZL = Math.min(...validL);
    const phiL = prFugacity(x, T, P, ZL, mixL.a_mix, mixL.b_mix, mixL.a_i, mixL.b_i, p);

    // Inner phi-phi loop: converge K and y at current T
    for(let inner=0; inner<20; inner++){
      const sumKx_inner = K.reduce((s,k,i) => s + k*x[i], 0);
      if(sumKx_inner <= 0 || !isFinite(sumKx_inner)) break;
      const y = K.map((k,i) => k*x[i]/sumKx_inner);
      // Vapor phase: EOS with composition y → largest valid Z
      const mixV = calcMixturePR(y, T, comps, p);
      const AV = mixV.a_mix*P/(R*R*T*T), BV = mixV.b_mix*P/(R*T);
      const rootsV = solveCubicEOS(1, -(1-BV), AV-3*BV*BV-2*BV, -(AV*BV-BV*BV-BV*BV*BV)) as number[];
      if(!rootsV || rootsV.length === 0) break;
      const validV = rootsV.filter(r => r > BV);
      if(validV.length === 0) break;
      const ZV = Math.max(...validV);
      const phiV = prFugacity(y, T, P, ZV, mixV.a_mix, mixV.b_mix, mixV.a_i, mixV.b_i, p);
      const Knew = phiL.map((pl,i) => pl / phiV[i]);
      if(Knew.some(k => !isFinite(k) || k <= 0)) break;
      const maxRelChange = Math.max(...K.map((k,i) => Math.abs((Knew[i]-k)/(k||1e-10))));
      K = Knew;
      if(maxRelChange < 1e-8) break;
    }

    const sumKx = K.reduce((s,k,i) => s + k*x[i], 0);
    const err = sumKx - 1;
    if(iter % 10 === 0 || Math.abs(err) < 1e-3) {
      if(_verboseBubbleT) console.log(`[solvePRBubbleT] iter ${iter}: T=${T.toFixed(2)}, err=${err.toExponential(3)}, K=[${K.map(v=>v.toFixed(4))}], ZL=${ZL.toFixed(4)}`);
    }
    if(!isFinite(err)){ console.warn(`[solvePRBubbleT] iter ${iter}: Non-finite error`); return null; }
    if(Math.abs(err) < 1e-5) {
      if(_verboseBubbleT) console.log(`[solvePRBubbleT] CONVERGED: T=${T.toFixed(2)} K, K=[${K.map(v=>v.toFixed(4))}]`);
      return {T, K};
    }

    // Secant step on T
    let dT: number;
    if (iter === 0) {
      dT = -Math.sign(err) * 5;
    } else {
      const denom = err - errprev;
      dT = denom === 0 ? -Math.sign(err) * 5 : -err * (T - Tprev) / denom;
    }
    if (!isFinite(dT)) dT = -Math.sign(err) * 2;
    Tprev = T;
    errprev = err;
    T += Math.max(-25, Math.min(25, dT));
    if(!isFinite(T) || T < 50 || T > 1500) { console.warn(`[solvePRBubbleT] iter ${iter}: T=${T} out of range`); return null; }
  }
  console.error(`[solvePRBubbleT] FAILED to converge for x=[${x.map(v=>v.toFixed(4))}]`);
  return null;
}
function odePr(x:number[],P:number,comps:CompoundData[],p:TernaryPrParams,Ti:number){const bub=solvePRBubbleT(normX_pr(x),P,comps,p,Ti);if(!bub)return null;const {T,K}=bub;const y=K.map((k,i)=>k*x[i]);const d=[0,1,2].map(i=>x[i]-y[i]);return{dxdxi:d,T};}
async function simulateODE_Pr(initX:number[],P:number,comps:CompoundData[],p:TernaryPrParams,step:number,maxSteps:number,Ti:number):Promise<ResidueCurve|null>{
  async function runDir(fwd:boolean){
    const curve:ResidueCurve=[];
    let x=[...initX];
    let T=Ti;
    let current_step = step;
    const min_step = 1e-5;
    const max_step = 0.05;
    const tol = 5e-4; // Increased tolerance for PR
    const safetyFactor = 0.9;

    for(let i=0;i<maxSteps;i++){
      const r_full = odePr(x, P, comps, p, T);
      if (!r_full) break;

      const x_full = x.map((v, j) => v + (fwd ? 1 : -1) * r_full.dxdxi[j] * current_step);

      const x_half = x.map((v, j) => v + (fwd ? 1 : -1) * r_full.dxdxi[j] * (current_step / 2));
      const r_half = odePr(x_half, P, comps, p, r_full.T);
      if (!r_half) break;
      
      const x_new_half = x_half.map((v, j) => v + (fwd ? 1 : -1) * r_half.dxdxi[j] * (current_step / 2));
      
      const error = Math.sqrt(x_new_half.reduce((sum, v, j) => sum + (v - x_full[j])**2, 0));

      if (error < tol) {
        T = r_half.T;
        x = x_new_half;
        curve.push({x: normX_pr(x), T_K:T, step:(fwd?1:-1)*i*current_step});
        if(x.some(v=>v<MIN_X_PR||v>1-MIN_X_PR)) break;
        current_step = Math.min(max_step, current_step * safetyFactor * Math.sqrt(tol / error));
      } else {
        current_step = Math.max(min_step, current_step * safetyFactor * Math.sqrt(tol / error));
        i--; // Redo this step
      }
    }
    return curve;
  }
  const fwd=await runDir(true);
  const bwd=await runDir(false);
  return [...bwd.reverse(),{x:normX_pr(initX),T_K:Ti},...fwd];
}

// =========================================================================
//  S R K  (Soave–Redlich–Kwong)      – trimmed maths
// =========================================================================
interface TernarySrkParams { k01:number;k10?:number;k02:number;k20?:number;k12:number;k21?:number; }
const MIN_X_SRK = 1e-9;
function normX_srk(x:number[]):number[]{const s=x.reduce((a,b)=>a+b,0)||1;const mx=1-MIN_X_SRK*(x.length-1);return x.map(v=>Math.max(MIN_X_SRK,Math.min(mx,v/s)));}
function getKijSrk(i:number,j:number,p:TernarySrkParams):number{
    // Handle asymmetric kij values first
    if (i === 1 && j === 0 && p.k10 !== undefined) return p.k10;
    if (i === 2 && j === 0 && p.k20 !== undefined) return p.k20;
    if (i === 2 && j === 1 && p.k21 !== undefined) return p.k21;

    const key = `${Math.min(i,j)}${Math.max(i,j)}`;
    if (i === j) return 0;

    switch(key){
        case '01': return p.k01;
        case '02': return p.k02;
        case '12': return p.k12;
    }
    return 0;
}
function calcPureSRK(T:number,srk:{Tc_K:number;Pc_Pa:number;omega:number}){const Tr=T/srk.Tc_K;const m=0.48+1.574*srk.omega-0.176*srk.omega*srk.omega;const alpha=(1+m*(1-Math.sqrt(Tr)))**2;const ac=0.42748*R*R*srk.Tc_K*srk.Tc_K/srk.Pc_Pa;const a=ac*alpha;const b=0.08664*R*srk.Tc_K/srk.Pc_Pa;return{a,b};}
function mixSRK(x:number[],T:number,comps:CompoundData[],p:TernarySrkParams){const a_i=comps.map(c=>calcPureSRK(T,c.srkParams!).a);const b_i=comps.map(c=>calcPureSRK(T,c.srkParams!).b);let a_mix=0;for(let i=0;i<3;i++){for(let j=0;j<3;j++){a_mix+=x[i]*x[j]*Math.sqrt(a_i[i]*a_i[j])*(1-getKijSrk(i,j,p));}}const b_mix=x.reduce((sum,xi,idx)=>sum+xi*b_i[idx],0);return{a_mix,b_mix,a_i,b_i};}
function srkFugacity(x:number[],T:number,P:number,Z:number,a_mix:number,b_mix:number,a_i:number[],b_i:number[],p:TernarySrkParams){const A=a_mix*P/(R*R*T*T);const B=b_mix*P/(R*T);const phi=[] as number[];for(let k=0;k<3;k++){let sum=0;for(let i=0;i<3;i++){sum+=x[i]*Math.sqrt(a_i[k]*a_i[i])*(1-getKijSrk(k,i,p));}const term1=b_i[k]/b_mix*(Z-1);const term2=-Math.log(Z-B);const term3=A/(B)*(b_i[k]/b_mix - (2*sum/a_mix) )*Math.log(1+B/Z);phi.push(Math.exp(term1+term2+term3));}return phi;}
function solveSRKBubbleT(x:number[],P:number,comps:CompoundData[],p:TernarySrkParams,Ti:number){
  if (comps.some(c => P > 0.95 * c.srkParams!.Pc_Pa)) return null;

  let T = Ti;
  let Tprev = T;
  let errprev = 0;
  // Initialize K from Raoult's law
  let K = comps.map(c => { const Psat = calculatePsat_Pa(c.antoine!, T); return Psat / P; });
  if(_verboseBubbleT) console.log(`[solveSRKBubbleT] START: x=[${x.map(v=>v.toFixed(4))}], P=${P.toFixed(0)}, Ti=${Ti.toFixed(2)}`);

  for(let iter=0;iter<80;iter++){
    // Liquid phase: EOS with composition x → smallest valid Z
    const mixL = mixSRK(x, T, comps, p);
    const AL = mixL.a_mix*P/(R*R*T*T), BL = mixL.b_mix*P/(R*T);
    const rootsL = solveCubicEOS(1, -1, AL-BL-BL*BL, -AL*BL) as number[];
    if(!rootsL || rootsL.length === 0) { console.error(`[solveSRKBubbleT] iter ${iter}: No liquid EOS roots`); return null; }
    const validL = rootsL.filter(r => r > BL);
    if(validL.length === 0) { console.error(`[solveSRKBubbleT] iter ${iter}: No valid liquid roots. roots=${JSON.stringify(rootsL)}, B=${BL}`); return null; }
    const ZL = Math.min(...validL);
    const phiL = srkFugacity(x, T, P, ZL, mixL.a_mix, mixL.b_mix, mixL.a_i, mixL.b_i, p);

    // Inner phi-phi loop: converge K and y at current T
    for(let inner=0; inner<20; inner++){
      const sumKx_inner = K.reduce((s,k,i) => s + k*x[i], 0);
      if(sumKx_inner <= 0 || !isFinite(sumKx_inner)) break;
      const y = K.map((k,i) => k*x[i]/sumKx_inner);
      // Vapor phase: EOS with composition y → largest valid Z
      const mixV = mixSRK(y, T, comps, p);
      const AV = mixV.a_mix*P/(R*R*T*T), BV = mixV.b_mix*P/(R*T);
      const rootsV = solveCubicEOS(1, -1, AV-BV-BV*BV, -AV*BV) as number[];
      if(!rootsV || rootsV.length === 0) break;
      const validV = rootsV.filter(r => r > BV);
      if(validV.length === 0) break;
      const ZV = Math.max(...validV);
      const phiV = srkFugacity(y, T, P, ZV, mixV.a_mix, mixV.b_mix, mixV.a_i, mixV.b_i, p);
      const Knew = phiL.map((pl,i) => pl / phiV[i]);
      if(Knew.some(k => !isFinite(k) || k <= 0)) break;
      const maxRelChange = Math.max(...K.map((k,i) => Math.abs((Knew[i]-k)/(k||1e-10))));
      K = Knew;
      if(maxRelChange < 1e-8) break;
    }

    const sumKx = K.reduce((s,k,i) => s + k*x[i], 0);
    const err = sumKx - 1;
    if(iter % 10 === 0 || Math.abs(err) < 1e-3) {
      if(_verboseBubbleT) console.log(`[solveSRKBubbleT] iter ${iter}: T=${T.toFixed(2)}, err=${err.toExponential(3)}, K=[${K.map(v=>v.toFixed(4))}], ZL=${ZL.toFixed(4)}`);
    }
    if(!isFinite(err)){ console.warn(`[solveSRKBubbleT] iter ${iter}: Non-finite error`); return null; }
    if(Math.abs(err) < 1e-5) {
      if(_verboseBubbleT) console.log(`[solveSRKBubbleT] CONVERGED: T=${T.toFixed(2)} K, K=[${K.map(v=>v.toFixed(4))}]`);
      return {T, K};
    }

    // Secant step on T
    let dT: number;
    if (iter === 0) {
      dT = -Math.sign(err) * 5;
    } else {
      const denom = err - errprev;
      dT = denom === 0 ? -Math.sign(err) * 5 : -err * (T - Tprev) / denom;
    }
    if (!isFinite(dT)) dT = -Math.sign(err) * 2;
    Tprev = T;
    errprev = err;
    T += Math.max(-25, Math.min(25, dT));
    if(!isFinite(T) || T < 50 || T > 1500) { console.warn(`[solveSRKBubbleT] iter ${iter}: T=${T} out of range`); return null; }
  }
  console.error(`[solveSRKBubbleT] FAILED to converge for x=[${x.map(v=>v.toFixed(4))}]`);
  return null;
}
function odeSrk(x:number[],P:number,comps:CompoundData[],p:TernarySrkParams,Ti:number){const bub=solveSRKBubbleT(normX_srk(x),P,comps,p,Ti);if(!bub)return null;const {T,K}=bub;const y=K.map((k,i)=>k*x[i]);const d=[0,1,2].map(i=>x[i]-y[i]);return{dxdxi:d,T};}
async function simulateODE_Srk(initX:number[],P:number,comps:CompoundData[],p:TernarySrkParams,step:number,maxSteps:number,Ti:number):Promise<ResidueCurve|null>{
  async function runDir(fwd:boolean){
    const curve:ResidueCurve=[];
    let x=[...initX];
    let T=Ti;
    let current_step = step;
    const min_step = 1e-5;
    const max_step = 0.05;
    const tol = 1e-4;
    const safetyFactor = 0.9;

    for(let i=0;i<maxSteps;i++){
      const r_full = odeSrk(x, P, comps, p, T);
      if (!r_full) break;

      const x_full = x.map((v, j) => v + (fwd ? 1 : -1) * r_full.dxdxi[j] * current_step);

      const x_half = x.map((v, j) => v + (fwd ? 1 : -1) * r_full.dxdxi[j] * (current_step / 2));
      const r_half = odeSrk(x_half, P, comps, p, r_full.T);
      if (!r_half) break;

      const x_new_half = x_half.map((v, j) => v + (fwd ? 1 : -1) * r_half.dxdxi[j] * (current_step / 2));
      
      const error = Math.sqrt(x_new_half.reduce((sum, v, j) => sum + (v - x_full[j])**2, 0));

      if (error < tol) {
        T = r_half.T;
        x = x_new_half;
        curve.push({x: normX_srk(x), T_K:T, step:(fwd?1:-1)*i*current_step});
        if(x.some(v=>v<MIN_X_SRK||v>1-MIN_X_SRK)) break;
        current_step = Math.min(max_step, current_step * safetyFactor * Math.sqrt(tol / error));
      } else {
        current_step = Math.max(min_step, current_step * safetyFactor * Math.sqrt(tol / error));
        i--; // Redo this step
      }
    }
    return curve;
  }
  const fwd=await runDir(true);
  const bwd=await runDir(false);
  return [...bwd.reverse(),{x:normX_srk(initX),T_K:Ti},...fwd];
}

// =========================================================================
//  U N I Q U A C   (trimmed)
// =========================================================================
interface TernaryUniquacParams {
  A01:number;B01:number;A10:number;B10:number;
  A02:number;B02:number;A20:number;B20:number;
  A12:number;B12:number;A21:number;B21:number;
}
const Z_UNIQ=10; const MIN_X_UNI=1e-9;
function normX_uni(x:number[]):number[]{const s=x.reduce((a,b)=>a+b,0)||1;const mx=1-MIN_X_UNI*(x.length-1);return x.map(v=>Math.max(MIN_X_UNI,Math.min(mx,v/s)));}
function tausUni(T:number,p:TernaryUniquacParams){
  // HYSYS UNIQUAC database: Aij, Bij in calories → use R_cal = 1.9872 cal·mol⁻¹·K⁻¹
  // Align to standard Aij/Bij convention: t_ij uses Aij/Bij.
  const R_cal = 1.9872; // cal·mol⁻¹·K⁻¹
  return {
    t01: Math.exp(-(p.A01/(R_cal*T) + p.B01)),
    t10: Math.exp(-(p.A10/(R_cal*T) + p.B10)),
    t02: Math.exp(-(p.A02/(R_cal*T) + p.B02)),
    t20: Math.exp(-(p.A20/(R_cal*T) + p.B20)),
    t12: Math.exp(-(p.A12/(R_cal*T) + p.B12)),
    t21: Math.exp(-(p.A21/(R_cal*T) + p.B21)),
  };
}
function gammaUni(x:number[],T:number,comps:CompoundData[],p:TernaryUniquacParams){if(comps.some(c=>!c.uniquacParams))return null;const r=comps.map(c=>c.uniquacParams!.r);const q=comps.map(c=>c.uniquacParams!.q);const tau=tausUni(T,p);const sum_r=x.reduce((s,xi,i)=>s+r[i]*xi,0);const phi=r.map((ri,i)=>ri*x[i]/sum_r);const sum_q=x.reduce((s,xi,i)=>s+q[i]*xi,0);const theta=q.map((qi,i)=>qi*x[i]/sum_q);const l=r.map((ri,i)=>(Z_UNIQ/2)*(ri-q[i])-(ri-1));const sum_xl=x.reduce((s,xi,i)=>s+xi*l[i],0);
const lnGammaC=phi.map((phi_i,i)=>Math.log(phi_i/x[i])+ (Z_UNIQ/2)*q[i]*Math.log(theta[i]/phi_i)+l[i]-(phi_i/x[i])*sum_xl);
// Residual part (only tau values needed):
const tauMat=[[1,tau.t01,tau.t02],[tau.t10,1,tau.t12],[tau.t20,tau.t21,1]];
// Precompute per-j denominators: sumThetaTau_j = Σ_k θ_k τ_{kj}
const sumThetaTau: number[] = [0,1,2].map(j => {
  let s = 0; for (let k = 0; k < 3; k++) s += theta[k] * tauMat[k][j]; return s;
});
const lnGammaR=[] as number[];
for(let i=0;i<3;i++){
  const sum_ji = sumThetaTau[i]; // Σ_j θ_j τ_{ji}
  let lastTerm = 0;
  for (let j = 0; j < 3; j++) { lastTerm += theta[j] * tauMat[i][j] / sumThetaTau[j]; }
  const ln = q[i] * (1 - Math.log(sum_ji) - lastTerm);
  lnGammaR.push(ln);
}
return lnGammaC.map((c,i)=>Math.exp(c+lnGammaR[i]));}
function bubbleQUni(x:number[],P:number,comps:CompoundData[],params:TernaryUniquacParams,Ti:number){
  let T=Ti;
  let Tprev = T;
  let errprev = 0;
  for(let iter=0;iter<80;iter++){
    const Ps=comps.map(c=>calculatePsat_Pa(c.antoine!,T));
    const g=gammaUni(x,T,comps,params);if(!g)return null;
    const Pcalc=x[0]*g[0]*Ps[0]+x[1]*g[1]*Ps[1]+x[2]*g[2]*Ps[2];
    const err=Pcalc-P;if(Math.abs(err/P)<1e-5){
      if(_verboseBubbleT) console.log(`[bubbleQUni] CONVERGED: T=${T.toFixed(2)}, gamma=[${g.map(v=>v.toFixed(4))}], Ps=[${Ps.map(v=>v.toFixed(0))}]`);
      return{T,g,Ps};
    }

    let dT: number;
    if (iter === 0) {
      dT = -Math.sign(err/P) * 5.0;
    } else {
      const denom = err - errprev;
      dT = denom === 0 ? -Math.sign(err/P) * 2.0 : -err * (T - Tprev) / denom;
    }
    
    Tprev = T;
    errprev = err;
    T += Math.max(-20, Math.min(20, dT));

    if(!isFinite(T) || T < 100) T = Math.max(100, Ti);
  }
  return null;
}
function odeUni(x:number[],P:number,comps:CompoundData[],params:TernaryUniquacParams,Ti:number){const bub=bubbleQUni(normX_uni(x),P,comps,params,Ti);if(!bub)return null;const {T,g,Ps}=bub;const y=[0,1,2].map(i=>x[i]*g[i]*Ps[i]/P);const d=[0,1,2].map(i=>x[i]-y[i]);return{dxdxi:d,T};}
async function simulateODE_Uniquac(initX:number[],P:number,comps:CompoundData[],params:TernaryUniquacParams,step:number,maxSteps:number,Ti:number):Promise<ResidueCurve|null>{
  async function runDir(forward:boolean){
    const curve:ResidueCurve=[];
    let x=[...initX];
    let T=Ti;
    for(let i=0;i<maxSteps;i++){
      const r=odeUni(x,P,comps,params,T);
      if(!r)break;
      T=r.T;
      curve.push({x:normX_uni(x),T_K:T,step:(forward?1:-1)*i*step});
      x=[0,1,2].map(j=>x[j]+(forward?1:-1)*r.dxdxi[j]*step);
      x=normX_uni(x);
      if(x.some(v=>v<MIN_X_UNI||v>1-MIN_X_UNI))break;
    }
    return curve;
  }
  const fwd=await runDir(true);
  const bwd=await runDir(false);
  return [...bwd.reverse(),{x:normX_uni(initX),T_K:Ti},...fwd];
}

// =========================================================================
//  G E N E R I C   A Z E O T R O P E   S O L V E R S  (x ≈ y)
// =========================================================================

function _solveAzeotropeActivity(
  x0:number[], P:number, comps:CompoundData[], Ti:number,
  bubbleFunc:(x:number[],P:number,comps:CompoundData[],params:any,Ti:number)=>any,
  params:any,
  normFunc:(x:number[])=>number[],
  maxIter:number=50, tol:number=1e-7
):AzeotropeResult|null {
  console.log(`[_solveAzeotropeActivity] START: x0=[${x0.map(v=>v.toFixed(4))}], Ti=${Ti.toFixed(2)}, P=${P.toFixed(0)}`);

  // --- Helper functions embedded to avoid scope issues ---
  const vecAdd = (v1: number[], v2: number[]) => v1.map((v, i) => v + v2[i]);
  const vecSub = (v1: number[], v2: number[]) => v1.map((v, i) => v - v2[i]);
  const vecDot = (v1: number[], v2: number[]) => v1.reduce((sum, v, i) => sum + v * v2[i], 0);
  const matVecMul = (m: number[][], v: number[]) => m.map(row => vecDot(row, v));
  const matAdd = (m1: number[][], m2: number[][]) => m1.map((row, i) => vecAdd(row, m2[i]));
  const vecMatMul = (v: number[], m: number[][]): number[] => {
    if(!m || m.length === 0 || !m[0] || m[0].length === 0) return [];
    const m_cols = m[0].length;
    const result = new Array(m_cols).fill(0);
    for (let j = 0; j < m_cols; j++) {
        for (let i = 0; i < v.length; i++) {
            result[j] += v[i] * m[i][j];
        }
    }
    return result;
  }
  function invertMatrix(m: number[][]): number[][] | null {
    const n = m.length;
    if (n === 1) return m[0][0] === 0 ? null : [[1 / m[0][0]]];
    if (n === 2) {
      const [a, b] = m[0];
      const [c, d] = m[1];
      const det = a * d - b * c;
      if (Math.abs(det) < 1e-12) return null;
      return [[d / det, -b / det], [-c / det, a / det]];
    }
        return null;
    }

  const n = comps.length;
  if (n <= 1) return null;
  const m = n - 1; // Number of independent variables

  let T = Ti;
  let x_full = [...x0];
  let X = x0.slice(0, m); // Use first m components as independent variables

  // Function to calculate the error vector F = K_i - 1 for i=0..m-1
  const getF = (currentX: number[], currentT: number): { F: number[], T_K: number, K: number[], x: number[] } | null => {
    const x_independent = [...currentX];
    const x_last = 1.0 - x_independent.reduce((a, b) => a + b, 0);

    if (x_independent.some(v => v < -tol) || x_last < -tol) return null; // Allow small numerical noise
    x_independent.forEach((v,i)=> { if(v < 0) x_independent[i]=0; });
    if(x_last < 0) {
      const norm = 1.0 / (1.0-x_last); // renormalize if last comp is negative
      for(let i=0; i<m; i++) x_independent[i] *= norm;
    }

    const x = [...x_independent, 1.0 - x_independent.reduce((a,b)=>a+b,0)];
    
    const bubble = bubbleFunc(x, P, comps, params, currentT);
    if (!bubble) return null;

    const { T: newT, g, Ps } = bubble;
    const K = g.map((gi: number, j: number) => gi * Ps[j] / P);
    const F = new Array(m);
    for (let i = 0; i < m; i++) F[i] = K[i] - 1;
    
    return { F, T_K: newT, K, x };
  };

  // Initial evaluation
  let evalResult = getF(X, T);
  if (!evalResult) return { x: x0, T_K: T, converged: false, errorNorm: Infinity };
  let { F, T_K, K, x: current_x_full } = evalResult;
  T = T_K;
  x_full = current_x_full;

  // Initial inverse Jacobian (B) via finite differences
  const h = 1e-5;
  const J = new Array(m).fill(0).map(() => new Array(m).fill(0));
  for (let j = 0; j < m; j++) {
    const X_pert = [...X];
    X_pert[j] += h;
    const pertResult = getF(X_pert, T);
    if (!pertResult) { // Fallback for singular points
      for (let i = 0; i < m; i++) J[i][i] = 1;
      break;
    }
    const F_pert = pertResult.F;
    for (let i = 0; i < m; i++) J[i][j] = (F_pert[i] - F[i]) / h;
  }

  let B = invertMatrix(J);
  if (!B) { // If singular, start with scaled identity matrix
    B = new Array(m).fill(0).map((_, i) => new Array(m).fill(0).map((__, j) => i === j ? -1 : 0));
  }

  for (let iter = 0; iter < maxIter; iter++) {
    const errorNorm = Math.sqrt(F.reduce((a, b) => a + b * b, 0));
    if (errorNorm < tol) {
      console.log(`[_solveAzeotropeActivity] CONVERGED iter=${iter}: x=[${normFunc(x_full).map(v=>v.toFixed(4))}], T=${T.toFixed(2)}K, err=${errorNorm.toExponential(3)}`);
      return { x: normFunc(x_full), T_K: T, errorNorm, converged: true };
    }
    if (iter % 10 === 0) {
      console.log(`[_solveAzeotropeActivity] iter ${iter}: x=[${x_full.map(v=>v.toFixed(4))}], T=${T.toFixed(2)}K, errNorm=${errorNorm.toExponential(3)}`);
    }

    // Update step: s_k = -B_k * F_k
    const s = matVecMul(B, F).map(v => -v);

    // Line search to keep compositions physical
    let alpha = 1.0;
    let nextX: number[] | undefined;
    while (alpha > 1e-8) {
        const currentNextX = vecAdd(X, s.map(v => v * alpha));
        const next_x_last = 1.0 - currentNextX.reduce((a, b) => a + b, 0);
        if (currentNextX.every(v => v >= 0) && next_x_last >= 0) {
          nextX = currentNextX;
          break;
        }
        alpha /= 2.0;
    }
    if (!nextX) return { x: normFunc(x_full), T_K: T, errorNorm, converged: false };
    
    const nextEval = getF(nextX, T);
    if (!nextEval) return { x: normFunc(x_full), T_K: T, errorNorm, converged: false };

    const { F: nextF, T_K: nextT, x: next_x_full } = nextEval;
    
    const dx = vecSub(nextX, X);
    const dy = vecSub(nextF, F);

    // Update B using Broyden's formula (Sherman-Morrison)
    const B_dy = matVecMul(B, dy);
    const dxT_B = vecMatMul(dx, B);
    const dxT_B_dy = vecDot(dxT_B, dy);
    
    if (Math.abs(dxT_B_dy) > 1e-12) {
      const term_vec = vecSub(dx, B_dy);
      const B_update = new Array(m).fill(0).map(() => new Array(m).fill(0));
      for (let i = 0; i < m; i++) {
        for (let j = 0; j < m; j++) {
          B_update[i][j] = (term_vec[i] * dxT_B[j]) / dxT_B_dy;
        }
      }
      B = matAdd(B, B_update);
    }
    
    X = nextX; T = nextT; F = nextF; x_full = next_x_full;
  }

  const finalEval = getF(X, T);
  if (finalEval) {
    const errorNorm = Math.sqrt(finalEval.F.reduce((a, b) => a + b * b, 0));
    if (errorNorm < tol * 10) {
      return { x: normFunc(finalEval.x), T_K: finalEval.T_K, errorNorm, converged: true };
    }
    return { x: normFunc(finalEval.x), T_K: finalEval.T_K, errorNorm, converged: false };
  }

  return { x: normFunc(x_full), T_K: T, converged: false, errorNorm: Infinity };
}

// Activity models wrappers
function azeotropeWilson(P:number, comps:CompoundData[], p:TernaryWilsonParams, xGuess:number[], Tguess:number):AzeotropeResult|null{
  _verboseBubbleT = true;
  const result = _solveAzeotropeActivity(xGuess,P,comps,Tguess,bubbleTWilson,p,normX);
  _verboseBubbleT = false;
  return result;
}
function azeotropeNrtlSolver(P:number, comps:CompoundData[], p:TernaryNrtlParams, xGuess:number[], Tguess:number):AzeotropeResult|null{
  _verboseBubbleT = true;
  const result = _solveAzeotropeActivity(xGuess,P,comps,Tguess,bubbleTNrtl,p,normX);
  _verboseBubbleT = false;
  return result;
}
function azeotropeUnifac(P:number, comps:CompoundData[], p:UnifacParameters, xGuess:number[], Tguess:number):AzeotropeResult|null{
  _verboseBubbleT = true;
  const result = _solveAzeotropeActivity(xGuess,P,comps,Tguess,bubbleTUnifac,p,normX_unifac);
  _verboseBubbleT = false;
  return result;
}
function azeotropeUniquacSolver(P:number, comps:CompoundData[], p:TernaryUniquacParams, xGuess:number[], Tguess:number):AzeotropeResult|null{
  _verboseBubbleT = true;
  const result = _solveAzeotropeActivity(xGuess,P,comps,Tguess,bubbleQUni,p,normX_uni);
  _verboseBubbleT = false;
  return result;
}
// Cubic EOS wrappers
function _solveAzeotropeCubic(
  x0: number[], P: number, comps: CompoundData[], Ti: number,
  bubbleFunc: (x: number[], P: number, comps: CompoundData[], params: any, Ti: number) => any,
  params: any,
  normFunc: (x: number[]) => number[],
  maxIter: number = 100, tol: number = 1e-7
): AzeotropeResult | null {
  console.log(`[_solveAzeotropeCubic] START: x0=[${x0.map(v=>v.toFixed(4))}], Ti=${Ti.toFixed(2)}, P=${P.toFixed(0)}`);

  // --- Linear Algebra Helpers (Scoped) ---
  const vecAdd = (v1: number[], v2: number[]) => v1.map((v, i) => v + v2[i]);
  const vecSub = (v1: number[], v2: number[]) => v1.map((v, i) => v - v2[i]);
  const vecDot = (v1: number[], v2: number[]) => v1.reduce((sum, v, i) => sum + v * v2[i], 0);
  const matVecMul = (m: number[][], v: number[]) => m.map(row => vecDot(row, v));
  const matAdd = (m1: number[][], m2: number[][]) => m1.map((row, i) => vecAdd(row, m2[i]));
  const vecMatMul = (v: number[], m: number[][]): number[] => {
    const m_cols = m[0].length;
    const result = new Array(m_cols).fill(0);
    for (let j = 0; j < m_cols; j++) {
      for (let i = 0; i < v.length; i++) {
        result[j] += v[i] * m[i][j];
      }
    }
    return result;
  };
  function invertMatrix(m: number[][]): number[][] | null {
    const n = m.length;
    if (n === 1) return m[0][0] === 0 ? null : [[1 / m[0][0]]];
    if (n === 2) {
      const [a, b] = m[0];
      const [c, d] = m[1];
      const det = a * d - b * c;
      if (Math.abs(det) < 1e-12) return null;
      return [[d / det, -b / det], [-c / det, a / det]];
    }
    return null; // Only 2x2 supported for ternary (2 independent vars)
  }

  const n = comps.length;
  if (n < 2) return null;
  const m = n - 1; // Number of independent variables

  let T = Ti;
  let x_full = [...x0];
  let X = x0.slice(0, m); // Independent vars

  // Error function F: Ki - 1 = 0
  const getF = (currentX: number[], currentT: number): { F: number[], T_K: number, K: number[], x: number[] } | null => {
    const x_independent = [...currentX];
    const x_last = 1.0 - x_independent.reduce((a, b) => a + b, 0);

    // Bounds check
    if (x_independent.some(v => v < -1e-5) || x_last < -1e-5) return null;
    
    // Normalize strictly for the bubble calculation
    let x_calc = [...x_independent, x_last];
    x_calc = x_calc.map(v => Math.max(0, v));
    const sum = x_calc.reduce((a, b) => a + b, 0);
    x_calc = x_calc.map(v => v / sum);

    const bubble = bubbleFunc(x_calc, P, comps, params, currentT);
    if (!bubble) return null;

    // EOS bubbleFunc returns {T, K} directly
    const { T: newT, K } = bubble;
    
    // F[i] = K[i] - 1
    const F = new Array(m);
    for (let i = 0; i < m; i++) F[i] = K[i] - 1;

    return { F, T_K: newT, K, x: x_calc };
  };

  // 1. Initial Evaluation
  let evalResult = getF(X, T);
  if (!evalResult) return null; // Bad initial guess
  let { F, T_K, x: current_x_full } = evalResult;
  T = T_K;
  x_full = current_x_full;

  // 2. Compute Initial Jacobian (Finite Difference)
  const h = 1e-5;
  const J = new Array(m).fill(0).map(() => new Array(m).fill(0));
  for (let j = 0; j < m; j++) {
    const X_pert = [...X];
    X_pert[j] += h;
    const pertResult = getF(X_pert, T);
    if (!pertResult) { // Fallback if perturbation fails
      for (let i = 0; i < m; i++) J[i][i] = 1; 
      break; 
    }
    const F_pert = pertResult.F;
    for (let i = 0; i < m; i++) J[i][j] = (F_pert[i] - F[i]) / h;
  }

  let B = invertMatrix(J);
  if (!B) {
    // If singular, fallback to Identity * damping
    B = new Array(m).fill(0).map((_, i) => new Array(m).fill(0).map((__, j) => i === j ? -0.1 : 0));
  }

  // 3. Main Newton-Broyden Loop
  for (let iter = 0; iter < maxIter; iter++) {
    const errorNorm = Math.sqrt(F.reduce((a, b) => a + b * b, 0));
    if (errorNorm < tol) {
      console.log(`[_solveAzeotropeCubic] CONVERGED iter=${iter}: x=[${normFunc(x_full).map(v=>v.toFixed(4))}], T=${T.toFixed(2)}K, err=${errorNorm.toExponential(3)}`);
      return { x: normFunc(x_full), T_K: T, errorNorm, converged: true };
    }
    if (iter % 10 === 0) {
      console.log(`[_solveAzeotropeCubic] iter ${iter}: x=[${x_full.map(v=>v.toFixed(4))}], T=${T.toFixed(2)}K, errNorm=${errorNorm.toExponential(3)}`);
    }

    // Newton step: s = -B * F
    const s = matVecMul(B, F).map(v => -v);

    // Line search (damped update)
    let alpha = 1.0;
    let nextX: number[] | undefined;
    let nextEval: ReturnType<typeof getF> | undefined;

    while (alpha > 0.05) {
      const currentNextX = vecAdd(X, s.map(v => v * alpha));
      // Quick bounds check before calling expensive EOS
      const next_last = 1.0 - currentNextX.reduce((a,b)=>a+b,0);
      if (currentNextX.every(v => v > -0.01) && next_last > -0.01) {
          const res = getF(currentNextX, T); // T updates inside here usually
          if (res) {
             // Check descent (optional, but good for stability)
             const nextErr = Math.sqrt(res.F.reduce((a,b)=>a+b*b,0));
             if (nextErr < errorNorm || alpha < 0.2) { // Accept if error decreases or step is small
                 nextX = currentNextX;
                 nextEval = res;
                 break;
             }
          }
      }
      alpha *= 0.5;
    }

    if (!nextX || !nextEval) {
      // Line search failed -> try simple steepest descent or just terminate
      return { x: normFunc(x_full), T_K: T, errorNorm, converged: false };
    }

    const { F: nextF, T_K: nextT, x: next_x_full } = nextEval;
    
    // Broyden Update
    const dx = vecSub(nextX, X);
    const dy = vecSub(nextF, F);
    const B_dy = matVecMul(B, dy);
    const dxT_B = vecMatMul(dx, B);
    const dxT_B_dy = vecDot(dxT_B, dy);

    if (Math.abs(dxT_B_dy) > 1e-12) {
      const term_vec = vecSub(dx, B_dy);
      // Outer product update
      for (let i = 0; i < m; i++) {
        for (let j = 0; j < m; j++) {
          B[i][j] += (term_vec[i] * dxT_B[j]) / dxT_B_dy;
        }
      }
    }

    X = nextX;
    T = nextT;
    F = nextF;
    x_full = next_x_full;
  }

  return { x: normFunc(x_full), T_K: T, errorNorm: Infinity, converged: false };
}
function azeotropePrSolver(P:number, comps:CompoundData[], p:TernaryPrParams, xGuess:number[], Tguess:number):AzeotropeResult|null{
  _verboseBubbleT = true;
  const result = _solveAzeotropeCubic(xGuess,P,comps,Tguess,solvePRBubbleT,p,normX_pr);
  _verboseBubbleT = false;
  return result;
}
function azeotropeSrkSolver(P:number, comps:CompoundData[], p:TernarySrkParams, xGuess:number[], Tguess:number):AzeotropeResult|null{
  _verboseBubbleT = true;
  const result = _solveAzeotropeCubic(xGuess,P,comps,Tguess,solveSRKBubbleT,p,normX_srk);
  _verboseBubbleT = false;
  return result;
}

// ========================  PUBLIC EXPORTS  ===============================

export function findAzeotropeWilson(P:number, comps:CompoundData[], params:TernaryWilsonParams, xGuess:number[], Tguess:number):AzeotropeResult|null{
  console.log(`[findAzeotropeWilson] xGuess=[${xGuess.map(v=>v.toFixed(4))}], Tguess=${Tguess.toFixed(2)}, P=${P.toFixed(0)}`);
  const result = azeotropeWilson(P,comps,params,xGuess,Tguess);
  console.log(`[findAzeotropeWilson] result: ${result ? `x=[${result.x.map(v=>v.toFixed(4))}], T=${result.T_K.toFixed(2)}K, conv=${result.converged}, err=${result.errorNorm?.toExponential(3)}` : 'null'}`);
  return result;
}
export function findAzeotropeUnifac(P:number, comps:CompoundData[], params:UnifacParameters, xGuess:number[], Tguess:number):AzeotropeResult|null{
  console.log(`[findAzeotropeUnifac] xGuess=[${xGuess.map(v=>v.toFixed(4))}], Tguess=${Tguess.toFixed(2)}, P=${P.toFixed(0)}`);
  const result = azeotropeUnifac(P,comps,params,xGuess,Tguess);
  console.log(`[findAzeotropeUnifac] result: ${result ? `x=[${result.x.map(v=>v.toFixed(4))}], T=${result.T_K.toFixed(2)}K, conv=${result.converged}, err=${result.errorNorm?.toExponential(3)}` : 'null'}`);
  return result;
}
export function findAzeotropeNRTL(P:number, comps:CompoundData[], params:TernaryNrtlParams, xGuess:number[], Tguess:number):AzeotropeResult|null{
  console.log(`[findAzeotropeNRTL] xGuess=[${xGuess.map(v=>v.toFixed(4))}], Tguess=${Tguess.toFixed(2)}, P=${P.toFixed(0)}`);
  const result = azeotropeNrtlSolver(P,comps,params,xGuess,Tguess);
  console.log(`[findAzeotropeNRTL] result: ${result ? `x=[${result.x.map(v=>v.toFixed(4))}], T=${result.T_K.toFixed(2)}K, conv=${result.converged}, err=${result.errorNorm?.toExponential(3)}` : 'null'}`);
  return result;
}
export function findAzeotropePr(P:number, comps:CompoundData[], params:TernaryPrParams, xGuess:number[], Tguess:number):AzeotropeResult|null{
  console.log(`[findAzeotropePr] xGuess=[${xGuess.map(v=>v.toFixed(4))}], Tguess=${Tguess.toFixed(2)}, P=${P.toFixed(0)}, kij={k01:${params.k01},k02:${params.k02},k12:${params.k12}}`);
  const result = azeotropePrSolver(P,comps,params,xGuess,Tguess);
  console.log(`[findAzeotropePr] result: ${result ? `x=[${result.x.map(v=>v.toFixed(4))}], T=${result.T_K.toFixed(2)}K, conv=${result.converged}, err=${result.errorNorm?.toExponential(3)}` : 'null'}`);
  return result;
}
export function findAzeotropeSrk(P:number, comps:CompoundData[], params:TernarySrkParams, xGuess:number[], Tguess:number):AzeotropeResult|null{
  console.log(`[findAzeotropeSrk] xGuess=[${xGuess.map(v=>v.toFixed(4))}], Tguess=${Tguess.toFixed(2)}, P=${P.toFixed(0)}, kij={k01:${params.k01},k02:${params.k02},k12:${params.k12}}`);
  const result = azeotropeSrkSolver(P,comps,params,xGuess,Tguess);
  console.log(`[findAzeotropeSrk] result: ${result ? `x=[${result.x.map(v=>v.toFixed(4))}], T=${result.T_K.toFixed(2)}K, conv=${result.converged}, err=${result.errorNorm?.toExponential(3)}` : 'null'}`);
  return result;
}
export function findAzeotropeUniquac(P:number, comps:CompoundData[], params:TernaryUniquacParams, xGuess:number[], Tguess:number):AzeotropeResult|null{
  console.log(`[findAzeotropeUniquac] xGuess=[${xGuess.map(v=>v.toFixed(4))}], Tguess=${Tguess.toFixed(2)}, P=${P.toFixed(0)}`);
  const result = azeotropeUniquacSolver(P,comps,params,xGuess,Tguess);
  console.log(`[findAzeotropeUniquac] result: ${result ? `x=[${result.x.map(v=>v.toFixed(4))}], T=${result.T_K.toFixed(2)}K, conv=${result.converged}, err=${result.errorNorm?.toExponential(3)}` : 'null'}`);
  return result;
}

// Re-export model parameter types for external use
export type { TernaryWilsonParams, TernaryNrtlParams, TernaryPrParams, TernarySrkParams, TernaryUniquacParams };

// ---------------------------------------------------------------------------
//  S Y S T E M A T I C   B I N A R Y   &   S M A R T   T E R N A R Y   S E A R C H
// ---------------------------------------------------------------------------

function clusterAzeotropes(results: AzeotropeResult[], tolerance: number = 0.08): AzeotropeResult[] {
  const allConverged = results.filter(r => r.converged);
  if (allConverged.length === 0) return [];

  let clusters: AzeotropeResult[][] = allConverged.map(r => [r]);

  let mergedInLastPass = true;
  while (mergedInLastPass) {
    mergedInLastPass = false;
    const nextClusters: AzeotropeResult[][] = [];
    const visited = new Array(clusters.length).fill(false);

    for (let i = 0; i < clusters.length; i++) {
      if (visited[i]) continue;
      
      const currentMergedCluster = [...clusters[i]];
      visited[i] = true;

      for (let j = i + 1; j < clusters.length; j++) {
        if (visited[j]) continue;

        const areAdjacent = currentMergedCluster.some(p1 =>
          clusters[j].some(p2 => {
            const distSq = (p1.x[0] - p2.x[0]) ** 2 + (p1.x[1] - p2.x[1]) ** 2 + (p1.x[2] - p2.x[2]) ** 2;
            return distSq < tolerance * tolerance;
          })
        );

        if (areAdjacent) {
          currentMergedCluster.push(...clusters[j]);
          visited[j] = true;
          mergedInLastPass = true;
        }
      }
      nextClusters.push(currentMergedCluster);
    }
    clusters = nextClusters;
  }
    
  return clusters.map(cluster => {
    if (cluster.length === 1) return cluster[0];

    const avgX = [0, 0, 0];
    let avgT = 0;
    let totalErr = 0;
    for (const member of cluster) {
      avgX[0] += member.x[0];
      avgX[1] += member.x[1];
      avgX[2] += member.x[2];
      avgT += member.T_K;
      totalErr += member.errorNorm || 0;
    }
    const n = cluster.length;
    const normFactor = avgX.reduce((s,v)=>s+v, 0) || 1;
    return {
      x: [avgX[0] / n / normFactor, avgX[1] / n / normFactor, avgX[2] / n / normFactor],
      T_K: avgT / n,
      converged: true,
      errorNorm: totalErr / n,
    };
  });
}

function _solveBinaryAzeotrope(
    i: number, j: number, x_guess: number, T_guess: number,
    evaluateKiKj: (x: number[], T: number, i: number, j: number) => { val: number; T: number } | null,
    buildX: (i: number, j: number, xi: number) => number[],
    maxIter: number = 30, tol: number = 1e-6
): AzeotropeResult | null {
    // ------------------------------------------------------------------
    // 1)  Initialise two points that (ideally) bracket the root.  We
    //     begin with a symmetric ±Δ perturbation from the supplied guess
    //     and expand if needed until a sign-change is detected or the
    //     ends hit the composition limits.
    // ------------------------------------------------------------------
    const expandUntilBracket = (xMid:number, fMid:number):{xl:number,xu:number,fl:number,fu:number}|null => {
        const MAX_EXPANSION=20;
        let delta=0.01; // initial step
        for(let n=0;n<MAX_EXPANSION;n++){
            let xl=Math.max(0,xMid-delta);
            let xu=Math.min(1,xMid+delta);
            const leftEval=evaluateKiKj(buildX(i,j,xl),T_guess,i,j);
            const rightEval=evaluateKiKj(buildX(i,j,xu),T_guess,i,j);
            if(!leftEval||!rightEval) { delta*=1.5; continue; }
            const fl=leftEval.val; const fu=rightEval.val;
            if(fl*fu<=0){
                // We have a bracket
                return {xl,xu,fl,fu};
            }
            // Otherwise expand search region
            delta*=1.5;
            if(xl===0&&xu===1) break;
        }
        return null;
    };

    const midEval=evaluateKiKj(buildX(i,j,x_guess),T_guess,i,j);
    if(!midEval) return null;
    let {xl,xu,fl,fu} = (()=>{
        const maybe=expandUntilBracket(x_guess,midEval.val);
        if(maybe) return maybe;
        // If we fail to bracket, fall back to secant with previous logic
        return {xl:Math.max(0,x_guess-0.05), xu:Math.min(1,x_guess+0.05),
                fl:midEval.val, fu:midEval.val};
    })();

    // ------------------------------------------------------------------
    // 2)  Bisection refinement – guarantees convergence once a bracket is
    //     available.  Perform a limited number of iterations to obtain a
    //     robust composition and temperature seed.
    // ------------------------------------------------------------------
    let xi = x_guess;
    let T = midEval.T;
    if(fl*fu<=0){
        for(let iter=0;iter<20;iter++){
            xi = 0.5*(xl+xu);
            const mid = evaluateKiKj(buildX(i,j,xi),T,i,j);
            if(!mid) break;
            const fm = mid.val; T = mid.T;
            if(Math.abs(fm) < tol) {
                return { x: buildX(i,j,xi), T_K: T, converged:true, errorNorm:Math.abs(fm) };
            }
            if(fl*fm<0){ xu=xi; fu=fm; }
            else { xl=xi; fl=fm; }
        }
    }

    // ------------------------------------------------------------------
    // 3)  Secant polish – start from the refined midpoint and one end of
    //     the bracket for a few quick super-linear iterations.
    // ------------------------------------------------------------------
    let x_prev = xl; // use lower bound as previous point
    let f_prev = fl;
    xi = xi; // current value from bisection
    let f_curr = evaluateKiKj(buildX(i,j,xi),T,i,j)?.val ?? 0;

    for(let iter=0; iter<maxIter; iter++){
        // Check convergence
        if(Math.abs(f_curr) < tol || Math.abs(xu-xl) < 1e-6){
            return { x: buildX(i,j,xi), T_K:T, converged:true, errorNorm:Math.abs(f_curr) };
        }

        // Secant candidate
        const denom = f_curr - f_prev;
        let xi_sec = Math.abs(denom) < 1e-12 ? xi : xi - f_curr * (xi - x_prev) / denom;

        // Ensure secant candidate lies within bracket; otherwise pick midpoint
        if(!(xi_sec > Math.min(xl,xu) && xi_sec < Math.max(xl,xu))) {
            xi_sec = 0.5*(xl+xu);
        }

        // Evaluate at candidate
        const candEval = evaluateKiKj(buildX(i,j,xi_sec),T,i,j);
        if(!candEval) break;
        const f_sec = candEval.val; const T_sec = candEval.T;

        // Update bracket
        if(fl*f_sec <= 0){ xu = xi_sec; fu = f_sec; }
        else { xl = xi_sec; fl = f_sec; }

        // Prepare next iteration (secant variables)
        x_prev = xi; f_prev = f_curr;
        xi = xi_sec; f_curr = f_sec; T = T_sec;
    }

    return null; // Convergence not achieved
}

/**
 * Scan each binary edge of the ternary composition triangle (0–1 mole fraction)
 * and locate points where the relative volatility between the two components
 * crosses unity ( K_i = K_j ).  A sign–change in ln(K_i/K_j) indicates a binary
 * azeotrope.  The routine brackets each crossing on a coarse grid then refines
 * it via bisection.
 *
 * The implementation supports all fluid-package options that already have a
 * bubble-temperature solver in this consolidated file.  The routine re-uses the
 * existing *bubbleT* helpers so no external dependencies are introduced.
 */
export function systematicBinaryAzeotropeSearch(
  pkg: FluidPackageResidue,
  P_system_Pa: number,
  comps: CompoundData[],
  pkgParams: any,
  initialT_K: number = 350,
  dx: number = 0.005 // Finer grid step for initial scan
): AzeotropeResult[] {
  console.log(`[systematicBinarySearch] pkg=${pkg}, P=${P_system_Pa.toFixed(0)}, comps=[${comps.map(c=>c.name)}], dx=${dx}`);
  _verboseBubbleT = true;
  const results: AzeotropeResult[] = [];
  const pairs: [number, number][] = [[0, 1], [0, 2], [1, 2]];
  const eps = 1e-5; // Slightly larger offset to minimize perturbation of Ki calculations

  // Builds a composition vector on a binary edge, with a small amount of the third component
  const buildX = (i: number, j: number, xi: number): number[] => {
    const k = 3 - i - j;
    const x = [0, 0, 0];
    x[i] = xi * (1 - eps);
    x[j] = (1.0 - xi) * (1 - eps);
    x[k] = eps;
    return x;
  };

  // Generic local solver call
  const runLocalSolver = (xGuess: number[], Tguess: number): AzeotropeResult | null => {
    switch (pkg) {
      case 'wilson': return findAzeotropeWilson(P_system_Pa, comps, pkgParams, xGuess, Tguess);
      case 'nrtl': return findAzeotropeNRTL(P_system_Pa, comps, pkgParams, xGuess, Tguess);
      case 'unifac': return findAzeotropeUnifac(P_system_Pa, comps, pkgParams, xGuess, Tguess);
      case 'pr': return findAzeotropePr(P_system_Pa, comps, pkgParams, xGuess, Tguess);
      case 'srk': return findAzeotropeSrk(P_system_Pa, comps, pkgParams, xGuess, Tguess);
      case 'uniquac': return findAzeotropeUniquac(P_system_Pa, comps, pkgParams, xGuess, Tguess);
      default: return null;
    }
  };

  // Evaluates K_i - K_j for a given binary pair
  const evaluateKiKj = (x: number[], Tguess: number, i: number, j: number): { val: number; T: number } | null => {
    let bubbleResult: { T: number, K: number[] } | null = null;
    
    const activityBubbleTFuncs = {
      'wilson': bubbleTWilson, 'nrtl': bubbleTNrtl, 'unifac': bubbleTUnifac, 'uniquac': bubbleQUni
    };
    const eosSolveBubbleTFuncs = { 
      'pr': solvePRBubbleT, 'srk': solveSRKBubbleT 
    };

    if (pkg in activityBubbleTFuncs) {
        const bubbleTFunc = activityBubbleTFuncs[pkg as keyof typeof activityBubbleTFuncs];
        const bubble = bubbleTFunc(x, P_system_Pa, comps, pkgParams, Tguess);
        if (!bubble) return null;
        const K = bubble.g.map((gi: number, idx: number) => gi * bubble.Ps[idx] / P_system_Pa);
        bubbleResult = { T: bubble.T, K };
    } else if (pkg in eosSolveBubbleTFuncs) {
        const solveBubbleTFunc = eosSolveBubbleTFuncs[pkg as keyof typeof eosSolveBubbleTFuncs];
        const bubble = solveBubbleTFunc(x, P_system_Pa, comps, pkgParams, Tguess);
        if (!bubble) return null;
        bubbleResult = { T: bubble.T, K: bubble.K };
    }

    if (bubbleResult) {
        const Ki = bubbleResult.K[i];
        const Kj = bubbleResult.K[j];
        if(Ki<=0||Kj<=0) return null;
        return { val: Math.log(Ki) - Math.log(Kj), T: bubbleResult.T };
    }
    return null;
  };

  for (const [i, j] of pairs) {
    let Tguess = initialT_K;
    const evaluatedPoints: {xi: number, val: number, T: number}[] = [];

    // 1. Grid Evaluation: Evaluate Ki-Kj across the entire edge
    for (let xi_step = 0; xi_step <= 1.0; xi_step += dx) {
        const xi = Math.max(0.001, Math.min(0.999, xi_step)); // Clamp xi to avoid pure components
        const x = buildX(i, j, xi);
        const currentEval = evaluateKiKj(x, Tguess, i, j);
        if (currentEval) {
            evaluatedPoints.push({ xi, val: currentEval.val, T: currentEval.T });
            Tguess = currentEval.T;
        }
    }

    // 2. Bracket Identification and Solving
    for (let k = 0; k < evaluatedPoints.length - 1; k++) {
        const p1 = evaluatedPoints[k];
        const p2 = evaluatedPoints[k+1];

        if (p1.val * p2.val < 0) { // A root is bracketed
            const x_guess = (p1.xi + p2.xi) / 2;
            const T_guess = (p1.T + p2.T) / 2;

            const sol = _solveBinaryAzeotrope(i, j, x_guess, T_guess, evaluateKiKj, buildX);
            
            if (sol) {
                // console.log(`[Debug] Binary search (${i}-${j}): Solver converged. x=[${sol.x.map(v=>v.toFixed(3)).join(', ')}]`);
                results.push(sol);
            }
        }
    }
  }

  _verboseBubbleT = false;
  const clustered = clusterAzeotropes(results);
  console.log(`[systematicBinarySearch] Found ${clustered.length} binary azeotropes:`, clustered.map(a=>`x=[${a.x.map(v=>v.toFixed(3))}] T=${a.T_K.toFixed(1)}K conv=${a.converged}`));
  return clustered;
}

/**
 * Given a list of converged binary azeotropes, generate a high-quality initial
 * guess for a potential ternary azeotrope by averaging the liquid composition
 * and temperature of the binary ones.
 */
export function estimateTernaryGuessFromBinary(binaryAzeos: AzeotropeResult[]): { x: number[]; T_K: number } | null {
  const good = binaryAzeos.filter(a => a.converged && a.x.length === 3);
  if (good.length === 0) return null;
  const xAvg = [0, 1, 2].map(idx => good.reduce((s, a) => s + a.x[idx], 0) / good.length);
  const sum = xAvg.reduce((s, v) => s + v, 0) || 1;
  const xNorm = xAvg.map(v => v / sum);
  const Tavg = good.reduce((s, a) => s + (a.T_K || 0), 0) / good.length;
  return { x: xNorm, T_K: Tavg };
}

/**
 * Systematic interior scan to bracket and refine potential ternary azeotropes.
 *
 * The algorithm samples a coarse barycentric grid (step `gridStep`) inside the
 * composition triangle.  At each node it evaluates two independent residual
 * functions that must simultaneously vanish at a ternary azeotrope:
 *   f1 = ln(K0/K2)
 *   f2 = ln(K1/K2)
 * (Any two of the three ln(Ki) differences are sufficient.)
 *
 * For every small triangular cell whose vertices exhibit a sign change in *both*
 * residuals, the routine assumes a root lies within the cell.  It takes the
 * centroid of the three compositions as an initial guess, estimates a bubble
 * temperature at that point, and then calls the existing package-specific
 * azeotrope solver to obtain a refined solution.
 *
 * The final list is deduplicated so that closely spaced solutions are reported
 * only once.  The function purposefully keeps the grid coarse – it is meant to
 * generate robust *initial guesses* rather than sub-Kelvin precision results.
 */
export function systematicTernaryAzeotropeSearch(
  pkg: FluidPackageResidue,
  P_system_Pa: number,
  comps: CompoundData[],
  pkgParams: any,
  initialT_K: number = 350,
  gridStep: number = 0.05 // finer step for better saddle-azeotrope detection
): AzeotropeResult[] {
  console.log(`[systematicTernarySearch] pkg=${pkg}, P=${P_system_Pa.toFixed(0)}, gridStep=${gridStep}`);
  _verboseBubbleT = true;
  const results: AzeotropeResult[] = [];

  const normComp = (x:number[]):number[]=>{const s=x.reduce((a,b)=>a+b,0)||1; return x.map(v=>v/s);} // simple normaliser

  // helper dedicated to ternary interior: returns f1 & f2 using local bubbleT evaluations
  const evalPoint = (x:number[], Tguess:number):{f1:number;f2:number;T:number}|null => {
    const attempt = (Ttrial:number) => {
      switch(pkg){
        case 'wilson':{
          const r=bubbleTWilson(x,P_system_Pa,comps,pkgParams as any,Ttrial); if(!r)return null;
          const K=r.g.map((g,i)=>g*r.Ps[i]/P_system_Pa); const lnK=K.map(Math.log); return {f1:lnK[0]-lnK[2], f2:lnK[1]-lnK[2], T:r.T};}
        case 'nrtl':{
          const r=bubbleTNrtl(x,P_system_Pa,comps,pkgParams as any,Ttrial); if(!r)return null;
          const K=r.g.map((g,i)=>g*r.Ps[i]/P_system_Pa); const lnK=K.map(Math.log); return {f1:lnK[0]-lnK[2], f2:lnK[1]-lnK[2], T:r.T};}
        case 'unifac':{
          const r=bubbleTUnifac(x,P_system_Pa,comps,pkgParams as any,Ttrial); if(!r)return null;
          const K=r.g.map((g,i)=>g*r.Ps[i]/P_system_Pa); const lnK=K.map(Math.log); return {f1:lnK[0]-lnK[2], f2:lnK[1]-lnK[2], T:r.T};}
        case 'uniquac':{
          const r=bubbleQUni(x,P_system_Pa,comps,pkgParams as any,Ttrial); if(!r)return null;
          const K=r.g.map((g,i)=>g*r.Ps[i]/P_system_Pa); const lnK=K.map(Math.log); return {f1:lnK[0]-lnK[2], f2:lnK[1]-lnK[2], T:r.T};}
        case 'pr':{
          const r=solvePRBubbleT(x,P_system_Pa,comps,pkgParams as any,Ttrial); if(!r)return null;
          const lnK=r.K.map(Math.log); return {f1:lnK[0]-lnK[2], f2:lnK[1]-lnK[2], T:r.T};}
        case 'srk':{
          const r=solveSRKBubbleT(x,P_system_Pa,comps,pkgParams as any,Ttrial); if(!r)return null;
          const lnK=r.K.map(Math.log); return {f1:lnK[0]-lnK[2], f2:lnK[1]-lnK[2], T:r.T};}
        default:return null;
      }
    };
    let out=attempt(Tguess); if(out) return out;
    for(const d of [-10,-5,5,10]){ out=attempt(Tguess+d); if(out) return out; }
    return null;
  };

  const n = Math.floor(1/gridStep);
  const grid:(ReturnType<typeof evalPoint>|null)[][] = Array.from({length:n+1},()=>[]);
  // Build grid and evaluate F at each point
  for(let i=0;i<=n;i++){
    for(let j=0;j<=n-i;j++){
      const x0=i*gridStep;
      const x1=j*gridStep;
      const x2=1-x0-x1;
      const x=[x0,x1,x2];
      const evalRes=evalPoint(x,initialT_K);
      grid[i][j]=evalRes;
    }
  }

  // Helper to call package-specific ternary solver
  const runLocalSolver = (xGuess:number[],Tguess:number):AzeotropeResult|null=>{
    switch(pkg){
      case 'wilson': return findAzeotropeWilson(P_system_Pa,comps,pkgParams,xGuess,Tguess);
      case 'nrtl':   return findAzeotropeNRTL(P_system_Pa,comps,pkgParams,xGuess,Tguess);
      case 'unifac': return findAzeotropeUnifac(P_system_Pa,comps,pkgParams,xGuess,Tguess);
      case 'pr':     return findAzeotropePr(P_system_Pa,comps,pkgParams,xGuess,Tguess);
      case 'srk':    return findAzeotropeSrk(P_system_Pa,comps,pkgParams,xGuess,Tguess);
      case 'uniquac':return findAzeotropeUniquac(P_system_Pa,comps,pkgParams,xGuess,Tguess);
      default: return null;
    }
  };

  const sign = (v:number)=> v>0?1: v<0? -1:0;

  for(let i=0;i<n;i++){
    for(let j=0;j<n-i;j++){
      const pA=grid[i][j];
      const pB=grid[i+1][j]; // along x0 axis
      const pC=grid[i][j+1]; // along x1 axis
      // first triangle ABC
      if(pA&&pB&&pC){
        const s11=[sign(pA.f1),sign(pB.f1),sign(pC.f1)];
        const s21=[sign(pA.f2),sign(pB.f2),sign(pC.f2)];
        if(new Set(s11).size>1 && new Set(s21).size>1){
          const xGuess=[(i+0.33)*gridStep,(j+0.33)*gridStep,1-((i+j)+0.66)*gridStep];
          const Tguess=(pA.T+pB.T+pC.T)/3;
          const sol=runLocalSolver(normComp(xGuess),Tguess);
          if(sol) results.push(sol);
        }
      }
      // second triangle BDC with D=(i+1,j+1)
      if(j+1<=n-i-1){
        const pD=grid[i+1][j+1];
        if(pB&&pC&&pD){
          const s12=[sign(pB.f1),sign(pD.f1),sign(pC.f1)];
          const s22=[sign(pB.f2),sign(pD.f2),sign(pC.f2)];
          if(new Set(s12).size>1 && new Set(s22).size>1){
            const xGuess=[(i+0.66)*gridStep,(j+0.66)*gridStep,1-((i+j)+1.32)*gridStep];
            const Tguess=(pB.T+pC.T+pD.T)/3;
            const sol=runLocalSolver(normComp(xGuess),Tguess);
            if(sol) results.push(sol);
          }
        }
      }
    }
  }

  // Deduplicate close results
  _verboseBubbleT = false;
  const clustered = clusterAzeotropes(results);
  console.log(`[systematicTernarySearch] Found ${clustered.length} ternary azeotropes:`, clustered.map(a=>`x=[${a.x.map(v=>v.toFixed(3))}] T=${a.T_K.toFixed(1)}K conv=${a.converged}`));
  return clustered;
}

// ---------------------------------------------------------------------------
//  E N D   O F   F I L E

// ------------------------------------------------------------------------- 

// ===========================================================================
//  Generic single-point ODE evaluator (exposed for Jacobian calculations)
// ===========================================================================
export function evaluateResidueODE(
  pkg: FluidPackageResidue,
  x: number[],
  T_K: number,
  P_system_Pa: number,
  comps: CompoundData[],
  pkgParams: any,
): ResidueODEPoint {
  switch (pkg) {
    case 'wilson': {
      const r = odeWilson(x, P_system_Pa, comps, pkgParams, T_K);
      if (!r) return null;
      return { d: r.dxdxi, T_K: r.T };
    }
    case 'nrtl': {
      const r = odeNrtl(x, P_system_Pa, comps, pkgParams, T_K);
      if (!r) return null;
      return { d: r.dxdxi, T_K: r.T };
    }
    case 'unifac': {
      const r = odeUnifac(x, P_system_Pa, comps, pkgParams, T_K);
      if (!r) return null;
      return { d: r.dxdxi, T_K: r.T };
    }
    case 'pr': {
      const r = odePr(x, P_system_Pa, comps, pkgParams, T_K);
      if (!r) return null;
      return { d: r.dxdxi, T_K: r.T };
    }
    case 'srk': {
      const r = odeSrk(x, P_system_Pa, comps, pkgParams, T_K);
      if (!r) return null;
      return { d: r.dxdxi, T_K: r.T };
    }
    case 'uniquac': {
      const r = odeUni(x, P_system_Pa, comps, pkgParams, T_K);
      if (!r) return null;
      return { d: r.dxdxi, T_K: r.T };
    }
    default:
      return null;
  }
}