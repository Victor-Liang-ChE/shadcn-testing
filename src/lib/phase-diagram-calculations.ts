"use server";

import { createClient } from '@supabase/supabase-js';
import type { SupabaseClient } from '@supabase/supabase-js';

const supabaseUrl = process.env.NEXT_PUBLIC_SUPABASE_URL;
const supabaseAnonKey = process.env.NEXT_PUBLIC_SUPABASE_ANON_KEY;

if (!supabaseUrl || !supabaseAnonKey) {
    throw new Error("Supabase environment variables not set for server-side calculations.");
}
const supabase = createClient(supabaseUrl, supabaseAnonKey);

// --- 1. SHARED TYPES & INTERFACES ---
export type FluidPackageType = 'unifac' | 'nrtl' | 'pr' | 'srk' | 'uniquac' | 'wilson';
export type DiagramType = 'txy' | 'pxy' | 'xy';

export interface AntoineParams { A: number; B: number; C: number; Tmin_K: number; Tmax_K: number; Units: string; EquationNo?: number | string; }
export interface PrPureComponentParams { Tc_K: number; Pc_Pa: number; omega: number; }
export interface SrkPureComponentParams { Tc_K: number; Pc_Pa: number; omega: number; }
export interface UniquacPureComponentParams { r: number; q: number; }
export interface WilsonPureComponentParams { V_L_m3mol: number; }
export interface UnifacGroupComposition { [subgroupId: number]: number; }

export interface CompoundData {
    name: string; cas_number?: string | null; antoine: AntoineParams | null; unifacGroups?: UnifacGroupComposition | null;
    prParams?: PrPureComponentParams | null; srkParams?: SrkPureComponentParams | null; uniquacParams?: UniquacPureComponentParams | null;
    wilsonParams?: WilsonPureComponentParams | null; r_i?: number; q_i?: number;
}
export interface BubbleDewResult { comp1_feed: number; comp1_equilibrium: number; T_K?: number; P_Pa?: number; error?: string; }
export interface VleChartData { x: number[]; y: number[]; t?: number[]; p?: number[]; }
interface PrInteractionParams { k_ij: number; }
interface SrkInteractionParams { k_ij: number; }
interface NrtlInteractionParams { A12: number; A21: number; alpha: number; }
interface UniquacInteractionParams { A12: number; A21: number; }
interface WilsonInteractionParams { a12_J_mol: number; a21_J_mol: number; }
interface UnifacParameters { Rk: { [id: number]: number }; Qk: { [id: number]: number }; mainGroupMap: { [id: number]: number }; a_mk: Map<string, number>; }

// --- 2. UTILITY FUNCTIONS ---
const R_gas_const_J_molK = 8.31446261815324, R_cal_molK = 1.98720425864083, JOULES_PER_CAL = 4.184;

function calculatePsat_Pa(params: AntoineParams | null, T_kelvin: number): number {
    if (!params || T_kelvin <= 0) return NaN;
    const conv = params.Units?.toLowerCase() === 'kpa' ? 1000 : 1;
    const logP = params.A - params.B / (T_kelvin + params.C);
    const P = (params.EquationNo === 1 || params.EquationNo === '1') ? Math.pow(10, logP) : Math.exp(logP);
    return isNaN(P) || P < 0 ? NaN : P * conv;
}

function solveCubicEOS(a: number, b: number, c: number, d: number): number[] | null {
    if (Math.abs(a) < 1e-12) return null;
    const p = b / a, q = c / a, r = d / a;
    const A = q - p * p / 3;
    const B_eos = 2 * p * p * p / 27 - p * q / 3 + r;
    let D = B_eos * B_eos / 4 + A * A * A / 27;

    if (Math.abs(D) < 1e-20) D = 0;

    let roots: number[];
    if (D > 0) {
        const root = Math.cbrt(-B_eos / 2 + Math.sqrt(D)) + Math.cbrt(-B_eos / 2 - Math.sqrt(D)) - p / 3;
        roots = [root];
    } else {
        const term1 = 2 * Math.sqrt(-A / 3);
        const term2_arg = -B_eos / (2 * Math.sqrt(-A * A * A / 27));
        const clamped_arg = Math.max(-1.0, Math.min(1.0, term2_arg));
        const angle = Math.acos(clamped_arg) / 3;
        const root1 = term1 * Math.cos(angle) - p / 3;
        const root2 = term1 * Math.cos(angle + 2 * Math.PI / 3) - p / 3;
        const root3 = term1 * Math.cos(angle + 4 * Math.PI / 3) - p / 3;
        roots = [root1, root2, root3];
    }

    const final = roots.filter(z => z > 1e-9).sort((x, y) => x - y);
    return final.length > 0 ? final : null;
}

function solveAntoineBoilingPoint(antoineParams: AntoineParams | null, P_target_Pa: number): number | null {
    if (!antoineParams) return null;
    try {
        const P_units = P_target_Pa / (antoineParams.Units?.toLowerCase() === 'kpa' ? 1000 : 1);
        const logP = (antoineParams.EquationNo === 1 || antoineParams.EquationNo === '1') ? Math.log10(P_units) : Math.log(P_units);
        if (antoineParams.A - logP === 0) return null;
        const T_K = antoineParams.B / (antoineParams.A - logP) - antoineParams.C;
        return (T_K > 0 && T_K < 1000) ? T_K : null;
    } catch { return null; }
}

// --- 3. MODEL-SPECIFIC CALCULATIONS ---
function calculatePrFugacityCoefficients(c: CompoundData[], phase_comp: number[], T: number, P: number, k: number, ph: 'liquid' | 'vapor'): [number, number] | null {
    if (!c[0].prParams || !c[1].prParams) return null;
    
    const p = c.map(d => {
        const { Tc_K, Pc_Pa, omega } = d.prParams!;
        const kap = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
        const alf = Math.pow(1 + kap * (1 - Math.sqrt(T / Tc_K)), 2);
        return {
            ac: 0.45724 * R_gas_const_J_molK ** 2 * Tc_K ** 2 / Pc_Pa * alf,
            b: 0.07780 * R_gas_const_J_molK * Tc_K / Pc_Pa
        };
    });

    const a12 = (1 - k) * Math.sqrt(p[0].ac * p[1].ac);
    
    const am = phase_comp[0] ** 2 * p[0].ac + phase_comp[1] ** 2 * p[1].ac + 2 * phase_comp[0] * phase_comp[1] * a12;
    const bm = phase_comp[0] * p[0].b + phase_comp[1] * p[1].b;
    
    const A = am * P / (R_gas_const_J_molK * T) ** 2;
    const B = bm * P / (R_gas_const_J_molK * T);
    
    const Z_r = solveCubicEOS(1, -(1 - B), A - 3 * B ** 2 - 2 * B, -(A * B - B ** 2 - B ** 3));
    if (!Z_r) return null;
    
    const Z = ph === 'liquid' ? Z_r[0] : Z_r[Z_r.length - 1];
    if (Z <= B) return null;

    if (bm === 0 || am === 0 || Z - B <= 0) return null;

    const term_ln_denom = (2 * Math.sqrt(2) * B);
    if (Math.abs(term_ln_denom) < 1e-12) return null;

    const term_ln_log_arg_num = Z + (1 + Math.sqrt(2)) * B;
    const term_ln_log_arg_den = Z + (1 - Math.sqrt(2)) * B;
    if (term_ln_log_arg_den <= 0) return null;

    const term_ln = A / term_ln_denom * Math.log(term_ln_log_arg_num / term_ln_log_arg_den);
    
    const t1 = (2 * (phase_comp[0] * p[0].ac + phase_comp[1] * a12) / am) - (p[0].b / bm);
    const t2 = (2 * (phase_comp[1] * p[1].ac + phase_comp[0] * a12) / am) - (p[1].b / bm);
    
    const phi1 = Math.exp((p[0].b / bm) * (Z - 1) - Math.log(Z - B) - term_ln * t1);
    const phi2 = Math.exp((p[1].b / bm) * (Z - 1) - Math.log(Z - B) - term_ln * t2);

    return [phi1, phi2];
}

function calculateSrkFugacityCoefficients(c: CompoundData[], phase_comp: number[], T: number, P: number, k: number, ph: 'liquid' | 'vapor'): [number, number] | null {
    if(!c[0].srkParams || !c[1].srkParams) return null;

    const p = c.map(d => {
        const {Tc_K, Pc_Pa, omega} = d.srkParams!;
        const m = 0.480 + 1.574 * omega - 0.176 * omega * omega;
        const a = Math.pow(1 + m * (1 - Math.sqrt(T/Tc_K)), 2);
        return {
            ac: 0.42748 * (R_gas_const_J_molK * Tc_K)**2 / Pc_Pa * a,
            b: 0.08664 * R_gas_const_J_molK * Tc_K / Pc_Pa
        };
    });

    const a12 = (1 - k) * Math.sqrt(p[0].ac * p[1].ac);
    const am = phase_comp[0]**2 * p[0].ac + phase_comp[1]**2 * p[1].ac + 2 * phase_comp[0] * phase_comp[1] * a12;
    const bm = phase_comp[0] * p[0].b + phase_comp[1] * p[1].b;

    const A = am * P / (R_gas_const_J_molK * T)**2;
    const B = bm * P / (R_gas_const_J_molK * T);
    
    const Z_r = solveCubicEOS(1, -1, A - B - B**2, -A * B);
    if(!Z_r) return null;

    const Z = ph === 'liquid' ? Z_r[0] : Z_r[Z_r.length - 1];
    if (Z <= B || Z <= 0) return null;
    if (bm === 0 || am === 0) return null;

    const term_A_div_B = A / B;
    const term_ln_ZB = Math.log(1 + B / Z);

    const term_sum_1 = 2 * (phase_comp[0] * p[0].ac + phase_comp[1] * a12) / am;
    const term_sum_2 = 2 * (phase_comp[1] * p[1].ac + phase_comp[0] * a12) / am;

    const phi1 = Math.exp((p[0].b / bm) * (Z - 1) - Math.log(Z - B) - term_A_div_B * (term_sum_1 - p[0].b / bm) * term_ln_ZB);
    const phi2 = Math.exp((p[1].b / bm) * (Z - 1) - Math.log(Z - B) - term_A_div_B * (term_sum_2 - p[1].b / bm) * term_ln_ZB);

    return [phi1, phi2];
}

function calculateActivityGamma(comps: CompoundData[], x: number[], T: number, pkg: FluidPackageType, p: any): [number, number]|null {
    switch(pkg){
        case 'nrtl': { const t12=p.A12/(R_cal_molK*T),t21=p.A21/(R_cal_molK*T),G12=Math.exp(-p.alpha*t12),G21=Math.exp(-p.alpha*t21);const d1=x[0]+x[1]*G21,d2=x[1]+x[0]*G12;if(d1===0||d2===0)return null;return[Math.exp(x[1]**2*(t21*(G21/d1)**2+t12*G12/d2**2)),Math.exp(x[0]**2*(t12*(G12/d2)**2+t21*G21/d1**2))]; }
        case 'wilson': { if(!comps[0].wilsonParams||!comps[1].wilsonParams)return null;const V1=comps[0].wilsonParams.V_L_m3mol,V2=comps[1].wilsonParams.V_L_m3mol;const L12=(V2/V1)*Math.exp(-p.a12_J_mol/(R_gas_const_J_molK*T)),L21=(V1/V2)*Math.exp(-p.a21_J_mol/(R_gas_const_J_molK*T));return[Math.exp(-Math.log(x[0]+L12*x[1])+x[1]*(L12/(x[0]+L12*x[1])-L21/(x[1]+L21*x[0]))),Math.exp(-Math.log(x[1]+L21*x[0])-x[0]*(L12/(x[0]+L12*x[1])-L21/(x[1]+L21*x[0])))];}
        case 'uniquac': { if(!comps[0].uniquacParams||!comps[1].uniquacParams)return null;const r=[comps[0].uniquacParams.r,comps[1].uniquacParams.r],q=[comps[0].uniquacParams.q,comps[1].uniquacParams.q];const l=r.map((ri,i)=>5*(ri-q[i])-(ri-1));const sxr=x[0]*r[0]+x[1]*r[1],sxq=x[0]*q[0]+x[1]*q[1];if(sxr<1e-9||sxq<1e-9)return null;const Phi=x.map((xi,i)=>xi*r[i]/sxr),Theta=x.map((xi,i)=>xi*q[i]/sxq);const sxl=x[0]*l[0]+x[1]*l[1];const lnGC=x.map((xi,i)=>(xi<1e-9||Phi[i]<1e-9||Theta[i]<1e-9)?0:Math.log(Phi[i]/xi)+5*q[i]*Math.log(Theta[i]/Phi[i])+l[i]-(r[i]/sxr)*sxl);const t12=Math.exp(-p.A12/T),t21=Math.exp(-p.A21/T);const lnGR=[q[0]*(1-Math.log(Theta[0]+Theta[1]*t21)-Theta[0]/(Theta[0]+Theta[1]*t21)-Theta[1]*t12/(Theta[1]+Theta[0]*t12)),q[1]*(1-Math.log(Theta[1]+Theta[0]*t12)-Theta[1]/(Theta[1]+Theta[0]*t12)-Theta[0]*t21/(Theta[0]+Theta[1]*t21))];return[Math.exp(lnGC[0]+lnGR[0]),Math.exp(lnGC[1]+lnGR[1])];}
        case 'unifac': { if(x[0]>=1||x[1]>=1)return[1,1];for(const c of comps){if(!c.r_i||!c.q_i){c.r_i=0;c.q_i=0;if(!c.unifacGroups)return null;for(const[s,k]of Object.entries(c.unifacGroups)){const sgId=parseInt(s);if(!p.Rk[sgId]||!p.Qk[sgId])return null;c.r_i+=k*p.Rk[sgId];c.q_i+=k*p.Qk[sgId];}}}const r=comps.map(c=>c.r_i!),q=comps.map(c=>c.q_i!),sxr=x[0]*r[0]+x[1]*r[1],sxq=x[0]*q[0]+x[1]*q[1];if(sxr<1e-9||sxq<1e-9)return null;const Phi=x.map((xi,i)=>xi*r[i]/sxr),Theta=x.map((xi,i)=>xi*q[i]/sxq),l=r.map((ri,i)=>5*(ri-q[i])-(ri-1)),sxl=x[0]*l[0]+x[1]*l[1];const lnGC=x.map((xi,i)=>(xi<1e-9||Phi[i]<1e-9||Theta[i]<1e-9)?0:Math.log(Phi[i]/xi)+5*q[i]*Math.log(Theta[i]/Phi[i])+l[i]-(r[i]/sxr)*sxl);const sg=Array.from(new Set(comps.flatMap(c=>Object.keys(c.unifacGroups||{})).map(Number))),v=comps.map(c=>sg.map(id=>c.unifacGroups?.[id]||0));const glg=(cx:number[])=>{let svk=cx.reduce((s,xi,i)=>s+v[i].reduce((ss,vk)=>ss+xi*vk,0),0);if(svk<1e-9)return Array(sg.length).fill(0);const Xm=sg.map((_,k)=>cx.reduce((s,xi,i)=>s+xi*v[i][k],0)/svk);let sXQ=Xm.reduce((s,Xk,k)=>s+Xk*p.Qk[sg[k]],0);if(sXQ<1e-9)return Array(sg.length).fill(0);const Tm=Xm.map((Xk,k)=>Xk*p.Qk[sg[k]]/sXQ),psi=sg.map((_,m)=>sg.map((_,n)=>Math.exp(-(p.a_mk.get(`${p.mainGroupMap[sg[m]]}-${p.mainGroupMap[sg[n]]}`)??0)/T)));return sg.map((_,k)=>{const s1=Tm.reduce((s,tm,m)=>s+tm*psi[m][k],0),s2=Tm.reduce((s,tm,m)=>{const s3=Tm.reduce((ss,tn,n)=>ss+tn*psi[n][m],0);return s+(s3>1e-9?(tm*psi[k][m])/s3:0)},0);if(s1<=0)return NaN;return p.Qk[sg[k]]*(1-Math.log(s1)-s2)})};const lgm=glg(x),lgp1=glg([1,0]),lgp2=glg([0,1]);if(lgm.some(isNaN)||lgp1.some(isNaN)||lgp2.some(isNaN))return null;let lnGR1=0,lnGR2=0;sg.forEach((_,k)=>{lnGR1+=v[0][k]*(lgm[k]-lgp1[k]);lnGR2+=v[1][k]*(lgm[k]-lgp2[k])});return[Math.exp(lnGC[0]+lnGR1),Math.exp(lnGC[1]+lnGR2)];}
        default: return null;
    }
}

// --- 4. BUBBLE POINT SOLVERS ---
function calculateBubbleTemperatureActivity(comps: CompoundData[], x1: number, P: number, initialT: number, gammaFunc: (x:number[], T:number) => [number,number]|null): BubbleDewResult {
    let T = initialT; const x = [x1, 1-x1];
    for (let i=0; i<50; i++) {
        const psat = [calculatePsat_Pa(comps[0].antoine, T), calculatePsat_Pa(comps[1].antoine, T)];
        if (psat.some(isNaN)) return { comp1_feed: x1, comp1_equilibrium: NaN, error: "Psat calc failed" };
        const gamma = gammaFunc(x, T);
        if (!gamma) return { comp1_feed: x1, comp1_equilibrium: NaN, error: "Gamma calc failed" };
        const P_calc = x[0]*gamma[0]*psat[0] + x[1]*gamma[1]*psat[1];
        if (isNaN(P_calc)) return { comp1_feed: x1, comp1_equilibrium: NaN, error: "P_calc is NaN"};
        if (Math.abs(P_calc - P) < 1e-3 * P) return { comp1_feed: x1, comp1_equilibrium: x[0]*gamma[0]*psat[0]/P_calc, T_K: T, P_Pa: P };
        T -= (P_calc - P) / (P*0.1);
        if (T <= 0) return { comp1_feed: x1, comp1_equilibrium: NaN, error: "T out of bounds" };
    }
    return { comp1_feed: x1, comp1_equilibrium: NaN, error: "Bubble-T max iterations" };
}

function calculateBubblePressureActivity(comps: CompoundData[], x1: number, T: number, gammaFunc: (x:number[], T:number) => [number,number]|null): BubbleDewResult {
    const x = [x1, 1-x1];
    const psat = [calculatePsat_Pa(comps[0].antoine, T), calculatePsat_Pa(comps[1].antoine, T)];
    if (psat.some(isNaN)) return { comp1_feed: x1, comp1_equilibrium: NaN, error: "Psat calc failed" };
    const gamma = gammaFunc(x, T);
    if (!gamma) return { comp1_feed: x1, comp1_equilibrium: NaN, error: "Gamma calc failed" };
    const P_calc = x[0]*gamma[0]*psat[0] + x[1]*gamma[1]*psat[1];
    if (isNaN(P_calc) || P_calc <= 0) return { comp1_feed: x1, comp1_equilibrium: NaN, error: "Invalid calculated P" };
    return { comp1_feed: x1, comp1_equilibrium: x[0]*gamma[0]*psat[0]/P_calc, T_K: T, P_Pa: P_calc };
}

function calculateBubbleTemperatureEos(comps: CompoundData[], x1: number, P: number, initialT: number, eosFugacityFunc: (c:CompoundData[],phase_comp:number[],T:number,P:number,p:any,ph:'liquid'|'vapor')=>[number,number]|null, params: any): BubbleDewResult {
    let T = initialT; const x = [x1, 1-x1]; let y = [...x];
    for (let i=0; i<100; i++) {
        const phi_L = eosFugacityFunc(comps, x, T, P, params, 'liquid');
        if (!phi_L) return { comp1_feed: x1, comp1_equilibrium: NaN, error: "Fugacity calculation failed for liquid phase." };
        const phi_V = eosFugacityFunc(comps, y, T, P, params, 'vapor');
        if (!phi_V) return { comp1_feed: x1, comp1_equilibrium: NaN, error: "Fugacity calculation failed for vapor phase." };
        const K = [phi_L[0]/phi_V[0], phi_L[1]/phi_V[1]];
        const sum_Kx = K[0]*x[0] + K[1]*x[1];
        if (Math.abs(sum_Kx - 1) < 1e-5) return { comp1_feed: x1, comp1_equilibrium: y[0], T_K: T, P_Pa: P };
        y = [K[0]*x[0]/sum_Kx, K[1]*x[1]/sum_Kx];
        T *= (1 + (1 - sum_Kx) * 0.1);
        if (T <= 0) return { comp1_feed: x1, comp1_equilibrium: NaN, error: "T out of bounds" };
    }
    return { comp1_feed: x1, comp1_equilibrium: NaN, error: "Bubble-T EOS max iterations" };
}

function calculateBubblePressureEos(comps: CompoundData[], x1: number, T: number, initialP: number, eosFugacityFunc: (c:CompoundData[],phase_comp:number[],T:number,P:number,p:any,ph:'liquid'|'vapor')=>[number,number]|null, params: any): BubbleDewResult {
    let P = initialP; const x = [x1, 1-x1]; let y = [...x];
    for(let i=0; i<50; i++) {
        const phi_L = eosFugacityFunc(comps, x, T, P, params, 'liquid');
        if (!phi_L) { P *= 0.9; continue; }
        const phi_V = eosFugacityFunc(comps, y, T, P, params, 'vapor');
        if (!phi_V) { P *= 1.1; continue; }
        const K = [phi_L[0]/phi_V[0], phi_L[1]/phi_V[1]];
        const sum_Kx = K[0]*x[0] + K[1]*x[1];
        y = [K[0]*x[0]/sum_Kx, K[1]*x[1]/sum_Kx];
        const P_new = P * sum_Kx;
        if (Math.abs(P_new - P) < 1e-4 * P) return { comp1_feed: x1, comp1_equilibrium: y[0], T_K: T, P_Pa: P_new };
        P = P_new;
        if (P <= 0) return { comp1_feed: x1, comp1_equilibrium: NaN, error: "P out of bounds" };
    }
    return { comp1_feed: x1, comp1_equilibrium: NaN, error: "Bubble-P EOS max iterations" };
}


// --- 5. DATA FETCHING ---
async function fetchInteractionParams(pkg: FluidPackageType, c1: CompoundData, c2: CompoundData): Promise<any> {
    if (pkg === 'unifac') {
        const allSgIds = Array.from(new Set([...Object.keys(c1.unifacGroups||{}), ...Object.keys(c2.unifacGroups||{})].map(Number)));
        if (allSgIds.length === 0) throw new Error("UNIFAC model requires UNIFAC groups for both components.");
        const { data: groupData, error: groupError } = await supabase.from('UNIFAC - Rk and Qk').select('"Subgroup #", "Main Group #", Rk, Qk').in('"Subgroup #"', allSgIds);
        if (groupError) throw groupError;
        if (!groupData || groupData.length < allSgIds.length) throw new Error("Missing UNIFAC Rk/Qk data for some subgroups.");
        const Rk: UnifacParameters['Rk'] = {}, Qk: UnifacParameters['Qk'] = {}, mainGroupMap: UnifacParameters['mainGroupMap'] = {};
        const mainGroupIds = new Set<number>();
        groupData.forEach(g => { const sgId=g["Subgroup #"]; const mgId=g["Main Group #"]; if (sgId!=null && mgId!=null && g.Rk!=null && g.Qk!=null) { Rk[sgId]=g.Rk; Qk[sgId]=g.Qk; mainGroupMap[sgId]=mgId; mainGroupIds.add(mgId); } });
        const { data: intData, error: intError } = await supabase.from('UNIFAC - a(ij)').select('i,j,"a(ij)"').in('i', Array.from(mainGroupIds)).in('j', Array.from(mainGroupIds));
        if (intError) throw intError;
        const a_mk = new Map<string, number>();
        intData?.forEach(i => { if(i.i != null && i.j != null && i["a(ij)"] != null) a_mk.set(`${i.i}-${i.j}`, i["a(ij)"]) });
        return { Rk, Qk, mainGroupMap, a_mk };
    }

    const cas1 = c1.cas_number, cas2 = c2.cas_number;
    if (!cas1 || !cas2 || cas1 === cas2) return (pkg === 'pr' || pkg === 'srk') ? { k_ij: 0 } : {};
    
    let table: string;
    switch(pkg) {
        case 'pr': table = 'peng-robinson parameters'; break;
        case 'srk': table = 'srk parameters'; break;
        case 'nrtl': table = 'nrtl parameters'; break;
        case 'uniquac': table = 'uniquac parameters'; break;
        case 'wilson': table = 'wilson parameters'; break;
        default: throw new Error(`Invalid package ${pkg}`);
    }
    const { data, error } = await supabase.from(table).select('*').or(`and(CASN1.eq.${cas1},CASN2.eq.${cas2}),and(CASN1.eq.${cas2},CASN2.eq.${cas1})`).limit(1);
    
    if (error) throw new Error(`Supabase query error for ${pkg}: ${error.message}`);
    if (!data || data.length === 0) return (pkg === 'pr' || pkg === 'srk') ? { k_ij: 0 } : {};

    const dbRow = data[0] as any, isSwapped = dbRow.CASN1 !== cas1;
    switch(pkg) {
        case 'pr': case 'srk': return { k_ij: dbRow.k12 ?? 0 };
        case 'nrtl': return { alpha: dbRow.alpha12 ?? 0.3, A12: (isSwapped ? dbRow.A21:dbRow.A12)??0, A21: (isSwapped ? dbRow.A12:dbRow.A21)??0 };
        case 'uniquac': return { A12: (isSwapped ? dbRow.A21:dbRow.A12)??0, A21: (isSwapped ? dbRow.A12:dbRow.A21)??0 };
        case 'wilson': return { a12_J_mol: ((isSwapped?dbRow.A21:dbRow.A12)??0)*JOULES_PER_CAL, a21_J_mol: ((isSwapped?dbRow.A12:dbRow.A21)??0)*JOULES_PER_CAL };
    }
}

export async function fetchCompoundData(compoundName: string): Promise<CompoundData> {
    const { data: compoundDbData, error: compoundError } = await supabase.from('compounds').select('id, name, cas_number').ilike('name', compoundName).limit(1).single();
    if (compoundError) throw new Error(`Compound '${compoundName}' not found.`);
    const { id: compoundId, name: foundName, cas_number: casNumber } = compoundDbData;
    const { data: propsData, error: propsError } = await supabase.from('compound_properties').select('properties').eq('compound_id', compoundId).limit(1).single();
    if (propsError) throw new Error(`No properties found for ${foundName}.`);
    
    const p = propsData.properties as any;
    const antoine = p.Antoine || p.AntoineVaporPressure;
    const Tc = p["Critical temperature"]?.value, Pc_val = p["Critical pressure"]?.value, omega = p["Acentric factor"]?.value;
    const Pc_units = p["Critical pressure"]?.units?.toLowerCase() ?? '', Pc_Pa = Pc_val * (Pc_units.includes('kpa') ? 1000 : Pc_units.includes('bar') ? 100000 : 1);
    const vL_obj = p["Liquid molar volume"] || p["Molar volume"] || p["Wilson volume"], vL_val = vL_obj?.value, vL_units = vL_obj?.units?.toLowerCase() ?? 'cm3/mol';
    let vL_m3mol = null;
    if (vL_val) { if (vL_units.includes('cm3')) vL_m3mol = vL_val * 1e-6; else if (vL_units.includes('m3/kmol')) vL_m3mol = vL_val / 1000; else if (vL_units.includes('m3')) vL_m3mol = vL_val; }
    return {
        name: foundName, cas_number: casNumber,
        antoine: antoine ? { A: antoine.A, B: antoine.B, C: antoine.C, Tmin_K: antoine.Tmin ?? 0, Tmax_K: antoine.Tmax ?? 10000, Units: antoine.units || 'Pa', EquationNo: antoine.eqno } : null,
        unifacGroups: p.elements_composition?.UNIFAC,
        prParams: (Tc&&Pc_val&&omega) ? { Tc_K: Tc, Pc_Pa, omega } : null,
        srkParams: (Tc&&Pc_val&&omega) ? { Tc_K: Tc, Pc_Pa, omega } : null,
        uniquacParams: (p["UNIQUAC r"]?.value && p["UNIQUAC q"]?.value) ? { r: p["UNIQUAC r"].value, q: p["UNIQUAC q"].value } : null,
        wilsonParams: vL_m3mol ? { V_L_m3mol: vL_m3mol } : null,
    };
}

// --- 6. MASTER CALCULATION FUNCTION ---
export async function calculateVleDiagramData(
    components: CompoundData[], diagramType: DiagramType,
    fluidPackage: FluidPackageType, fixedConditionValue: number, isTempFixed: boolean
): Promise<VleChartData> {
    const params = await fetchInteractionParams(fluidPackage, components[0], components[1]);
    
    let finalResults: BubbleDewResult[] = [];

    if (!isTempFixed) { // This is a T-x-y diagram
        const P = fixedConditionValue * 1e5;

        if (fluidPackage === 'pr' || fluidPackage === 'srk') {
            const eosFugacityFunc = fluidPackage === 'pr' ? calculatePrFugacityCoefficients : calculateSrkFugacityCoefficients;
            const k_ij = params.k_ij;

            const solveEosBoilingPoint = (comp: CompoundData, fallback_comp: CompoundData): number | null => {
                const antoine_T = solveAntoineBoilingPoint(comp.antoine, P);
                if (!antoine_T) return 373.15;

                let T_low = antoine_T - 100, T_high = antoine_T + 100;
                
                const calc_K_minus_1 = (T: number): number | null => {
                    if (T <= 0) return null;
                    const phi_L = eosFugacityFunc([comp, comp], [1, 0], T, P, 0, 'liquid');
                    const phi_V = eosFugacityFunc([comp, comp], [1, 0], T, P, 0, 'vapor');
                    if (!phi_L || !phi_V) return null;
                    return (phi_L[0] / phi_V[0]) - 1;
                };

                let f_low = calc_K_minus_1(T_low);
                let f_high = calc_K_minus_1(T_high);

                if (f_low === null || f_high === null || f_low * f_high > 0) {
                     const result = calculateBubbleTemperatureEos([comp, fallback_comp], 0.999999, P, antoine_T, eosFugacityFunc, k_ij);
                     return result.T_K ?? antoine_T;
                }
                
                for (let i = 0; i < 30; i++) {
                    const T_mid = (T_low + T_high) / 2;
                    if(T_high - T_low < 0.01) return T_mid;
                    const f_mid = calc_K_minus_1(T_mid);
                    if (f_mid === null) return T_mid;
                    if (Math.abs(f_mid) < 1e-5) return T_mid;
                    if (f_mid * f_low < 0) { T_high = T_mid; } 
                    else { T_low = T_mid; f_low = f_mid; }
                }
                return (T_low + T_high) / 2;
            };

            const Tbp1 = solveEosBoilingPoint(components[0], components[1]);
            const Tbp2 = solveEosBoilingPoint(components[1], components[0]);

            const calculationRunner = (x1: number) => {
                const initialT = (Tbp1 && Tbp2) ? x1 * Tbp1 + (1 - x1) * Tbp2 : 373.15;
                return calculateBubbleTemperatureEos(components, x1, P, initialT, eosFugacityFunc, k_ij);
            };

            const x_values_inner = Array.from({ length: 49 }, (_, i) => parseFloat(((i + 1) * 0.02).toFixed(4)));
            const innerResults = x_values_inner.map(calculationRunner);
            
            if (Tbp2) { finalResults.push({ comp1_feed: 0, comp1_equilibrium: 0, T_K: Tbp2, P_Pa: P }); }
            finalResults.push(...innerResults);
            if (Tbp1) { finalResults.push({ comp1_feed: 1, comp1_equilibrium: 1, T_K: Tbp1, P_Pa: P }); }
        
        } else { // Activity coefficient models
            const Tbp1 = solveAntoineBoilingPoint(components[0].antoine, P);
            const Tbp2 = solveAntoineBoilingPoint(components[1].antoine, P);
            
            const calculationRunner = (x1: number) => {
                const initialT = (Tbp1 && Tbp2) ? x1*Tbp1 + (1-x1)*Tbp2 : 373.15;
                return calculateBubbleTemperatureActivity(components, x1, P, initialT, (x, T) => calculateActivityGamma(components,x,T,fluidPackage,params));
            };
            const x_values_inner = Array.from({ length: 49 }, (_, i) => parseFloat(((i + 1) * 0.02).toFixed(4)));
            const innerResults = x_values_inner.map(calculationRunner);
            
            if (Tbp2) { finalResults.push({ comp1_feed: 0, comp1_equilibrium: 0, T_K: Tbp2, P_Pa: P }); }
            finalResults.push(...innerResults);
            if (Tbp1) { finalResults.push({ comp1_feed: 1, comp1_equilibrium: 1, T_K: Tbp1, P_Pa: P }); }
        }

    } else { // This is a P-x-y or isothermal x-y diagram
        const T = fixedConditionValue + 273.15;
        const x_values = Array.from({ length: 51 }, (_, i) => i * 0.02);
        const calculationRunner = (x1: number) => {
            const initialP = calculatePsat_Pa(components[0].antoine, T) * x1 + calculatePsat_Pa(components[1].antoine, T) * (1-x1) || 101325;
            switch (fluidPackage) {
                case 'pr': return calculateBubblePressureEos(components, x1, T, initialP, calculatePrFugacityCoefficients, params.k_ij);
                case 'srk': return calculateBubblePressureEos(components, x1, T, initialP, calculateSrkFugacityCoefficients, params.k_ij);
                default: return calculateBubblePressureActivity(components, x1, T, (x, T) => calculateActivityGamma(components,x,T,fluidPackage,params));
            }
        };
        finalResults = x_values.map(calculationRunner);
    }

    const validResults = finalResults.filter(r => r && !r.error);
    const chartData: VleChartData = {
        x: validResults.map(p => p.comp1_feed),
        y: validResults.map(p => p.comp1_equilibrium),
        t: diagramType === 'txy' ? validResults.map(p => p.T_K).filter((v): v is number => v !== undefined) : undefined,
        p: diagramType === 'pxy' ? validResults.map(p => p.P_Pa).filter((v): v is number => v !== undefined) : undefined,
    };
    return chartData;
}